"""
Marker sets and the cell-type call.

A marker set is a cell type and a list of genes. Sets come from the optional
`markers` input (markers.csv, long `cell_type,gene` or wide `cell_type,genes`
with `;`-separated genes) or from the built-in table (builtin_markers.json).
Gene symbols are matched to the panel case-insensitively, so a human table
matches mouse symbols (CD3E ~ Cd3e).

Each set is scored per observation with scanpy.tl.score_genes (mean of the
set's genes minus the mean of a matched control set, Satija et al. 2015), the
score is z-scored across all observations, and the call is the set with the
highest z-score: per cluster (the mean z over the cluster's observations) or
per cell. Below MIN_SCORE the call is `Unassigned`.
"""
from __future__ import annotations

import csv
import json
import logging
import re
from dataclasses import dataclass, field
from pathlib import Path

log = logging.getLogger("st-celltype")

BUILTIN_PATH = Path(__file__).resolve().parent / "builtin_markers.json"
UNASSIGNED = "Unassigned"
SAFE_RE = re.compile(r"[^A-Za-z0-9_+./ -]+")


@dataclass
class MarkerSet:
    cell_type: str
    genes: list[str]                              # as given
    present: list[str] = field(default_factory=list)   # var_names actually used
    missing: list[str] = field(default_factory=list)


def _clean_type(text: str) -> str:
    text = SAFE_RE.sub("", (text or "").strip()).strip()
    return re.sub(r"\s+", " ", text)


def _split_genes(text: str) -> list[str]:
    return [g.strip() for g in re.split(r"[;,|\s]+", text or "") if g.strip()]


def builtin_sets() -> tuple[list[MarkerSet], str]:
    with open(BUILTIN_PATH) as fh:
        data = json.load(fh)
    return [MarkerSet(s["cell_type"], list(s["genes"])) for s in data["sets"]], data.get("description", "")


def read_markers_csv(path: Path) -> list[MarkerSet]:
    """markers.csv in either shape: long (cell_type, gene; one gene per row) or wide
    (cell_type, genes; genes separated by ';'). Column names are matched case-insensitively;
    a set's genes are deduplicated in order; empty cell types or genes are ignored."""
    with open(path, newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        if not reader.fieldnames:
            raise ValueError(f"{path.name}: empty file")
        cols = {c.strip().lower(): c for c in reader.fieldnames if c}
        type_col = cols.get("cell_type") or cols.get("celltype") or cols.get("label") or cols.get("cell type")
        gene_col = cols.get("gene") or cols.get("genes") or cols.get("marker") or cols.get("markers") or cols.get("symbol")
        if type_col is None or gene_col is None:
            raise ValueError(f"{path.name}: needs a cell_type column and a gene (long) or genes (wide) column; "
                             f"found {list(cols)}")
        order: list[str] = []
        genes: dict[str, list[str]] = {}
        for row in reader:
            ct = _clean_type(row.get(type_col) or "")
            if not ct:
                continue
            if ct not in genes:
                order.append(ct)
                genes[ct] = []
            for g in _split_genes(row.get(gene_col) or ""):
                if g not in genes[ct]:
                    genes[ct].append(g)
    sets = [MarkerSet(ct, genes[ct]) for ct in order if genes[ct]]
    if not sets:
        raise ValueError(f"{path.name}: no marker sets with genes")
    return sets


def load_sets(marker_set: str, markers_csv: Path | None) -> tuple[list[MarkerSet], str]:
    """The marker sets to score and where they came from: a markers.csv when one is linked
    (it wins over MARKER_SET=builtin), the built-in table, or none (clusters only)."""
    if marker_set == "none":
        return [], "none"
    if markers_csv is not None:
        return read_markers_csv(markers_csv), f"markers.csv ({markers_csv.name})"
    if marker_set == "builtin":
        sets, _ = builtin_sets()
        return sets, "builtin"
    raise ValueError(f"MARKER_SET must be builtin or none, got {marker_set!r}")


def resolve(sets: list[MarkerSet], var_names, min_genes_present: int) -> tuple[list[MarkerSet], list[MarkerSet]]:
    """Match each set's genes to the panel case-insensitively (first panel symbol per upper-cased
    key). Returns (sets with >= min_genes_present present genes, the skipped ones)."""
    lookup: dict[str, str] = {}
    for name in var_names:
        key = str(name).upper()
        lookup.setdefault(key, str(name))
    kept, skipped = [], []
    for s in sets:
        present, missing, seen = [], [], set()
        for g in s.genes:
            hit = lookup.get(g.upper())
            if hit is None:
                missing.append(g)
            elif hit not in seen:
                seen.add(hit)
                present.append(hit)
        s.present, s.missing = present, missing
        (kept if len(present) >= min_genes_present else skipped).append(s)
    return kept, skipped


def score_sets(adata, sets: list[MarkerSet], seed: int):
    """Raw scanpy.tl.score_genes scores, one column per set (DataFrame indexed like adata.obs)."""
    import numpy as np
    import pandas as pd
    import scanpy as sc

    n_genes = adata.n_vars
    # score_genes draws ctrl_size control genes per expression bin from the whole panel; a small
    # panel can't give 50 per bin, so the control set is bounded by what the panel has.
    ctrl_size = int(min(50, max(5, n_genes // 10)))
    n_bins = int(min(25, max(2, n_genes // 10)))
    scores = pd.DataFrame(index=adata.obs_names)
    for s in sets:
        name = f"score_{s.cell_type}"
        sc.tl.score_genes(adata, s.present, score_name=name, ctrl_size=ctrl_size, n_bins=n_bins,
                          random_state=seed, use_raw=False)
        scores[s.cell_type] = np.asarray(adata.obs[name], dtype=np.float64)
        del adata.obs[name]
    return scores


def zscore(scores):
    """Each column standardised over all observations (a constant column stays 0)."""
    import numpy as np

    z = scores.copy()
    for c in z.columns:
        v = z[c].to_numpy(dtype=np.float64)
        sd = float(v.std())
        z[c] = (v - v.mean()) / sd if sd > 0 else 0.0
    return z


def _best_two(values, names):
    """(best name, best value, second name, second value) for a 1-D array of set values."""
    import numpy as np

    v = np.asarray(values, dtype=np.float64)
    if v.size == 0:
        return UNASSIGNED, float("nan"), "", float("nan")
    order = np.argsort(-v, kind="stable")
    best = (names[order[0]], float(v[order[0]]))
    second = (names[order[1]], float(v[order[1]])) if v.size > 1 else ("", float("nan"))
    return best[0], best[1], second[0], second[1]


def annotate_by_cluster(z, clusters, min_score: float):
    """One call per cluster from the mean z-score of its observations.
    Returns (cell_type per observation (list), annotation rows per cluster)."""
    import numpy as np
    import pandas as pd

    names = list(z.columns)
    clusters = pd.Series(np.asarray(clusters, dtype=object), index=z.index)
    rows, call_of = [], {}
    for cl in sorted(clusters.unique()):
        mask = (clusters == cl).to_numpy()
        means = z.loc[mask].mean(axis=0).to_numpy() if names else np.array([])
        best, best_v, second, second_v = _best_two(means, names)
        call = best if (names and best_v >= min_score) else UNASSIGNED
        call_of[cl] = call
        rows.append({"cluster": cl, "cell_type": call, "score": best_v, "n_obs": int(mask.sum()),
                     "best_set": best if names else "", "second_best": second, "second_score": second_v})
    return [call_of[c] for c in clusters], rows


def annotate_by_cell(z, clusters, min_score: float):
    """One call per observation; the cluster table then reports each cluster's majority call."""
    import numpy as np
    import pandas as pd

    names = list(z.columns)
    n = len(z.index)
    if not names:
        calls = np.array([UNASSIGNED] * n, dtype=object)
        best_v = np.full(n, np.nan)
    else:
        values = z.to_numpy(dtype=np.float64)
        idx = values.argmax(axis=1)
        best_v = values[np.arange(n), idx]
        calls = np.array(names, dtype=object)[idx]
        calls[best_v < min_score] = UNASSIGNED
    clusters = pd.Series(np.asarray(clusters, dtype=object), index=z.index)
    rows = []
    for cl in sorted(clusters.unique()):
        mask = (clusters == cl).to_numpy()
        counts = pd.Series(calls[mask]).value_counts()
        top = counts.index[0]
        second = counts.index[1] if len(counts) > 1 else ""
        means = z.loc[mask].mean(axis=0) if names else None
        rows.append({"cluster": cl, "cell_type": top, "score": float(means[top]) if (means is not None and top in means) else float("nan"),
                     "n_obs": int(mask.sum()), "best_set": top if top != UNASSIGNED else "",
                     "second_best": second, "second_score": float(counts.iloc[1] / mask.sum()) if second else float("nan"),
                     "fraction": float(counts.iloc[0] / mask.sum())})
    return list(calls), rows
