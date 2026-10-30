"""
Fixtures that follow the workflow contract (spatial/<sample_id>.h5ad as st-qc
writes them, metadata/samples.csv).

  python -m tests.make_fixture data/fixture [--layout both|linked|flat] [--species human|mouse]
      synthetic: two samples, ~2000 cells each, 60 genes, 7 planted cell types with known
      markers laid out in spatial patches (truth in obs["annotation"] and truth.csv)
  python -m tests.make_fixture data/demo --from-prep DIR [--max-obs 20000]
      st-prep output (raw counts in X) normalized the way st-qc does it (a stand-in for st-qc)
  python -m tests.make_fixture data/demo --from-xenium DIR [--max-obs 20000]
      raw Xenium folders listed in DIR/samples.csv (sample_id, path) read, filtered, normalized
      (a stand-in for st-prep + st-qc while those apps are being built)
"""
from __future__ import annotations

import argparse
import csv
import json
import shutil
from pathlib import Path

import numpy as np

# Planted types use genes that belong to exactly one built-in marker set, so the expected call
# is unambiguous (CD3E/CD3D would also score the CD8 / CD4 sets).
PLANTED = {
    "T cell": ["CD3G", "CD2", "TRAC"],
    "B cell": ["MS4A1", "CD79A", "CD19"],
    "Plasma cell": ["MZB1", "JCHAIN", "XBP1"],
    "Macrophage": ["CD68", "CD163", "MRC1"],
    "Mast cell": ["TPSAB1", "CPA3", "KIT"],
    "Fibroblast": ["COL1A1", "DCN", "LUM"],
    "Endothelial": ["PECAM1", "VWF", "CDH5"],
}
N_FILLER = 60 - sum(len(g) for g in PLANTED.values())
SAMPLES = ["s1", "s2"]
TARGET_SUM = 100.0
FIELD = 2000.0   # microns


def mouse_symbol(gene: str) -> str:
    return gene[0] + gene[1:].lower()


def _planted_cells(rng: np.random.Generator, n: int):
    """Types laid out in a 4 x 2 grid of patches (one per type, the last patch a mixture)."""
    types = list(PLANTED)
    x = rng.uniform(0, FIELD, n)
    y = rng.uniform(0, FIELD, n)
    col = np.minimum((x / (FIELD / 4)).astype(int), 3)
    row = np.minimum((y / (FIELD / 2)).astype(int), 1)
    patch = row * 4 + col
    labels = np.empty(n, dtype=object)
    for i in range(n):
        labels[i] = types[patch[i]] if patch[i] < len(types) else types[rng.integers(len(types))]
    return x, y, labels


def _counts(rng: np.random.Generator, labels, genes: list[str], depth: float):
    n = len(labels)
    base = rng.lognormal(mean=-0.5, sigma=0.8, size=len(genes))   # filler genes, same in every type
    mean = np.tile(base, (n, 1))
    for t, markers in PLANTED.items():
        rows = labels == t
        for g in markers:
            j = genes.index(g)
            mean[:, j] = 0.15
            mean[rows, j] = 6.0
    size = rng.lognormal(0, 0.3, n)[:, None] * depth
    return rng.poisson(mean * size).astype(np.int32)


def write_sample(path: Path, sample_id: str, counts, genes: list[str], x, y, labels, seed: int) -> None:
    import anndata as ad
    import pandas as pd
    import scanpy as sc
    import scipy.sparse as sp

    n = counts.shape[0]
    obs_names = [f"{sample_id}:{i + 1}" for i in range(n)]
    obs = pd.DataFrame({
        "sample_id": pd.Categorical([sample_id] * n), "x": x, "y": y, "in_tissue": np.ones(n, dtype=np.int32),
        "total_counts": counts.sum(axis=1).astype(np.float32), "n_genes_by_counts": (counts > 0).sum(axis=1).astype(np.int32),
        "pct_counts_mt": np.zeros(n, dtype=np.float32), "cell_area": np.full(n, 80.0, dtype=np.float32),
        "control_counts": np.zeros(n, dtype=np.float32), "annotation": pd.Categorical(labels), "condition": "demo",
    }, index=obs_names)
    var = pd.DataFrame({"gene_ids": genes, "feature_types": "Gene Expression", "mt": False,
                        "highly_variable": True}, index=genes)
    a = ad.AnnData(X=sp.csr_matrix(counts.astype(np.float32)), obs=obs, var=var)
    a.layers["counts"] = sp.csr_matrix(counts.astype(np.int32))
    sc.pp.normalize_total(a, target_sum=TARGET_SUM)
    sc.pp.log1p(a)
    a.X = sp.csr_matrix(a.X, dtype=np.float32)
    a.obsm["spatial"] = np.column_stack([x, y]).astype(np.float64)
    sc.pp.pca(a, n_comps=20, random_state=seed)   # per sample, as st-qc does: not a joint embedding
    a.uns["pennsieve"] = {"platform": "xenium", "kind": "cell", "coordinate_unit": "micron", "spot_diameter": 10.0,
                          "grid": "none", "source": sample_id, "prep_version": "fixture", "normalized": True,
                          "target_sum": TARGET_SUM}
    path.parent.mkdir(parents=True, exist_ok=True)
    a.write_h5ad(path)


def make_synthetic(root: Path, n_cells: int = 2000, seed: int = 7, species: str = "human") -> Path:
    rng = np.random.default_rng(seed)
    genes = [g for gs in PLANTED.values() for g in gs] + [f"G{i + 1:02d}" for i in range(N_FILLER)]
    rng.shuffle(genes)
    out_genes = [mouse_symbol(g) if species == "mouse" else g for g in genes]
    rows, truth = [], []
    for k, sid in enumerate(SAMPLES):
        x, y, labels = _planted_cells(rng, n_cells)
        counts = _counts(rng, labels, genes, depth=1.0 if k == 0 else 0.6)   # sample 2 is shallower (a batch effect)
        write_sample(root / "spatial" / f"{sid}.h5ad", sid, counts, out_genes, x, y, labels, seed + k)
        rows.append({"sample_id": sid, "path": sid, "platform": "xenium", "kind": "cell", "n_obs": n_cells,
                     "n_genes": len(genes), "n_obs_in_tissue": n_cells, "condition": "demo"})
        truth.extend({"sample_id": sid, "obs_id": f"{sid}:{i + 1}", "cell_type": labels[i]} for i in range(n_cells))
    meta = root / "metadata"
    meta.mkdir(parents=True, exist_ok=True)
    with open(meta / "samples.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    with open(meta / "prep.json", "w") as fh:
        json.dump({"samples": {r["sample_id"]: {"platform": "xenium", "files": [], "warnings": []} for r in rows}}, fh, indent=2)
    with open(root / "truth.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["sample_id", "obs_id", "cell_type"])
        w.writeheader()
        w.writerows(truth)
    return root


def flatten(linked: Path, flat: Path) -> Path:
    """The same files side by side at the root (no links used)."""
    flat.mkdir(parents=True, exist_ok=True)
    for p in sorted(linked.rglob("*")):
        if p.is_file():
            shutil.copy2(p, flat / p.name)
    return flat


def make_fixture(root: Path, layout: str = "both", species: str = "human", n_cells: int = 2000) -> dict[str, Path]:
    root = Path(root)
    out = {}
    linked = make_synthetic(root / "linked", n_cells=n_cells, species=species)
    if layout in ("linked", "both"):
        out["linked"] = linked
    if layout in ("flat", "both"):
        out["flat"] = flatten(linked, root / "flat")
        if layout == "flat":
            shutil.rmtree(linked)
    return out


def write_markers_csv(path: Path, shape: str, sets: dict[str, list[str]] | None = None) -> Path:
    """A markers.csv in the long (cell_type, gene) or wide (cell_type, genes) shape."""
    sets = sets or PLANTED
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        if shape == "long":
            w.writerow(["cell_type", "gene"])
            for ct, genes in sets.items():
                for g in genes:
                    w.writerow([ct, g])
        else:
            w.writerow(["cell_type", "genes"])
            for ct, genes in sets.items():
                w.writerow([ct, ";".join(genes)])
    return path


# --- stand-ins for the upstream apps on real data -----------------------------------------

def _qc_like(a, sample_id: str, max_obs: int, seed: int, min_counts: int = 10):
    """Filter, subsample, normalize and PCA one sample with raw counts in X, as st-qc would."""
    import pandas as pd
    import scanpy as sc
    import scipy.sparse as sp

    a.var_names_make_unique()
    a.obs["sample_id"] = pd.Categorical([sample_id] * a.n_obs)
    keep = np.asarray(a.X.sum(axis=1)).ravel() >= min_counts
    a = a[keep].copy()
    if max_obs and a.n_obs > max_obs:
        idx = np.sort(np.random.default_rng(seed).choice(a.n_obs, size=max_obs, replace=False))
        a = a[idx].copy()
    a.layers["counts"] = sp.csr_matrix(a.X, dtype=np.int32)
    a.X = sp.csr_matrix(a.X, dtype=np.float32)
    sc.pp.normalize_total(a, target_sum=TARGET_SUM)
    sc.pp.log1p(a)
    a.X = sp.csr_matrix(a.X, dtype=np.float32)
    if a.n_vars > 1000:
        sc.pp.highly_variable_genes(a, n_top_genes=2000, flavor="seurat")
    else:
        a.var["highly_variable"] = True
    sc.pp.pca(a, n_comps=min(30, a.n_vars - 1), random_state=seed)
    pen = dict(a.uns.get("pennsieve", {}))
    pen.update({"normalized": True, "target_sum": TARGET_SUM})
    a.uns["pennsieve"] = pen
    return a


def read_xenium(folder: Path, sample_id: str):
    """cell_feature_matrix.h5 + cells.parquet (or cells.csv.gz) -> AnnData with raw counts in X."""
    import anndata as ad
    import h5py
    import pandas as pd
    import scipy.sparse as sp

    with h5py.File(folder / "cell_feature_matrix.h5", "r") as f:
        m = f["matrix"]
        shape = tuple(int(v) for v in m["shape"][:])     # genes x cells
        X = sp.csc_matrix((m["data"][:], m["indices"][:], m["indptr"][:]), shape=shape).T.tocsr()
        barcodes = [b.decode() for b in m["barcodes"][:]]
        names = [n.decode() for n in m["features"]["name"][:]]
        ids = [n.decode() for n in m["features"]["id"][:]]
        types = [t.decode() for t in m["features"]["feature_type"][:]]
    is_gene = np.array([t == "Gene Expression" for t in types])
    control = np.asarray(X[:, ~is_gene].sum(axis=1)).ravel()
    X = X[:, is_gene]
    cells = pd.read_parquet(folder / "cells.parquet") if (folder / "cells.parquet").exists() \
        else pd.read_csv(folder / "cells.csv.gz")
    cells["cell_id"] = cells["cell_id"].astype(str)
    cells = cells.set_index("cell_id").reindex(barcodes)
    obs = pd.DataFrame({
        "x": cells["x_centroid"].to_numpy(float), "y": cells["y_centroid"].to_numpy(float), "in_tissue": 1,
        "total_counts": np.asarray(X.sum(axis=1)).ravel().astype(np.float32),
        "n_genes_by_counts": np.asarray((X > 0).sum(axis=1)).ravel().astype(np.int32), "pct_counts_mt": np.float32(0.0),
        "cell_area": cells["cell_area"].to_numpy(float), "nucleus_area": cells["nucleus_area"].to_numpy(float),
        "control_counts": control.astype(np.float32),
    }, index=[f"{sample_id}:{b}" for b in barcodes])
    var = pd.DataFrame({"gene_ids": [i for i, g in zip(ids, is_gene) if g], "feature_types": "Gene Expression"},
                       index=[n for n, g in zip(names, is_gene) if g])
    var["mt"] = var.index.str.upper().str.startswith("MT-")
    a = ad.AnnData(X=sp.csr_matrix(X, dtype=np.float32), obs=obs, var=var)
    a.obsm["spatial"] = obs[["x", "y"]].to_numpy(np.float64)
    a.uns["pennsieve"] = {"platform": "xenium", "kind": "cell", "coordinate_unit": "micron", "spot_diameter": 10.0,
                          "grid": "none", "source": str(folder), "prep_version": "fixture-harness"}
    return a


def from_real(root: Path, source: Path, mode: str, max_obs: int, seed: int = 42) -> Path:
    import anndata as ad

    root = Path(root)
    (root / "spatial").mkdir(parents=True, exist_ok=True)
    (root / "metadata").mkdir(parents=True, exist_ok=True)
    rows = []
    if mode == "prep":
        with open(source / "metadata" / "samples.csv", newline="") as fh:
            sheet = list(csv.DictReader(fh))
        for r in sheet:
            a = ad.read_h5ad(next(source.rglob(f"{r['sample_id']}.h5ad")))
            a = _qc_like(a, r["sample_id"], max_obs, seed)
            a.write_h5ad(root / "spatial" / f"{r['sample_id']}.h5ad")
            rows.append({**r, "n_obs": a.n_obs, "n_genes": a.n_vars, "n_obs_in_tissue": a.n_obs})
            print(f"{r['sample_id']}: {a.n_obs} x {a.n_vars}")
    else:
        with open(source / "samples.csv", newline="") as fh:
            sheet = list(csv.DictReader(fh))
        for r in sheet:
            a = read_xenium(source / r["path"], r["sample_id"])
            n_raw = a.n_obs
            a = _qc_like(a, r["sample_id"], max_obs, seed)
            a.write_h5ad(root / "spatial" / f"{r['sample_id']}.h5ad")
            rows.append({"sample_id": r["sample_id"], "path": r["path"], "platform": "xenium", "kind": "cell",
                         "n_obs": a.n_obs, "n_genes": a.n_vars, "n_obs_in_tissue": a.n_obs,
                         **{k: v for k, v in r.items() if k not in ("sample_id", "path")}})
            print(f"{r['sample_id']}: {n_raw} cells read, {a.n_obs} kept x {a.n_vars} genes")
    with open(root / "metadata" / "samples.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    return root


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out")
    ap.add_argument("--layout", choices=["both", "linked", "flat"], default="both")
    ap.add_argument("--species", choices=["human", "mouse"], default="human")
    ap.add_argument("--cells", type=int, default=2000, help="cells per synthetic sample")
    ap.add_argument("--from-prep", metavar="DIR", help="st-prep output folder (spatial/, metadata/samples.csv)")
    ap.add_argument("--from-xenium", metavar="DIR", help="folder with samples.csv (sample_id, path) and raw Xenium folders")
    ap.add_argument("--max-obs", type=int, default=0, help="seeded subsample per sample (0 = all)")
    args = ap.parse_args()
    out = Path(args.out)
    if out.exists():
        shutil.rmtree(out)
    if args.from_prep:
        from_real(out, Path(args.from_prep), "prep", args.max_obs)
    elif args.from_xenium:
        from_real(out, Path(args.from_xenium), "xenium", args.max_obs)
    else:
        for layout, path in make_fixture(out, args.layout, args.species, args.cells).items():
            print(f"{layout}: {path}")


if __name__ == "__main__":
    main()
