"""
st-celltype: clusters and cell types for imaging-based spatial
transcriptomics (Xenium, MERSCOPE, CosMx) and Visium HD bins.

Reads the qc'd <sample_id>.h5ad files (port `spatial`), samples.csv (port
`metadata`) and an optional markers.csv (port `markers`); clusters all samples
jointly (PCA, Harmony by sample, neighbours, Leiden), scores marker sets per
observation and calls a cell type per cluster (or per cell), and writes

  labels/   labels.csv (sample_id, obs_id, cluster, cell_type), cluster_markers.csv,
            cell_type_summary.csv, cluster_summary.csv, cluster_annotation.csv, scores.csv,
            agreement.csv (when obs["annotation"] exists), summary.json, provenance.json
  figures/  umap_cluster.png, umap_cell_type.png, <sample_id>_cell_type.png,
            dotplot_marker_sets.png, cell_type_by_sample.png, summary.json, provenance.json

See processor/cluster.py and processor/markers.py for the method.
"""
from __future__ import annotations

import csv
import json
import logging
import math
import os
import platform
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

from processor import threads as T

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
log = logging.getLogger("st-celltype")
log.setLevel(logging.INFO)
for noisy in ("matplotlib", "PIL", "fontTools", "numba", "harmonypy"):
    logging.getLogger(noisy).setLevel(logging.WARNING)

APP_NAME = "st-celltype"
REPO_ROOT = Path(__file__).resolve().parent.parent

DEFAULTS = {
    "RESOLUTION": "0.5",
    "N_NEIGHBORS": "15",
    "N_PCS": "30",
    "INTEGRATION": "harmony",
    "MARKER_SET": "builtin",
    "ANNOTATE": "cluster",
    "MIN_SCORE": "0.1",
    "MIN_GENES_PRESENT": "2",
    "MAX_OBS_FOR_UMAP": "50000",
    "SEED": "42",
    "THREADS": "0",
}
PARAMETERS = tuple(DEFAULTS)
CHOICES = {"INTEGRATION": ("harmony", "none"), "MARKER_SET": ("builtin", "none"), "ANNOTATE": ("cluster", "cell")}


class ParameterError(ValueError):
    pass


def read_parameters(env=os.environ) -> dict:
    raw = {k: (env.get(k) if env.get(k) not in (None, "") else v) for k, v in DEFAULTS.items()}
    p: dict = {}
    for name, choices in CHOICES.items():
        value = raw[name].strip().lower()
        if value not in choices:
            raise ParameterError(f"{name} must be one of {', '.join(choices)}, got {raw[name]!r}")
        p[name] = value
    try:
        p["RESOLUTION"] = float(raw["RESOLUTION"])
        p["MIN_SCORE"] = float(raw["MIN_SCORE"])
        for name in ("N_NEIGHBORS", "N_PCS", "MIN_GENES_PRESENT", "MAX_OBS_FOR_UMAP", "SEED", "THREADS"):
            p[name] = int(raw[name])
    except ValueError as exc:
        raise ParameterError(f"bad numeric parameter: {exc}") from exc
    if p["RESOLUTION"] <= 0 or p["N_NEIGHBORS"] < 2 or p["N_PCS"] < 2 or p["MIN_GENES_PRESENT"] < 1 \
            or p["MAX_OBS_FOR_UMAP"] < 10 or p["THREADS"] < 0:
        raise ParameterError("RESOLUTION > 0, N_NEIGHBORS >= 2, N_PCS >= 2, MIN_GENES_PRESENT >= 1, "
                             "MAX_OBS_FOR_UMAP >= 10, THREADS >= 0")
    return p


def _rel(path: Path | None, root: Path) -> str | None:
    if path is None:
        return None
    try:
        return str(Path(path).relative_to(root))
    except ValueError:
        return str(path)


def _fmt(v):
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, float):
        return "" if math.isnan(v) else f"{v:.6g}"
    return v


def write_csv(path: Path, rows, columns: list[str]) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(columns)
        for r in rows:
            w.writerow([_fmt(r.get(c, "")) for c in columns])


def _clean_json(obj):
    if isinstance(obj, float):
        return None if math.isnan(obj) or math.isinf(obj) else obj
    if isinstance(obj, dict):
        return {str(k): _clean_json(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_clean_json(v) for v in obj]
    if hasattr(obj, "item"):   # numpy scalars
        return _clean_json(obj.item())
    return obj


def write_json(path: Path, data) -> None:
    with open(path, "w") as fh:
        json.dump(_clean_json(data), fh, indent=2)
        fh.write("\n")


def tool_versions() -> dict:
    from importlib.metadata import version

    out = {"python": platform.python_version()}
    for pkg in ("scanpy", "anndata", "harmonypy", "igraph", "umap-learn", "scikit-learn", "numpy", "scipy",
                "pandas", "matplotlib", "h5py"):
        try:
            out[pkg] = version(pkg)
        except Exception:
            out[pkg] = None
    return out


def peak_rss_mb() -> float | None:
    try:
        import resource

        peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        return round(peak / (1024 * 1024 if sys.platform == "darwin" else 1024), 1)
    except (ImportError, OSError):
        return None


def cgroup_peak_mb() -> float | None:
    try:
        return round(int(Path("/sys/fs/cgroup/memory.peak").read_text()) / (1024 * 1024), 1)
    except (OSError, ValueError):
        return None


def app_version() -> str:
    try:
        return (REPO_ROOT / "VERSION").read_text().strip()
    except OSError:
        return "unknown"


def main() -> int:
    started = datetime.now(timezone.utc)
    t0 = time.monotonic()
    input_dir, output_dir = os.environ.get("INPUT_DIR"), os.environ.get("OUTPUT_DIR")
    if not input_dir or not output_dir:
        log.error("INPUT_DIR and OUTPUT_DIR must be set")
        return 2
    input_dir, output_dir = Path(input_dir), Path(output_dir)
    if not input_dir.is_dir():
        log.error("INPUT_DIR %s is not a folder", input_dir)
        return 2
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
    os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba")
    try:
        params = read_parameters()
    except ParameterError as exc:
        log.error("%s", exc)
        return 2
    n_threads, thread_source = T.configure(params["THREADS"])
    log.info("threads: %d (%s)", n_threads, thread_source)

    # numpy and friends only after the thread variables are set
    import numpy as np
    import pandas as pd

    from processor import cluster as C
    from processor import figures as F
    from processor import markers as M
    from processor.inputs import discover

    try:
        found = discover(input_dir)
    except (OSError, ValueError) as exc:
        log.error("cannot read the inputs: %s", exc)
        return 1
    for name, reason in found.skipped.items():
        log.warning("%s: skipped, %s", name, reason)
    if not found.samples:
        log.error("no .h5ad samples found under %s", input_dir)
        return 1
    try:
        sets, marker_source = M.load_sets(params["MARKER_SET"], found.markers_csv)
    except ValueError as exc:
        log.error("marker sets: %s", exc)
        return 2
    log.info("%d samples, marker sets: %s (%d sets)", len(found.samples), marker_source, len(sets))

    # --- read, one sample at a time ---------------------------------------------------------
    t_read = time.monotonic()
    adatas, statuses = [], []
    for s in found.samples:
        t = time.monotonic()
        status = {"sample_id": s.sample_id, "file": _rel(s.path, input_dir)}
        try:
            a = C.read_sample(s.path, s.sample_id)
            adatas.append(a)
            status.update(status="ok", n_obs=int(a.n_obs), n_genes=int(a.n_vars),
                          has_pca="X_pca" in a.obsm, seconds=round(time.monotonic() - t, 3))
            log.info("%s: %d observations x %d genes%s, %.2fs", s.sample_id, a.n_obs, a.n_vars,
                     "" if "X_pca" in a.obsm else " (no X_pca)", time.monotonic() - t)
        except Exception as exc:
            status.update(status="failed", error=f"{type(exc).__name__}: {exc}", seconds=round(time.monotonic() - t, 3))
            log.error("%s: failed: %s", s.sample_id, exc)
        statuses.append(status)
    if not adatas:
        log.error("every sample failed")
        return 1
    adata, joint_pca = C.concatenate(adatas)
    del adatas
    t_read = time.monotonic() - t_read
    rss_read = peak_rss_mb()
    log.info("joint matrix: %d observations x %d genes, %.1fs, peak RSS %s MB", adata.n_obs, adata.n_vars, t_read, rss_read)
    if adata.n_obs < 10:
        log.error("only %d observations in total; nothing to cluster", adata.n_obs)
        return 1

    # --- embedding and clustering -----------------------------------------------------------
    t_cluster = time.monotonic()
    rep, embed_info = C.embed(adata, params["N_PCS"], params["SEED"], params["INTEGRATION"], n_threads)
    if joint_pca:
        embed_info["pca"] = "reused from st-qc (same loadings in every sample)"
    log.info("PCA %s (%d PCs, %.1fs); integration: %s", embed_info["pca"], embed_info["n_pcs"],
             embed_info["pca_seconds"], embed_info["integration"])
    labels, cluster_info = C.cluster(adata, rep, params["N_NEIGHBORS"], params["RESOLUTION"], params["SEED"])
    adata.obs["cluster"] = pd.Categorical(labels)
    t_cluster = time.monotonic() - t_cluster
    log.info("%d clusters (neighbours %.1fs, leiden %.1fs)", cluster_info["n_clusters"],
             cluster_info["neighbors_seconds"], cluster_info["leiden_seconds"])

    # --- marker scores and cell types --------------------------------------------------------
    t_score = time.monotonic()
    used, skipped_sets = M.resolve(sets, adata.var_names, params["MIN_GENES_PRESENT"])
    for s in skipped_sets:
        log.info("marker set %s skipped: %d of %d genes in the panel (%s)", s.cell_type, len(s.present),
                 len(s.genes), ", ".join(s.present) or "none")
    for s in used:
        log.info("marker set %s: %d of %d genes present%s", s.cell_type, len(s.present), len(s.genes),
                 f", missing {', '.join(s.missing)}" if s.missing else "")
    if sets and not used:
        log.warning("no marker set has %d genes in the panel: every observation is Unassigned", params["MIN_GENES_PRESENT"])
    scores = M.score_sets(adata, used, params["SEED"]) if used else pd.DataFrame(index=adata.obs_names)
    z = M.zscore(scores)
    if params["ANNOTATE"] == "cluster":
        cell_types, annotation_rows = M.annotate_by_cluster(z, labels, params["MIN_SCORE"])
    else:
        cell_types, annotation_rows = M.annotate_by_cell(z, labels, params["MIN_SCORE"])
    adata.obs["cell_type"] = pd.Categorical(cell_types)
    t_score = time.monotonic() - t_score
    for r in annotation_rows:
        log.info("%s: %d obs -> %s (%.2f%s)", r["cluster"], r["n_obs"], r["cell_type"],
                 r["score"] if not math.isnan(r["score"]) else float("nan"),
                 f", then {r['second_best']} {r['second_score']:.2f}" if r["second_best"] else "")

    # --- cluster markers --------------------------------------------------------------------
    t_markers = time.monotonic()
    marker_rows = []
    try:
        marker_rows = C.cluster_markers(adata, "cluster", params["SEED"])
    except Exception as exc:
        log.exception("cluster markers failed: %s", exc)
    t_markers = time.monotonic() - t_markers
    log.info("cluster markers: %d rows, %.1fs", len(marker_rows), t_markers)

    # --- tables -----------------------------------------------------------------------------
    labels_dir = output_dir / "labels"
    figures_dir = output_dir / "figures"
    labels_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)
    sample_ids = adata.obs["sample_id"].astype(str).to_numpy()
    obs_ids = adata.obs_names.to_numpy()
    samples_in_order = [s.sample_id for s in found.samples if s.sample_id in set(sample_ids)]
    cluster_order = sorted(set(labels))
    type_order = [s.cell_type for s in used if s.cell_type in set(cell_types)]
    if M.UNASSIGNED in set(cell_types):
        type_order.append(M.UNASSIGNED)

    frame = pd.DataFrame({"sample_id": sample_ids, "obs_id": obs_ids, "cluster": labels, "cell_type": cell_types})
    frame.to_csv(labels_dir / "labels.csv", index=False)
    write_csv(labels_dir / "cluster_markers.csv", marker_rows,
              ["cluster", "gene", "logfoldchange", "pval_adj", "pct_in", "pct_out", "rank", "score"])

    def summarise(column: str, order: list[str]) -> list[dict]:
        rows = []
        for sid in samples_in_order:
            sub = frame[frame["sample_id"] == sid]
            counts = sub[column].value_counts()
            for lab in order:
                n = int(counts.get(lab, 0))
                rows.append({"sample_id": sid, column: lab, "n_obs": n, "fraction": n / len(sub) if len(sub) else 0.0})
        return rows

    type_rows = summarise("cell_type", type_order)
    cluster_rows = summarise("cluster", cluster_order)
    write_csv(labels_dir / "cell_type_summary.csv", type_rows, ["sample_id", "cell_type", "n_obs", "fraction"])
    write_csv(labels_dir / "cluster_summary.csv", cluster_rows, ["sample_id", "cluster", "n_obs", "fraction"])
    write_csv(labels_dir / "cluster_annotation.csv", annotation_rows,
              ["cluster", "cell_type", "score", "n_obs", "second_best", "second_score"])
    scores_out = z.copy()
    scores_out.insert(0, "obs_id", obs_ids)
    scores_out.insert(0, "sample_id", sample_ids)
    scores_out.to_csv(labels_dir / "scores.csv", index=False, float_format="%.5g")

    agreement = None
    if "annotation" in adata.obs:
        agreement = compute_agreement(frame, adata.obs["annotation"], samples_in_order)
        write_csv(labels_dir / "agreement.csv", agreement["rows"], ["sample_id", "label", "ari", "nmi", "n_obs"])
        write_csv(labels_dir / "confusion.csv", agreement["confusion"], ["label", "value", "annotation", "n_obs"])
        for r in agreement["rows"]:
            if r["sample_id"] == "all":
                log.info("agreement with obs['annotation'], %s: ARI %.3f, NMI %.3f", r["label"], r["ari"], r["nmi"])

    # --- figures ----------------------------------------------------------------------------
    t_fig = time.monotonic()
    figure_paths = []
    try:
        figure_paths = make_figures(adata, rep, frame, z, used, cluster_order, type_order, samples_in_order,
                                    params, figures_dir, C, F)
    except Exception as exc:
        log.exception("figures failed: %s", exc)
    t_fig = time.monotonic() - t_fig

    # --- summary and provenance -------------------------------------------------------------
    overall = frame["cell_type"].value_counts()
    summary = {
        "n_samples": len(samples_in_order), "n_obs": int(adata.n_obs), "n_genes": int(adata.n_vars),
        "n_clusters": cluster_info["n_clusters"],
        "n_cell_types": int(sum(1 for t in overall.index if t != M.UNASSIGNED)),
        "fraction_unassigned": float(overall.get(M.UNASSIGNED, 0) / adata.n_obs),
        "cell_types": [{"cell_type": t, "n_obs": int(overall[t]), "fraction": float(overall[t] / adata.n_obs)}
                       for t in type_order],
        "cell_types_by_sample": {sid: {r["cell_type"]: r["fraction"] for r in type_rows if r["sample_id"] == sid}
                                 for sid in samples_in_order},
        "clusters": [{"cluster": r["cluster"], "cell_type": r["cell_type"], "n_obs": r["n_obs"], "score": r["score"],
                      "second_best": r["second_best"], "second_score": r["second_score"]} for r in annotation_rows],
        "marker_sets": {"source": marker_source, "n_given": len(sets), "n_used": len(used),
                        "used": {s.cell_type: {"genes_used": s.present, "genes_missing": s.missing} for s in used},
                        "skipped": {s.cell_type: {"genes_present": s.present, "genes_missing": s.missing} for s in skipped_sets}},
        "embedding": embed_info, "clustering": cluster_info,
        "agreement": agreement["rows"] if agreement else None,
        "samples": {s["sample_id"]: {k: s[k] for k in ("status", "n_obs", "n_genes") if k in s} for s in statuses},
        "parameters": params, "figures": [p.name for p in figure_paths],
    }
    write_json(labels_dir / "summary.json", summary)
    write_json(figures_dir / "summary.json", {"figures": [p.name for p in figure_paths], "n_samples": len(samples_in_order),
                                              "n_clusters": cluster_info["n_clusters"], "n_cell_types": summary["n_cell_types"]})
    n_failed = sum(1 for s in statuses if s["status"] != "ok")
    provenance = {
        "app": APP_NAME, "version": app_version(), "tools": tool_versions(), "parameters": params,
        "threads": {"n": n_threads, "source": thread_source},
        "method": {
            "pca": embed_info["pca"], "n_pcs": embed_info["n_pcs"], "integration": embed_info["integration"],
            "neighbors": f"scanpy.pp.neighbors(n_neighbors={params['N_NEIGHBORS']}, use_rep={rep})",
            "clustering": f"scanpy.tl.leiden(resolution={params['RESOLUTION']}, flavor=igraph, n_iterations=2)",
            "scores": "scanpy.tl.score_genes per marker set on the log-normalized X, z-scored over all observations",
            "annotation": f"{params['ANNOTATE']}: highest mean z-score, Unassigned below MIN_SCORE={params['MIN_SCORE']}",
            "cluster_markers": f"scanpy.tl.rank_genes_groups(wilcoxon, vs rest, BH), top {C.N_TOP_MARKERS}, "
                               f"on at most {C.MAX_OBS_FOR_MARKERS} seeded observations",
            "umap": f"scanpy.tl.umap on at most {params['MAX_OBS_FOR_UMAP']} seeded observations",
            "marker_sets": marker_source,
        },
        "inputs": {"spatial": [_rel(p, input_dir) for p in found.h5ad_seen], "samples_csv": _rel(found.samples_csv, input_dir),
                   "markers_csv": _rel(found.markers_csv, input_dir), "skipped": found.skipped},
        "started": started.isoformat(), "finished": datetime.now(timezone.utc).isoformat(),
        "timings_seconds": {"read": round(t_read, 2), "cluster": round(t_cluster, 2), "score": round(t_score, 2),
                            "markers": round(t_markers, 2), "figures": round(t_fig, 2), "total": round(time.monotonic() - t0, 2)},
        "peak_rss_mb": {"after_reading": rss_read, "total": peak_rss_mb(), "cgroup_peak": cgroup_peak_mb()},
        "samples": statuses, "n_ok": len(statuses) - n_failed, "n_failed": n_failed,
    }
    write_json(labels_dir / "provenance.json", {**provenance, "outputs": sorted(p.name for p in labels_dir.iterdir())})
    write_json(figures_dir / "provenance.json", {**provenance, "outputs": sorted(p.name for p in figures_dir.iterdir())})
    log.info("done: %d samples (%d failed), %d observations, %d clusters, %d cell types (%.1f%% unassigned), "
             "%.1fs, peak RSS %s MB -> %s", len(statuses) - n_failed, n_failed, adata.n_obs, cluster_info["n_clusters"],
             summary["n_cell_types"], 100 * summary["fraction_unassigned"], time.monotonic() - t0, peak_rss_mb(), output_dir)
    return 0


def compute_agreement(frame, annotation, samples_in_order: list[str]) -> dict:
    """ARI / NMI of cluster and cell_type against obs['annotation'], per sample and overall."""
    import pandas as pd
    from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

    ann = pd.Series(annotation.astype(str).to_numpy(), index=frame.index)
    valid = ~ann.isin(["", "nan", "NaN", "None"])
    rows, confusion = [], []
    for label in ("cluster", "cell_type"):
        for sid in samples_in_order + ["all"]:
            m = valid & ((frame["sample_id"] == sid) if sid != "all" else True)
            if int(m.sum()) < 2:
                continue
            rows.append({"sample_id": sid, "label": label, "ari": float(adjusted_rand_score(ann[m], frame.loc[m, label])),
                         "nmi": float(normalized_mutual_info_score(ann[m], frame.loc[m, label])), "n_obs": int(m.sum())})
        ct = pd.crosstab(frame.loc[valid, label], ann[valid])
        for value in ct.index:
            for a in ct.columns:
                n = int(ct.loc[value, a])
                if n:
                    confusion.append({"label": label, "value": value, "annotation": a, "n_obs": n})
    return {"rows": rows, "confusion": confusion}


def make_figures(adata, rep, frame, z, used, cluster_order, type_order, samples_in_order, params, out, C, F):
    import numpy as np

    paths = []
    t = time.monotonic()
    idx, coords = C.umap_subsample(adata, rep, params["MAX_OBS_FOR_UMAP"], params["N_NEIGHBORS"], params["SEED"])
    log.info("UMAP of %d observations, %.1fs", len(idx), time.monotonic() - t)
    sub = frame.iloc[idx]
    n_note = f" ({len(idx):,} of {adata.n_obs:,} observations)" if len(idx) < adata.n_obs else f" ({adata.n_obs:,} observations)"
    paths.append(F.scatter_categorical(coords, sub["cluster"].to_numpy(), cluster_order, "UMAP by cluster" + n_note,
                                       out / "umap_cluster.png"))
    paths.append(F.scatter_categorical(coords, sub["cell_type"].to_numpy(), type_order, "UMAP by cell type" + n_note,
                                       out / "umap_cell_type.png"))
    spatial = np.asarray(adata.obsm["spatial"], dtype=np.float64)
    uns_spatial = adata.uns.get("spatial", {}) if isinstance(adata.uns.get("spatial", {}), dict) else {}
    unit = ""
    pen = adata.uns.get("pennsieve") if isinstance(adata.uns.get("pennsieve"), dict) else None
    for sid in samples_in_order:
        m = (frame["sample_id"] == sid).to_numpy()
        xy = spatial[m]
        image, scale = F.lowres_image(uns_spatial, sid)
        if image is not None:
            xy = xy * scale
        paths.append(F.scatter_categorical(xy, frame.loc[m, "cell_type"].to_numpy(), type_order,
                                           f"{sid}: cell types in space ({int(m.sum()):,} observations)",
                                           out / f"{F.safe_name(sid)}_cell_type.png", xlabel="x" + unit, ylabel="y" + unit,
                                           invert_y=True, image=image))
    if used:
        sets = [s.cell_type for s in used]
        clusters = frame["cluster"].to_numpy()
        means = np.zeros((len(cluster_order), len(sets)))
        frac = np.zeros_like(means)
        zv = z[sets].to_numpy(dtype=float)
        for i, cl in enumerate(cluster_order):
            m = clusters == cl
            means[i] = zv[m].mean(axis=0)
            frac[i] = (zv[m] > 0).mean(axis=0)
        paths.append(F.dotplot(means, frac, cluster_order, sets, "Marker-set scores by cluster", out / "dotplot_marker_sets.png"))
    table = {}
    for sid in samples_in_order:
        sub = frame[frame["sample_id"] == sid]
        counts = sub["cell_type"].value_counts()
        table[sid] = {lab: float(counts.get(lab, 0) / len(sub)) for lab in type_order} if len(sub) else {}
    paths.append(F.stacked_bar(table, samples_in_order, type_order, "Cell types by sample", out / "cell_type_by_sample.png"))
    return paths


if __name__ == "__main__":
    sys.exit(main())
