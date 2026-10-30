"""
Joint clustering of all samples: read each qc'd .h5ad (X, obs, var, obsm only;
layers["counts"] stays on disk), concatenate on the shared genes, PCA (reused
from st-qc when it is one joint embedding, else computed here), Harmony by
sample, scanpy neighbours, Leiden (igraph flavour), a UMAP of a subsample for
the figures, and Wilcoxon cluster markers.
"""
from __future__ import annotations

import logging
import time
from pathlib import Path

log = logging.getLogger("st-celltype")

MAX_OBS_FOR_MARKERS = 100_000   # rank_genes_groups runs on at most this many cells (seeded)
N_TOP_MARKERS = 25
SMALL_PANEL = 1000              # at most this many genes: PCA on all genes, else on the HVGs


def read_sample(path: Path, sample_id: str):
    """X (log-normalized), obs, var, obsm["spatial"/"X_pca"], varm["PCs"], uns["spatial"] and
    uns["pennsieve"] of one qc'd file, read element by element so the counts layer is never loaded."""
    import anndata as ad
    import h5py
    import numpy as np
    import scipy.sparse as sp
    from anndata.io import read_elem

    with h5py.File(path, "r") as f:
        X = read_elem(f["X"])
        obs = read_elem(f["obs"])
        var = read_elem(f["var"])
        obsm = {}
        if "obsm" in f:
            for key in ("spatial", "X_pca"):
                if key in f["obsm"]:
                    obsm[key] = np.asarray(read_elem(f["obsm"][key]))
        varm = {}
        if "varm" in f and "PCs" in f["varm"]:
            varm["PCs"] = np.asarray(read_elem(f["varm"]["PCs"]))
        uns = {}
        if "uns" in f:
            for key in ("pennsieve", "spatial"):
                if key in f["uns"]:
                    uns[key] = read_elem(f["uns"][key])
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    X = sp.csr_matrix(X, dtype=np.float32)
    a = ad.AnnData(X=X, obs=obs, var=var, obsm=obsm, varm=varm, uns=uns)
    if "sample_id" not in a.obs:
        a.obs["sample_id"] = sample_id
    a.obs["sample_id"] = a.obs["sample_id"].astype(str)
    if "spatial" not in a.obsm:
        if {"x", "y"} <= set(a.obs.columns):
            a.obsm["spatial"] = a.obs[["x", "y"]].to_numpy(dtype=np.float64)
        else:
            raise ValueError("no obsm['spatial'] and no obs x/y columns")
    pen = a.uns.get("pennsieve", {})
    if not (isinstance(pen, dict) and pen.get("normalized")):
        log.warning("%s: uns['pennsieve']['normalized'] is not set; X is taken as log-normalized anyway", sample_id)
    a.var_names_make_unique()
    a.obs_names_make_unique()
    return a


def concatenate(adatas: list):
    """Inner join on genes; per-sample HVG flags are OR-ed; per-sample PCA is reused only when
    every sample carries the same loadings (one joint embedding), else dropped."""
    import anndata as ad
    import numpy as np

    genes = set(adatas[0].var_names)
    for a in adatas[1:]:
        genes &= set(a.var_names)
    n_lost = max(a.n_vars for a in adatas) - len(genes)
    if n_lost:
        log.warning("%d genes are not in every sample and are dropped from the joint analysis", n_lost)
    hv = None
    for a in adatas:
        if "highly_variable" in a.var:
            flags = a.var["highly_variable"].reindex(sorted(genes)).fillna(False).astype(bool)
            hv = flags if hv is None else (hv | flags)
    joint_pca = all("PCs" in a.varm and "X_pca" in a.obsm for a in adatas)
    if joint_pca and len(adatas) > 1:
        ref = adatas[0]
        order = sorted(genes)
        p0 = ref[:, order].varm["PCs"]
        for a in adatas[1:]:
            p = a[:, order].varm["PCs"]
            if p.shape != p0.shape or not np.allclose(p, p0, atol=1e-5):
                joint_pca = False
                break
    for a in adatas:
        if "PCs" in a.varm:
            del a.varm["PCs"]
        if not joint_pca and "X_pca" in a.obsm:
            del a.obsm["X_pca"]
        for key in list(a.uns):
            if key != "spatial":
                del a.uns[key]
    merged = ad.concat(adatas, join="inner", merge="first", uns_merge="first", index_unique=None)
    if hv is not None:
        merged.var["highly_variable"] = hv.reindex(merged.var_names).fillna(False).astype(bool).to_numpy()
    merged.obs["sample_id"] = merged.obs["sample_id"].astype("category")
    return merged, joint_pca


def embed(adata, n_pcs: int, seed: int, integration: str, threads: int) -> tuple[str, dict]:
    """obsm["X_pca"] (reused or computed) and, with Harmony, obsm["X_harmony"]; returns the
    representation to build neighbours on and what was done."""
    import numpy as np
    import scanpy as sc

    info: dict = {}
    t = time.monotonic()
    if "X_pca" in adata.obsm and adata.obsm["X_pca"].shape[1] >= min(n_pcs, 2):
        info["pca"] = "reused from st-qc"
    else:
        if "X_pca" in adata.obsm:
            del adata.obsm["X_pca"]
        n_comps = int(min(n_pcs, adata.n_vars - 1, adata.n_obs - 1))
        use_hvg = adata.n_vars > SMALL_PANEL and "highly_variable" in adata.var and int(adata.var["highly_variable"].sum()) >= n_comps + 1
        sc.pp.pca(adata, n_comps=n_comps, mask_var="highly_variable" if use_hvg else None, random_state=seed)
        info["pca"] = f"computed on {'the %d highly variable genes' % int(adata.var['highly_variable'].sum()) if use_hvg else 'all %d genes' % adata.n_vars}"
    X = np.asarray(adata.obsm["X_pca"], dtype=np.float32)[:, :n_pcs]
    adata.obsm["X_pca"] = X
    info["n_pcs"] = int(X.shape[1])
    info["pca_seconds"] = round(time.monotonic() - t, 2)
    n_samples = adata.obs["sample_id"].nunique()
    rep = "X_pca"
    if integration == "harmony" and n_samples > 1:
        import harmonypy

        t = time.monotonic()
        ho = harmonypy.run_harmony(X, adata.obs[["sample_id"]], "sample_id", random_state=seed, verbose=False,
                                   ncores=max(1, threads))
        Z = np.asarray(ho.Z_corr)
        if Z.shape[0] != adata.n_obs:
            Z = Z.T
        adata.obsm["X_harmony"] = np.ascontiguousarray(Z, dtype=np.float32)
        rep = "X_harmony"
        info["integration"] = "harmony by sample_id"
        info["harmony_seconds"] = round(time.monotonic() - t, 2)
    else:
        info["integration"] = "none" if integration == "none" else "none (one sample)"
    return rep, info


def cluster(adata, rep: str, n_neighbors: int, resolution: float, seed: int) -> tuple[list[str], dict]:
    """Neighbours on `rep`, Leiden (igraph, 2 iterations); labels C01.. by decreasing size."""
    import numpy as np
    import scanpy as sc

    info: dict = {}
    t = time.monotonic()
    n_neighbors = int(min(n_neighbors, max(2, adata.n_obs - 1)))
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep=rep, random_state=seed)
    info["neighbors_seconds"] = round(time.monotonic() - t, 2)
    t = time.monotonic()
    sc.tl.leiden(adata, resolution=resolution, flavor="igraph", n_iterations=2, directed=False,
                 random_state=seed, key_added="leiden")
    info["leiden_seconds"] = round(time.monotonic() - t, 2)
    raw = adata.obs["leiden"].astype(str).to_numpy()
    sizes = {}
    for c in raw:
        sizes[c] = sizes.get(c, 0) + 1
    order = sorted(sizes, key=lambda c: (-sizes[c], c))
    width = max(2, len(str(len(order))))
    rename = {c: f"C{i + 1:0{width}d}" for i, c in enumerate(order)}
    labels = [rename[c] for c in raw]
    del adata.obs["leiden"]
    info["n_clusters"] = len(order)
    return labels, info


def umap_subsample(adata, rep: str, max_obs: int, n_neighbors: int, seed: int):
    """UMAP of at most max_obs seeded observations (their own neighbour graph). Returns
    (indices into adata, (n, 2) coordinates)."""
    import numpy as np
    import scanpy as sc

    rng = np.random.default_rng(seed)
    n = adata.n_obs
    if n > max_obs:
        idx = np.sort(rng.choice(n, size=max_obs, replace=False))
        sub = adata[idx].copy()
    else:
        idx = np.arange(n)
        sub = adata
    if sub.n_obs < 5:
        return idx, np.zeros((sub.n_obs, 2), dtype=np.float32)
    if sub is not adata:
        sc.pp.neighbors(sub, n_neighbors=int(min(n_neighbors, sub.n_obs - 1)), use_rep=rep, random_state=seed)
    sc.tl.umap(sub, random_state=seed)
    coords = np.asarray(sub.obsm["X_umap"], dtype=np.float32)
    if sub is adata:
        del adata.obsm["X_umap"]
    return idx, coords


def cluster_markers(adata, cluster_key: str, seed: int, max_obs: int = MAX_OBS_FOR_MARKERS,
                    n_top: int = N_TOP_MARKERS) -> list[dict]:
    """Top genes per cluster (Wilcoxon, each cluster vs the rest, BH-adjusted), with the fraction of
    cells expressing the gene inside and outside the cluster."""
    import numpy as np
    import scanpy as sc

    n = adata.n_obs
    if n > max_obs:
        idx = np.sort(np.random.default_rng(seed).choice(n, size=max_obs, replace=False))
        sub = adata[idx].copy()
    else:
        sub = adata
    counts = sub.obs[cluster_key].value_counts()
    groups = [c for c in sorted(counts.index) if counts[c] >= 3]
    if len(groups) < 2:
        return []
    sub.obs[cluster_key] = sub.obs[cluster_key].astype("category")
    sc.tl.rank_genes_groups(sub, cluster_key, groups=groups, reference="rest", method="wilcoxon", pts=True,
                            use_raw=False, key_added="rgg")
    rows = []
    for g in groups:
        df = sc.get.rank_genes_groups_df(sub, group=g, key="rgg").head(n_top)
        for rank, r in enumerate(df.itertuples(index=False), start=1):
            rows.append({"cluster": g, "gene": r.names, "logfoldchange": float(r.logfoldchanges),
                         "pval_adj": float(r.pvals_adj), "score": float(r.scores),
                         "pct_in": float(getattr(r, "pct_nz_group", float("nan"))),
                         "pct_out": float(getattr(r, "pct_nz_reference", float("nan"))), "rank": rank})
    if sub is adata:
        del adata.uns["rgg"]
    return rows
