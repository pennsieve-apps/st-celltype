# st-celltype

Clusters and cell types for **single-cell spatial transcriptomics**: 10x Xenium, Vizgen
MERSCOPE and NanoString CosMx (one observation = one segmented cell), and also Visium HD bins.
It is the third step of the *Single-cell Spatial Profiling* workflow
(`st-prep → st-qc → st-celltype → {st-niches, st-ccc} → st-report`).

**For a biologist.** Imaging panels measure a few hundred to a few thousand genes per cell. This
app groups cells with similar expression into clusters (`C01`, `C02`, ... largest first) across
all your slides at once, so the same cluster means the same thing on every slide, and names each
cluster after the cell type whose canonical marker genes it expresses most (T cell, fibroblast,
endothelial, ...). You can give your own marker table; otherwise a built-in table of 37 cell types
taken from PanglaoDB is used, keeping only the cell types with at least two marker genes in your
panel. Each cell then carries a `cluster` and a `cell_type` that the niche and cell-cell
communication apps use, with maps of the cell types in tissue space.

**For a computational analyst.** All samples are concatenated on their shared genes (`X` read
element by element; `layers["counts"]` stays on disk); a joint PCA is computed on the
log-normalized `X` (all genes for panels ≤ 1,000 genes, else the union of st-qc's highly variable
genes; st-qc's per-sample PCA is reused only when every sample carries the same loadings), Harmony
corrects the PCs for `sample_id` (when > 1 sample), then `scanpy.pp.neighbors` and
`scanpy.tl.leiden(flavor="igraph", n_iterations=2)`. Marker sets are scored per observation with
`scanpy.tl.score_genes` (mean expression of the set minus a matched random control set, the Seurat
approach of Satija et al. 2015); each set's score is z-scored over all observations, and the call
is the set with the highest mean z per cluster (`ANNOTATE=cluster`) or per cell
(`ANNOTATE=cell`), `Unassigned` below `MIN_SCORE`. Cluster markers come from a Wilcoxon test of each
cluster against the rest.

## Method

1. **Inputs.** `spatial/<sample_id>.h5ad` from st-qc and `metadata/samples.csv` (samples and
   their order); an `.h5ad` not listed in samples.csv, or a listed sample without a file, is
   logged and skipped. A sample whose file cannot be read is skipped; the run fails (exit 1) only
   when no sample is usable.
2. **Joint PCA** (`N_PCS`, seeded), **Harmony** (`harmonypy.run_harmony` on `sample_id`,
   `INTEGRATION=harmony`), **neighbours** (`N_NEIGHBORS`), **Leiden** (`RESOLUTION`; igraph
   implementation, 2 iterations). Clusters are renamed `C01..` by decreasing size.
3. **Marker sets.** From `markers/markers.csv` when linked (it wins over `MARKER_SET`), else the
   built-in table (`MARKER_SET=builtin`), or none (`MARKER_SET=none`: clusters only). Symbols are
   matched to `var_names` case-insensitively (`CD3E` matches mouse `Cd3e`); a set with fewer than
   `MIN_GENES_PRESENT` genes in the panel is skipped and listed in `summary.json`.
4. **Scores.** `scanpy.tl.score_genes` per set, the control set sized to the panel
   (`ctrl_size = min(50, max(5, n_genes / 10))`, as many expression bins), then z-scored per set.
5. **Call.** `cluster`: the set with the highest mean z over the cluster's cells, with the
   runner-up in `cluster_annotation.csv`; `cell`: the best set per cell (the cluster table then
   gives each cluster's majority call). Below `MIN_SCORE`: `Unassigned`.
6. **Cluster markers.** `scanpy.tl.rank_genes_groups(method="wilcoxon", reference="rest")` on at
   most 100,000 seeded cells, top 25 per cluster with the fraction expressing inside/outside.
7. **Agreement.** When the input has `obs["annotation"]`, ARI and NMI of `cluster` and
   `cell_type` against it (per sample and overall) and a confusion table.
8. **Figures.** UMAP (on a seeded subsample of at most `MAX_OBS_FOR_UMAP` cells with its own
   neighbour graph) by cluster and by cell type; per sample, cell types in space (`obsm["spatial"]`,
   y inverted; on the lowres H&E for Visium HD when present); a dotplot of the mean z-scored
   marker-set score per cluster; cell-type composition per sample.

### The built-in marker table

`processor/builtin_markers.json`: 37 cell types (epithelial, luminal hormone-sensing,
basal/myoepithelial, T / CD8 T / CD4 T / regulatory T, B, plasma, NK, macrophage, monocyte, cDC1,
cDC2, plasmacytoid DC, mast, neutrophil, erythrocyte, fibroblast, smooth muscle/pericyte,
endothelial, lymphatic endothelial, adipocyte, and neural, liver, gut, lung, kidney, heart and
skin types), 3–8 genes each, 188 genes in all. Every gene was picked by hand among the **human
markers PanglaoDB lists for that cell type** (Franzén et al. 2019; marker table of 27 March 2020),
preferring canonical, specific (low ubiquitousness) genes that targeted panels tend to include;
each set records its PanglaoDB cell type(s). `scripts/verify_builtin_markers.py` re-checks every
gene against the PanglaoDB table (read through OmniPath, Türei et al. 2021, since panglaodb.se
refuses scripted downloads): 188 of 188 found. Proliferation ("cycling") is a state, not a
lineage, and is not in the table; add it in a markers.csv if you want it.

## Ports

| Port | Direction | Content |
|---|---|---|
| `spatial` | input, required | st-qc's `<sample_id>.h5ad` (log-normalized `X`, `layers["counts"]`, `obsm["spatial"]`, `uns["pennsieve"]`) |
| `metadata` | input, required | st-prep's `samples.csv` (`prep.json` ignored) |
| `markers` | input, optional | `markers.csv`: long (`cell_type,gene`) or wide (`cell_type,genes` with `;`-separated genes) |
| `labels` | output | `labels/` tables below |
| `figures` | output | `figures/` PNGs below |

All inputs use `placement: subfolder`; without links the files are found by name anywhere under
`INPUT_DIR` (`*.h5ad`, `samples.csv`, `markers.csv`).

## Parameters

| Name | Default | What it does |
|---|---|---|
| `RESOLUTION` | `0.5` | Leiden resolution. Higher: more, smaller clusters (1–2 splits immune and stromal subtypes; see results). |
| `N_NEIGHBORS` | `15` | Neighbours per cell in the expression graph. |
| `N_PCS` | `30` | Principal components for Harmony and the graph (capped by genes and cells). |
| `INTEGRATION` | `harmony` | `harmony`: correct the PCs for `sample_id` (only with > 1 sample). `none`: no correction. |
| `MARKER_SET` | `builtin` | `builtin` table or `none` (clusters only). A linked `markers.csv` always wins. |
| `ANNOTATE` | `cluster` | `cluster`: one call per cluster. `cell`: one call per cell (noisier). |
| `MIN_SCORE` | `0.1` | Best mean z-score below this gives `Unassigned`. |
| `MIN_GENES_PRESENT` | `2` | Marker genes a set needs in the panel to be scored. |
| `MAX_OBS_FOR_UMAP` | `50000` | Cells in the UMAP figures (seeded subsample). |
| `SEED` | `42` | Seed for PCA, Harmony, Leiden, UMAP and subsamples. |
| `THREADS` | `0` | BLAS / numba / Harmony threads; 0 = the container's CPU quota (4 without a quota). |

Bad values exit with code 2.

## Outputs

`labels/` (single-run app: at the port root):

| File | Columns / content |
|---|---|
| `labels.csv` | `sample_id, obs_id, cluster, cell_type`: every observation of every sample |
| `cluster_markers.csv` | `cluster, gene, logfoldchange, pval_adj, pct_in, pct_out, rank, score` (Wilcoxon, BH-adjusted) |
| `cell_type_summary.csv` | `sample_id, cell_type, n_obs, fraction` |
| `cluster_summary.csv` | `sample_id, cluster, n_obs, fraction` |
| `cluster_annotation.csv` | `cluster, cell_type, score, n_obs, second_best, second_score` (`score` = mean z of the call; with `ANNOTATE=cell`, `second_*` is the runner-up call and its fraction) |
| `scores.csv` | `sample_id, obs_id`, one z-scored column per marker set scored |
| `agreement.csv`, `confusion.csv` | Only when `obs["annotation"]` exists: `sample_id, label, ari, nmi, n_obs`; `label, value, annotation, n_obs` |
| `summary.json` | counts, fraction unassigned, cell types overall and by sample, cluster calls, marker sets used / skipped with genes found and missing, embedding and clustering details |
| `provenance.json` | app version, tool versions, parameters, method, inputs seen, timings, per-sample status, peak RSS |

`figures/`: `umap_cluster.png`, `umap_cell_type.png`, `<sample_id>_cell_type.png`,
`dotplot_marker_sets.png`, `cell_type_by_sample.png`, `summary.json`, `provenance.json`.

## Results on the demo data

10x *FFPE Human Breast Cancer* Xenium, replicates 1 and 2 (313 genes), st-prep → st-qc output
(`data/work/xenium/`): 277,475 cells after QC. Defaults: Harmony by slide, 10 clusters, 8 cell
types, 0 % Unassigned:

| Cluster | Cells | Call (mean z) | Runner-up |
|---|---|---|---|
| C01, C05, C07 | 54,722 / 24,797 / 20,412 | Epithelial (1.07 / 1.04 / 1.23) | Luminal hormone-sensing (0.78–0.88) |
| C02 | 51,500 | T cell (1.74) | CD4 T cell (1.40) |
| C03 | 42,614 | Fibroblast (1.59) | Smooth muscle/pericyte |
| C04 | 28,394 | Macrophage (1.99) | Monocyte |
| C06 | 22,462 | Endothelial (2.59) | Smooth muscle/pericyte |
| C08 | 17,117 | Basal/myoepithelial (2.69) | Smooth muscle/pericyte |
| C09 | 8,289 | Plasma cell (3.54) | B cell |
| C10 | 7,168 | B cell (4.44) | Plasmacytoid DC |

The spatial maps show what the tissue should: epithelial (tumour) nests and ducts rimmed by
basal/myoepithelial cells (the DCIS ducts keep their myoepithelial layer), fibroblast-rich stroma
with endothelium, and dense T/B/plasma-cell infiltrates with macrophages; composition is similar
on both replicates (epithelial 40 % / 31 %, T cells 18 % / 19 %, fibroblasts 13 % / 18 %). At
`RESOLUTION=1.5` (stand-in input, 282 k cells, 22 clusters) the same data also gives CD8 T cells
(17 k), a cDC2-scored myeloid cluster, pericytes/smooth muscle and a luminal hormone-sensing
epithelial cluster. Rare populations (mast cells, pDCs, ~1 % or fewer) do not form a cluster of their
own at either resolution and are absorbed by their neighbours; `ANNOTATE=cell` names them but
also over-calls them (per-cell scores on 2–4 panel genes are noisy). `data/work/xenium-20k/`
(20 k cells per slide before QC, 38,858 after): 9 clusters, the same types except that the
B cells do not form a cluster of their own.

On the synthetic fixture (7 planted types, two slides with a depth batch effect) the ARI of both
`cluster` and `cell_type` against the planted truth is 0.996 (Docker run of `make run`).

## Runtime

`computeTypes: [standard]`, 4 vCPU / 8 GiB, 4 h timeout. Clustering is joint across samples, so
all samples are held at once (log-normalized `X` only, CSR float32, plus the PCs).

Measured in the linux/amd64 image under emulation on an M-series Mac (Docker, `--memory 4g`,
4 threads), full Xenium breast demo after st-qc (2 slides, 277,475 cells × 313 genes): **peak
RSS 3.0 GB (cgroup `memory.peak` 3.04 GB), wall time 209 s** (read 3 s, PCA + Harmony +
neighbours + Leiden 155 s, of which the neighbour graph is ~130 s, scores 3 s, cluster markers
10 s, figures 29 s). The 20 k-per-slide subsample: 1.37 GB, 172 s (mostly numba compilation
under emulation). Memory grows with cells × detected genes: a 5,000-gene Xenium panel or a Visium
HD 8 µm sample (~500 k × 18 k) holds several times more non-zeros, hence 8 GiB; give such runs
16 GiB. Image: 950 MB on disk (227 MB compressed).

## Running it locally

```sh
python3.12 -m venv .venv && .venv/bin/pip install -r requirements-dev.txt
make test PYTHON=.venv/bin/python       # pytest on the synthetic fixture, no Docker
make run PYTHON=.venv/bin/python        # Docker over data/fixture/linked -> data/output/{labels,figures}
make validate PYTHON=.venv/bin/python   # app.yml and CITATION.cff
```

By hand, on st-qc / st-prep output:

```sh
docker build --platform linux/amd64 -t st-celltype:dev .
docker run --rm --platform linux/amd64 --memory 4g --network none \
  -e INPUT_DIR=/data/input -e OUTPUT_DIR=/data/output -e RESOLUTION=1.0 \
  -v /path/to/qc/spatial:/data/input/spatial:ro -v /path/to/prep/metadata:/data/input/metadata:ro \
  -v "$PWD/data/output:/data/output" st-celltype:dev
```

## Limitations

- **Marker scoring is not a reference mapping.** A call says which marker set a cluster expresses
  most, relative to the other clusters of *this* dataset (scores are z-scored over all cells). It
  does not tell tumour from normal epithelium, activation states, or subtypes the table lacks;
  every cluster gets its best set even when none fits well (check `score` and `second_best`).
  Use a reference-based method or a curated markers.csv when you need finer types.
- **The panel decides what can be called.** Sets with fewer than `MIN_GENES_PRESENT` panel genes
  are skipped (on the 313-gene breast panel 16 of 37: cDC1, lymphatic endothelium and the
  non-breast tissue types). Two-gene sets (cDC2, pDC, neutrophil, erythrocyte, adipocyte) are
  fragile.
- **Rare populations** (< ~1–2 % of cells) usually merge into a larger cluster at the default
  resolution; raise `RESOLUTION` or use `ANNOTATE=cell` (noisier).
- **Segmentation errors** (transcripts from neighbouring cells assigned to a cell) blend signatures,
  e.g. epithelial genes in tumour-infiltrating macrophages; this app does not correct for them.
- **Harmony removes slide-level differences**, including real biology that differs systematically
  between slides; use `INTEGRATION=none` when slides are expected to differ in composition.
- **Overlapping sets** (T cell vs CD4/CD8 T; macrophage vs monocyte; myoepithelial vs smooth
  muscle) are separated only as well as the panel's specific genes allow.
- Visium HD 2 µm bins are out of scope; 8 µm bins are not single cells, so calls are of the
  dominant type in a bin.

## Licences

Our code: Apache-2.0 (`LICENSE`). The image bundles and imports **igraph** (GPL-2.0-or-later)
and **harmonypy** (GPL-3.0-or-later since 2.0), so the image as a whole is distributed under the
GPL-3.0; everything else is BSD/MIT/Apache/MPL/PSF. Details in `THIRD_PARTY.md`. The built-in
markers are gene symbols selected from PanglaoDB (credited below).

## Citations

All verified on Crossref and PubMed (title, authors, year) on 2026-10-30; full author lists in
`CITATION.cff`.

- Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis.
  *Genome Biol* 2018;19:15. doi:10.1186/s13059-017-1382-0. PMID 29409532.
- Virshup I, Rybakov S, Theis FJ, Angerer P, Wolf FA. anndata: Access and store annotated data
  matrices. *J Open Source Softw* 2024;9(101):4371. doi:10.21105/joss.04371.
- Korsunsky I, Millard N, Fan J, et al. Fast, sensitive and accurate integration of single-cell
  data with Harmony. *Nat Methods* 2019;16:1289-1296. doi:10.1038/s41592-019-0619-0. PMID 31740819.
- Traag VA, Waltman L, van Eck NJ. From Louvain to Leiden: guaranteeing well-connected
  communities. *Sci Rep* 2019;9:5233. doi:10.1038/s41598-019-41695-z. PMID 30914743.
- McInnes L, Healy J, Saul N, Großberger L. UMAP: Uniform Manifold Approximation and Projection.
  *J Open Source Softw* 2018;3(29):861. doi:10.21105/joss.00861.
- Satija R, Farrell JA, Gennert D, Schier AF, Regev A. Spatial reconstruction of single-cell gene
  expression data. *Nat Biotechnol* 2015;33:495-502. doi:10.1038/nbt.3192. PMID 25867923.
  (The gene-set score reproduced by `scanpy.tl.score_genes`.)
- Franzén O, Gan LM, Björkegren JLM. PanglaoDB: a web server for exploration of mouse and human
  single-cell RNA sequencing data. *Database* 2019;2019:baz046. doi:10.1093/database/baz046.
  PMID 30951143. (Source of the built-in marker table.)
- Türei D, Valdeolivas A, Gul L, et al. Integrated intra- and intercellular signaling knowledge
  for multicellular omics analysis. *Mol Syst Biol* 2021;17:e9923. doi:10.15252/msb.20209923.
  PMID 33749993. (OmniPath, used to check the marker table.)
- igraph (python-igraph): https://github.com/igraph/python-igraph (Leiden implementation).
- Demo data: 10x Genomics, *FFPE Human Breast Cancer* Xenium replicates 1 and 2 (CC BY 4.0);
  Janesick A, Shelansky R, Gottscho AD, et al. High resolution mapping of the tumor
  microenvironment using integrated single-cell, spatial and in situ analysis. *Nat Commun*
  2023;14:8353. doi:10.1038/s41467-023-43458-x. PMID 38114474.
