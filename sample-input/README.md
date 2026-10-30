# sample-input

This app reads the output of `st-qc` (port `spatial`: `<sample_id>.h5ad` with log-normalized
`X`, `layers["counts"]`, `obsm["spatial"]`, `uns["pennsieve"]`), `st-prep`'s
`metadata/samples.csv`, and optionally a `markers/markers.csv`. The contract is in
`workflow-spatial-transcriptomics/docs/conventions.md`.

## Synthetic fixture (no download)

```sh
make fixture PYTHON=.venv/bin/python   # data/fixture/{linked,flat}/
make run PYTHON=.venv/bin/python       # Docker over data/fixture/linked/ -> data/output/{labels,figures}/
```

Two samples (`s1`, `s2`, the second sequenced 40 % shallower as a batch effect), 2,000 cells
each, a 60-gene panel. Seven cell types (T, B, plasma, macrophage, mast, fibroblast,
endothelial) are planted in spatial patches, each expressing three of its canonical markers; the
truth is in `obs["annotation"]` (so the app reports ARI/NMI) and in `truth.csv`.
`python -m tests.make_fixture DIR --species mouse` writes the same with mouse-style symbols
(`Cd3g`, `Ms4a1`, ...).

## A markers.csv

Either shape works (column names are matched case-insensitively):

```csv
cell_type,gene
T cell,CD3E
T cell,CD3D
B cell,MS4A1
```

```csv
cell_type,genes
T cell,CD3E;CD3D;CD2
B cell,MS4A1;CD79A
```

## Real data

Run `st-prep` and `st-qc` on an imaging-based dataset, for example the workflow's demo data
(`scripts/fetch-sample-data.sh` in the workflow repository):

- 10x Genomics *FFPE Human Breast Cancer* Xenium replicates 1 and 2 (Janesick et al. 2023,
  https://doi.org/10.1038/s41467-023-43458-x; CC BY 4.0): 313-gene panel, ~167 k and ~118 k cells.

Then point this app at their outputs (one folder per input port):

```sh
mkdir -p in && ln -s /path/to/st-qc/output/spatial in/spatial && ln -s /path/to/st-prep/output/metadata in/metadata
docker run --rm --platform linux/amd64 --memory 4g --network none \
  -e INPUT_DIR=/data/input -e OUTPUT_DIR=/data/output \
  -v "$PWD/in:/data/input:ro" -v "$PWD/data/output:/data/output" st-celltype:dev
```

(Docker does not follow symlinks that point outside a mount; mount the two folders directly, e.g.
`-v /path/to/st-qc/output/spatial:/data/input/spatial:ro -v /path/to/st-prep/output/metadata:/data/input/metadata:ro`.)

While st-prep / st-qc are not at hand, `python -m tests.make_fixture data/xenium --from-xenium
<folder with samples.csv and Rep1/ Rep2/>` reads the raw Xenium output, filters cells with < 10
counts, normalizes to 100 counts + log1p and computes a PCA, as a stand-in for the two apps.
