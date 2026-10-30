# Third-party software

This app's own code is Apache-2.0 (see `LICENSE`). The image it builds contains the following
third-party components, installed from PyPI at the pinned versions in `requirements.txt`; their
licences apply to the image.

**GPL notice.** The image includes, and the app imports in-process, two GPL libraries:
**igraph** (python-igraph, GPL-2.0-or-later; the Leiden implementation behind
`scanpy.tl.leiden(flavor="igraph")`) and **harmonypy** (GPL-3.0-or-later since version 2.0, with a
compiled C++ core; Harmony integration). The image as a whole is therefore distributed under the
terms of the GPL-3.0 (the combination of GPL-2.0-or-later and GPL-3.0-or-later code). Our source
code stays Apache-2.0, which is GPL-3.0-compatible; no GPL code is copied into it, the libraries
are only called. leidenalg is not installed.

| Component | What | Licence | Source |
|---|---|---|---|
| scanpy 1.12.5 | `pca`, `neighbors`, `leiden`, `umap`, `score_genes`, `rank_genes_groups` | BSD-3-Clause | https://github.com/scverse/scanpy |
| anndata 0.13.4 | Reads the `.h5ad` files (element by element, so `layers["counts"]` is never loaded). | BSD-3-Clause | https://github.com/scverse/anndata |
| igraph 1.0.0 (python-igraph), texttable 1.7.2 | Leiden community detection (`community_leiden`). | GPL-2.0-or-later; MIT | https://github.com/igraph/python-igraph |
| harmonypy 2.1.0 | Harmony integration of the PCA by sample. | GPL-3.0-or-later | https://github.com/slowkow/harmonypy |
| umap-learn 0.5.12, pynndescent 0.7.0 | UMAP of the figures; approximate nearest neighbours. | BSD-3-Clause; BSD-2-Clause | https://github.com/lmcinnes/umap |
| scikit-learn 1.9.1 | PCA solver (through scanpy), ARI / NMI against `obs["annotation"]`. | BSD-3-Clause | https://scikit-learn.org |
| NumPy 2.5.4, SciPy 1.18.1 | Arrays, CSR sparse matrices, the Wilcoxon test. | BSD-3-Clause (NumPy also bundles 0BSD/MIT/Zlib/CC0 parts) | https://numpy.org, https://scipy.org |
| pandas 3.1.0 | Tables. | BSD-3-Clause | https://pandas.pydata.org |
| h5py 3.16.0 | HDF5 backend of anndata. | BSD-3-Clause | https://www.h5py.org |
| Matplotlib 3.11.2, Pillow 13.0.0 | PNG figures. | Matplotlib licence (PSF-based, BSD-compatible); MIT-CMU | https://matplotlib.org, https://python-pillow.org |
| numba 0.68.1, llvmlite 0.50.0 | Compiled kernels used by scanpy and UMAP. | BSD-2-Clause (llvmlite also Apache-2.0 with LLVM exception) | https://numba.pydata.org |
| statsmodels 0.15.0, patsy 1.0.3, formulaic, narwhals, interface_meta, seaborn 0.13.2, networkx 3.7, joblib 1.6.0, cloudpickle, threadpoolctl 3.7.0, natsort 8.4.0 | Dependencies of scanpy / scikit-learn. | BSD-3-Clause / BSD-2-Clause (patsy) / MIT (formulaic, narwhals, interface_meta, natsort) | PyPI |
| fast-array-utils 1.5.1, legacy-api-wrap 1.5, session-info2 0.4.2, scverse-misc 0.1.8 | scverse helpers pulled in by scanpy/anndata. | MPL-2.0 (scverse-misc BSD-3-Clause) | https://github.com/scverse |
| zarr 3.4.2, numcodecs 0.17.1, donfig 0.8.1, google-crc32c 1.9.0, msgspec, PyYAML, pydantic, pydantic-settings, pydantic-core, python-dotenv, annotated-types, typing-inspection | Dependencies of anndata (zarr backend, unused here) and of those. | MIT / BSD-3-Clause (msgspec, python-dotenv) / Apache-2.0 (google-crc32c) | PyPI |
| array-api-compat, certifi, contourpy, cycler, fonttools, kiwisolver, packaging, pyparsing, python-dateutil, six, tqdm, typing_extensions, tzdata, wrapt | Dependencies of the above (versions in `requirements.txt`). | MIT / BSD / Apache-2.0 or BSD / MPL-2.0 (certifi) / MPL-2.0+MIT (tqdm) / PSF-2.0 | PyPI |
| python:3.12.16-slim-bookworm | Base image (Debian 12 + CPython 3.12). | PSF (Python); Debian packages under their own free licences | https://hub.docker.com/_/python |

## Data bundled in the image

- `processor/builtin_markers.json`: 37 marker gene sets selected by hand from the PanglaoDB
  marker table (27 March 2020 release; Franzén O, Gan LM, Björkegren JLM. *Database* 2019;
  baz046, doi:10.1093/database/baz046). Only gene symbols and cell-type names are used (facts,
  not the database itself); PanglaoDB is credited in the README and `CITATION.cff`.

## Not bundled

- **Demo data** (not in the image or the repository): 10x Genomics *FFPE Human Breast Cancer*
  Xenium replicates 1 and 2 (CC BY 4.0; Janesick et al. 2023, doi:10.1038/s41467-023-43458-x).
