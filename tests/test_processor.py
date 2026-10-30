"""The synthetic fixture through the processor (both input layouts), the marker table, the
annotation rules and parameter handling. No Docker."""
from __future__ import annotations

import csv
import json
import shutil
from pathlib import Path

import numpy as np
import pytest

from processor import main as app
from processor import markers as M
from tests.make_fixture import PLANTED, make_fixture, write_markers_csv

LABEL_FILES = {"labels.csv", "cluster_markers.csv", "cell_type_summary.csv", "cluster_summary.csv",
               "cluster_annotation.csv", "scores.csv", "agreement.csv", "confusion.csv", "summary.json", "provenance.json"}
FIGURE_FILES = {"umap_cluster.png", "umap_cell_type.png", "s1_cell_type.png", "s2_cell_type.png",
                "dotplot_marker_sets.png", "cell_type_by_sample.png", "summary.json", "provenance.json"}


@pytest.fixture(scope="session")
def fixture_root(tmp_path_factory) -> dict[str, Path]:
    return make_fixture(tmp_path_factory.mktemp("fixture"), "both")


@pytest.fixture(scope="session")
def mouse_root(tmp_path_factory) -> Path:
    return make_fixture(tmp_path_factory.mktemp("mouse"), "linked", species="mouse")["linked"]


def run(input_dir: Path, output_dir: Path, monkeypatch, **params) -> int:
    monkeypatch.setenv("INPUT_DIR", str(input_dir))
    monkeypatch.setenv("OUTPUT_DIR", str(output_dir))
    monkeypatch.setenv("MPLCONFIGDIR", "/tmp/matplotlib")
    for name in app.PARAMETERS:
        monkeypatch.delenv(name, raising=False)
    monkeypatch.setenv("THREADS", "2")
    for name, value in params.items():
        monkeypatch.setenv(name, str(value))
    return app.main()


def rows(path: Path) -> list[dict[str, str]]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def truth_of(root: Path) -> dict[str, str]:
    truth = root / "truth.csv" if (root / "truth.csv").exists() else root.parent / "linked" / "truth.csv"
    return {r["obs_id"]: r["cell_type"] for r in rows(truth)}


def agreement(labels: list[dict], truth: dict[str, str]) -> tuple[float, float]:
    """(fraction of cells whose cell_type is the planted type, ARI of cluster vs planted type)."""
    from sklearn.metrics import adjusted_rand_score

    planted = [truth[r["obs_id"]] for r in labels]
    hit = np.mean([r["cell_type"] == t for r, t in zip(labels, planted)])
    ari = adjusted_rand_score(planted, [r["cluster"] for r in labels])
    return float(hit), float(ari)


@pytest.mark.parametrize("layout", ["linked", "flat"])
def test_planted_types_are_recovered(fixture_root, tmp_path, monkeypatch, layout):
    out = tmp_path / "out"
    assert run(fixture_root[layout], out, monkeypatch) == 0
    assert sorted(p.name for p in out.iterdir()) == ["figures", "labels"]
    assert LABEL_FILES <= {p.name for p in (out / "labels").iterdir()}
    assert FIGURE_FILES <= {p.name for p in (out / "figures").iterdir()}

    labels = rows(out / "labels" / "labels.csv")
    assert list(labels[0]) == ["sample_id", "obs_id", "cluster", "cell_type"]
    assert len(labels) == 4000
    assert all(r["cluster"].startswith("C") and r["cluster"][1:].isdigit() for r in labels)
    hit, ari = agreement(labels, truth_of(fixture_root[layout]))
    assert hit >= 0.95, f"only {hit:.3f} of the cells got their planted type"
    assert ari >= 0.6, f"cluster ARI {ari:.3f}"

    summary = json.loads((out / "labels" / "summary.json").read_text())
    assert summary["n_samples"] == 2 and summary["n_obs"] == 4000 and summary["n_genes"] == 60
    assert summary["marker_sets"]["source"] == "builtin"
    assert set(PLANTED) <= set(summary["marker_sets"]["used"])
    # Sets whose genes are not in the panel are skipped, not scored.
    assert "Neuron" in summary["marker_sets"]["skipped"]
    assert summary["embedding"]["integration"] == "harmony by sample_id"
    assert summary["embedding"]["pca"].startswith("computed on all 60 genes")   # per-sample PCA is not reused
    assert summary["fraction_unassigned"] < 0.05
    found = {c["cell_type"] for c in summary["cell_types"]}
    assert set(PLANTED) <= found

    # The planted truth is in obs["annotation"], so the app reports its own agreement.
    agree = {(r["sample_id"], r["label"]): float(r["ari"]) for r in rows(out / "labels" / "agreement.csv")}
    assert agree[("all", "cell_type")] >= 0.9
    assert agree[("all", "cluster")] == pytest.approx(ari, abs=1e-6)

    ann = rows(out / "labels" / "cluster_annotation.csv")
    assert list(ann[0]) == ["cluster", "cell_type", "score", "n_obs", "second_best", "second_score"]
    assert sum(int(r["n_obs"]) for r in ann) == 4000
    markers = rows(out / "labels" / "cluster_markers.csv")
    assert list(markers[0])[:7] == ["cluster", "gene", "logfoldchange", "pval_adj", "pct_in", "pct_out", "rank"]
    # The top marker of a cluster called T cell is one of the planted T-cell genes.
    t_clusters = {r["cluster"] for r in ann if r["cell_type"] == "T cell"}
    tops = {r["gene"] for r in markers if r["cluster"] in t_clusters and r["rank"] == "1"}
    assert tops and tops <= set(PLANTED["T cell"])

    scores = rows(out / "labels" / "scores.csv")
    assert list(scores[0])[:2] == ["sample_id", "obs_id"] and set(PLANTED) <= set(scores[0])
    per_sample = rows(out / "labels" / "cell_type_summary.csv")
    assert list(per_sample[0]) == ["sample_id", "cell_type", "n_obs", "fraction"]
    assert abs(sum(float(r["fraction"]) for r in per_sample if r["sample_id"] == "s1") - 1.0) < 1e-6

    prov = json.loads((out / "labels" / "provenance.json").read_text())
    assert prov["app"] == "st-celltype" and prov["version"] == "0.1.0"
    assert prov["n_ok"] == 2 and prov["n_failed"] == 0
    assert {"scanpy", "anndata", "harmonypy"} <= set(prov["tools"])
    assert (out / "figures" / "provenance.json").exists()


def test_builtin_markers_resolve_mouse_symbols(mouse_root, tmp_path, monkeypatch):
    out = tmp_path / "out"
    assert run(mouse_root, out, monkeypatch, INTEGRATION="none") == 0
    labels = rows(out / "labels" / "labels.csv")
    hit, _ = agreement(labels, truth_of(mouse_root))
    assert hit >= 0.95
    summary = json.loads((out / "labels" / "summary.json").read_text())
    assert summary["marker_sets"]["used"]["T cell"]["genes_used"] == ["Cd3g", "Trac", "Cd2"]
    assert summary["embedding"]["integration"] == "none"


def test_resolve_is_case_insensitive():
    sets = [M.MarkerSet("T cell", ["CD3E", "CD3D", "CD2"]), M.MarkerSet("B cell", ["MS4A1", "CD79A"]),
            M.MarkerSet("Neuron", ["SNAP25", "SYT1"])]
    kept, skipped = M.resolve(sets, ["Cd3e", "cd3d", "Ms4a1", "Snap25", "Other"], min_genes_present=2)
    assert [s.cell_type for s in kept] == ["T cell", "B cell"] or [s.cell_type for s in kept] == ["T cell"]
    assert kept[0].present == ["Cd3e", "cd3d"] and kept[0].missing == ["CD2"]
    assert "Neuron" in {s.cell_type for s in skipped}
    kept2, _ = M.resolve([M.MarkerSet("B cell", ["MS4A1", "CD79A"])], ["Ms4a1"], min_genes_present=1)
    assert kept2[0].present == ["Ms4a1"]


@pytest.mark.parametrize("shape", ["long", "wide"])
def test_markers_csv_shapes(tmp_path, shape):
    path = write_markers_csv(tmp_path / "markers.csv", shape, {"My T": ["CD3E", "CD3D", "CD3E"], "My B": ["MS4A1", "CD79A"]})
    sets = M.read_markers_csv(path)
    assert [(s.cell_type, s.genes) for s in sets] == [("My T", ["CD3E", "CD3D"]), ("My B", ["MS4A1", "CD79A"])]


def test_markers_csv_header_variants(tmp_path):
    p = tmp_path / "markers.csv"
    p.write_text("Cell_Type,Genes\nFoo,\"A; B ,C\"\n,X\nBar,\n")
    sets = M.read_markers_csv(p)
    assert [(s.cell_type, s.genes) for s in sets] == [("Foo", ["A", "B", "C"])]
    p.write_text("name,stuff\nFoo,A\n")
    with pytest.raises(ValueError):
        M.read_markers_csv(p)


@pytest.mark.parametrize("shape", ["long", "wide"])
def test_custom_markers_port(fixture_root, tmp_path, monkeypatch, shape):
    # A markers port with custom names, lacking the fibroblast set: fibroblasts become Unassigned.
    src = tmp_path / "in"
    shutil.copytree(fixture_root["linked"], src)
    custom = {f"my {t}": g for t, g in PLANTED.items() if t != "Fibroblast"}
    write_markers_csv(src / "markers" / "markers.csv", shape, custom)
    out = tmp_path / "out"
    assert run(src, out, monkeypatch) == 0
    labels = rows(out / "labels" / "labels.csv")
    truth = truth_of(fixture_root["linked"])
    expected = {t: f"my {t}" for t in custom and PLANTED if t != "Fibroblast"}
    expected["Fibroblast"] = "Unassigned"
    hit = np.mean([r["cell_type"] == expected[truth[r["obs_id"]]] for r in labels])
    assert hit >= 0.95, f"{hit:.3f}"
    summary = json.loads((out / "labels" / "summary.json").read_text())
    assert summary["marker_sets"]["source"].startswith("markers.csv")
    assert set(summary["marker_sets"]["used"]) == set(custom)
    assert summary["fraction_unassigned"] == pytest.approx(
        np.mean([truth[r["obs_id"]] == "Fibroblast" for r in labels]), abs=0.05)


def test_unassigned_when_nothing_scores(fixture_root, tmp_path, monkeypatch):
    out = tmp_path / "out"
    assert run(fixture_root["linked"], out, monkeypatch, MIN_SCORE="100") == 0
    labels = rows(out / "labels" / "labels.csv")
    assert {r["cell_type"] for r in labels} == {"Unassigned"}
    ann = rows(out / "labels" / "cluster_annotation.csv")
    assert all(r["cell_type"] == "Unassigned" and r["second_best"] for r in ann)   # the scores are still reported
    summary = json.loads((out / "labels" / "summary.json").read_text())
    assert summary["n_cell_types"] == 0 and summary["fraction_unassigned"] == 1.0


def test_marker_set_none_gives_clusters_only(fixture_root, tmp_path, monkeypatch):
    out = tmp_path / "out"
    assert run(fixture_root["linked"], out, monkeypatch, MARKER_SET="none", MAX_OBS_FOR_UMAP="500") == 0
    labels = rows(out / "labels" / "labels.csv")
    assert {r["cell_type"] for r in labels} == {"Unassigned"}
    assert list(rows(out / "labels" / "scores.csv")[0]) == ["sample_id", "obs_id"]
    assert not (out / "figures" / "dotplot_marker_sets.png").exists()
    _, ari = agreement(labels, truth_of(fixture_root["linked"]))
    assert ari >= 0.6


def test_annotate_per_cell(fixture_root, tmp_path, monkeypatch):
    out = tmp_path / "out"
    assert run(fixture_root["linked"], out, monkeypatch, ANNOTATE="cell", MIN_GENES_PRESENT="3") == 0
    labels = rows(out / "labels" / "labels.csv")
    truth = truth_of(fixture_root["linked"])
    summary = json.loads((out / "labels" / "summary.json").read_text())
    assert summary["parameters"]["ANNOTATE"] == "cell"
    # With 3 genes needed, the sets with only 2 panel genes drop out: the built-in plasma-cell set
    # has only MZB1 and JCHAIN in the fixture panel, so planted plasma cells cannot be called.
    assert all(len(v["genes_used"]) >= 3 for v in summary["marker_sets"]["used"].values())
    assert "Plasma cell" in summary["marker_sets"]["skipped"]
    kept = [r for r in labels if truth[r["obs_id"]] != "Plasma cell"]
    hit, _ = agreement(kept, truth)
    assert hit >= 0.9
    assert not any(r["cell_type"] == "Plasma cell" for r in labels)


def test_annotation_rules():
    import pandas as pd

    z = pd.DataFrame({"A": [2.0, 2.0, -1.0, -1.0, 0.05], "B": [-1.0, -1.0, 1.5, 1.5, 0.0]},
                     index=[f"c{i}" for i in range(5)])
    clusters = ["C1", "C1", "C2", "C2", "C3"]
    calls, table = M.annotate_by_cluster(z, clusters, 0.1)
    assert calls == ["A", "A", "B", "B", "Unassigned"]
    assert [(r["cluster"], r["cell_type"], r["second_best"]) for r in table] == [("C1", "A", "B"), ("C2", "B", "A"), ("C3", "Unassigned", "B")]
    assert table[0]["score"] == pytest.approx(2.0) and table[0]["second_score"] == pytest.approx(-1.0)
    calls, table = M.annotate_by_cell(z, clusters, 0.1)
    assert calls == ["A", "A", "B", "B", "Unassigned"]
    assert table[2]["cell_type"] == "Unassigned"
    empty = pd.DataFrame(index=z.index)
    assert M.annotate_by_cluster(empty, clusters, 0.1)[0] == ["Unassigned"] * 5
    assert M.annotate_by_cell(empty, clusters, 0.1)[0] == ["Unassigned"] * 5


def test_joint_pca_is_reused(fixture_root):
    from processor import cluster as C

    a = C.read_sample(fixture_root["linked"] / "spatial" / "s1.h5ad", "s1")
    b = C.read_sample(fixture_root["linked"] / "spatial" / "s2.h5ad", "s2")
    assert "X_pca" in a.obsm and "PCs" in a.varm
    merged, joint = C.concatenate([a.copy(), b.copy()])
    assert not joint and "X_pca" not in merged.obsm
    b.varm["PCs"] = a.varm["PCs"].copy()   # the same loadings in every sample: one joint embedding
    merged, joint = C.concatenate([a, b])
    assert joint and merged.obsm["X_pca"].shape == (4000, 20)
    assert merged.n_vars == 60 and bool(merged.var["highly_variable"].all())


@pytest.mark.parametrize("params", [{"INTEGRATION": "scvi"}, {"ANNOTATE": "both"}, {"RESOLUTION": "0"},
                                    {"N_NEIGHBORS": "x"}, {"MARKER_SET": "panglao"}])
def test_bad_parameters(fixture_root, tmp_path, monkeypatch, params):
    assert run(fixture_root["linked"], tmp_path / "out", monkeypatch, **params) == 2


def test_no_input_fails(tmp_path, monkeypatch):
    (tmp_path / "in").mkdir()
    assert run(tmp_path / "in", tmp_path / "out", monkeypatch) == 1


def test_a_bad_file_is_skipped(fixture_root, tmp_path, monkeypatch):
    src = tmp_path / "in"
    shutil.copytree(fixture_root["linked"], src)
    (src / "spatial" / "s2.h5ad").write_bytes(b"not an h5ad")
    out = tmp_path / "out"
    assert run(src, out, monkeypatch, MAX_OBS_FOR_UMAP="500") == 0
    prov = json.loads((out / "labels" / "provenance.json").read_text())
    assert prov["n_failed"] == 1 and prov["n_ok"] == 1
    assert len(rows(out / "labels" / "labels.csv")) == 2000
