"""
Finding the inputs.

Linked layout (placement: subfolder): INPUT_DIR/spatial/<sample_id>.h5ad,
INPUT_DIR/metadata/samples.csv and, optionally, INPUT_DIR/markers/markers.csv.
Without links the same files are found by name anywhere under INPUT_DIR
(symlinks followed; INPUT_DIR entries are symlinks).
"""
from __future__ import annotations

import logging
import os
from dataclasses import dataclass, field
from pathlib import Path

log = logging.getLogger("st-celltype")


@dataclass
class Sample:
    sample_id: str
    path: Path
    row: dict[str, str] = field(default_factory=dict)


@dataclass
class Inputs:
    samples: list[Sample]
    samples_csv: Path | None
    markers_csv: Path | None
    h5ad_seen: list[Path]
    skipped: dict[str, str]        # sample_id or file -> reason


def walk(root: Path) -> list[Path]:
    """Every file under root, following symlinks, sorted."""
    found: list[Path] = []
    for dirpath, dirnames, filenames in os.walk(root, followlinks=True):
        dirnames[:] = sorted(d for d in dirnames if not d.startswith("."))
        found.extend(Path(dirpath) / f for f in sorted(filenames) if not f.startswith("."))
    return found


def find(root: Path, port: str, predicate) -> list[Path]:
    """Files of a port: its own folder when it exists, else anywhere under INPUT_DIR."""
    folder = root / port
    base = folder if folder.is_dir() else root
    return [p for p in walk(base) if predicate(p)]


def one(root: Path, port: str, name: str) -> Path | None:
    hits = find(root, port, lambda p: p.name.lower() == name)
    if len(hits) > 1:
        log.warning("%d files named %s; using %s", len(hits), name, hits[0])
    return hits[0] if hits else None


def read_samples_csv(path: Path) -> list[dict[str, str]]:
    import csv

    with open(path, newline="", encoding="utf-8-sig") as fh:
        rows = list(csv.DictReader(fh))
    if not rows or "sample_id" not in rows[0]:
        raise ValueError(f"{path.name}: no sample_id column")
    return [{k.strip(): (v or "").strip() for k, v in r.items() if k} for r in rows]


def discover(input_dir: Path) -> Inputs:
    h5ads = find(input_dir, "spatial", lambda p: p.suffix.lower() == ".h5ad")
    samples_csv = one(input_dir, "metadata", "samples.csv")
    markers_csv = one(input_dir, "markers", "markers.csv")
    skipped: dict[str, str] = {}
    by_stem: dict[str, Path] = {}
    for p in h5ads:
        if p.stem in by_stem:
            log.warning("two .h5ad files for sample %s; using %s", p.stem, by_stem[p.stem])
            continue
        by_stem[p.stem] = p

    samples: list[Sample] = []
    if samples_csv is not None:
        listed = set()
        for row in read_samples_csv(samples_csv):
            sid = row["sample_id"]
            listed.add(sid)
            if sid not in by_stem:
                skipped[sid] = "listed in samples.csv but no .h5ad found"
                continue
            samples.append(Sample(sid, by_stem[sid], row))
        for stem in by_stem:
            if stem not in listed:
                skipped[stem] = ".h5ad not listed in samples.csv"
    else:
        log.warning("no samples.csv: every .h5ad found is a sample")
        samples = [Sample(stem, path, {"sample_id": stem}) for stem, path in sorted(by_stem.items())]
    return Inputs(samples, samples_csv, markers_csv, h5ads, skipped)
