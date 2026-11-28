#!/usr/bin/env python3
"""Manifest loading, validation, graph resolution, and QC merging for BPC assemblies.

Implemented using only the Python 3.10+ standard library so it can be imported
directly by Snakefile or invoked from the CLI without extra dependencies.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import secrets
import sys
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

VALID_FORMAT_MODES: tuple[str, ...] = ("intact", "already_placed", "ragtag_preserve", "ragtag_cut")
DEFAULT_CHROMOSOMES: list[str] = [str(i) for i in range(1, 30)] + ["X", "Y", "MT"]
SAFE_ID_RE = re.compile(r"^[A-Za-z0-9._-]+$")
HEX_SAFE_LETTERS: frozenset[str] = frozenset("ABCDF")

DEFAULT_QC_THRESHOLDS: dict[str, float] = {
    "min_genome_bp": 2_300_000_000.0,
    "max_genome_bp": 3_300_000_000.0,
    "min_ng50": 1_000_000.0,
    "min_autosome_cov": 0.80,
    "min_busco_complete": 90.0,
    "x_cov_female_min": 0.70,
    "y_cov_male_min": 0.15,
}

_TRUTHY = frozenset({"y", "yes", "true", "1", "t"})
_FALSY_OR_PASS = frozenset({"", "n", "no", "false", "0", "f", "pass", "ok", "na", "none"})
_SEX_MAP = {
    "": "",
    "m": "M",
    "male": "M",
    "f": "F",
    "female": "F",
    "u": "U",
    "unk": "U",
    "unknown": "U",
    "na": "U",
}

CANONICAL_TSV_COLUMNS: list[str] = [
    "animal_id",
    "haplotype",
    "sex",
    "species",
    "subspecies",
    "breed",
    "source",
    "format_mode",
    "chrom_map",
    "checksum",
    "qc_fail",
    "exclude_reason",
    "graphs",
    "tags",
    "reference_annotation",
]


@dataclass(frozen=True)
class AssemblyRecord:
    animal_id: str
    haplotype: str
    sample: str
    sex: str = ""
    species: str = ""
    subspecies: str = ""
    breed: str = ""
    source: str = ""
    format_mode: str = "ragtag_cut"
    chrom_map: str = ""
    checksum: str = ""
    qc_fail: bool = False
    exclude_reason: str = ""
    graphs: tuple[str, ...] = ()
    tags: tuple[str, ...] = ()
    reference_annotation: bool = False
    extra: dict[str, str] = field(default_factory=dict)

    def to_metadata_row(self) -> dict[str, str]:
        """Return a dict compatible with legacy Snakefile Polars column names."""
        return {
            "Animal ID": self.animal_id,
            "Haplotype": self.haplotype,
            "Filename": self.sample,
            "Sex": self.sex,
            "Species": self.species,
            "Subspecies": self.subspecies,
            "Breed": self.breed,
            "Source": self.source,
            "Format mode": self.format_mode,
            "Chrom map": self.chrom_map,
            "Checksum": self.checksum,
            "Reference annotation": "Y" if (self.reference_annotation or self.format_mode == "ragtag_preserve") else "",
            "QC Fail": "Y" if self.qc_fail else "",
            "Exclude reason": self.exclude_reason,
            "Graphs": ";".join(self.graphs),
            "Tags": ";".join(self.tags),
        }


def _norm_key(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", name.strip().lower()).strip("_")


def _split_list(val: str) -> tuple[str, ...]:
    if not val:
        return ()
    parts = [p.strip() for p in re.split(r"[;,]", val) if p.strip()]
    seen: list[str] = []
    for p in parts:
        if p not in seen:
            seen.append(p)
    return tuple(seen)


def _normalize_sex(raw: str, row_num: int) -> str:
    key = raw.strip().lower()
    if key not in _SEX_MAP:
        raise ValueError(
            f"Row {row_num}: invalid Sex value {raw!r}; expected M, F, U, or empty"
        )
    return _SEX_MAP[key]


def _get_field(row_norm: dict[str, str], *aliases: str) -> str:
    for alias in aliases:
        if alias in row_norm and row_norm[alias] != "":
            return row_norm[alias]
    return ""


def detect_delimiter(path: Path) -> str:
    with path.open("r", encoding="utf-8-sig", newline="") as fin:
        first_line = next((line for line in fin if line.strip()), "")
    if "\t" in first_line:
        return "\t"
    return ","


def load_raw_rows(path: Path) -> tuple[list[str], list[dict[str, str]], str]:
    """Read CSV/TSV stripping UTF-8 BOM, CRLF, and surrounding cell whitespace."""
    if not path.exists():
        raise FileNotFoundError(f"Manifest file not found: {path}")
    delim = detect_delimiter(path)
    with path.open("r", encoding="utf-8-sig", newline="") as fin:
        reader = csv.reader(fin, delimiter=delim)
        raw_headers = next(reader, None)
        if not raw_headers:
            raise ValueError(f"Manifest {path} is empty")
        headers = [h.strip() for h in raw_headers]
        rows: list[dict[str, str]] = []
        for raw_row in reader:
            if not any(cell.strip() for cell in raw_row):
                continue
            padded = list(raw_row) + [""] * max(0, len(headers) - len(raw_row))
            row_dict = {headers[i]: padded[i].strip() for i in range(len(headers))}
            rows.append(row_dict)
    return headers, rows, delim


def validate_rows(
    headers: Sequence[str],
    rows: Sequence[Mapping[str, str]],
    source_name: str = "<manifest>",
    default_format_mode: str = "intact",
) -> list[AssemblyRecord]:
    if default_format_mode not in VALID_FORMAT_MODES:
        raise ValueError(
            f"Invalid default_format_mode {default_format_mode!r}; expected one of {VALID_FORMAT_MODES}"
        )
    norm_headers = {_norm_key(h) for h in headers if h.strip()}
    if "animal_id" not in norm_headers or "haplotype" not in norm_headers:
        raise ValueError(
            f"{source_name}: manifest must contain 'Animal ID' (or 'animal_id') "
            f"and 'Haplotype' (or 'haplotype') columns; found {list(headers)}"
        )
    if not rows:
        raise ValueError(f"{source_name}: manifest contains no assembly rows")

    records: list[AssemblyRecord] = []
    seen_samples: dict[str, int] = {}

    known_norm_keys = {
        "animal_id",
        "haplotype",
        "sex",
        "species",
        "subspecies",
        "breed",
        "source",
        "raw_path",
        "public_download",
        "format_mode",
        "chrom_map",
        "checksum",
        "qc_fail",
        "exclude_reason",
        "graphs",
        "tags",
        "reference_annotation",
    }

    for idx, row in enumerate(rows, start=2):
        row_norm = {_norm_key(k): (v.strip() if v is not None else "") for k, v in row.items() if k}

        animal_id = _get_field(row_norm, "animal_id")
        haplotype = _get_field(row_norm, "haplotype")
        if not animal_id:
            raise ValueError(f"{source_name} row {idx}: missing Animal ID")
        if not haplotype:
            raise ValueError(f"{source_name} row {idx}: missing Haplotype")

        if not SAFE_ID_RE.match(animal_id):
            raise ValueError(
                f"{source_name} row {idx}: invalid Animal ID {animal_id!r}; "
                "must match ^[A-Za-z0-9._-]+$ (no '#', ':', spaces, or slashes)"
            )
        if not SAFE_ID_RE.match(haplotype):
            raise ValueError(
                f"{source_name} row {idx}: invalid Haplotype {haplotype!r}; "
                "must match ^[A-Za-z0-9._-]+$"
            )

        sample = f"{animal_id}_{haplotype}"
        if sample in seen_samples:
            raise ValueError(
                f"{source_name} row {idx}: duplicate sample '{sample}' "
                f"(first seen on row {seen_samples[sample]})"
            )
        seen_samples[sample] = idx

        sex = _normalize_sex(_get_field(row_norm, "sex"), idx)
        species = _get_field(row_norm, "species")
        subspecies = _get_field(row_norm, "subspecies")
        breed = _get_field(row_norm, "breed")
        source = _get_field(row_norm, "source", "raw_path", "public_download")

        ref_ann_raw = _get_field(row_norm, "reference_annotation")
        reference_annotation = ref_ann_raw.lower() in _TRUTHY

        format_mode_raw = _get_field(row_norm, "format_mode")
        if not format_mode_raw:
            format_mode = "ragtag_preserve" if reference_annotation else default_format_mode
        else:
            format_mode = format_mode_raw.lower()
            if format_mode not in VALID_FORMAT_MODES:
                raise ValueError(
                    f"{source_name} row {idx} ({sample}): invalid format_mode {format_mode_raw!r}; "
                    f"expected one of {VALID_FORMAT_MODES}"
                )

        chrom_map = _get_field(row_norm, "chrom_map")
        checksum = _get_field(row_norm, "checksum")
        exclude_reason = _get_field(row_norm, "exclude_reason")
        qc_fail_raw = _get_field(row_norm, "qc_fail")
        qc_fail = (qc_fail_raw.lower() not in _FALSY_OR_PASS) or bool(exclude_reason)

        graphs = _split_list(_get_field(row_norm, "graphs"))
        tags = _split_list(_get_field(row_norm, "tags"))

        extra = {
            k: v
            for k, v in row.items()
            if k and _norm_key(k) not in known_norm_keys and v != ""
        }

        records.append(
            AssemblyRecord(
                animal_id=animal_id,
                haplotype=haplotype,
                sample=sample,
                sex=sex,
                species=species,
                subspecies=subspecies,
                breed=breed,
                source=source,
                format_mode=format_mode,
                chrom_map=chrom_map,
                checksum=checksum,
                qc_fail=qc_fail,
                exclude_reason=exclude_reason,
                graphs=graphs,
                tags=tags,
                reference_annotation=reference_annotation,
                extra=extra,
            )
        )

    return records


def load_manifest(path: str | Path, config: Mapping[str, Any] | None = None) -> list[AssemblyRecord]:
    """Load and validate an assembly manifest CSV/TSV."""
    p = Path(path)
    headers, rows, _ = load_raw_rows(p)
    default_mode = str((config or {}).get("default_format_mode", "intact")).strip().lower()
    return validate_rows(headers, rows, source_name=str(p), default_format_mode=default_mode)


def cactus_seqfile_name(
    record: AssemblyRecord,
    reference_ids: Iterable[str] = (),
) -> str:
    """Return Cactus seqfile genome identifier, avoiding '.0' suffix on reference haplotypes."""
    ref_set = {str(r).strip() for r in reference_ids if str(r).strip()}
    if record.haplotype == "0" and record.animal_id in ref_set:
        return record.animal_id
    return f"{record.animal_id}.{record.haplotype}"


def resolve_graph_choices(
    records: Sequence[AssemblyRecord],
    config: Mapping[str, Any] | None = None,
) -> list[str]:
    """Return deduplicated list of graph names from config['graphs'] and manifest rows."""
    choices: list[str] = ["All"]
    cfg_graphs = (config or {}).get("graphs") or {}
    if isinstance(cfg_graphs, Mapping):
        for g in cfg_graphs:
            g_str = str(g).strip()
            if g_str and g_str not in choices:
                choices.append(g_str)
    for rec in records:
        for g in rec.graphs:
            if g not in choices:
                choices.append(g)
    return choices


def _haplotype_rank(hap: str) -> tuple[int, str]:
    order = {"0": 0, "1": 1, "2": 2}
    return (order.get(hap, 3), hap)


def _safe_float(val: Any) -> float:
    if val is None:
        return math.nan
    s = str(val).strip()
    if not s or s.upper() in ("NA", "NAN", "NONE", "."):
        return math.nan
    try:
        return float(s)
    except ValueError:
        return math.nan


def resolve_graph_samples(
    records: Sequence[AssemblyRecord],
    graph: str,
    config: Mapping[str, Any] | None = None,
    qc_rows: Mapping[str, Mapping[str, Any]] | None = None,
) -> list[str]:
    """Resolve sample IDs belonging to `graph` using config['graphs'] or legacy tags."""
    cfg_graphs = (config or {}).get("graphs") or {}
    if isinstance(cfg_graphs, Mapping) and graph in cfg_graphs:
        spec = cfg_graphs[graph] or {}
        if not isinstance(spec, Mapping):
            raise ValueError(f"config['graphs'][{graph!r}] must be a mapping, got {type(spec).__name__}")

        exclude_qc_fail = bool(spec.get("exclude_qc_fail", graph != "All"))
        include_set = {str(s).strip() for s in (spec.get("include") or []) if str(s).strip()}
        exclude_set = {str(s).strip() for s in (spec.get("exclude") or []) if str(s).strip()}

        def _lower_set(key: str) -> set[str]:
            vals = spec.get(key) or []
            if isinstance(vals, str):
                vals = [vals]
            return {str(v).strip().lower() for v in vals if str(v).strip()}

        species_filter = _lower_set("species")
        subspecies_filter = _lower_set("subspecies")
        breed_filter = _lower_set("breed")
        sex_filter = _lower_set("sex")
        mode_filter = _lower_set("format_modes")
        tags_filter = _lower_set("tags")

        min_ng50 = _safe_float(spec.get("min_ng50"))
        min_auto_cov = _safe_float(spec.get("min_autosome_cov"))
        min_gsize = _safe_float(spec.get("min_genome_size"))
        max_gsize = _safe_float(spec.get("max_genome_size"))
        one_per_animal = bool(spec.get("one_per_animal", False))

        filtered: list[AssemblyRecord] = []
        for rec in records:
            if rec.sample in exclude_set:
                continue
            force_inc = rec.sample in include_set
            if not force_inc:
                if exclude_qc_fail and rec.qc_fail:
                    continue
                if species_filter and rec.species.lower() not in species_filter:
                    continue
                if subspecies_filter and rec.subspecies.lower() not in subspecies_filter:
                    continue
                if breed_filter and rec.breed.lower() not in breed_filter:
                    continue
                if sex_filter and rec.sex.lower() not in sex_filter:
                    continue
                if mode_filter and rec.format_mode.lower() not in mode_filter:
                    continue
                if tags_filter:
                    rec_tags = {t.lower() for t in (rec.tags + rec.graphs)}
                    if not (rec_tags & tags_filter):
                        continue
                if qc_rows and rec.sample in qc_rows:
                    q = qc_rows[rec.sample]
                    ng50 = _safe_float(q.get("NG50"))
                    acov = _safe_float(q.get("autosomes covered"))
                    gsize = _safe_float(q.get("genome size"))
                    if not math.isnan(min_ng50) and not math.isnan(ng50) and ng50 < min_ng50:
                        continue
                    if not math.isnan(min_auto_cov) and not math.isnan(acov) and acov < min_auto_cov:
                        continue
                    if not math.isnan(min_gsize) and not math.isnan(gsize) and gsize < min_gsize:
                        continue
                    if not math.isnan(max_gsize) and not math.isnan(gsize) and gsize > max_gsize:
                        continue
            filtered.append(rec)

        if one_per_animal:
            by_animal: dict[str, list[AssemblyRecord]] = {}
            for rec in filtered:
                by_animal.setdefault(rec.animal_id, []).append(rec)
            selected: list[AssemblyRecord] = []
            for recs in by_animal.values():
                if len(recs) == 1:
                    selected.append(recs[0])
                    continue
                if qc_rows:
                    def _score(r: AssemblyRecord) -> tuple[int, float, tuple[int, str]]:
                        q = qc_rows.get(r.sample, {})
                        ng = _safe_float(q.get("NG50"))
                        has_ng = 1 if not math.isnan(ng) else 0
                        ng_val = ng if not math.isnan(ng) else -1.0
                        # Negate haplotype rank comparison by sorting descending on score
                        hr = _haplotype_rank(r.haplotype)
                        return (has_ng, ng_val, (-hr[0], r.haplotype))

                    best = max(recs, key=_score)
                else:
                    best = min(recs, key=lambda r: _haplotype_rank(r.haplotype))
                selected.append(best)
            return [r.sample for r in selected]

        return [r.sample for r in filtered]

    if graph == "All":
        return [r.sample for r in records]

    return [r.sample for r in records if graph in r.graphs and not r.qc_fail]


def _parse_busco_complete(raw: str) -> float:
    """Extract total complete BUSCO % (Single + Duplicated) from summary field."""
    s = (raw or "").strip()
    if not s:
        return math.nan
    m_c = re.search(r"C:(\d+(?:\.\d+)?)%", s)
    if m_c:
        return float(m_c.group(1))
    m_s = re.search(r"S:(\d+(?:\.\d+)?)%", s)
    m_d = re.search(r"D:(\d+(?:\.\d+)?)%", s)
    if m_s and m_d:
        return float(m_s.group(1)) + float(m_d.group(1))
    parts = [p.strip().rstrip("%") for p in s.split(";") if p.strip()]
    nums: list[float] = []
    for p in parts[:2]:
        val = _safe_float(p.split(":")[-1])
        if math.isnan(val):
            return math.nan
        nums.append(val)
    if len(nums) == 2:
        return sum(nums)
    if len(nums) == 1:
        return nums[0]
    return math.nan


def infer_sex_and_qc_flags(
    record: AssemblyRecord,
    qc_row: Mapping[str, Any] | None,
    thresholds: Mapping[str, float] | None = None,
) -> tuple[str, list[str]]:
    """Infer biological sex from X/Y coverage and compute automated QC flags."""
    if not qc_row:
        return ("", [])
    th = dict(DEFAULT_QC_THRESHOLDS)
    if thresholds:
        for k, v in thresholds.items():
            th[k] = float(v)

    flags: list[str] = []
    gsize = _safe_float(qc_row.get("genome size"))
    ng50 = _safe_float(qc_row.get("NG50"))
    acov = _safe_float(qc_row.get("autosomes covered"))
    xcov = _safe_float(qc_row.get("X covered"))
    ycov = _safe_float(qc_row.get("Y covered"))
    busco_comp = _parse_busco_complete(str(qc_row.get("BUSCOs", "")))

    if not math.isnan(gsize):
        if gsize < th["min_genome_bp"]:
            flags.append("genome_size_low")
        elif gsize > th["max_genome_bp"]:
            flags.append("genome_size_high")

    if not math.isnan(ng50) and ng50 < th["min_ng50"]:
        flags.append("low_ng50")

    if not math.isnan(acov) and acov < th["min_autosome_cov"]:
        flags.append("low_autosome_cov")

    if not math.isnan(busco_comp) and busco_comp < th["min_busco_complete"]:
        flags.append("low_busco")

    inferred_sex = "U"
    if not math.isnan(ycov) and ycov >= th["y_cov_male_min"]:
        inferred_sex = "M"
    elif not math.isnan(xcov) and xcov >= th["x_cov_female_min"] and (
        math.isnan(ycov) or ycov < th["y_cov_male_min"]
    ):
        inferred_sex = "F"

    if record.sex in ("M", "F") and inferred_sex in ("M", "F") and record.sex != inferred_sex:
        flags.append("sex_mismatch")

    return inferred_sex, flags


def load_qc_summary(qc_path: Path) -> dict[str, dict[str, str]]:
    """Load analyses/QC_summary.{graph}.csv into sample -> row dict."""
    if not qc_path.exists():
        raise FileNotFoundError(f"QC summary file not found: {qc_path}")
    delim = detect_delimiter(qc_path)
    qc_map: dict[str, dict[str, str]] = {}
    with qc_path.open("r", encoding="utf-8-sig", newline="") as fin:
        reader = csv.DictReader(fin, delimiter=delim)
        for row in reader:
            cleaned = {k.strip(): (v.strip() if v is not None else "") for k, v in row.items() if k}
            sid = cleaned.get("sample") or cleaned.get("Filename") or ""
            if not sid and cleaned.get("Animal ID") and cleaned.get("Haplotype"):
                sid = f"{cleaned['Animal ID']}_{cleaned['Haplotype']}"
            if sid:
                qc_map[sid] = cleaned
    return qc_map


def merge_manifest_and_qc(
    manifest_path: str | Path,
    qc_csv_path: str | Path,
    output_tsv_path: str | Path,
    thresholds: Mapping[str, float] | None = None,
) -> Path:
    """Join input manifest with QC_summary CSV and write a clean UTF-8 TSV report."""
    records = load_manifest(manifest_path)
    qc_map = load_qc_summary(Path(qc_csv_path))

    qc_cols = [
        "genome size",
        "N contigs",
        "NG50",
        "BUSCOs",
        "SNPs",
        "InDels",
        "autosomes covered",
        "X covered",
        "Y covered",
        "MT covered",
    ]
    out_cols = [
        "sample",
        "animal_id",
        "haplotype",
        "sex",
        "inferred_sex",
        "species",
        "subspecies",
        "breed",
        "format_mode",
        "qc_fail",
        "exclude_reason",
        "auto_qc_flags",
        *qc_cols,
    ]

    out_path = Path(output_tsv_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as fout:
        writer = csv.DictWriter(fout, fieldnames=out_cols, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for rec in records:
            q = qc_map.get(rec.sample)
            inferred_sex, flags = infer_sex_and_qc_flags(rec, q, thresholds=thresholds)
            row_out: dict[str, str] = {
                "sample": rec.sample,
                "animal_id": rec.animal_id,
                "haplotype": rec.haplotype,
                "sex": rec.sex,
                "inferred_sex": inferred_sex,
                "species": rec.species,
                "subspecies": rec.subspecies,
                "breed": rec.breed,
                "format_mode": rec.format_mode,
                "qc_fail": "Y" if rec.qc_fail else "N",
                "exclude_reason": rec.exclude_reason,
                "auto_qc_flags": ";".join(flags) if flags else "PASS" if q else "",
            }
            for col in qc_cols:
                row_out[col] = q.get(col, "") if q else ""
            writer.writerow(row_out)

    return out_path


def generate_unique_hex_id(existing_ids: Iterable[str], length: int = 8) -> str:
    """Generate an uppercase hex ID guaranteed not to collide or coerce to a number.

    Requiring at least one letter from {A, B, C, D, F} avoids Excel/R/Polars
    interpreting all-digit strings as integers or `1234E567` as scientific notation.
    """
    if length < 2:
        raise ValueError("Hex ID length must be >= 2")
    used = {str(x).strip().upper() for x in existing_ids if str(x).strip()}
    n_bytes = (length + 1) // 2
    for _ in range(10_000):
        candidate = secrets.token_hex(n_bytes).upper()[:length]
        if any(ch in HEX_SAFE_LETTERS for ch in candidate) and candidate not in used:
            return candidate
    raise RuntimeError("Failed to generate a unique hex ID after 10,000 attempts")


def load_id_map(path: str | Path) -> tuple[dict[str, str], set[str]]:
    """Load private external_id -> animal_id crosswalk TSV/CSV."""
    p = Path(path)
    if not p.exists():
        return {}, set()
    delim = detect_delimiter(p)
    ext_to_animal: dict[str, str] = {}
    used_animals: set[str] = set()
    with p.open("r", encoding="utf-8-sig", newline="") as fin:
        reader = csv.DictReader(fin, delimiter=delim)
        if not reader.fieldnames:
            return {}, set()
        for row in reader:
            norm = {_norm_key(k): (v.strip() if v else "") for k, v in row.items() if k}
            aid = norm.get("animal_id", "")
            ext = norm.get("external_id", "")
            if aid:
                used_animals.add(aid.upper())
                if ext and ext not in ext_to_animal:
                    ext_to_animal[ext] = aid
    return ext_to_animal, used_animals


def append_id_map(path: str | Path, entries: Sequence[tuple[str, str]]) -> Path:
    """Append (animal_id, external_id) mappings with UTC timestamp to crosswalk TSV."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    needs_header = (not p.exists()) or p.stat().st_size == 0
    ts = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    with p.open("a", encoding="utf-8", newline="") as fout:
        writer = csv.writer(fout, delimiter="\t", lineterminator="\n")
        if needs_header:
            writer.writerow(["animal_id", "external_id", "added_utc"])
        for animal_id, external_id in entries:
            writer.writerow([animal_id, external_id, ts])
    return p


def _load_config_file(path: Path) -> dict[str, Any]:
    text = path.read_text(encoding="utf-8")
    if path.suffix.lower() == ".json":
        data = json.loads(text)
        return dict(data) if isinstance(data, dict) else {}
    try:
        import yaml  # type: ignore

        data = yaml.safe_load(text)
        return dict(data) if isinstance(data, dict) else {}
    except ImportError:
        return {}


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Validate and inspect BPC assembly manifests.")
    sub = parser.add_subparsers(dest="command", required=True)

    p_val = sub.add_parser("validate", help="Validate manifest schema, IDs, and format modes.")
    p_val.add_argument("--manifest", required=True, type=Path, help="Input CSV/TSV manifest.")
    p_val.add_argument("--config", type=Path, help="Optional config YAML/JSON.")

    p_graphs = sub.add_parser("list-graphs", help="List available graph subsets.")
    p_graphs.add_argument("--manifest", required=True, type=Path, help="Input CSV/TSV manifest.")
    p_graphs.add_argument("--config", type=Path, help="Optional config YAML/JSON.")

    p_merge = sub.add_parser("merge-qc", help="Join manifest with QC_summary CSV into TSV.")
    p_merge.add_argument("--manifest", required=True, type=Path, help="Input CSV/TSV manifest.")
    p_merge.add_argument("--qc", required=True, type=Path, help="QC_summary CSV file.")
    p_merge.add_argument("--output", required=True, type=Path, help="Output merged TSV file.")
    p_merge.add_argument("--thresholds-json", type=str, help="JSON dict overriding QC thresholds.")

    p_new = sub.add_parser("new-id", help="Generate unique Excel-safe hex animal ID(s).")
    p_new.add_argument("--manifest", type=Path, help="Existing manifest CSV/TSV to avoid collisions.")
    p_new.add_argument("--id-map", type=Path, help="Optional private ID crosswalk TSV.")
    p_new.add_argument("--count", type=int, default=1, help="Number of IDs to generate.")
    p_new.add_argument("--length", type=int, default=8, help="Hex ID length (default: 8).")

    args = parser.parse_args(argv)

    if args.command == "validate":
        cfg = _load_config_file(args.config) if args.config else None
        records = load_manifest(args.manifest, config=cfg)
        graphs = resolve_graph_choices(records, config=cfg)
        print(
            f"Valid manifest: {len(records)} assemblies, "
            f"graphs={','.join(graphs)}"
        )
        return 0

    if args.command == "list-graphs":
        cfg = _load_config_file(args.config) if args.config else None
        records = load_manifest(args.manifest, config=cfg)
        for g in resolve_graph_choices(records, config=cfg):
            samples = resolve_graph_samples(records, g, config=cfg)
            print(f"{g}\t{len(samples)}\t{','.join(samples)}")
        return 0

    if args.command == "merge-qc":
        thresholds = json.loads(args.thresholds_json) if args.thresholds_json else None
        out = merge_manifest_and_qc(
            args.manifest, args.qc, args.output, thresholds=thresholds
        )
        print(f"Wrote {out}")
        return 0

    if args.command == "new-id":
        used: set[str] = set()
        if args.manifest and args.manifest.exists():
            used.update(r.animal_id for r in load_manifest(args.manifest))
        if args.id_map:
            _, map_used = load_id_map(args.id_map)
            used.update(map_used)
        for _ in range(max(1, args.count)):
            new_id = generate_unique_hex_id(used, length=args.length)
            used.add(new_id)
            print(new_id)
        return 0

    return 1


if __name__ == "__main__":
    sys.exit(main())
