#!/usr/bin/env python3
"""Single-command CLI helper to register new assemblies into the BPC manifest.

Supports single-assembly, paired diploid (--hap1 / --hap2), and batch TSV/CSV
intake, optional automatic Excel-safe hex animal ID generation, private ID
crosswalk mapping (--id-map), pre-flight validation, and SHA256 calculation.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import re
import shutil
import subprocess
import sys
from pathlib import Path

SCRIPTS_DIR = Path(__file__).resolve().parent
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from manifest import (  # noqa: E402
    CANONICAL_TSV_COLUMNS,
    VALID_FORMAT_MODES,
    append_id_map,
    generate_unique_hex_id,
    load_id_map,
    load_raw_rows,
    validate_rows,
)

NCBI_ACC_RE = re.compile(r"^GC[AF]_\d{9}(?:\.\d+)?$")
URL_RE = re.compile(r"^(?:https?|ftp)://", re.IGNORECASE)
AGC_RE = re.compile(r"^agc://", re.IGNORECASE)


def _is_remote_or_archive(source: str) -> bool:
    return bool(NCBI_ACC_RE.match(source) or URL_RE.match(source) or AGC_RE.match(source))


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fin:
        while chunk := fin.read(1024 * 1024):
            h.update(chunk)
    return f"sha256:{h.hexdigest()}"


def _peek_fasta_header(path: Path) -> str:
    with path.open("rb") as probe:
        magic = probe.read(4)
    if magic.startswith(b"\x1f\x8b"):
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fin:
            for line in fin:
                s = line.strip()
                if s:
                    return s
    elif magic == b"\x28\xb5\x2f\xfd":
        if not shutil.which("zstd"):
            return ">"  # cannot inspect without zstd, defer to stage_assembly
        proc = subprocess.Popen(["zstd", "-dc", str(path)], stdout=subprocess.PIPE)
        assert proc.stdout is not None
        try:
            for raw in proc.stdout:
                s = raw.decode("utf-8", errors="replace").strip()
                if s:
                    return s
        finally:
            proc.terminate()
            proc.wait()
    else:
        with path.open("rt", encoding="utf-8", errors="replace") as fin:
            for line in fin:
                s = line.strip()
                if s:
                    return s
    return ""


def _inspect_local_source(source: str, checksum: str, compute_checksum: bool) -> str:
    if not source or _is_remote_or_archive(source):
        return checksum
    p = Path(source)
    if not p.exists():
        raise FileNotFoundError(f"Local source assembly does not exist: {p}")
    first_line = _peek_fasta_header(p)
    if not first_line.startswith(">"):
        raise ValueError(
            f"Source {p} does not look like a FASTA file (first non-empty line: {first_line!r})"
        )
    if "_RagTag" in first_line:
        print(
            f"Warning: {p} contains '_RagTag' in header {first_line!r}; "
            "stage_assembly will strip '_RagTag' suffixes.",
            file=sys.stderr,
        )
    if not checksum and compute_checksum:
        return _sha256_file(p)
    return checksum


def _map_to_existing_headers(
    headers: list[str],
    canonical_row: dict[str, str],
) -> dict[str, str]:
    """Map canonical keys onto an existing manifest's header names."""
    alias_lookup = {
        "animal_id": ("animal_id", "Animal ID"),
        "haplotype": ("haplotype", "Haplotype"),
        "sex": ("sex", "Sex"),
        "species": ("species", "Species"),
        "subspecies": ("subspecies", "Subspecies"),
        "breed": ("breed", "Breed"),
        "source": ("source", "Source", "raw_path", "Raw path", "Public download"),
        "format_mode": ("format_mode", "Format mode"),
        "chrom_map": ("chrom_map", "Chrom map"),
        "checksum": ("checksum", "Checksum"),
        "qc_fail": ("qc_fail", "QC Fail"),
        "exclude_reason": ("exclude_reason", "Exclude reason"),
        "graphs": ("graphs", "Graphs"),
        "tags": ("tags", "Tags"),
        "reference_annotation": ("reference_annotation", "Reference annotation"),
    }
    norm_to_header = {
        re.sub(r"[^a-z0-9]+", "_", h.strip().lower()).strip("_"): h for h in headers
    }
    out = {h: "" for h in headers}
    for key, val in canonical_row.items():
        if key in norm_to_header:
            out[norm_to_header[key]] = val
            continue
        for alias in alias_lookup.get(key, ()):
            norm_a = re.sub(r"[^a-z0-9]+", "_", alias.strip().lower()).strip("_")
            if norm_a in norm_to_header:
                out[norm_to_header[norm_a]] = val
                break
    return out


def add_assemblies(
    manifest_path: Path,
    incoming_rows: list[dict[str, str]],
    compute_checksum: bool = True,
    dry_run: bool = False,
    id_map_path: Path | None = None,
) -> tuple[list[str], Path]:
    if manifest_path.exists():
        headers, existing_rows, delim = load_raw_rows(manifest_path)
        existing_records = (
            validate_rows(headers, existing_rows, source_name=str(manifest_path))
            if existing_rows
            else []
        )
    else:
        headers = list(CANONICAL_TSV_COLUMNS)
        existing_rows = []
        existing_records = []
        delim = "," if manifest_path.suffix.lower() == ".csv" else "\t"

    used_ids: set[str] = {r.animal_id.upper() for r in existing_records}
    ext_map: dict[str, str] = {}
    if id_map_path is not None:
        ext_map, map_used = load_id_map(id_map_path)
        used_ids.update(map_used)

    group_assigned: dict[str, str] = {}
    new_id_map_entries: list[tuple[str, str]] = []

    enriched_incoming: list[dict[str, str]] = []
    for raw_in in incoming_rows:
        row = dict(raw_in)
        ext_id = (row.pop("external_id", "") or "").strip()
        group_key = (row.pop("_group_key", "") or "").strip()
        aid = (row.get("animal_id") or "").strip()

        if not aid or aid.lower() == "auto":
            if ext_id:
                if ext_id in ext_map:
                    aid = ext_map[ext_id]
                elif ext_id in group_assigned:
                    aid = group_assigned[ext_id]
                else:
                    aid = generate_unique_hex_id(used_ids)
                    used_ids.add(aid)
                    group_assigned[ext_id] = aid
                    new_id_map_entries.append((aid, ext_id))
            elif group_key:
                if group_key in group_assigned:
                    aid = group_assigned[group_key]
                else:
                    aid = generate_unique_hex_id(used_ids)
                    used_ids.add(aid)
                    group_assigned[group_key] = aid
            else:
                aid = generate_unique_hex_id(used_ids)
                used_ids.add(aid)
            row["animal_id"] = aid
        else:
            used_ids.add(aid.upper())
            if ext_id and ext_id not in ext_map and ext_id not in group_assigned:
                group_assigned[ext_id] = aid
                new_id_map_entries.append((aid, ext_id))

        row["checksum"] = _inspect_local_source(
            row.get("source", ""), row.get("checksum", ""), compute_checksum=compute_checksum
        )
        enriched_incoming.append(row)

    for key in ("source", "format_mode", "chrom_map", "checksum", "tags"):
        if any(r.get(key) for r in enriched_incoming):
            norm_existing = {
                re.sub(r"[^a-z0-9]+", "_", h.strip().lower()).strip("_") for h in headers
            }
            if key not in norm_existing and not (
                key == "source" and ("raw_path" in norm_existing or "public_download" in norm_existing)
            ):
                headers.append(key)

    mapped_new = [_map_to_existing_headers(headers, r) for r in enriched_incoming]
    combined = list(existing_rows) + mapped_new
    records = validate_rows(headers, combined, source_name=str(manifest_path))
    added_samples = [r.sample for r in records[len(existing_rows) :]]

    if not dry_run:
        manifest_path.parent.mkdir(parents=True, exist_ok=True)
        with manifest_path.open("w", encoding="utf-8", newline="") as fout:
            writer = csv.DictWriter(
                fout, fieldnames=headers, delimiter=delim, lineterminator="\n"
            )
            writer.writeheader()
            for r in combined:
                writer.writerow({h: r.get(h, "") for h in headers})
        if id_map_path is not None and new_id_map_entries:
            append_id_map(id_map_path, new_id_map_entries)

    return added_samples, manifest_path


def _load_batch_rows(batch_path: Path) -> list[dict[str, str]]:
    _, raw_rows, _ = load_raw_rows(batch_path)
    out: list[dict[str, str]] = []
    for row in raw_rows:
        norm = {
            re.sub(r"[^a-z0-9]+", "_", k.strip().lower()).strip("_"): v.strip()
            for k, v in row.items()
            if k
        }
        out.append(
            {
                "animal_id": norm.get("animal_id", ""),
                "external_id": norm.get("external_id") or norm.get("original_id") or "",
                "haplotype": norm.get("haplotype", ""),
                "sex": norm.get("sex", ""),
                "species": norm.get("species", ""),
                "subspecies": norm.get("subspecies", ""),
                "breed": norm.get("breed", ""),
                "source": norm.get("source") or norm.get("raw_path") or norm.get("public_download") or "",
                "format_mode": norm.get("format_mode", ""),
                "chrom_map": norm.get("chrom_map", ""),
                "checksum": norm.get("checksum", ""),
                "qc_fail": norm.get("qc_fail", ""),
                "exclude_reason": norm.get("exclude_reason", ""),
                "graphs": norm.get("graphs", ""),
                "tags": norm.get("tags", ""),
                "reference_annotation": norm.get("reference_annotation", ""),
            }
        )
    return out


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Add one or more assemblies to the BPC manifest with pre-flight checks."
    )
    parser.add_argument(
        "--manifest",
        required=True,
        type=Path,
        help="Target manifest TSV/CSV (created if missing).",
    )
    parser.add_argument("--batch", type=Path, help="Batch TSV/CSV file of new assemblies.")
    parser.add_argument("--source", default="", help="Local path, URL, GCA_/GCF_, or agc://...")
    parser.add_argument("--hap1", default="", help="Haplotype 1 source for paired diploid intake.")
    parser.add_argument("--hap2", default="", help="Haplotype 2 source for paired diploid intake.")
    parser.add_argument(
        "--animal-id",
        default="auto",
        help="Animal ID (e.g. BOV0123), or 'auto' (default) to mint a unique hex ID.",
    )
    parser.add_argument("--haplotype", default="", help="Haplotype (e.g. 0, 1, 2).")
    parser.add_argument(
        "--external-id",
        default="",
        help="Optional private external/sample ID stored only in --id-map, never in manifest.",
    )
    parser.add_argument(
        "--id-map",
        type=Path,
        default=None,
        help="Optional private TSV mapping external_id -> generated hex animal_id.",
    )
    parser.add_argument("--sex", default="", help="Sex (M, F, U).")
    parser.add_argument("--species", default="", help="Species (e.g. 'Bos taurus').")
    parser.add_argument("--subspecies", default="", help="Subspecies.")
    parser.add_argument("--breed", default="", help="Breed.")
    parser.add_argument(
        "--format-mode",
        default="",
        choices=("", *VALID_FORMAT_MODES),
        help="Formatting mode (ragtag_cut, ragtag_preserve, already_placed).",
    )
    parser.add_argument("--chrom-map", default="", help="Optional chromosome mapping file.")
    parser.add_argument("--checksum", default="", help="Expected sha256:<hex> or md5:<hex>.")
    parser.add_argument("--checksum1", default="", help="Expected checksum for --hap1.")
    parser.add_argument("--checksum2", default="", help="Expected checksum for --hap2.")
    parser.add_argument("--tags", default="", help="Optional semicolon-separated tags.")
    parser.add_argument("--graphs", default="", help="Optional semicolon-separated graph names.")
    parser.add_argument(
        "--reference-annotation",
        default="",
        help="Set 'Y' if assembly carries reference annotation.",
    )
    parser.add_argument(
        "--no-checksum",
        action="store_true",
        help="Skip computing sha256 checksum for local files.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Validate without modifying the manifest on disk.",
    )

    args = parser.parse_args(argv)

    if args.batch:
        if args.hap1 or args.hap2:
            parser.error("Cannot combine --batch with --hap1/--hap2.")
        incoming = _load_batch_rows(args.batch)
    elif args.hap1 or args.hap2:
        if args.source or args.haplotype:
            parser.error("Do not combine --hap1/--hap2 with --source or --haplotype.")
        shared_group = args.external_id or "__diploid_pair__"
        incoming = []
        for hap_num, hap_src, hap_ck in (
            ("1", args.hap1, args.checksum1),
            ("2", args.hap2, args.checksum2),
        ):
            if not hap_src:
                continue
            incoming.append(
                {
                    "animal_id": args.animal_id,
                    "external_id": args.external_id,
                    "_group_key": shared_group,
                    "haplotype": hap_num,
                    "sex": args.sex,
                    "species": args.species,
                    "subspecies": args.subspecies,
                    "breed": args.breed,
                    "source": hap_src,
                    "format_mode": args.format_mode,
                    "chrom_map": args.chrom_map,
                    "checksum": hap_ck,
                    "qc_fail": "",
                    "exclude_reason": "",
                    "graphs": args.graphs,
                    "tags": args.tags,
                    "reference_annotation": args.reference_annotation,
                }
            )
    else:
        if not args.haplotype:
            parser.error("Provide --haplotype (or --hap1/--hap2, or --batch).")
        incoming = [
            {
                "animal_id": args.animal_id,
                "external_id": args.external_id,
                "haplotype": args.haplotype,
                "sex": args.sex,
                "species": args.species,
                "subspecies": args.subspecies,
                "breed": args.breed,
                "source": args.source,
                "format_mode": args.format_mode,
                "chrom_map": args.chrom_map,
                "checksum": args.checksum,
                "qc_fail": "",
                "exclude_reason": "",
                "graphs": args.graphs,
                "tags": args.tags,
                "reference_annotation": args.reference_annotation,
            }
        ]

    added_samples, target = add_assemblies(
        args.manifest,
        incoming,
        compute_checksum=not args.no_checksum,
        dry_run=args.dry_run,
        id_map_path=args.id_map,
    )

    action = "Validated (dry-run)" if args.dry_run else "Added"
    print(f"{action} {len(added_samples)} assembly(ies) -> {target}: {', '.join(added_samples)}")
    qc_targets = " ".join(f"analyses/QC/summary/{s}.csv" for s in added_samples)
    print("Next step (ingest & QC new sample(s)):")
    print(f"  snakemake -s Snakefile --configfile config.yaml {qc_targets}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
