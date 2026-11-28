#!/usr/bin/env python3
"""Safe PanSN FASTA renaming and single-pass per-chromosome splitting.

Replaces fragile awk/grep pipelines:
- Full IUPAC reverse complement (preserving N and ambiguity bases on '-' strand)
- Low-memory natural chromosome ordering without buffering the whole genome in `sort`
- Deterministic unplaced contig naming and coordinate validation
- Single-pass PanSN splitting across all chromosomes (including X/Y/MT and absent chroms)
"""

from __future__ import annotations

import argparse
import gzip
import json
import os
import re
import shutil
import subprocess
import sys
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

DEFAULT_CHROMOSOMES: tuple[str, ...] = tuple(str(i) for i in range(1, 30)) + ("X", "Y", "MT")

_REVCOMP_TRANS = str.maketrans(
    "ACGTRYKMSWBDHVNacgtrykmswbdhvn",
    "TGCAYRMKSWVHDBNtgcayrmkswvhdbn",
)
_CHR_PREFIX_RE = re.compile(r"^(?:chr|bta)[_-]?", re.IGNORECASE)
_PANSN_UNDERSCORE_COORD_RE = re.compile(r"^(.+?)_(\d+-\d+)$")


@dataclass(frozen=True)
class OutputPlacement:
    header: str
    chrom: str
    is_placed: bool
    comp_id: str
    comp_start: int  # 1-based inclusive
    comp_end: int  # 1-based inclusive
    orient: str
    chrom_rank: tuple[int, int, int, int, str]


@dataclass(frozen=True)
class ContigLocation:
    length: int
    offset: int = -1
    byte_len: int = -1


def reverse_complement(seq: str) -> str:
    """Return full IUPAC-preserving reverse complement."""
    return seq.translate(_REVCOMP_TRANS)[::-1]


def format_pansn_header(
    animal_id: str,
    haplotype: str,
    chrom: str,
    start: int,
    end: int,
    schema: str = "cactus",
) -> str:
    prefix = f"{animal_id}#{haplotype}#{chrom}"
    if schema == "cactus":
        return f"{prefix}:{start}-{end}"
    if schema == "vg":
        return f"{prefix}[{start - 1}-{end}]"
    if schema == "_":
        return f"{prefix}_{start}-{end}"
    raise ValueError(f"Unsupported PanSN schema {schema!r}; expected cactus, vg, or _")


def canonicalize_chrom_name(raw_name: str, valid_chroms: set[str]) -> str | None:
    """Normalize common chromosome labels (chr1, chrM, MT, etc.) against valid_chroms."""
    name = raw_name.strip()
    if name.endswith("_RagTag"):
        name = name[: -len("_RagTag")]
    if name in valid_chroms:
        return name
    stripped = _CHR_PREFIX_RE.sub("", name)
    upper = stripped.upper()
    if upper in ("M", "MT") and "MT" in valid_chroms:
        return "MT"
    if upper in valid_chroms:
        return upper
    if stripped in valid_chroms:
        return stripped
    return None


def load_chrom_map(path: Path, valid_chroms: set[str]) -> dict[str, str]:
    """Load 2-column TSV/CSV, NCBI assembly_report.txt, or UCSC chromAlias.txt."""
    if not path.exists():
        raise FileNotFoundError(f"Chromosome map file not found: {path}")

    text_lines = path.read_text(encoding="utf-8-sig").splitlines()
    mapping: dict[str, str] = {}

    # Detect NCBI assembly_report.txt (usually has comment header with Sequence-Name / Assigned-Molecule)
    is_ncbi_report = any(
        "Sequence-Name" in ln and "Assigned-Molecule" in ln
        for ln in text_lines[:40]
        if ln.startswith("#")
    )

    for raw in text_lines:
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if "\t" in line:
            fields = [f.strip() for f in line.split("\t")]
        elif "," in line:
            fields = [f.strip() for f in line.split(",")]
        else:
            fields = line.split()

        if is_ncbi_report and len(fields) >= 7:
            seq_role = fields[1].lower()
            assigned = fields[2]
            if seq_role not in ("assembled-molecule", "assembled_molecule") and assigned.upper() not in ("MT", "M"):
                continue
            canon = canonicalize_chrom_name(assigned, valid_chroms)
            if not canon:
                continue
            for col_idx in (0, 4, 6, 9):
                if col_idx < len(fields) and fields[col_idx] and fields[col_idx].lower() != "na":
                    mapping[fields[col_idx]] = canon
            continue

        if len(fields) == 2:
            raw_id, target = fields[0], fields[1]
            canon = canonicalize_chrom_name(target, valid_chroms) or target
            if raw_id:
                mapping[raw_id] = canon
            continue

        # UCSC chromAlias.txt or multi-column alias table: find a field matching valid_chroms
        if len(fields) > 2:
            canon: str | None = None
            for f in fields:
                canon = canonicalize_chrom_name(f, valid_chroms)
                if canon:
                    break
            if canon:
                for f in fields:
                    if f:
                        mapping[f] = canon

    return mapping


def _is_gzip(path: Path) -> bool:
    with path.open("rb") as probe:
        return probe.read(2) == b"\x1f\x8b"


def _index_plain_fasta(path: Path) -> dict[str, ContigLocation]:
    """Scan plain FASTA once to record lengths and exact byte offsets."""
    locs: dict[str, ContigLocation] = {}
    current_id: str | None = None
    seq_start = 0
    seq_len = 0
    with path.open("rb") as fin:
        while True:
            pos = fin.tell()
            line = fin.readline()
            if not line:
                break
            if line.startswith(b">"):
                if current_id is not None:
                    locs[current_id] = ContigLocation(
                        length=seq_len, offset=seq_start, byte_len=pos - seq_start
                    )
                cid = line[1:].strip().split(None, 1)[0].decode("utf-8", errors="strict")
                if cid in locs:
                    raise ValueError(f"Duplicate contig '{cid}' in {path}")
                current_id = cid
                seq_start = fin.tell()
                seq_len = 0
            else:
                seq_len += len(line.strip())
        if current_id is not None:
            end_pos = fin.tell()
            locs[current_id] = ContigLocation(
                length=seq_len, offset=seq_start, byte_len=end_pos - seq_start
            )
    return locs


def _read_fai(fai_path: Path) -> dict[str, ContigLocation]:
    locs: dict[str, ContigLocation] = {}
    with fai_path.open("r", encoding="utf-8") as fin:
        for line in fin:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            locs[fields[0]] = ContigLocation(length=int(fields[1]))
    return locs


class SequenceFetcher:
    """Memory-efficient random-access sequence reader for plain or bgzip+fai FASTA."""

    def __init__(self, fasta_path: Path) -> None:
        self.fasta_path = fasta_path
        self.is_gz = _is_gzip(fasta_path)
        self._plain_fh = None
        self._mem_seqs: dict[str, str] | None = None

        if not self.is_gz:
            self.locations = _index_plain_fasta(fasta_path)
            self._plain_fh = fasta_path.open("rb")
        else:
            fai_path = Path(f"{fasta_path}.fai")
            if fai_path.exists() and shutil.which("samtools"):
                self.locations = _read_fai(fai_path)
            else:
                # Non-bgzip/unindexed gzip fallback
                seqs: dict[str, str] = {}
                locs: dict[str, ContigLocation] = {}
                cur_id: str | None = None
                chunks: list[str] = []
                with gzip.open(fasta_path, "rt", encoding="utf-8") as fin:
                    for line in fin:
                        s = line.strip()
                        if not s:
                            continue
                        if s.startswith(">"):
                            if cur_id is not None:
                                joined = "".join(chunks).upper()
                                seqs[cur_id] = joined
                                locs[cur_id] = ContigLocation(length=len(joined))
                            cur_id = s[1:].split(None, 1)[0]
                            chunks = []
                        else:
                            chunks.append(s)
                    if cur_id is not None:
                        joined = "".join(chunks).upper()
                        seqs[cur_id] = joined
                        locs[cur_id] = ContigLocation(length=len(joined))
                self._mem_seqs = seqs
                self.locations = locs

    def close(self) -> None:
        if self._plain_fh is not None:
            self._plain_fh.close()
            self._plain_fh = None

    def fetch(self, contig_id: str, start_1based: int | None = None, end_1based: int | None = None) -> str:
        loc = self.locations.get(contig_id)
        if loc is None:
            raise KeyError(f"Contig '{contig_id}' not found in {self.fasta_path}")

        if self._mem_seqs is not None:
            seq = self._mem_seqs[contig_id]
        elif self._plain_fh is not None:
            self._plain_fh.seek(loc.offset)
            raw = self._plain_fh.read(loc.byte_len)
            seq = "".join(raw.decode("ascii").split()).upper()
        else:
            region = (
                f"{contig_id}:{start_1based}-{end_1based}"
                if start_1based is not None and end_1based is not None
                else contig_id
            )
            res = subprocess.run(
                ["samtools", "faidx", str(self.fasta_path), region],
                capture_output=True,
                text=True,
                check=True,
            )
            lines = [ln.strip() for ln in res.stdout.splitlines() if ln and not ln.startswith(">")]
            return "".join(lines).upper()

        if start_1based is not None and end_1based is not None:
            return seq[start_1based - 1 : end_1based]
        return seq


def _build_placements_from_agp(
    agp_path: Path,
    contig_locs: dict[str, ContigLocation],
    animal_id: str,
    haplotype: str,
    chroms: Sequence[str],
    schema: str,
) -> list[OutputPlacement]:
    chrom_order = {c: i for i, c in enumerate(chroms)}
    valid_chroms = set(chroms)

    placements: list[OutputPlacement] = []
    covered_contigs: set[str] = set()
    unplaced_obj_ids: dict[str, int] = {}

    with agp_path.open("r", encoding="utf-8") as fin:
        for line_num, raw in enumerate(fin, start=1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 9:
                raise ValueError(f"{agp_path} line {line_num}: expected >= 9 AGP columns")
            comp_type = fields[4]
            if comp_type != "W":
                continue

            obj_id = fields[0]
            obj_start = int(fields[1])
            obj_end = int(fields[2])
            comp_id = fields[5]
            comp_start = int(fields[6])
            comp_end = int(fields[7])
            orient = fields[8] if fields[8] in ("+", "-") else "+"

            if comp_id not in contig_locs:
                raise ValueError(
                    f"{agp_path} line {line_num}: component '{comp_id}' missing from input FASTA"
                )
            covered_contigs.add(comp_id)

            stripped_obj = obj_id[: -len("_RagTag")] if obj_id.endswith("_RagTag") else obj_id
            if obj_id.endswith("_RagTag") and stripped_obj in valid_chroms:
                chrom = stripped_obj
                is_placed = True
                rank = (0, chrom_order[chrom], obj_start, obj_end, comp_id)
            elif stripped_obj in valid_chroms and obj_id in valid_chroms:
                chrom = stripped_obj
                is_placed = True
                rank = (0, chrom_order[chrom], obj_start, obj_end, comp_id)
            else:
                if obj_id not in unplaced_obj_ids:
                    unplaced_obj_ids[obj_id] = len(unplaced_obj_ids) + 1
                uid = unplaced_obj_ids[obj_id]
                chrom = f"unplaced_{uid}"
                is_placed = False
                rank = (1, uid, obj_start, obj_end, comp_id)

            header = format_pansn_header(
                animal_id, haplotype, chrom, obj_start, obj_end, schema=schema
            )
            placements.append(
                OutputPlacement(
                    header=header,
                    chrom=chrom,
                    is_placed=is_placed,
                    comp_id=comp_id,
                    comp_start=comp_start,
                    comp_end=comp_end,
                    orient=orient,
                    chrom_rank=rank,
                )
            )

    # Guardrail: emit any input contig absent from AGP (should not occur, prevents silent loss)
    next_unplaced = len(unplaced_obj_ids) + 1
    for cid, loc in contig_locs.items():
        if cid not in covered_contigs:
            chrom = f"unplaced_{next_unplaced}"
            header = format_pansn_header(
                animal_id, haplotype, chrom, 1, loc.length, schema=schema
            )
            placements.append(
                OutputPlacement(
                    header=header,
                    chrom=chrom,
                    is_placed=False,
                    comp_id=cid,
                    comp_start=1,
                    comp_end=loc.length,
                    orient="+",
                    chrom_rank=(1, next_unplaced, 1, loc.length, cid),
                )
            )
            next_unplaced += 1

    placements.sort(key=lambda p: p.chrom_rank)
    return placements


def _build_placements_already_placed(
    contig_locs: dict[str, ContigLocation],
    animal_id: str,
    haplotype: str,
    chroms: Sequence[str],
    schema: str,
    chrom_map: dict[str, str],
) -> list[OutputPlacement]:
    chrom_order = {c: i for i, c in enumerate(chroms)}
    valid_chroms = set(chroms)

    placements: list[OutputPlacement] = []
    seen_placed: set[str] = set()
    unplaced_counter = 1

    for cid, loc in contig_locs.items():
        mapped = chrom_map.get(cid)
        if mapped and mapped in valid_chroms:
            canon = mapped
        else:
            canon = canonicalize_chrom_name(cid, valid_chroms)

        if canon and canon in valid_chroms and canon not in seen_placed:
            seen_placed.add(canon)
            chrom = canon
            is_placed = True
            rank = (0, chrom_order[chrom], 1, loc.length, cid)
        else:
            chrom = f"unplaced_{unplaced_counter}"
            unplaced_counter += 1
            is_placed = False
            rank = (1, unplaced_counter, 1, loc.length, cid)

        header = format_pansn_header(
            animal_id, haplotype, chrom, 1, loc.length, schema=schema
        )
        placements.append(
            OutputPlacement(
                header=header,
                chrom=chrom,
                is_placed=is_placed,
                comp_id=cid,
                comp_start=1,
                comp_end=loc.length,
                orient="+",
                chrom_rank=rank,
            )
        )

    placements.sort(key=lambda p: p.chrom_rank)
    return placements


def _build_placements_intact(
    contig_locs: dict[str, ContigLocation],
    animal_id: str,
    haplotype: str,
    chroms: Sequence[str],
    chrom_map: dict[str, str],
) -> list[OutputPlacement]:
    """Keep all scaffolds intact, mapping canonical chroms while preserving other IDs."""
    chrom_order = {c: i for i, c in enumerate(chroms)}
    valid_chroms = set(chroms)

    placements: list[OutputPlacement] = []
    seen_placed: set[str] = set()
    seen_names: set[str] = set()

    for cid, loc in contig_locs.items():
        mapped = chrom_map.get(cid)
        canon = mapped if (mapped and mapped in valid_chroms) else canonicalize_chrom_name(cid, valid_chroms)

        if canon and canon in valid_chroms and canon not in seen_placed:
            seen_placed.add(canon)
            chrom = canon
            is_placed = True
            rank = (0, chrom_order[chrom], 0, 0, cid)
        else:
            sanitized = re.sub(r"[#\s]+", "_", cid).strip("_") or "scaffold"
            candidate = sanitized
            dup_idx = 2
            while candidate in seen_names:
                candidate = f"{sanitized}_dup{dup_idx}"
                dup_idx += 1
            chrom = candidate
            is_placed = False
            rank = (1, 0, -loc.length, 0, cid)

        seen_names.add(chrom)
        header = f"{animal_id}#{haplotype}#{chrom}"
        placements.append(
            OutputPlacement(
                header=header,
                chrom=chrom,
                is_placed=is_placed,
                comp_id=cid,
                comp_start=1,
                comp_end=loc.length,
                orient="+",
                chrom_rank=rank,
            )
        )

    placements.sort(key=lambda p: p.chrom_rank)
    return placements


def run_format(
    fasta_path: Path,
    animal_id: str,
    haplotype: str,
    output_path: Path,
    agp_path: Path | None = None,
    already_placed: bool = False,
    intact: bool = False,
    chrom_map_path: Path | None = None,
    optional_chrom_map_path: Path | None = None,
    schema: str = "cactus",
    chroms: Sequence[str] = DEFAULT_CHROMOSOMES,
    threads: int = 2,
    line_width: int = 0,
    report_path: Path | None = None,
) -> dict[str, object]:
    mode_count = int(agp_path is not None) + int(already_placed) + int(intact)
    if mode_count != 1:
        raise ValueError("Specify exactly one of --agp, --already-placed, or --intact")
    if not shutil.which("bgzip"):
        raise RuntimeError("'bgzip' not found on PATH")
    if not shutil.which("samtools"):
        raise RuntimeError("'samtools' not found on PATH")

    valid_chroms = set(chroms)
    chrom_map: dict[str, str] = {}
    if chrom_map_path is not None:
        chrom_map = load_chrom_map(chrom_map_path, valid_chroms)
    elif optional_chrom_map_path is not None and optional_chrom_map_path.exists() and optional_chrom_map_path.stat().st_size > 0:
        chrom_map = load_chrom_map(optional_chrom_map_path, valid_chroms)

    fetcher = SequenceFetcher(fasta_path)
    try:
        if intact:
            placements = _build_placements_intact(
                fetcher.locations, animal_id, haplotype, chroms, chrom_map
            )
        elif already_placed:
            placements = _build_placements_already_placed(
                fetcher.locations, animal_id, haplotype, chroms, schema, chrom_map
            )
        else:
            assert agp_path is not None
            placements = _build_placements_from_agp(
                agp_path, fetcher.locations, animal_id, haplotype, chroms, schema
            )

        output_path.parent.mkdir(parents=True, exist_ok=True)
        tmp_out = output_path.parent / f".{output_path.name}.tmp.{os.getpid()}"
        bgzip_cmd = ["bgzip", "-@", str(max(1, threads)), "-c"]

        placed_contigs = 0
        unplaced_contigs = 0
        placed_bp = 0
        unplaced_bp = 0
        revcomp_contigs = 0
        chroms_present: set[str] = set()

        try:
            with tmp_out.open("wb") as raw_out:
                bgzip_proc = subprocess.Popen(
                    bgzip_cmd, stdin=subprocess.PIPE, stdout=raw_out
                )
                assert bgzip_proc.stdin is not None

                for p in placements:
                    seq = fetcher.fetch(p.comp_id, p.comp_start, p.comp_end)
                    if p.orient == "-":
                        seq = reverse_complement(seq)
                        revcomp_contigs += 1

                    if p.is_placed:
                        placed_contigs += 1
                        placed_bp += len(seq)
                        chroms_present.add(p.chrom)
                    else:
                        unplaced_contigs += 1
                        unplaced_bp += len(seq)

                    bgzip_proc.stdin.write(f">{p.header}\n".encode("utf-8"))
                    if line_width > 0:
                        for i in range(0, len(seq), line_width):
                            bgzip_proc.stdin.write(seq[i : i + line_width].encode("ascii") + b"\n")
                    else:
                        bgzip_proc.stdin.write(seq.encode("ascii") + b"\n")

                bgzip_proc.stdin.close()
                rc = bgzip_proc.wait()
                if rc != 0:
                    raise subprocess.CalledProcessError(rc, bgzip_cmd)

            os.replace(tmp_out, output_path)
        except Exception:
            if tmp_out.exists():
                tmp_out.unlink()
            raise
    finally:
        fetcher.close()

    subprocess.run(["samtools", "faidx", str(output_path)], check=True)

    present_ordered = [c for c in chroms if c in chroms_present]
    missing_ordered = [c for c in chroms if c not in chroms_present]
    report: dict[str, object] = {
        "sample": f"{animal_id}_{haplotype}",
        "schema": schema,
        "placed_contigs": placed_contigs,
        "unplaced_contigs": unplaced_contigs,
        "placed_bp": placed_bp,
        "unplaced_bp": unplaced_bp,
        "revcomp_contigs": revcomp_contigs,
        "chroms_present": present_ordered,
        "chroms_missing": missing_ordered,
    }
    if report_path is not None:
        report_path.parent.mkdir(parents=True, exist_ok=True)
        report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    return report


def parse_pansn_chrom(header: str) -> str:
    """Extract chromosome token from cactus, vg, or underscore PanSN headers."""
    token = header.lstrip(">").strip().split(None, 1)[0]
    parts = token.split("#")
    third = parts[2] if len(parts) >= 3 else parts[-1]
    if ":" in third:
        return third.split(":", 1)[0]
    if "[" in third:
        return third.split("[", 1)[0]
    m_under = _PANSN_UNDERSCORE_COORD_RE.match(third)
    if m_under:
        return m_under.group(1)
    return third


def _write_empty_indexed_bgzf(fa_gz_path: Path) -> None:
    """Create a valid empty BGZF archive with empty .fai and valid 8-byte .gzi."""
    res = subprocess.run(
        ["bgzip", "-c"],
        input=b"",
        capture_output=True,
        check=True,
    )
    fa_gz_path.write_bytes(res.stdout)
    Path(f"{fa_gz_path}.fai").write_bytes(b"")
    Path(f"{fa_gz_path}.gzi").write_bytes(b"\x00" * 8)


def run_split(
    fasta_path: Path,
    out_dir: Path,
    sample: str,
    chroms: Sequence[str] = DEFAULT_CHROMOSOMES,
    threads: int = 1,
) -> dict[str, int]:
    """Split a PanSN FASTA into per-chromosome bgzip+faidx files in a single streaming pass."""
    if not shutil.which("bgzip"):
        raise RuntimeError("'bgzip' not found on PATH")
    if not shutil.which("samtools"):
        raise RuntimeError("'samtools' not found on PATH")

    out_dir.mkdir(parents=True, exist_ok=True)
    valid_chroms = set(chroms)
    tmp_handles: dict[str, Path] = {}
    counts: dict[str, int] = {c: 0 for c in chroms}

    open_fn = gzip.open if _is_gzip(fasta_path) else Path.open
    fhs = {}
    try:
        for c in chroms:
            tmp_fa = out_dir / f".{sample}.{c}.tmp.{os.getpid()}.fa"
            tmp_handles[c] = tmp_fa
            fhs[c] = tmp_fa.open("w", encoding="utf-8")

        active_fh = None
        with open_fn(fasta_path, "rt", encoding="utf-8") as fin:
            for line in fin:
                if line.startswith(">"):
                    chrom = parse_pansn_chrom(line)
                    if chrom in valid_chroms:
                        active_fh = fhs[chrom]
                        counts[chrom] += 1
                        active_fh.write(line)
                    else:
                        active_fh = None
                elif active_fh is not None:
                    active_fh.write(line)
    finally:
        for fh in fhs.values():
            fh.close()

    try:
        for c in chroms:
            dest = out_dir / f"{sample}.{c}.fa.gz"
            tmp_fa = tmp_handles[c]
            if counts[c] > 0:
                tmp_gz = out_dir / f".{sample}.{c}.tmp.{os.getpid()}.fa.gz"
                with tmp_fa.open("rb") as fin, tmp_gz.open("wb") as fout:
                    subprocess.run(
                        ["bgzip", "-@", str(max(1, threads)), "-c"],
                        stdin=fin,
                        stdout=fout,
                        check=True,
                    )
                os.replace(tmp_gz, dest)
                subprocess.run(["samtools", "faidx", str(dest)], check=True)
            else:
                _write_empty_indexed_bgzf(dest)
    finally:
        for tmp_fa in tmp_handles.values():
            if tmp_fa.exists():
                tmp_fa.unlink()

    return counts


def run_strip_pansn(
    fasta_path: Path,
    output_path: Path,
    threads: int = 2,
) -> int:
    """Strip leading '{sample}#{hap}#' prefix from PanSN FASTA headers for Cactus."""
    if not shutil.which("bgzip"):
        raise RuntimeError("'bgzip' not found on PATH")
    if not shutil.which("samtools"):
        raise RuntimeError("'samtools' not found on PATH")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_out = output_path.parent / f".{output_path.name}.tmp.{os.getpid()}"
    open_fn = gzip.open if _is_gzip(fasta_path) else Path.open
    bgzip_cmd = ["bgzip", "-@", str(max(1, threads)), "-c"]
    n_seqs = 0

    try:
        with tmp_out.open("wb") as raw_out:
            bgzip_proc = subprocess.Popen(
                bgzip_cmd, stdin=subprocess.PIPE, stdout=raw_out
            )
            assert bgzip_proc.stdin is not None
            try:
                with open_fn(fasta_path, "rt", encoding="utf-8") as fin:
                    for line in fin:
                        if line.startswith(">"):
                            token = line[1:].strip().split(None, 1)[0]
                            parts = token.split("#", 2)
                            stripped = parts[2] if len(parts) == 3 else token
                            n_seqs += 1
                            bgzip_proc.stdin.write(f">{stripped}\n".encode("utf-8"))
                        else:
                            bgzip_proc.stdin.write(line.encode("ascii"))
                bgzip_proc.stdin.close()
                rc = bgzip_proc.wait()
                if rc != 0:
                    raise subprocess.CalledProcessError(rc, bgzip_cmd)
            except Exception:
                if not bgzip_proc.stdin.closed:
                    bgzip_proc.stdin.close()
                bgzip_proc.terminate()
                bgzip_proc.wait()
                raise
        os.replace(tmp_out, output_path)
    except Exception:
        if tmp_out.exists():
            tmp_out.unlink()
        raise

    subprocess.run(["samtools", "faidx", str(output_path)], check=True)
    return n_seqs


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Format assemblies to PanSN and split by chromosome.")
    sub = parser.add_subparsers(dest="command", required=True)

    p_fmt = sub.add_parser("format", help="Rename and orient assembly sequences into PanSN format.")
    p_fmt.add_argument("--fasta", required=True, type=Path, help="Input FASTA (.fa or .fa.gz).")
    p_fmt.add_argument("--animal-id", required=True, help="Animal ID.")
    p_fmt.add_argument("--haplotype", required=True, help="Haplotype ID.")
    p_fmt.add_argument("--output", required=True, type=Path, help="Output bgzipped FASTA.")
    p_fmt.add_argument("--agp", type=Path, help="RagTag scaffold AGP file.")
    p_fmt.add_argument(
        "--already-placed",
        action="store_true",
        help="Skip AGP and map existing chromosome headers with PanSN coordinates.",
    )
    p_fmt.add_argument(
        "--intact",
        action="store_true",
        help="Keep scaffolds intact without RagTag and emit Animal#Hap#Contig headers.",
    )
    p_fmt.add_argument("--chrom-map", type=Path, help="Required chromosome alias map.")
    p_fmt.add_argument(
        "--optional-chrom-map",
        type=Path,
        help="Optional chromosome alias map (ignored if missing or empty).",
    )
    p_fmt.add_argument("--schema", choices=("cactus", "vg", "_"), default="cactus")
    p_fmt.add_argument(
        "--chroms",
        default=",".join(DEFAULT_CHROMOSOMES),
        help="Comma-separated canonical chromosomes in sort order.",
    )
    p_fmt.add_argument("--threads", type=int, default=2, help="bgzip compression threads.")
    p_fmt.add_argument("--line-width", type=int, default=0, help="Wrap width (0 = unwrapped).")
    p_fmt.add_argument("--report", type=Path, help="Optional JSON report output.")

    p_spl = sub.add_parser("split", help="Split a PanSN FASTA into per-chromosome bgzipped FASTAs.")
    p_spl.add_argument("--fasta", required=True, type=Path, help="Input PanSN .fa.gz.")
    p_spl.add_argument("--out-dir", required=True, type=Path, help="Output directory for splits.")
    p_spl.add_argument("--sample", required=True, help="Sample prefix (<AnimalID>_<Haplotype>).")
    p_spl.add_argument(
        "--chroms",
        default=",".join(DEFAULT_CHROMOSOMES),
        help="Comma-separated canonical chromosomes.",
    )
    p_spl.add_argument("--threads", type=int, default=1, help="bgzip compression threads.")

    p_strip = sub.add_parser(
        "strip-pansn",
        help="Strip '{sample}#{hap}#' prefix from PanSN FASTA headers for Cactus.",
    )
    p_strip.add_argument("--fasta", required=True, type=Path, help="Input PanSN .fa.gz.")
    p_strip.add_argument("--output", required=True, type=Path, help="Output bgzipped FASTA.")
    p_strip.add_argument("--threads", type=int, default=2, help="bgzip compression threads.")

    args = parser.parse_args(argv)

    if args.command == "format":
        chroms = [c.strip() for c in args.chroms.split(",") if c.strip()]
        run_format(
            fasta_path=args.fasta,
            animal_id=args.animal_id,
            haplotype=args.haplotype,
            output_path=args.output,
            agp_path=args.agp,
            already_placed=args.already_placed,
            intact=args.intact,
            chrom_map_path=args.chrom_map,
            optional_chrom_map_path=args.optional_chrom_map,
            schema=args.schema,
            chroms=chroms,
            threads=args.threads,
            line_width=args.line_width,
            report_path=args.report,
        )
        return 0

    if args.command == "split":
        chroms = [c.strip() for c in args.chroms.split(",") if c.strip()]
        run_split(
            fasta_path=args.fasta,
            out_dir=args.out_dir,
            sample=args.sample,
            chroms=chroms,
            threads=args.threads,
        )
        return 0

    if args.command == "strip-pansn":
        run_strip_pansn(
            fasta_path=args.fasta,
            output_path=args.output,
            threads=args.threads,
        )
        return 0

    return 1


if __name__ == "__main__":
    sys.exit(main())
