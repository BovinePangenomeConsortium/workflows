#!/usr/bin/env python3
"""Summarize an SRF BED file by satellite family/subfamily and chromosome.

Works both on pre-existing per-sample `.srf-aln.bed` files (when `--membership`
from `cluster_srf_motifs.py` is supplied) and on unified-catalog `.bed` files
directly. Normalizes satellite coverage against exact chromosome and total
assembly lengths parsed from the assembly `.fai`.
"""

from __future__ import annotations

import argparse
import collections
import csv
import re
from typing import Dict, Optional, Tuple

_PLACED_CHROMS = set(map(str, range(1, 30))) | {"X", "Y", "MT"}
_FAM_PREFIX_RE = re.compile(r"^(C\d+)_S\d+")
_SUBFAM_PREFIX_RE = re.compile(r"^(C\d+_S\d+)")


def parse_chrom(contig_name: str) -> Tuple[str, bool]:
    """Extract canonical chromosome label and placed flag from panSN or raw name."""
    parts = contig_name.split("#")
    suffix = parts[2] if len(parts) >= 3 else parts[-1]
    base = suffix.split(":")[0].split("[")[0]
    if base in _PLACED_CHROMS:
        return base, True
    return "unplaced", False


def load_fai(fai_path: str) -> Tuple[int, Dict[str, int]]:
    """Return (assembly_bp, chrom_bp_by_chrom) from a FASTA .fai index."""
    assembly_bp = 0
    chrom_bp: Dict[str, int] = collections.defaultdict(int)
    with open(fai_path) as fin:
        for line in fin:
            line = line.strip()
            if not line:
                continue
            fields = line.split("\t")
            length = int(fields[1])
            chrom, _ = parse_chrom(fields[0])
            assembly_bp += length
            chrom_bp[chrom] += length
    return assembly_bp, dict(chrom_bp)


def load_membership(
    membership_path: str,
) -> Dict[str, Tuple[str, str, str, str]]:
    """Map original_id -> (family_id, subfamily_id, catalog_id, ref_name)."""
    mapping: Dict[str, Tuple[str, str, str, str]] = {}
    with open(membership_path, newline="") as fin:
        reader = csv.DictReader(fin, delimiter="\t")
        for row in reader:
            mapping[row["original_id"]] = (
                row["family_id"],
                row["subfamily_id"],
                row["catalog_id"],
                row["ref_name"],
            )
    return mapping


def load_catalog_refs(catalog_fa_path: str) -> Dict[str, str]:
    """Extract catalog_id -> ref_name from catalog FASTA headers if provided."""
    ref_map: Dict[str, str] = {}
    with open(catalog_fa_path) as fin:
        for line in fin:
            if not line.startswith(">"):
                continue
            parts = line[1:].strip().split()
            cid = parts[0]
            ref_name = "unannotated"
            for tok in parts[1:]:
                if tok.startswith("ref="):
                    val = tok[4:]
                    ref_name = val.split("(")[0]
                    break
            ref_map[cid] = ref_name
    return ref_map


def extract_bed_div(fields: list[str], interval_len: int) -> Optional[float]:
    """Extract numeric divergence from srfutils.js paf2bed extra columns if present."""
    if len(fields) <= 4:
        return None
    candidates: list[Tuple[int, float]] = []
    for idx in range(4, len(fields)):
        val_str = fields[idx]
        if val_str in ("+", "-", ".", "NA", "NaN"):
            continue
        try:
            val = float(val_str)
        except ValueError:
            continue
        candidates.append((idx, val))
    if not candidates:
        return None
    # If col 4 is an integer equal to interval_len and col 5 is a float, prefer col 5
    if (
        len(candidates) >= 2
        and candidates[0][0] == 4
        and abs(candidates[0][1] - interval_len) < 1e-6
    ):
        return candidates[1][1]
    return candidates[0][1]


def summarize_bed(
    bed_path: str,
    fai_path: str,
    sample: str,
    output_path: str,
    membership_path: Optional[str] = None,
    catalog_fa_path: Optional[str] = None,
) -> None:
    assembly_bp, chrom_bp_map = load_fai(fai_path)
    membership = load_membership(membership_path) if membership_path else {}
    catalog_refs = load_catalog_refs(catalog_fa_path) if catalog_fa_path else {}

    # Key: (family_id, subfamily_id, catalog_id, ref_name, chrom, is_placed)
    agg: Dict[Tuple[str, str, str, str, str, bool], list] = {}

    with open(bed_path) as fin:
        for line in fin:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 4:
                continue
            ctg = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            motif_id = fields[3]
            length = end - start
            if length <= 0:
                continue

            chrom, is_placed = parse_chrom(ctg)

            if motif_id in membership:
                family_id, subfamily_id, catalog_id, ref_name = membership[motif_id]
            else:
                catalog_id = motif_id
                m_sub = _SUBFAM_PREFIX_RE.match(motif_id)
                m_fam = _FAM_PREFIX_RE.match(motif_id)
                subfamily_id = m_sub.group(1) if m_sub else motif_id
                family_id = m_fam.group(1) if m_fam else motif_id
                ref_name = catalog_refs.get(catalog_id, "unannotated")

            div = extract_bed_div(fields, length)
            key = (family_id, subfamily_id, catalog_id, ref_name, chrom, is_placed)
            entry = agg.get(key)
            if entry is None:
                # [total_bp, n_intervals, div_weighted_sum, div_bp_sum]
                entry = [0, 0, 0.0, 0]
                agg[key] = entry
            entry[0] += length
            entry[1] += 1
            if div is not None:
                entry[2] += div * length
                entry[3] += length

    def chrom_sort_key(item: Tuple[Tuple[str, str, str, str, str, bool], list]) -> Tuple[str, str, int, str]:
        fam, subfam, cid, _, chrom, is_placed = item[0]
        if chrom.isdigit():
            c_ord = (0, int(chrom), chrom)
        elif chrom == "X":
            c_ord = (1, 0, chrom)
        elif chrom == "Y":
            c_ord = (1, 1, chrom)
        elif chrom == "MT":
            c_ord = (1, 2, chrom)
        else:
            c_ord = (2, 0, chrom)
        return (fam, subfam, c_ord[0] * 100 + c_ord[1], cid)

    with open(output_path, "w", newline="") as fout:
        writer = csv.writer(fout, delimiter="\t")
        writer.writerow(
            [
                "sample",
                "family_id",
                "subfamily_id",
                "catalog_id",
                "ref_name",
                "chrom",
                "is_placed",
                "bp",
                "n_intervals",
                "weighted_mean_div",
                "chrom_bp",
                "assembly_bp",
                "chrom_frac",
                "genome_frac",
            ]
        )
        for key, stats in sorted(agg.items(), key=chrom_sort_key):
            family_id, subfamily_id, catalog_id, ref_name, chrom, is_placed = key
            bp, n_int, div_wsum, div_bp = stats
            mean_div = f"{div_wsum / div_bp:.6f}" if div_bp > 0 else "-1"
            c_bp = chrom_bp_map.get(chrom, 0)
            c_frac = f"{bp / c_bp:.8f}" if c_bp > 0 else "0"
            g_frac = f"{bp / assembly_bp:.8f}" if assembly_bp > 0 else "0"
            writer.writerow(
                [
                    sample,
                    family_id,
                    subfamily_id,
                    catalog_id,
                    ref_name,
                    chrom,
                    "TRUE" if is_placed else "FALSE",
                    str(bp),
                    str(n_int),
                    mean_div,
                    str(c_bp),
                    str(assembly_bp),
                    c_frac,
                    g_frac,
                ]
            )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bed", required=True, help="Input .srf-aln.bed file")
    parser.add_argument("--fai", required=True, help="Input assembly .fai index")
    parser.add_argument("--sample", required=True, help="Sample identifier")
    parser.add_argument("--output", required=True, help="Output TSV file")
    parser.add_argument(
        "--membership",
        help="Optional <graph>.membership.tsv to map per-sample motif IDs to canonical families",
    )
    parser.add_argument(
        "--catalog",
        help="Optional <graph>.catalog.fa to annotate ref_name when --membership is not used",
    )
    args = parser.parse_args()
    summarize_bed(
        bed_path=args.bed,
        fai_path=args.fai,
        sample=args.sample,
        output_path=args.output,
        membership_path=args.membership,
        catalog_fa_path=args.catalog,
    )


if __name__ == "__main__":
    main()
