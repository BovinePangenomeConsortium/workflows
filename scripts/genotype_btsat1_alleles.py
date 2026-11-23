#!/usr/bin/env python3
"""Direct assembly genotyping of BTSAT1 (1.715) 1400 vs 1401 bp haplotypes.

Scans an assembly FASTA (.fa or .fa.gz/.bgz) in a single regex pass for the
two-locus diagnostic haplotypes spanning the pos 1019 homopolymer (6A vs 7A)
and pos 1034 transversion (T vs G), including the recombinant/error controls
(6A+G and 7A+T). Reports per-chromosome and genome-wide allele counts,
transition purity, homopolymer error rate, and empirical physical monomer
spacings (pos[i+1] - pos[i]).
"""

from __future__ import annotations

import argparse
import collections
import csv
import gzip
import re
from typing import Dict, Iterable, List, Tuple

_PLACED_CHROMS = set(map(str, range(1, 30))) | {"X", "Y", "MT"}
_COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}

# Matches forward GG(A{4,9})CCCCAGGTTCCAAA([ACGT])ACAG or its reverse complement
_DIAG_RE = re.compile(
    r"GG(A{4,9})CCCCAGGTTCCAAA([ACGT])ACAG"
    r"|CTGT([ACGT])TTTGGAACCTGGGG(T{4,9})CC"
)


def parse_chrom(contig_name: str) -> Tuple[str, bool]:
    parts = contig_name.split("#")
    suffix = parts[2] if len(parts) >= 3 else parts[-1]
    base = suffix.split(":")[0].split("[")[0]
    if base in _PLACED_CHROMS:
        return base, True
    return "unplaced", False


def iter_fasta(path: str) -> Iterable[Tuple[str, str]]:
    with open(path, "rb") as probe:
        magic = probe.read(2)
    opener = gzip.open if magic == b"\x1f\x8b" else open

    hdr: str | None = None
    chunks: List[str] = []
    with opener(path, "rt") as fin:
        for line in fin:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if hdr is not None:
                    yield hdr, "".join(chunks).upper()
                hdr = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        if hdr is not None:
            yield hdr, "".join(chunks).upper()


def classify_match(m: re.Match[str]) -> Tuple[str, str, int]:
    """Return (strand, state_label, anchor_coord) for a diagnostic regex hit."""
    if m.group(1) is not None:
        n_a = len(m.group(1))
        snp = m.group(2)
        strand = "+"
        anchor = m.start()
    else:
        n_a = len(m.group(4))
        snp = _COMP[m.group(3)]
        strand = "-"
        anchor = m.end()

    if n_a == 6 and snp == "T":
        state = "6A_T"
    elif n_a == 7 and snp == "G":
        state = "7A_G"
    elif n_a == 6 and snp == "G":
        state = "6A_G"
    elif n_a == 7 and snp == "T":
        state = "7A_T"
    else:
        state = "other"
    return strand, state, anchor


def make_empty_bucket() -> Dict[str, object]:
    return {
        "counts": collections.Counter(),
        "trans": collections.Counter(),
        "spacings_6AT": collections.Counter(),
        "spacings_7AG": collections.Counter(),
    }


def format_allele_row(
    sample: str, chrom: str, is_placed_str: str, bucket: Dict[str, object]
) -> List[str]:
    counts: collections.Counter[str] = bucket["counts"]  # type: ignore[assignment]
    trans: collections.Counter[str] = bucket["trans"]  # type: ignore[assignment]
    sp_6at: collections.Counter[int] = bucket["spacings_6AT"]  # type: ignore[assignment]
    sp_7ag: collections.Counter[int] = bucket["spacings_7AG"]  # type: ignore[assignment]

    n_6at = counts["6A_T"]
    n_7ag = counts["7A_G"]
    n_6ag = counts["6A_G"]
    n_7at = counts["7A_T"]
    n_other = counts["other"]

    n_canon = n_6at + n_7ag
    n_ctrl = n_6ag + n_7at
    frac_7ag = f"{n_7ag / n_canon:.6f}" if n_canon > 0 else "NA"
    ctrl_rate = f"{n_ctrl / (n_canon + n_ctrl):.6f}" if (n_canon + n_ctrl) > 0 else "NA"

    t_6at = trans["6A_T->6A_T"]
    t_7ag = trans["7A_G->7A_G"]
    t_mix = trans["mixed_6AT_7AG"]
    t_tot = t_6at + t_7ag + t_mix
    purity = f"{(t_6at + t_7ag) / t_tot:.6f}" if t_tot > 0 else "NA"

    mode_6at = str(sp_6at.most_common(1)[0][0]) if sp_6at else "NA"
    mode_7ag = str(sp_7ag.most_common(1)[0][0]) if sp_7ag else "NA"

    return [
        sample,
        chrom,
        is_placed_str,
        str(n_6at),
        str(n_7ag),
        str(n_6ag),
        str(n_7at),
        str(n_other),
        frac_7ag,
        ctrl_rate,
        str(t_6at),
        str(t_7ag),
        str(t_mix),
        purity,
        str(sp_6at[1400]),
        str(sp_6at[1401]),
        mode_6at,
        str(sp_7ag[1400]),
        str(sp_7ag[1401]),
        mode_7ag,
    ]


def run_genotype(
    fasta_path: str,
    sample: str,
    out_alleles: str,
    out_spacing: str,
    min_spacing: int,
    max_spacing: int,
) -> None:
    chrom_buckets: Dict[Tuple[str, bool], Dict[str, object]] = {}
    genome_bucket = make_empty_bucket()
    spacing_hist: Dict[str, collections.Counter[int]] = collections.defaultdict(
        collections.Counter
    )

    for contig_name, seq in iter_fasta(fasta_path):
        chrom, is_placed = parse_chrom(contig_name)
        key = (chrom, is_placed)
        bucket = chrom_buckets.get(key)
        if bucket is None:
            bucket = make_empty_bucket()
            chrom_buckets[key] = bucket

        prev_by_strand: Dict[str, Tuple[int, str]] = {}
        for m in _DIAG_RE.finditer(seq):
            strand, state, anchor = classify_match(m)
            bucket["counts"][state] += 1  # type: ignore[index]
            genome_bucket["counts"][state] += 1  # type: ignore[index]

            prev = prev_by_strand.get(strand)
            if prev is not None:
                prev_anchor, prev_state = prev
                dist = anchor - prev_anchor
                if min_spacing <= dist <= max_spacing:
                    if prev_state == "6A_T" and state == "6A_T":
                        pair_type = "6A_T->6A_T"
                        bucket["trans"][pair_type] += 1  # type: ignore[index]
                        genome_bucket["trans"][pair_type] += 1  # type: ignore[index]
                        bucket["spacings_6AT"][dist] += 1  # type: ignore[index]
                        genome_bucket["spacings_6AT"][dist] += 1  # type: ignore[index]
                    elif prev_state == "7A_G" and state == "7A_G":
                        pair_type = "7A_G->7A_G"
                        bucket["trans"][pair_type] += 1  # type: ignore[index]
                        genome_bucket["trans"][pair_type] += 1  # type: ignore[index]
                        bucket["spacings_7AG"][dist] += 1  # type: ignore[index]
                        genome_bucket["spacings_7AG"][dist] += 1  # type: ignore[index]
                    elif {prev_state, state} == {"6A_T", "7A_G"}:
                        pair_type = "mixed_6AT_7AG"
                        bucket["trans"][pair_type] += 1  # type: ignore[index]
                        genome_bucket["trans"][pair_type] += 1  # type: ignore[index]
                    else:
                        pair_type = "other"
                    spacing_hist[pair_type][dist] += 1
            prev_by_strand[strand] = (anchor, state)

    def chrom_key(item: Tuple[Tuple[str, bool], Dict[str, object]]) -> Tuple[int, str]:
        c = item[0][0]
        if c.isdigit():
            return (0, f"{int(c):03d}")
        if c == "X":
            return (1, "0")
        if c == "Y":
            return (1, "1")
        if c == "MT":
            return (1, "2")
        return (2, c)

    with open(out_alleles, "w", newline="") as fout:
        writer = csv.writer(fout, delimiter="\t")
        writer.writerow(
            [
                "sample",
                "chrom",
                "is_placed",
                "n_6A_T",
                "n_7A_G",
                "n_6A_G",
                "n_7A_T",
                "n_other",
                "frac_7A_G",
                "control_rate",
                "trans_6AT_6AT",
                "trans_7AG_7AG",
                "trans_mixed",
                "purity_frac",
                "n_pairs_6AT_1400",
                "n_pairs_6AT_1401",
                "mode_spacing_6AT",
                "n_pairs_7AG_1400",
                "n_pairs_7AG_1401",
                "mode_spacing_7AG",
            ]
        )
        for (chrom, is_placed), bucket in sorted(chrom_buckets.items(), key=chrom_key):
            # Only emit chromosome rows that actually have at least one diagnostic hit
            if sum(bucket["counts"].values()) == 0:  # type: ignore[union-attr]
                continue
            writer.writerow(
                format_allele_row(
                    sample, chrom, "TRUE" if is_placed else "FALSE", bucket
                )
            )
        writer.writerow(format_allele_row(sample, "GENOME_TOTAL", "ALL", genome_bucket))

    with open(out_spacing, "w", newline="") as fout:
        writer = csv.writer(fout, delimiter="\t")
        writer.writerow(["sample", "pair_type", "spacing_bp", "count"])
        for pair_type in sorted(spacing_hist):
            for dist, cnt in sorted(spacing_hist[pair_type].items()):
                writer.writerow([sample, pair_type, str(dist), str(cnt)])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", required=True, help="Assembly FASTA (.fa or .fa.gz)")
    parser.add_argument("--sample", required=True, help="Sample identifier")
    parser.add_argument("--out-alleles", required=True, help="Per-chrom & genome allele TSV")
    parser.add_argument("--out-spacing", required=True, help="Monomer spacing histogram TSV")
    parser.add_argument("--min-spacing", type=int, default=1200, help="Min consecutive spacing (bp)")
    parser.add_argument("--max-spacing", type=int, default=1600, help="Max consecutive spacing (bp)")
    args = parser.parse_args()
    run_genotype(
        fasta_path=args.fasta,
        sample=args.sample,
        out_alleles=args.out_alleles,
        out_spacing=args.out_spacing,
        min_spacing=args.min_spacing,
        max_spacing=args.max_spacing,
    )


if __name__ == "__main__":
    main()
