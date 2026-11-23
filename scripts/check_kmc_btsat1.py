#!/usr/bin/env python3
"""Query a KMC k-mer dump (.kmc.txt) for central-offset BTSAT1 diagnostic alleles.

Because KMC stores canonical k-mers, selecting a single deterministic motif
offset pair `(off_fwd, off_rc)` where `off_rc = k_len - off_fwd - motif_len`
ensures each assembled BTSAT1 monomer copy maps to at most one canonical
k-mer in `.kmc.txt`, avoiding ~120x double-counting across sliding windows.
"""

from __future__ import annotations

import argparse
import csv
import gzip
from typing import Dict, List, Tuple

_RC_TRANS = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq: str) -> str:
    return seq.translate(_RC_TRANS)[::-1]


_MOTIFS: List[Tuple[str, str]] = [
    ("6A_T", "GGAAAAAACCCCAGGTTCCAAATACAG"),
    ("7A_G", "GGAAAAAAACCCCAGGTTCCAAAGACAG"),
    ("6A_G", "GGAAAAAACCCCAGGTTCCAAAGACAG"),
    ("7A_T", "GGAAAAAAACCCCAGGTTCCAAATACAG"),
]
_CORE_FWD = "CCCCAGGTTCCAAA"
_CORE_REV = revcomp(_CORE_FWD)


def check_kmc(kmc_path: str, sample: str, output_path: str) -> None:
    stats: Dict[str, List[int]] = {
        name: [0, 0, 0] for name, _ in _MOTIFS  # [n_distinct, total_count, max_count]
    }
    rc_motifs = {name: (seq, revcomp(seq), len(seq)) for name, seq in _MOTIFS}

    with open(kmc_path, "rb") as probe:
        magic = probe.read(2)
    opener = gzip.open if magic == b"\x1f\x8b" else open

    with opener(kmc_path, "rt") as fin:
        for line in fin:
            if _CORE_FWD not in line and _CORE_REV not in line:
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            kmer = parts[0].upper()
            count = int(parts[1])
            k_len = len(kmer)

            for name, (mf, mr, m_len) in rc_motifs.items():
                off_fwd = (k_len - m_len) // 2
                off_rc = k_len - off_fwd - m_len
                if k_len < m_len:
                    continue
                if (
                    kmer[off_fwd : off_fwd + m_len] == mf
                    or kmer[off_rc : off_rc + m_len] == mr
                ):
                    st = stats[name]
                    st[0] += 1
                    st[1] += count
                    if count > st[2]:
                        st[2] = count
                    break

    tot_6at = stats["6A_T"][1]
    tot_7ag = stats["7A_G"][1]
    max_6at = stats["6A_T"][2]
    max_7ag = stats["7A_G"][2]

    frac_total = (
        f"{tot_7ag / (tot_6at + tot_7ag):.6f}" if (tot_6at + tot_7ag) > 0 else "NA"
    )
    frac_max = (
        f"{max_7ag / (max_6at + max_7ag):.6f}" if (max_6at + max_7ag) > 0 else "NA"
    )

    with open(output_path, "w", newline="") as fout:
        writer = csv.writer(fout, delimiter="\t")
        writer.writerow(
            [
                "sample",
                "state",
                "n_distinct_kmers",
                "total_kmer_count",
                "max_kmer_count",
                "frac_7A_G_by_total",
                "frac_7A_G_by_max",
            ]
        )
        for name, _ in _MOTIFS:
            n_dist, tot_cnt, max_cnt = stats[name]
            writer.writerow(
                [
                    sample,
                    name,
                    str(n_dist),
                    str(tot_cnt),
                    str(max_cnt),
                    frac_total,
                    frac_max,
                ]
            )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--kmc", required=True, help="Input KMC dump .kmc.txt")
    parser.add_argument("--sample", required=True, help="Sample identifier")
    parser.add_argument("--output", required=True, help="Output TSV file")
    args = parser.parse_args()
    check_kmc(args.kmc, args.sample, args.output)


if __name__ == "__main__":
    main()
