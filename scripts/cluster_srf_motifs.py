#!/usr/bin/env python3
"""Rotation- and strand-invariant clustering of SRF satellite motifs.

Clusters per-sample SRF consensus motifs into canonical satellite families using
circular canonical k-mer containment, orients and rotates all motifs onto a shared
strand/origin (anchored to verified published reference sequences when available),
splits families into sequence-identity subfamilies, and emits a unified catalog.
"""

from __future__ import annotations

import argparse
import collections
import csv
import difflib
import os
import re
import sys
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Set, Tuple

_RC_TRANS = str.maketrans("ACGTNacgtn", "TGCANtgcan")
_SRF_HDR_RE = re.compile(
    r"^>(?P<orig_id>(?P<sample>\S+?)#circ(?P<contig>\d+)-(?P<length>\d+))"
    r"(?:\s+min=(?P<min>\d+),max=(?P<max>\d+),avg=(?P<avg>\d+))?"
)
_ABUN_RE = re.compile(
    r"^(?P<orig_id>(?P<sample>\S+?)#circ(?P<contig>\d+)-(?P<length>\d+))\t"
    r"(?P<bp>\d+)\t(?P<mean_div>\S+)\t(?P<filtered>\S+)\t(?P<unfiltered>\S+)"
)


@dataclass
class MotifRecord:
    original_id: str
    sample: str
    contig: str
    length: int
    seq: str
    bp: int
    mean_divergence: float
    k_min: int = 0
    k_max: int = 0
    k_avg: int = 0


@dataclass
class RefRecord:
    short_name: str
    accession: str
    species: str
    seq: str
    k15: Set[str]
    k13: Set[str]


def revcomp(seq: str) -> str:
    return seq.translate(_RC_TRANS)[::-1]


def circular_canonical_kmers(seq: str, k: int) -> Set[str]:
    n = len(seq)
    if n == 0:
        return set()
    k_eff = min(k, n)
    ext = seq + seq[: k_eff - 1]
    out: Set[str] = set()
    for i in range(n):
        km = ext[i : i + k_eff]
        rkm = revcomp(km)
        out.add(km if km < rkm else rkm)
    return out


def orient_and_rotate(
    seq: str, ref: str, k: int = 15, max_kmer_occ: int = 25
) -> Tuple[str, str, int, int]:
    """Orient and circularly rotate `seq` to match `ref` via modal k-mer offset voting.

    Returns (rotated_seq, strand, shift, best_hits).
    """
    n_seq = len(seq)
    n_ref = len(ref)
    if n_seq == 0 or n_ref == 0:
        return seq, "+", 0, 0
    k_eff = min(k, n_seq, n_ref)
    ref_ext = ref + ref[: k_eff - 1]
    ref_pos: Dict[str, List[int]] = collections.defaultdict(list)
    for i in range(n_ref):
        ref_pos[ref_ext[i : i + k_eff]].append(i)

    best_strand = "+"
    best_hits = -1
    best_shift = 0
    for strand, candidate in (("+", seq), ("-", revcomp(seq))):
        c_ext = candidate + candidate[: k_eff - 1]
        shifts: collections.Counter[int] = collections.Counter()
        hits = 0
        for j in range(n_seq):
            km = c_ext[j : j + k_eff]
            positions = ref_pos.get(km)
            if not positions:
                continue
            hits += 1
            if len(positions) <= max_kmer_occ:
                for i in positions:
                    shifts[(j - i) % n_seq] += 1
        if hits > best_hits:
            best_hits = hits
            best_strand = strand
            best_shift = shifts.most_common(1)[0][0] if shifts else 0

    oriented = seq if best_strand == "+" else revcomp(seq)
    rotated = oriented[best_shift:] + oriented[:best_shift]
    return rotated, best_strand, best_shift, best_hits


def parse_fasta_iter(path: str) -> Iterable[Tuple[str, str]]:
    hdr: Optional[str] = None
    chunks: List[str] = []
    with open(path) as fin:
        for line in fin:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if hdr is not None:
                    yield hdr, "".join(chunks).upper()
                hdr = line
                chunks = []
            else:
                chunks.append(line)
        if hdr is not None:
            yield hdr, "".join(chunks).upper()


def load_refs(refs_path: Optional[str]) -> List[RefRecord]:
    if not refs_path:
        return []
    refs: List[RefRecord] = []
    for hdr, seq in parse_fasta_iter(refs_path):
        token = hdr[1:].split()[0]
        parts = token.split("|")
        short_name = parts[0]
        acc = parts[1] if len(parts) > 1 else short_name
        species = parts[2] if len(parts) > 2 else "unknown"
        refs.append(
            RefRecord(
                short_name=short_name,
                accession=acc,
                species=species,
                seq=seq,
                k15=circular_canonical_kmers(seq, 15),
                k13=circular_canonical_kmers(seq, 13),
            )
        )
    return refs


def load_motifs_from_csv(csv_path: str) -> List[MotifRecord]:
    records: List[MotifRecord] = []
    with open(csv_path, newline="") as fin:
        reader = csv.DictReader(fin)
        for row in reader:
            seq = row["fasta"].strip().upper()
            sample = row["sample"].strip()
            contig = row["contig"].strip()
            length = int(row["length"])
            div_raw = row.get("mean_divergence", row.get("mean_identity", "-1"))
            try:
                mean_div = float(div_raw)
            except ValueError:
                mean_div = -1.0
            records.append(
                MotifRecord(
                    original_id=f"{sample}#circ{contig}-{length}",
                    sample=sample,
                    contig=contig,
                    length=length,
                    seq=seq,
                    bp=int(row["bp"]),
                    mean_divergence=mean_div,
                    k_min=int(row.get("min", 0) or 0),
                    k_max=int(row.get("max", 0) or 0),
                    k_avg=int(row.get("avg", 0) or 0),
                )
            )
    return records


def load_motifs_from_srf(
    srf_paths: List[str], abun_paths: Optional[List[str]]
) -> List[MotifRecord]:
    abun_map: Dict[str, Tuple[int, float]] = {}
    if abun_paths:
        for ap in abun_paths:
            with open(ap) as fin:
                for line in fin:
                    m = _ABUN_RE.match(line.rstrip("\n"))
                    if not m:
                        continue
                    div_str = m.group("mean_div")
                    mean_div = -1.0 if div_str == "NaN" else float(div_str)
                    abun_map[m.group("orig_id")] = (int(m.group("bp")), mean_div)

    records: List[MotifRecord] = []
    for sp in srf_paths:
        for hdr, seq in parse_fasta_iter(sp):
            m = _SRF_HDR_RE.match(hdr)
            if m:
                orig_id = m.group("orig_id")
                sample = m.group("sample")
                contig = m.group("contig")
                length = int(m.group("length"))
                k_min = int(m.group("min") or 0)
                k_max = int(m.group("max") or 0)
                k_avg = int(m.group("avg") or 0)
            else:
                orig_id = hdr[1:].split()[0]
                sample = os.path.basename(sp).split(".")[0]
                contig = orig_id
                length = len(seq)
                k_min = k_max = k_avg = 0

            bp, mean_div = abun_map.get(orig_id, (k_avg * length, -1.0))
            records.append(
                MotifRecord(
                    original_id=orig_id,
                    sample=sample,
                    contig=contig,
                    length=length,
                    seq=seq,
                    bp=bp,
                    mean_divergence=mean_div,
                    k_min=k_min,
                    k_max=k_max,
                    k_avg=k_avg,
                )
            )
    return records


def match_reference(
    seq: str, refs: List[RefRecord], k_orient: int
) -> Optional[Tuple[RefRecord, float, float, bool, str, int]]:
    """Find best matching published reference for `seq`.

    Returns (ref, report_identity, global_ratio, is_full_length, strand, shift) or None.
    """
    if not refs or not seq:
        return None
    use_k13 = len(seq) < 100
    qk = circular_canonical_kmers(seq, 13 if use_k13 else 15)
    if not qk:
        return None

    best_full: Optional[Tuple[float, float, RefRecord, str, int]] = None
    best_partial: Optional[Tuple[float, float, RefRecord, str, int]] = None

    for ref in refs:
        rk = ref.k13 if (use_k13 or len(ref.seq) < 100) else ref.k15
        if not rk:
            continue
        inter = len(qk & rk)
        if inter == 0:
            continue
        cont = inter / min(len(qk), len(rk))
        if cont < 0.15:
            continue

        len_ratio = min(len(seq), len(ref.seq)) / max(len(seq), len(ref.seq))
        if len_ratio >= 0.70:
            rot_q, strand, shift, _ = orient_and_rotate(seq, ref.seq, k=min(k_orient, 13))
            global_id = seq_identity(ref.seq, rot_q)
            matches = global_id * (len(seq) + len(ref.seq)) / 2.0
            overlap_id = min(1.0, matches / min(len(seq), len(ref.seq)))
        elif len(seq) >= len(ref.seq):
            rot_q, strand, shift, _ = orient_and_rotate(seq, ref.seq, k=min(k_orient, 13))
            win = rot_q[: min(len(rot_q), len(ref.seq) + 20)]
            win_id = seq_identity(ref.seq, win)
            matches = min(float(len(ref.seq)), win_id * (len(ref.seq) + len(win)) / 2.0)
            overlap_id = matches / len(ref.seq)
            global_id = (2.0 * matches) / (len(seq) + len(ref.seq))
        else:
            rot_r, strand, shift, _ = orient_and_rotate(ref.seq, seq, k=min(k_orient, 13))
            win = rot_r[: min(len(rot_r), len(seq) + 20)]
            win_id = seq_identity(seq, win)
            matches = min(float(len(seq)), win_id * (len(seq) + len(win)) / 2.0)
            overlap_id = matches / len(seq)
            global_id = (2.0 * matches) / (len(seq) + len(ref.seq))

        if len_ratio >= 0.70 and global_id >= 0.75:
            cand = (global_id, overlap_id, ref, strand, shift)
            if best_full is None or cand[:2] > best_full[:2]:
                best_full = cand
        elif (
            len_ratio < 0.70
            and min(len(seq), len(ref.seq)) >= 150
            and overlap_id >= 0.80
            and cont >= 0.40
        ):
            cand = (overlap_id, global_id, ref, strand, shift)
            if best_partial is None or cand[:2] > best_partial[:2]:
                best_partial = cand

    if best_full is not None:
        gid, _, ref, strand, shift = best_full
        return ref, gid, gid, True, strand, shift
    if best_partial is not None:
        oid, gid, ref, strand, shift = best_partial
        return ref, oid, gid, False, strand, shift
    return None


def cluster_families(
    records: List[MotifRecord], k_cluster: int, min_containment: float
) -> List[List[int]]:
    ksets = [circular_canonical_kmers(r.seq, k_cluster) for r in records]
    order = sorted(range(len(records)), key=lambda i: (records[i].bp, records[i].length), reverse=True)

    clusters: List[Tuple[Set[str], List[int]]] = []
    for idx in order:
        ks = ksets[idx]
        if not ks:
            clusters.append((ks, [idx]))
            continue
        best_c: Optional[int] = None
        best_ov = 0.0
        for ci, (rep_ks, _) in enumerate(clusters):
            if not rep_ks:
                continue
            inter = len(ks & rep_ks)
            if inter == 0:
                continue
            cont = inter / min(len(ks), len(rep_ks))
            if cont > best_ov:
                best_ov = cont
                best_c = ci
        if best_c is not None and best_ov >= min_containment:
            clusters[best_c][1].append(idx)
        else:
            clusters.append((ks, [idx]))

    clusters.sort(key=lambda c: sum(records[i].bp for i in c[1]), reverse=True)
    return [members for _, members in clusters]


def seq_identity(a: str, b: str) -> float:
    """Compute global sequence identity via anchor-chained banded alignment."""
    if a == b:
        return 1.0
    la, lb = len(a), len(b)
    if la == 0 or lb == 0:
        return 0.0
    if max(la, lb) <= 250:
        return difflib.SequenceMatcher(None, a, b, autojunk=False).ratio()

    anchor_k = 15
    b_pos_raw: Dict[str, List[int]] = collections.defaultdict(list)
    for j in range(lb - anchor_k + 1):
        b_pos_raw[b[j : j + anchor_k]].append(j)
    b_pos = {km: pos for km, pos in b_pos_raw.items() if len(pos) <= 6}

    band = max(120, abs(la - lb) + 80)
    blocks: List[Tuple[int, int, int]] = []
    last_i = 0
    last_j = 0
    i = 0
    while i <= la - anchor_k:
        positions = b_pos.get(a[i : i + anchor_k])
        if positions:
            exp_diag = (blocks[-1][1] - blocks[-1][0]) if blocks else 0
            best_j: Optional[int] = None
            best_d = band + 1
            for j in positions:
                if j < last_j:
                    continue
                d = abs((j - i) - exp_diag)
                if d <= band and d < best_d:
                    best_d = d
                    best_j = j
            if best_j is not None:
                sz = anchor_k
                while i + sz < la and best_j + sz < lb and a[i + sz] == b[best_j + sz]:
                    sz += 1
                blocks.append((i, best_j, sz))
                last_i = i + sz
                last_j = best_j + sz
                i = last_i
                continue
        i += 1

    if not blocks:
        ka = circular_canonical_kmers(a, 11)
        kb = circular_canonical_kmers(b, 11)
        return (2.0 * len(ka & kb)) / (len(ka) + len(kb)) if (ka and kb) else 0.0

    matches = 0
    pi, pj = 0, 0
    for bi, bj, bsz in blocks + [(la, lb, 0)]:
        matches += bsz
        ga = bi - pi
        gb = bj - pj
        if ga > 0 and gb > 0:
            if max(ga, gb) <= 250:
                sm = difflib.SequenceMatcher(None, a[pi:bi], b[pj:bj], autojunk=False)
                matches += sum(mb.size for mb in sm.get_matching_blocks())
            else:
                matches += sum(1 for ca, cb in zip(a[pi:bi], b[pj:bj]) if ca == cb)
        pi = bi + bsz
        pj = bj + bsz

    return (2.0 * matches) / (la + lb)


def run_pipeline(args: argparse.Namespace) -> None:
    if args.csv:
        records = load_motifs_from_csv(args.csv)
    else:
        records = load_motifs_from_srf(args.srf, args.abun)

    if not records:
        raise ValueError("No SRF motifs found in input.")

    refs = load_refs(args.refs)
    families = cluster_families(records, args.k_cluster, args.min_containment)
    total_all_bp = sum(r.bp for r in records) or 1

    out_dir = os.path.dirname(os.path.abspath(args.out_prefix))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    catalog_fa_path = f"{args.out_prefix}.catalog.fa"
    membership_path = f"{args.out_prefix}.membership.tsv"
    summary_path = f"{args.out_prefix}.catalog_summary.tsv"

    membership_rows: List[List[str]] = []
    summary_rows: List[List[str]] = []
    catalog_entries: List[Tuple[str, str, str]] = []

    for fam_num, members in enumerate(families, start=1):
        family_id = f"C{fam_num:02d}"
        # Rep is already the highest-bp member due to greedy sort order
        fam_rep_idx = max(members, key=lambda i: (records[i].bp, records[i].length))
        fam_rep_rec = records[fam_rep_idx]
        fam_ref_seq = fam_rep_rec.seq

        # Orient/rotate family representative onto matching full-length published ref if available
        fam_ref_match = match_reference(fam_ref_seq, refs, args.k_orient)
        if fam_ref_match is not None and fam_ref_match[3]:  # is_full_length
            ref_rec, _, _, _, r_strand, r_shift = fam_ref_match
            base = fam_ref_seq if r_strand == "+" else revcomp(fam_ref_seq)
            fam_ref_seq = base[r_shift:] + base[:r_shift]

        # Orient and rotate every member against family reference sequence
        rotated_by_idx: Dict[int, Tuple[str, str, int, float]] = {}
        k15_by_idx: Dict[int, Set[str]] = {}
        for idx in members:
            rec = records[idx]
            rot, strand, shift, _ = orient_and_rotate(rec.seq, fam_ref_seq, args.k_orient)
            sim_fam = seq_identity(fam_ref_seq, rot)
            rotated_by_idx[idx] = (rot, strand, shift, sim_fam)
            k15_by_idx[idx] = circular_canonical_kmers(rot, min(15, len(rot)))

        # Greedy subfamily clustering within family (ordered by bp descending)
        members_sorted = sorted(
            members, key=lambda i: (records[i].bp, records[i].length), reverse=True
        )
        subfamilies: List[Tuple[int, str, List[Tuple[int, float]]]] = []
        for idx in members_sorted:
            rot = rotated_by_idx[idx][0]
            qk = k15_by_idx[idx]
            assigned = False
            best_sf: Optional[int] = None
            best_sim = -1.0
            for si, (sf_rep_idx, sf_rep_rot, _) in enumerate(subfamilies):
                len_ratio = min(len(rot), len(sf_rep_rot)) / max(len(rot), len(sf_rep_rot))
                if len_ratio < args.subfamily_len_ratio:
                    continue
                rk = k15_by_idx[sf_rep_idx]
                if qk and rk:
                    k_cont = len(qk & rk) / min(len(qk), len(rk))
                    if k_cont < max(0.20, args.subfamily_id - 0.35):
                        continue
                sim = seq_identity(sf_rep_rot, rot)
                if sim >= args.subfamily_id and sim > best_sim:
                    best_sim = sim
                    best_sf = si
            if best_sf is not None:
                subfamilies[best_sf][2].append((idx, best_sim))
                assigned = True
            if not assigned:
                subfamilies.append((idx, rot, [(idx, 1.0)]))

        # Sort subfamilies by total bp descending
        subfamilies.sort(
            key=lambda sf: sum(records[idx].bp for idx, _ in sf[2]), reverse=True
        )

        for sf_num, (sf_rep_idx, sf_rep_rot, sf_members) in enumerate(subfamilies, start=1):
            subfamily_id = f"{family_id}_S{sf_num:02d}"
            sf_rep_rec = records[sf_rep_idx]
            catalog_id = f"{subfamily_id}_L{sf_rep_rec.length}"

            sf_ref_match = match_reference(sf_rep_rot, refs, args.k_orient)
            if sf_ref_match is not None:
                ref_obj, ref_id, _, is_full, _, _ = sf_ref_match
                ref_name = ref_obj.short_name if is_full else f"{ref_obj.short_name}(partial)"
                ref_acc = ref_obj.accession
                ref_id_str = f"{ref_id:.4f}"
            else:
                ref_name = "unannotated"
                ref_acc = "NA"
                ref_id_str = "NA"

            sf_total_bp = sum(records[idx].bp for idx, _ in sf_members)
            sf_samples = len({records[idx].sample for idx, _ in sf_members})
            in_catalog = sf_total_bp >= args.min_catalog_bp
            len_counter = collections.Counter(records[idx].length for idx, _ in sf_members)
            len_str = ",".join(f"{l}:{c}" for l, c in len_counter.most_common(8))

            summary_rows.append(
                [
                    catalog_id,
                    family_id,
                    subfamily_id,
                    sf_rep_rec.original_id,
                    str(sf_rep_rec.length),
                    str(len(sf_members)),
                    str(sf_samples),
                    str(sf_total_bp),
                    f"{100.0 * sf_total_bp / total_all_bp:.4f}",
                    ref_name,
                    ref_acc,
                    ref_id_str,
                    len_str,
                    "TRUE" if in_catalog else "FALSE",
                ]
            )

            if in_catalog:
                hdr_meta = (
                    f"family={family_id} ref={ref_name}({ref_acc}:{ref_id_str}) "
                    f"rep={sf_rep_rec.original_id} n_motifs={len(sf_members)} "
                    f"n_samples={sf_samples} total_bp={sf_total_bp}"
                )
                catalog_entries.append((catalog_id, hdr_meta, sf_rep_rot))

            for idx, sim_sf in sf_members:
                rec = records[idx]
                _, strand, shift, sim_fam = rotated_by_idx[idx]
                membership_rows.append(
                    [
                        rec.original_id,
                        rec.sample,
                        rec.contig,
                        str(rec.length),
                        family_id,
                        subfamily_id,
                        catalog_id,
                        ref_name,
                        ref_acc,
                        ref_id_str,
                        strand,
                        str(shift),
                        f"{sim_fam:.4f}",
                        f"{sim_sf:.4f}",
                        str(rec.bp),
                        f"{rec.mean_divergence:.6f}" if rec.mean_divergence >= 0 else "-1",
                    ]
                )

    with open(catalog_fa_path, "w") as fout:
        for cid, meta, seq in catalog_entries:
            fout.write(f">{cid} {meta}\n")
            for i in range(0, len(seq), 80):
                fout.write(seq[i : i + 80] + "\n")

    with open(membership_path, "w", newline="") as fout:
        writer = csv.writer(fout, delimiter="\t")
        writer.writerow(
            [
                "original_id",
                "sample",
                "contig",
                "length",
                "family_id",
                "subfamily_id",
                "catalog_id",
                "ref_name",
                "ref_acc",
                "ref_identity",
                "strand_vs_rep",
                "shift_vs_rep",
                "sim_to_family_rep",
                "sim_to_subfamily_rep",
                "bp",
                "mean_divergence",
            ]
        )
        writer.writerows(membership_rows)

    with open(summary_path, "w", newline="") as fout:
        writer = csv.writer(fout, delimiter="\t")
        writer.writerow(
            [
                "catalog_id",
                "family_id",
                "subfamily_id",
                "rep_original_id",
                "rep_length",
                "n_rows",
                "n_samples",
                "total_bp",
                "pct_all_bp",
                "ref_name",
                "ref_acc",
                "ref_identity",
                "length_counts",
                "in_catalog",
            ]
        )
        writer.writerows(summary_rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    src = parser.add_mutually_exclusive_group(required=True)
    src.add_argument("--csv", help="Gathered SRF CSV from srf_gather (e.g., All.csv)")
    src.add_argument("--srf", nargs="+", help="Per-sample .srf.fa files")
    parser.add_argument("--abun", nargs="+", help="Optional .srf-aln.abun files matching --srf")
    parser.add_argument("--refs", help="Optional curated reference FASTA of known satellites")
    parser.add_argument("--out-prefix", required=True, help="Output prefix")
    parser.add_argument("--k-cluster", type=int, default=17, help="Canonical k-mer length for family clustering")
    parser.add_argument("--min-containment", type=float, default=0.20, help="Min circular k-mer containment for family clustering")
    parser.add_argument("--k-orient", type=int, default=15, help="Forward k-mer length for strand/offset voting")
    parser.add_argument("--subfamily-id", type=float, default=0.97, help="Rotated sequence identity threshold for subfamily merging")
    parser.add_argument("--subfamily-len-ratio", type=float, default=0.90, help="Min length ratio for subfamily merging")
    parser.add_argument("--min-catalog-bp", type=int, default=500000, help="Min total bp across samples to include subfamily in catalog FASTA")
    args = parser.parse_args()
    run_pipeline(args)


if __name__ == "__main__":
    main()
