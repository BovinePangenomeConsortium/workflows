#!/usr/bin/env python3
"""Stage, validate, filter, checksum-verify, and bgzip-index raw genome assemblies.

Supports local FASTA files (.fa/.fasta/.fna, plain, gzip/bgzip, or zstd),
HTTP/HTTPS/FTP URLs, NCBI GenBank/RefSeq accessions (GCA_*/GCF_* via `datasets`),
and AGC archives (`agc://archive.agc#sample`).
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import urllib.request
import zipfile
from pathlib import Path

NCBI_ACC_RE = re.compile(r"^GC[AF]_\d{9}(?:\.\d+)?$")
URL_RE = re.compile(r"^(?:https?|ftp)://", re.IGNORECASE)
AGC_RE = re.compile(r"^agc://(.+?)#(.+)$")
MT_ID_RE = re.compile(r"^(?:chr)?(?:mt|m)$", re.IGNORECASE)
VALID_IUPAC_BYTES = frozenset(b"ACGTURYKMSWBDHVN")
_TRANSLATE_UPPER = bytes.maketrans(b"abcdefghijklmnopqrstuvwxyz", b"ABCDEFGHIJKLMNOPQRSTUVWXYZ")
_N_SPLIT_RE = re.compile(rb"N+|[^N]+")


def _compute_file_hash(path: Path, alg: str) -> str:
    h = hashlib.new(alg)
    with path.open("rb") as fin:
        while chunk := fin.read(1024 * 1024):
            h.update(chunk)
    return h.hexdigest()


def verify_checksum(path: Path, expected_spec: str) -> str:
    """Verify checksum (`sha256:<hex>`, `md5:<hex>`, or raw hex) and return sha256."""
    actual_sha256 = _compute_file_hash(path, "sha256")
    spec = expected_spec.strip()
    if not spec:
        return actual_sha256

    if ":" in spec:
        alg, expected_hex = spec.split(":", 1)
        alg = alg.strip().lower()
    else:
        expected_hex = spec
        if len(expected_hex) == 64:
            alg = "sha256"
        elif len(expected_hex) == 32:
            alg = "md5"
        else:
            raise ValueError(
                f"Unrecognized checksum format {expected_spec!r}; use 'sha256:<hex>' or 'md5:<hex>'"
            )

    expected_hex = expected_hex.strip().lower()
    if alg not in ("sha256", "md5"):
        raise ValueError(f"Unsupported checksum algorithm {alg!r}; expected sha256 or md5")

    actual = actual_sha256 if alg == "sha256" else _compute_file_hash(path, "md5")
    if actual.lower() != expected_hex:
        raise ValueError(
            f"Checksum mismatch for {path} ({alg}): expected {expected_hex}, got {actual}"
        )
    return actual_sha256


def _extract_ncbi_chrom_map(seq_report_path: Path, out_path: Path) -> None:
    """Parse NCBI sequence_report.jsonl into 2-column accession -> canonical chrom TSV."""
    rows: list[tuple[str, str]] = []
    with seq_report_path.open("r", encoding="utf-8") as fin:
        for line in fin:
            line = line.strip()
            if not line:
                continue
            rec = json.loads(line)
            role = (rec.get("role") or "").lower()
            chr_name = (rec.get("chrName") or "").strip()
            if not chr_name:
                continue
            if role == "assembled-molecule" or rec.get("assignedMoleculeLocationType") == "Mitochondrion":
                canon = "MT" if chr_name.upper() in ("MT", "M", "CHRM", "CHRMT") else chr_name
                for key in ("genbankAccession", "refseqAccession", "ucscStyleName", "sequenceName"):
                    val = (rec.get(key) or "").strip()
                    if val:
                        rows.append((val, canon))
    if rows:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with out_path.open("w", encoding="utf-8") as fout:
            for acc, canon in rows:
                fout.write(f"{acc}\t{canon}\n")


def _resolve_source_to_local(
    source: str,
    tmpdir: Path,
    chrom_map_out: Path | None = None,
) -> Path:
    src = source.strip()
    if not src:
        raise ValueError("Empty --source provided")

    if NCBI_ACC_RE.match(src):
        if not shutil.which("datasets"):
            raise RuntimeError(
                f"Cannot download NCBI accession {src}: 'datasets' CLI not found on PATH"
            )
        pkg_zip = tmpdir / f"{src}.zip"
        subprocess.run(
            [
                "datasets",
                "download",
                "genome",
                "accession",
                src,
                "--include",
                "genome,seq-report",
                "--filename",
                str(pkg_zip),
            ],
            check=True,
        )
        extract_dir = tmpdir / "ncbi_pkg"
        with zipfile.ZipFile(pkg_zip, "r") as zf:
            zf.extractall(extract_dir)
        fna_files = sorted(extract_dir.glob("ncbi_dataset/data/*/*_genomic.fna"))
        if not fna_files:
            fna_files = sorted(extract_dir.rglob("*.fna"))
        if not fna_files:
            raise FileNotFoundError(f"No genomic FASTA found in NCBI package for {src}")
        if chrom_map_out is not None:
            seq_reports = sorted(extract_dir.rglob("sequence_report.jsonl"))
            if seq_reports:
                _extract_ncbi_chrom_map(seq_reports[0], chrom_map_out)
        return fna_files[0]

    if URL_RE.match(src):
        suffix = ".fa.gz" if src.lower().endswith(".gz") else ".fa"
        dest = tmpdir / f"downloaded{suffix}"
        with urllib.request.urlopen(src) as resp, dest.open("wb") as fout:
            shutil.copyfileobj(resp, fout, length=1024 * 1024)
        return dest

    m_agc = AGC_RE.match(src)
    if m_agc:
        archive, sample_name = m_agc.group(1), m_agc.group(2)
        if not shutil.which("agc"):
            raise RuntimeError(f"Cannot extract {src}: 'agc' CLI not found on PATH")
        dest = tmpdir / "agc_extracted.fa"
        with dest.open("wb") as fout:
            subprocess.run(["agc", "getset", archive, sample_name], stdout=fout, check=True)
        return dest

    local_path = Path(src)
    if not local_path.exists():
        raise FileNotFoundError(f"Local source assembly not found: {local_path}")
    return local_path


def _open_fasta_binary(path: Path) -> tuple[io.BufferedIOBase, subprocess.Popen[bytes] | None]:
    with path.open("rb") as probe:
        magic = probe.read(4)
    if magic.startswith(b"\x1f\x8b"):
        return gzip.open(path, "rb"), None  # type: ignore[return-value]
    if magic == b"\x28\xb5\x2f\xfd":
        if not shutil.which("zstd"):
            raise RuntimeError(f"Cannot decompress {path}: 'zstd' CLI not found on PATH")
        proc = subprocess.Popen(["zstd", "-dc", str(path)], stdout=subprocess.PIPE)
        assert proc.stdout is not None
        return proc.stdout, proc  # type: ignore[return-value]
    return path.open("rb"), None


def _compute_nx(lengths: list[int], target_bp: int) -> int:
    if not lengths or target_bp <= 0:
        return 0
    acc = 0
    for length in sorted(lengths, reverse=True):
        acc += length
        if acc >= target_bp:
            return length
    return 0


def stage_assembly(
    sample: str,
    source: str,
    output: Path,
    report_path: Path | None = None,
    chrom_map_out: Path | None = None,
    checksum: str = "",
    threads: int = 2,
    line_width: int = 80,
    min_contig_len: int = 0,
    keep_mt: bool = True,
    genome_size: int = 3_000_000_000,
    gap_len: int = 10,
) -> dict[str, object]:
    """Resolve, verify, validate, filter, normalize, bgzip-compress, and faidx-index."""
    if not shutil.which("bgzip"):
        raise RuntimeError("'bgzip' not found on PATH")
    if not shutil.which("samtools"):
        raise RuntimeError("'samtools' not found on PATH")
    if line_width < 1:
        raise ValueError("--line-width must be >= 1")
    if min_contig_len < 0:
        raise ValueError("--min-contig-len must be >= 0")
    if gap_len < 1:
        raise ValueError("--gap-len must be >= 1")

    output.parent.mkdir(parents=True, exist_ok=True)
    if chrom_map_out is not None:
        chrom_map_out.parent.mkdir(parents=True, exist_ok=True)
        chrom_map_out.write_text("", encoding="utf-8")

    with tempfile.TemporaryDirectory(prefix=f"stage_{sample}_") as tmp:
        tmpdir = Path(tmp)
        local_src = _resolve_source_to_local(source, tmpdir, chrom_map_out=chrom_map_out)
        source_sha256 = verify_checksum(local_src, checksum)

        tmp_out = output.parent / f".{output.name}.tmp.{os.getpid()}"
        fin, zstd_proc = _open_fasta_binary(local_src)

        n_contigs = 0
        total_bp = 0
        n_bases = 0
        gc_bases = 0
        n_gap_runs = 0
        min_len: int | None = None
        max_len = 0
        dropped_contigs = 0
        dropped_bp = 0
        scaffold_lens: list[int] = []
        contig_lens: list[int] = []

        seen_ids: set[str] = set()
        current_id: str | None = None
        current_len = 0
        current_n = 0
        current_gc = 0
        current_gap_runs = 0
        in_gap = False
        pending_n = 0
        piece_len = 0
        cur_contig_pieces: list[int] = []
        cur_exempt = True
        cur_flushed = False
        contig_buf = bytearray()
        wrap_buf = bytearray()

        bgzip_cmd = ["bgzip", "-@", str(max(1, threads)), "-c"]
        try:
            with tmp_out.open("wb") as raw_out:
                bgzip_proc = subprocess.Popen(
                    bgzip_cmd, stdin=subprocess.PIPE, stdout=raw_out
                )
                assert bgzip_proc.stdin is not None
                try:
                    def _flush_wrap(final: bool = False) -> None:
                        idx = 0
                        buf_len = len(wrap_buf)
                        while buf_len - idx >= line_width:
                            bgzip_proc.stdin.write(wrap_buf[idx : idx + line_width])
                            bgzip_proc.stdin.write(b"\n")
                            idx += line_width
                        if idx:
                            del wrap_buf[:idx]
                        if final and wrap_buf:
                            bgzip_proc.stdin.write(wrap_buf)
                            bgzip_proc.stdin.write(b"\n")
                            wrap_buf.clear()

                    def _begin_output_contig() -> None:
                        nonlocal cur_flushed
                        assert current_id is not None
                        cur_flushed = True
                        bgzip_proc.stdin.write(f">{current_id}\n".encode("utf-8"))
                        if contig_buf:
                            wrap_buf.extend(contig_buf)
                            contig_buf.clear()
                            _flush_wrap(final=False)

                    def _update_piece_stats(seq: bytes, n_count: int) -> None:
                        nonlocal in_gap, current_gap_runs, pending_n, piece_len
                        seq_len = len(seq)
                        if n_count == 0:
                            in_gap = False
                            if pending_n > 0:
                                if pending_n >= gap_len:
                                    if piece_len > 0:
                                        cur_contig_pieces.append(piece_len)
                                    piece_len = seq_len
                                else:
                                    piece_len += pending_n + seq_len
                                pending_n = 0
                            else:
                                piece_len += seq_len
                        elif n_count == seq_len:
                            if not in_gap:
                                current_gap_runs += 1
                                in_gap = True
                            pending_n += seq_len
                        else:
                            for m in _N_SPLIT_RE.finditer(seq):
                                tok = m.group(0)
                                tlen = len(tok)
                                if tok[0] == 78:  # 'N'
                                    if not in_gap:
                                        current_gap_runs += 1
                                        in_gap = True
                                    pending_n += tlen
                                else:
                                    in_gap = False
                                    if pending_n > 0:
                                        if pending_n >= gap_len:
                                            if piece_len > 0:
                                                cur_contig_pieces.append(piece_len)
                                            piece_len = tlen
                                        else:
                                            piece_len += pending_n + tlen
                                        pending_n = 0
                                    else:
                                        piece_len += tlen

                    def _finish_contig() -> None:
                        nonlocal n_contigs, total_bp, n_bases, gc_bases, n_gap_runs
                        nonlocal min_len, max_len, dropped_contigs, dropped_bp
                        nonlocal piece_len, pending_n
                        if current_id is None:
                            return
                        if current_len == 0:
                            raise ValueError(
                                f"Assembly {sample}: contig '{current_id}' has 0 bases"
                            )
                        if not cur_flushed:
                            if cur_exempt or current_len >= min_contig_len:
                                _begin_output_contig()
                            else:
                                dropped_contigs += 1
                                dropped_bp += current_len
                                contig_buf.clear()
                                cur_contig_pieces.clear()
                                piece_len = 0
                                pending_n = 0
                                return

                        _flush_wrap(final=True)
                        if pending_n < gap_len:
                            piece_len += pending_n
                        if piece_len > 0:
                            cur_contig_pieces.append(piece_len)
                        contig_lens.extend(cur_contig_pieces)
                        cur_contig_pieces.clear()
                        piece_len = 0
                        pending_n = 0

                        n_contigs += 1
                        total_bp += current_len
                        n_bases += current_n
                        gc_bases += current_gc
                        n_gap_runs += current_gap_runs
                        scaffold_lens.append(current_len)
                        min_len = current_len if min_len is None else min(min_len, current_len)
                        if current_len > max_len:
                            max_len = current_len

                    for raw_line in fin:
                        line = raw_line.strip()
                        if not line:
                            continue
                        if line.startswith(b">"):
                            _finish_contig()
                            header_bytes = line[1:].split(None, 1)[0]
                            contig_id = header_bytes.decode("utf-8", errors="strict")
                            if contig_id.endswith("_RagTag"):
                                contig_id = contig_id[: -len("_RagTag")]
                            if not contig_id:
                                raise ValueError(f"Assembly {sample}: empty FASTA header encountered")
                            if "#" in contig_id:
                                raise ValueError(
                                    f"Assembly {sample}: contig ID '{contig_id}' contains '#' "
                                    "which conflicts with PanSN delimiters"
                                )
                            if contig_id in seen_ids:
                                raise ValueError(
                                    f"Assembly {sample}: duplicate contig ID '{contig_id}'"
                                )
                            seen_ids.add(contig_id)
                            current_id = contig_id
                            current_len = 0
                            current_n = 0
                            current_gc = 0
                            current_gap_runs = 0
                            in_gap = False
                            pending_n = 0
                            piece_len = 0
                            cur_contig_pieces.clear()
                            contig_buf.clear()
                            cur_exempt = (min_contig_len == 0) or (
                                keep_mt and bool(MT_ID_RE.match(contig_id))
                            )
                            cur_flushed = False
                            if cur_exempt:
                                _begin_output_contig()
                        else:
                            if current_id is None:
                                raise ValueError(
                                    f"Assembly {sample}: sequence data found before first '>' header"
                                )
                            seq = line.translate(_TRANSLATE_UPPER, b" \t\r")
                            if not seq:
                                continue
                            invalid = set(seq) - VALID_IUPAC_BYTES
                            if invalid:
                                bad_chars = "".join(chr(c) for c in sorted(invalid))
                                raise ValueError(
                                    f"Assembly {sample} contig '{current_id}': "
                                    f"invalid non-IUPAC nucleotide character(s): {bad_chars!r}"
                                )
                            seq_len = len(seq)
                            n_count = seq.count(b"N")
                            current_len += seq_len
                            current_n += n_count
                            current_gc += seq.count(b"G") + seq.count(b"C")
                            _update_piece_stats(seq, n_count)

                            if cur_flushed:
                                wrap_buf.extend(seq)
                                _flush_wrap(final=False)
                            else:
                                contig_buf.extend(seq)
                                if cur_exempt or current_len >= min_contig_len:
                                    _begin_output_contig()

                    _finish_contig()
                    if n_contigs == 0:
                        raise ValueError(
                            f"Assembly {sample}: no contigs retained after staging "
                            f"(dropped_contigs={dropped_contigs}, min_contig_len={min_contig_len})"
                        )

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

            if zstd_proc is not None:
                zrc = zstd_proc.wait()
                if zrc != 0:
                    raise subprocess.CalledProcessError(zrc, ["zstd", "-dc"])

            os.replace(tmp_out, output)
        except Exception:
            if tmp_out.exists():
                tmp_out.unlink()
            raise
        finally:
            fin.close()

    subprocess.run(["samtools", "faidx", str(output)], check=True)

    acgt_bases = total_bp - n_bases
    gc_frac = (gc_bases / acgt_bases) if acgt_bases > 0 else 0.0
    scaffold_n50 = _compute_nx(scaffold_lens, (total_bp + 1) // 2)
    scaffold_ng50 = _compute_nx(scaffold_lens, (genome_size + 1) // 2)
    contig_total_bp = sum(contig_lens)
    contig_n50 = _compute_nx(contig_lens, (contig_total_bp + 1) // 2)
    contig_ng50 = _compute_nx(contig_lens, (genome_size + 1) // 2)

    report: dict[str, object] = {
        "sample": sample,
        "source": source,
        "output": str(output),
        "source_sha256": source_sha256,
        "n_contigs": n_contigs,
        "total_bp": total_bp,
        "n_bases": n_bases,
        "n_gap_runs": n_gap_runs,
        "min_len": min_len or 0,
        "max_len": max_len,
        "gc_frac": round(gc_frac, 6),
        "dropped_contigs": dropped_contigs,
        "dropped_bp": dropped_bp,
        "scaffold_n50": scaffold_n50,
        "scaffold_ng50": scaffold_ng50,
        "contig_n50": contig_n50,
        "contig_ng50": contig_ng50,
    }

    if report_path is not None:
        report_path.parent.mkdir(parents=True, exist_ok=True)
        report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    return report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Stage, validate, filter, and bgzip-index a raw assembly FASTA."
    )
    parser.add_argument("--sample", required=True, help="Sample identifier (AnimalID_Haplotype).")
    parser.add_argument(
        "--source",
        required=True,
        help="Local path, URL, NCBI accession (GCA_*/GCF_*), or agc://archive#sample.",
    )
    parser.add_argument("--output", required=True, type=Path, help="Output .fa.gz path.")
    parser.add_argument("--report", type=Path, help="Optional JSON staging report path.")
    parser.add_argument(
        "--chrom-map-out",
        type=Path,
        help="Optional output TSV for NCBI sequence_report chromosome mapping (created empty if non-NCBI).",
    )
    parser.add_argument(
        "--checksum",
        default="",
        help="Optional expected checksum (sha256:<hex>, md5:<hex>, or raw hex).",
    )
    parser.add_argument("--threads", type=int, default=2, help="bgzip compression threads.")
    parser.add_argument("--line-width", type=int, default=80, help="FASTA wrap width.")
    parser.add_argument(
        "--min-contig-len",
        type=int,
        default=0,
        help="Drop contigs shorter than N bp (default: 0 = keep all).",
    )
    parser.add_argument(
        "--keep-mt",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Exempt MT/chrM contigs from --min-contig-len filtering (default: true).",
    )
    parser.add_argument(
        "--genome-size",
        type=int,
        default=3_000_000_000,
        help="Target genome size for scaffold/contig NG50 (default: 3000000000).",
    )
    parser.add_argument(
        "--gap-len",
        type=int,
        default=10,
        help="Minimum N-run length separating contigs for contig N50/NG50 (default: 10).",
    )

    args = parser.parse_args(argv)
    stage_assembly(
        sample=args.sample,
        source=args.source,
        output=args.output,
        report_path=args.report,
        chrom_map_out=args.chrom_map_out,
        checksum=args.checksum,
        threads=args.threads,
        line_width=args.line_width,
        min_contig_len=args.min_contig_len,
        keep_mt=args.keep_mt,
        genome_size=args.genome_size,
        gap_len=args.gap_len,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
