#!/usr/bin/env python3
"""Unit and integration tests for BPC assembly ingestion and PanSN formatting."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPTS_DIR = REPO_ROOT / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

import add_assembly  # noqa: E402
import format_pansn  # noqa: E402
import manifest  # noqa: E402
import stage_assembly  # noqa: E402


def _read_fasta_dict(path: Path) -> dict[str, str]:
    opener = gzip.open if path.suffix == ".gz" else Path.open
    records: dict[str, str] = {}
    hdr: str | None = None
    chunks: list[str] = []
    with opener(path, "rt", encoding="utf-8") as fin:
        for line in fin:
            s = line.strip()
            if not s:
                continue
            if s.startswith(">"):
                if hdr is not None:
                    records[hdr] = "".join(chunks)
                hdr = s[1:]
                chunks = []
            else:
                chunks.append(s)
        if hdr is not None:
            records[hdr] = "".join(chunks)
    return records


class ManifestTests(unittest.TestCase):
    def test_bom_crlf_and_legacy_schema(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / "legacy.csv"
            content = (
                "\ufeffAnimal ID , Haplotype , Sex , Species , Reference annotation , QC Fail , Graphs \r\n"
                "A1,1,male,Bos taurus,Y,,All;Bovinae\r\n"
                "A1,2,F,Bos taurus,,Y,All;Bovinae\r\n"
                "B2,0,unknown,Bos indicus,,,All;Bovinae\r\n"
            )
            p.write_bytes(content.encode("utf-8"))
            records = manifest.load_manifest(p)
            self.assertEqual(len(records), 3)
            self.assertEqual(records[0].sample, "A1_1")
            self.assertEqual(records[0].sex, "M")
            self.assertEqual(records[0].format_mode, "ragtag_preserve")
            self.assertFalse(records[0].qc_fail)
            self.assertEqual(records[1].sample, "A1_2")
            self.assertEqual(records[1].format_mode, "intact")
            self.assertTrue(records[1].qc_fail)
            self.assertEqual(records[2].sex, "U")

            ref_rec = manifest.AssemblyRecord(
                animal_id="ARS_UCD2.0", haplotype="0", sample="ARS_UCD2.0_0"
            )
            self.assertEqual(
                manifest.cactus_seqfile_name(ref_rec, ["ARS_UCD2.0"]),
                "ARS_UCD2.0",
            )
            self.assertEqual(
                manifest.cactus_seqfile_name(records[0], ["ARS_UCD2.0"]),
                "A1.1",
            )

            cfg = {
                "graphs": {
                    "All": {"exclude_qc_fail": False},
                    "Taurus_n1": {
                        "exclude_qc_fail": True,
                        "species": ["Bos taurus"],
                        "one_per_animal": True,
                    },
                }
            }
            self.assertEqual(
                manifest.resolve_graph_choices(records, cfg),
                ["All", "Taurus_n1", "Bovinae"],
            )
            self.assertEqual(
                manifest.resolve_graph_samples(records, "All", cfg),
                ["A1_1", "A1_2", "B2_0"],
            )
            self.assertEqual(
                manifest.resolve_graph_samples(records, "Taurus_n1", cfg),
                ["A1_1"],
            )
            self.assertEqual(
                manifest.resolve_graph_samples(records, "Bovinae", cfg),
                ["A1_1", "B2_0"],
            )

    def test_validation_rejects_bad_ids_duplicates_and_modes(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / "bad.tsv"
            p.write_text("animal_id\thaplotype\nBAD#ID\t1\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "invalid Animal ID"):
                manifest.load_manifest(p)

            p.write_text("animal_id\thaplotype\nA1\t1\nA1\t1\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "duplicate sample"):
                manifest.load_manifest(p)

            p.write_text("animal_id\thaplotype\tformat_mode\nA1\t1\tbogus\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "invalid format_mode"):
                manifest.load_manifest(p)

    def test_merge_qc_and_sex_mismatch_flags(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_dir = Path(tmp)
            m_path = tmp_dir / "manifest.tsv"
            m_path.write_text(
                "animal_id\thaplotype\tsex\tspecies\tformat_mode\n"
                "AN1\t1\tF\tBos taurus\tragtag_cut\n"
                "AN2\t0\tF\tBos taurus\talready_placed\n",
                encoding="utf-8",
            )
            qc_path = tmp_dir / "QC_summary.All.csv"
            qc_path.write_text(
                "sample,genome size,N contigs,NG50,BUSCOs,SNPs,InDels,autosomes covered,X covered,Y covered,MT covered\n"
                "AN1_1,2700000000,120,25000000,S:97.5%;D:1.0%;F:0.5%;I:0.1%;M:0.9%,1000,100,0.96,0.92,0.55,1.0\n"
                "AN2_0,1500000000,5000,200000,S:70.0%;D:0.5%;F:5.0%;I:1.0%;M:23.5%,500,50,0.65,0.88,0.01,1.0\n",
                encoding="utf-8",
            )
            out_tsv = tmp_dir / "merged.tsv"
            manifest.merge_manifest_and_qc(m_path, qc_path, out_tsv)
            with out_tsv.open("r", encoding="utf-8") as fin:
                rows = list(csv.DictReader(fin, delimiter="\t"))
            self.assertEqual(rows[0]["inferred_sex"], "M")
            self.assertIn("sex_mismatch", rows[0]["auto_qc_flags"])
            self.assertEqual(rows[1]["inferred_sex"], "F")
            for flag in ("genome_size_low", "low_ng50", "low_autosome_cov", "low_busco"):
                self.assertIn(flag, rows[1]["auto_qc_flags"])


class StageAssemblyTests(unittest.TestCase):
    def test_stage_plain_and_gzip_with_checksum_and_validation(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_dir = Path(tmp)
            raw_fa = tmp_dir / "input.fa"
            raw_fa.write_text(
                ">ctg1_RagTag extra description\nacgtNNNNacgt\n>ctg2\nGCRYT\n",
                encoding="utf-8",
            )
            sha = hashlib.sha256(raw_fa.read_bytes()).hexdigest()

            out_gz = tmp_dir / "staged.fa.gz"
            report_json = tmp_dir / "stage.json"
            rep = stage_assembly.stage_assembly(
                sample="TEST_1",
                source=str(raw_fa),
                output=out_gz,
                report_path=report_json,
                checksum=f"sha256:{sha}",
                threads=1,
                line_width=6,
            )

            self.assertTrue(out_gz.exists())
            self.assertTrue(Path(f"{out_gz}.fai").exists())
            self.assertTrue(Path(f"{out_gz}.gzi").exists())
            subprocess.run(["bgzip", "-t", str(out_gz)], check=True)

            seqs = _read_fasta_dict(out_gz)
            self.assertEqual(seqs, {"ctg1": "ACGTNNNNACGT", "ctg2": "GCRYT"})
            self.assertEqual(rep["n_contigs"], 2)
            self.assertEqual(rep["total_bp"], 17)
            self.assertEqual(rep["n_bases"], 4)
            self.assertEqual(rep["n_gap_runs"], 1)
            loaded_rep = json.loads(report_json.read_text(encoding="utf-8"))
            self.assertEqual(loaded_rep["source_sha256"], sha)

            # Regular gzip input + checksum mismatch failure
            raw_gz = tmp_dir / "input.fasta.gz"
            with gzip.open(raw_gz, "wb") as fout:
                fout.write(b">c1\nACGTACGT\n")
            with self.assertRaisesRegex(ValueError, "Checksum mismatch"):
                stage_assembly.stage_assembly(
                    sample="TEST_2",
                    source=str(raw_gz),
                    output=tmp_dir / "fail.fa.gz",
                    checksum="sha256:" + ("0" * 64),
                )

            # Duplicate contig rejection
            dup_fa = tmp_dir / "dup.fa"
            dup_fa.write_text(">c1\nACGT\n>c1\nTTTT\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "duplicate contig"):
                stage_assembly.stage_assembly(
                    sample="TEST_DUP",
                    source=str(dup_fa),
                    output=tmp_dir / "dup.fa.gz",
                )

            # Invalid base rejection
            bad_fa = tmp_dir / "bad_base.fa"
            bad_fa.write_text(">c1\nACGTZ\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "invalid non-IUPAC"):
                stage_assembly.stage_assembly(
                    sample="TEST_BAD",
                    source=str(bad_fa),
                    output=tmp_dir / "bad.fa.gz",
                )

    def test_stage_min_contig_len_mt_exempt_contiguity_and_empty_chrom_map(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_dir = Path(tmp)
            raw_fa = tmp_dir / "draft.fa"
            # scaf1 has two 100bp contigs separated by 12 Ns -> contig lens [100, 100], scaf len 212
            # tiny is 40bp (< 50bp -> dropped), MT is 30bp (< 50bp -> kept due to MT exemption)
            raw_fa.write_text(
                f">scaf1\n{'A' * 100}{'N' * 12}{'C' * 100}\n"
                f">tiny\n{'G' * 40}\n"
                f">MT\n{'T' * 30}\n",
                encoding="utf-8",
            )
            out_gz = tmp_dir / "staged.fa.gz"
            cmap_out = tmp_dir / "ncbi.tsv"
            rep = stage_assembly.stage_assembly(
                sample="DRAFT_1",
                source=str(raw_fa),
                output=out_gz,
                chrom_map_out=cmap_out,
                min_contig_len=50,
                keep_mt=True,
                genome_size=240,
                gap_len=10,
                threads=1,
            )
            self.assertTrue(cmap_out.exists())
            self.assertEqual(cmap_out.read_text(encoding="utf-8"), "")
            seqs = _read_fasta_dict(out_gz)
            self.assertEqual(set(seqs.keys()), {"scaf1", "MT"})
            self.assertEqual(rep["dropped_contigs"], 1)
            self.assertEqual(rep["dropped_bp"], 40)
            self.assertEqual(rep["scaffold_n50"], 212)
            self.assertEqual(rep["scaffold_ng50"], 212)
            self.assertEqual(rep["contig_n50"], 100)
            self.assertEqual(rep["contig_ng50"], 100)


class FormatPanSNTests(unittest.TestCase):
    def test_agp_revcomp_preserves_n_and_iupac_and_splits_xy_mt(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_dir = Path(tmp)
            staged_in = tmp_dir / "staged.fa"
            # ctg_neg has N and IUPAC bases R, Y, K, M, S, W, B, D, H, V
            seq_neg = "AANNNRYKMSWBDHVTT"
            staged_in.write_text(
                f">ctg_x\nACGTACGT\n>ctg_neg\n{seq_neg}\n>ctg_mt\nGGGGAAAA\n>ctg_unp\nCCCC\n",
                encoding="utf-8",
            )
            agp = tmp_dir / "ragtag.agp"
            agp.write_text(
                "X_RagTag\t1\t8\t1\tW\tctg_x\t1\t8\t+\n"
                "1_RagTag\t1\t17\t1\tW\tctg_neg\t1\t17\t-\n"
                "MT_RagTag\t1\t8\t1\tW\tctg_mt\t1\t8\t+\n"
                "scaf_unplaced\t1\t4\t1\tW\tctg_unp\t1\t4\t+\n",
                encoding="utf-8",
            )

            out_gz = tmp_dir / "AN1_1.fa.gz"
            report_json = tmp_dir / "pansn.json"
            chroms = ["1", "2", "X", "Y", "MT"]
            rep = format_pansn.run_format(
                fasta_path=staged_in,
                animal_id="AN1",
                haplotype="1",
                output_path=out_gz,
                agp_path=agp,
                schema="cactus",
                chroms=chroms,
                threads=1,
                report_path=report_json,
            )

            seqs = _read_fasta_dict(out_gz)
            headers = list(seqs.keys())
            self.assertEqual(
                headers,
                [
                    "AN1#1#1:1-17",
                    "AN1#1#X:1-8",
                    "AN1#1#MT:1-8",
                    "AN1#1#unplaced_1:1-4",
                ],
            )
            expected_rev = format_pansn.reverse_complement(seq_neg)
            self.assertEqual(len(seqs["AN1#1#1:1-17"]), len(seq_neg))
            self.assertEqual(seqs["AN1#1#1:1-17"], expected_rev)
            self.assertIn("NNN", seqs["AN1#1#1:1-17"])
            self.assertEqual(rep["chroms_present"], ["1", "X", "MT"])
            self.assertEqual(rep["chroms_missing"], ["2", "Y"])
            self.assertEqual(rep["revcomp_contigs"], 1)

            # Split across 1, 2, X, Y, MT (2 and Y are absent)
            split_dir = tmp_dir / "chromosomes"
            counts = format_pansn.run_split(
                fasta_path=out_gz,
                out_dir=split_dir,
                sample="AN1_1",
                chroms=chroms,
            )
            self.assertEqual(counts, {"1": 1, "2": 0, "X": 1, "Y": 0, "MT": 1})
            self.assertEqual(_read_fasta_dict(split_dir / "AN1_1.X.fa.gz"), {"AN1#1#X:1-8": "ACGTACGT"})
            self.assertEqual(_read_fasta_dict(split_dir / "AN1_1.MT.fa.gz"), {"AN1#1#MT:1-8": "GGGGAAAA"})

            # Verify concatenating an empty split + non-empty split indexes cleanly with samtools faidx
            gathered = tmp_dir / "gathered_Y.fa.gz"
            with gathered.open("wb") as fout:
                fout.write((split_dir / "AN1_1.Y.fa.gz").read_bytes())
                fout.write((split_dir / "AN1_1.X.fa.gz").read_bytes())
            subprocess.run(["samtools", "faidx", str(gathered)], check=True)
            self.assertTrue(Path(f"{gathered}.fai").exists())

    def test_already_placed_with_ncbi_assembly_report(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_dir = Path(tmp)
            raw_fa = tmp_dir / "ncbi.fa"
            raw_fa.write_text(
                ">NC_037328.1\nACGT\n>NC_006853.1\nTTTT\n>NW_020190.1\nGGGG\n",
                encoding="utf-8",
            )
            staged_gz = tmp_dir / "ncbi.staged.fa.gz"
            stage_assembly.stage_assembly(
                sample="REF_0", source=str(raw_fa), output=staged_gz, threads=1
            )

            report_txt = tmp_dir / "assembly_report.txt"
            report_txt.write_text(
                "# Sequence-Name\tSequence-Role\tAssigned-Molecule\tAssigned-Molecule-Location/Type\tGenBank-Accn\tRelationship\tRefSeq-Accn\tAssembly-Unit\tSequence-Length\tUCSC-style-name\n"
                "1\tassembled-molecule\t1\tChromosome\tCM008168.2\t=\tNC_037328.1\tPrimary Assembly\t4\tchr1\n"
                "MT\tassembled-molecule\tMT\tMitochondrion\tAY526085.1\t=\tNC_006853.1\tnon-nuclear\t4\tchrM\n"
                "unplaced1\tunplaced-scaffold\tna\tna\tNKLS020001\t=\tNW_020190.1\tPrimary Assembly\t4\tna\n",
                encoding="utf-8",
            )

            out_gz = tmp_dir / "REF_0.fa.gz"
            format_pansn.run_format(
                fasta_path=staged_gz,
                animal_id="REF",
                haplotype="0",
                output_path=out_gz,
                already_placed=True,
                chrom_map_path=report_txt,
                schema="vg",
                chroms=["1", "X", "MT"],
                threads=1,
            )
            seqs = _read_fasta_dict(out_gz)
            self.assertEqual(
                list(seqs.keys()),
                ["REF#0#1[0-4]", "REF#0#MT[0-4]", "REF#0#unplaced_1[0-4]"],
            )

    def test_vg_schema_agp_uses_zero_based_half_open_coords(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_dir = Path(tmp)
            fa = tmp_dir / "in.fa"
            fa.write_text(">c1\nACGTACGT\n>c2\nGGGGAA\n", encoding="utf-8")
            agp = tmp_dir / "in.agp"
            agp.write_text(
                "1_RagTag\t1\t8\t1\tW\tc1\t1\t8\t+\n"
                "1_RagTag\t109\t114\t3\tW\tc2\t1\t6\t+\n",
                encoding="utf-8",
            )
            out_gz = tmp_dir / "VG_1.fa.gz"
            format_pansn.run_format(
                fasta_path=fa,
                animal_id="VG",
                haplotype="1",
                output_path=out_gz,
                agp_path=agp,
                schema="vg",
                chroms=["1"],
                threads=1,
            )
            seqs = _read_fasta_dict(out_gz)
            self.assertEqual(
                list(seqs.keys()),
                ["VG#1#1[0-8]", "VG#1#1[108-114]"],
            )
            for hdr, seq in seqs.items():
                coords = hdr.split("[", 1)[1].rstrip("]")
                s, e = map(int, coords.split("-"))
                self.assertEqual(e - s, len(seq))

    def test_intact_format_optional_chrom_map_and_strip_pansn(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_dir = Path(tmp)
            fa = tmp_dir / "scaffolds.fa"
            fa.write_text(
                ">scaf_short\nACGT\n>chr2\nGGGGAAAA\n>chrM\nTTTT\n>scaf_long\nCCCCCCCCCC\n",
                encoding="utf-8",
            )
            empty_cmap = tmp_dir / "empty.tsv"
            empty_cmap.write_text("", encoding="utf-8")
            out_gz = tmp_dir / "INT_1.fa.gz"
            format_pansn.run_format(
                fasta_path=fa,
                animal_id="INT",
                haplotype="1",
                output_path=out_gz,
                intact=True,
                optional_chrom_map_path=empty_cmap,
                chroms=["1", "2", "X", "Y", "MT"],
                threads=1,
            )
            seqs = _read_fasta_dict(out_gz)
            self.assertEqual(
                list(seqs.keys()),
                [
                    "INT#1#2",
                    "INT#1#MT",
                    "INT#1#scaf_long",
                    "INT#1#scaf_short",
                ],
            )

            stripped_gz = tmp_dir / "INT_1.cactus.fa.gz"
            n_seqs = format_pansn.run_strip_pansn(out_gz, stripped_gz, threads=1)
            self.assertEqual(n_seqs, 4)
            stripped_seqs = _read_fasta_dict(stripped_gz)
            self.assertEqual(
                list(stripped_seqs.keys()),
                ["2", "MT", "scaf_long", "scaf_short"],
            )
            self.assertTrue(Path(f"{stripped_gz}.fai").exists())


class AddAssemblyTests(unittest.TestCase):
    def test_hex_id_generator_is_excel_safe_and_unique(self) -> None:
        existing = {"85D7A68E", "aaaaaaaa"}
        seen: set[str] = set()
        for _ in range(50):
            hid = manifest.generate_unique_hex_id(existing | seen)
            self.assertRegex(hid, r"^[0-9A-F]{8}$")
            self.assertTrue(any(c in "ABCDF" for c in hid))
            self.assertNotIn(hid.upper(), {x.upper() for x in existing | seen})
            seen.add(hid)

    def test_auto_hex_diploid_pairing_and_id_map(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_dir = Path(tmp)
            h1 = tmp_dir / "hap1.fa"
            h2 = tmp_dir / "hap2.fa"
            h3 = tmp_dir / "hap2_retry.fa"
            h1.write_text(">c1\nACGTACGT\n", encoding="utf-8")
            h2.write_text(">c1\nTTTTGGGG\n", encoding="utf-8")
            h3.write_text(">c1\nCCCCAAAA\n", encoding="utf-8")
            manifest_tsv = tmp_dir / "assemblies.tsv"
            id_map_tsv = tmp_dir / "private_id_map.tsv"

            rc = add_assembly.main(
                [
                    "--manifest",
                    str(manifest_tsv),
                    "--hap1",
                    str(h1),
                    "--hap2",
                    str(h2),
                    "--external-id",
                    "LAB_COW_42",
                    "--id-map",
                    str(id_map_tsv),
                    "--sex",
                    "M",
                    "--species",
                    "Bos taurus",
                ]
            )
            self.assertEqual(rc, 0)
            recs = manifest.load_manifest(manifest_tsv)
            self.assertEqual(len(recs), 2)
            self.assertEqual(recs[0].animal_id, recs[1].animal_id)
            self.assertRegex(recs[0].animal_id, r"^[0-9A-F]{8}$")
            self.assertEqual([r.haplotype for r in recs], ["1", "2"])
            self.assertNotIn("LAB_COW_42", manifest_tsv.read_text(encoding="utf-8"))
            self.assertIn("LAB_COW_42", id_map_tsv.read_text(encoding="utf-8"))

            # Batch mode reusing external_id for another haplotype/sample & minting new ID for a 2nd animal
            batch_tsv = tmp_dir / "batch.tsv"
            batch_tsv.write_text(
                "external_id\thaplotype\tsex\tspecies\tsource\n"
                f"LAB_COW_42\t0\tM\tBos taurus\t{h3}\n"
                f"LAB_COW_99\t1\tF\tBos taurus\t{h1}\n"
                f"LAB_COW_99\t2\tF\tBos taurus\t{h2}\n",
                encoding="utf-8",
            )
            rc2 = add_assembly.main(
                [
                    "--manifest",
                    str(manifest_tsv),
                    "--batch",
                    str(batch_tsv),
                    "--id-map",
                    str(id_map_tsv),
                ]
            )
            self.assertEqual(rc2, 0)
            recs2 = manifest.load_manifest(manifest_tsv)
            self.assertEqual(len(recs2), 5)
            self.assertEqual(recs2[2].animal_id, recs[0].animal_id)
            self.assertEqual(recs2[2].haplotype, "0")
            self.assertEqual(recs2[3].animal_id, recs2[4].animal_id)
            self.assertNotEqual(recs2[3].animal_id, recs[0].animal_id)
            self.assertRegex(recs2[3].animal_id, r"^[0-9A-F]{8}$")

    def test_add_single_and_batch_and_duplicate_rejection(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_dir = Path(tmp)
            fa1 = tmp_dir / "s1.fa"
            fa1.write_text(">c1\nACGTACGT\n", encoding="utf-8")
            fa2 = tmp_dir / "s2.fa"
            fa2.write_text(">c1\nTTTTGGGG\n", encoding="utf-8")
            manifest_tsv = tmp_dir / "assemblies.tsv"

            rc = add_assembly.main(
                [
                    "--manifest",
                    str(manifest_tsv),
                    "--source",
                    str(fa1),
                    "--animal-id",
                    "BOV1",
                    "--haplotype",
                    "1",
                    "--sex",
                    "M",
                    "--species",
                    "Bos taurus",
                ]
            )
            self.assertEqual(rc, 0)
            recs = manifest.load_manifest(manifest_tsv)
            self.assertEqual(len(recs), 1)
            self.assertTrue(recs[0].checksum.startswith("sha256:"))

            batch_tsv = tmp_dir / "batch.tsv"
            batch_tsv.write_text(
                f"animal_id\thaplotype\tsex\tspecies\tsource\tformat_mode\n"
                f"BOV2\t0\tF\tBos indicus\t{fa2}\talready_placed\n",
                encoding="utf-8",
            )
            rc2 = add_assembly.main(
                ["--manifest", str(manifest_tsv), "--batch", str(batch_tsv)]
            )
            self.assertEqual(rc2, 0)
            recs2 = manifest.load_manifest(manifest_tsv)
            self.assertEqual([r.sample for r in recs2], ["BOV1_1", "BOV2_0"])
            self.assertEqual(recs2[1].format_mode, "already_placed")

            with self.assertRaisesRegex(ValueError, "duplicate sample"):
                add_assembly.add_assemblies(
                    manifest_tsv,
                    [{"animal_id": "BOV1", "haplotype": "1", "source": str(fa1)}],
                )


class TestFixturesAndSnakefilesTests(unittest.TestCase):
    def test_repo_test_fixtures_stage_format_and_split(self) -> None:
        test_dir = REPO_ROOT / ".test"
        recs = manifest.load_manifest(test_dir / "assemblies.tsv")
        self.assertEqual([r.sample for r in recs], ["AAAA0000_1", "AAAA0000_2", "BBBB1111_0"])
        self.assertTrue(all(not r.qc_fail for r in recs))
        self.assertEqual(
            [r.format_mode for r in recs],
            ["intact", "ragtag_preserve", "already_placed"],
        )
        self.assertEqual(
            manifest.resolve_graph_choices(recs),
            ["All", "Bovinae", "BosTaurus_n1"],
        )

        # Validate legacy CSV still parses cleanly too
        legacy = manifest.load_manifest(test_dir / "BPC_metadata.csv")
        self.assertEqual(len(legacy), 3)

        # Stage and PanSN-format the already_placed test fixture BBBB1111_0 end-to-end
        with tempfile.TemporaryDirectory() as tmp:
            tmp_dir = Path(tmp)
            b_rec = recs[2]
            staged = tmp_dir / "BBBB1111_0.staged.fa.gz"
            stage_assembly.stage_assembly(
                sample=b_rec.sample,
                source=str(test_dir / b_rec.source),
                output=staged,
                threads=1,
            )
            formatted = tmp_dir / "BBBB1111_0.fa.gz"
            format_pansn.run_format(
                fasta_path=staged,
                animal_id=b_rec.animal_id,
                haplotype=b_rec.haplotype,
                output_path=formatted,
                already_placed=True,
                chrom_map_path=test_dir / b_rec.chrom_map,
                chroms=["1", "2", "X", "Y", "MT"],
                threads=1,
            )
            seqs = _read_fasta_dict(formatted)
            self.assertIn("BBBB1111#0#1:1-120", seqs)
            self.assertIn("BBBB1111#0#X:1-120", seqs)
            self.assertIn("BBBB1111#0#MT:1-120", seqs)


if __name__ == "__main__":
    unittest.main()
