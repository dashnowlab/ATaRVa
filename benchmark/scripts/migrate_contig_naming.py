#!/usr/bin/env python3
"""One-off migration: rename pseudo-contigs from "<chrom>_<start>_<end>" to
"<chrom>s<start>e<end>" across all already-built files in benchmark/data/, in place --
no need to re-fetch anything from remote. Not part of the regular pipeline; run once if
you're carrying forward a bundle built before this naming change (see window_name() in
02_build_subset.py for why: ATaRVa v0.7.1's own region-file validation requires the CHROM
column to be str.isalnum(), which rejects underscores).
"""
import glob
import gzip
import os
import re
import subprocess
import sys


def run(cmd, **kwargs):
    print(f"> {' '.join(cmd)}", flush=True)
    return subprocess.run(cmd, check=True, **kwargs)


def rename(name):
    return re.sub(r"^(chr[0-9XYM]+)_(\d+)_(\d+)$", r"\1s\2e\3", name)


def migrate_plain_text(path, rename_first_column=True):
    """Rewrite a bgzipped, tab-delimited file's first column in place."""
    tmp_path = path[:-3]
    with gzip.open(path, "rt") as src, open(tmp_path, "w") as dst:
        for line in src:
            if line.startswith("#"):
                dst.write(re.sub(r"(chr[0-9XYM]+)_(\d+)_(\d+)", r"\1s\2e\3", line))
                continue
            fields = line.rstrip("\n").split("\t")
            fields[0] = rename(fields[0])
            dst.write("\t".join(fields) + "\n")
    run(["bash", "-c", f"bgzip -f {tmp_path}"])


def migrate_fasta(ref_fasta_path):
    tmp_path = ref_fasta_path + ".renamed"
    with open(ref_fasta_path) as src, open(tmp_path, "w") as dst:
        for line in src:
            if line.startswith(">"):
                dst.write(">" + rename(line[1:].strip()) + "\n")
            else:
                dst.write(line)
    os.replace(tmp_path, ref_fasta_path)
    for ext in (".fai",):
        if os.path.exists(ref_fasta_path + ext):
            os.remove(ref_fasta_path + ext)
    run(["samtools", "faidx", ref_fasta_path])


def migrate_bam(bam_path):
    header_sam = subprocess.run(["samtools", "view", "-H", bam_path],
                                check=True, capture_output=True, text=True).stdout
    new_header_lines = []
    for line in header_sam.splitlines():
        if line.startswith("@SQ"):
            new_header_lines.append(re.sub(r"(chr[0-9XYM]+)_(\d+)_(\d+)", r"\1s\2e\3", line))
        else:
            new_header_lines.append(line)
    new_header_path = bam_path + ".newheader.sam"
    with open(new_header_path, "w") as f:
        f.write("\n".join(new_header_lines) + "\n")

    tmp_bam_path = bam_path + ".renamed.bam"
    run(["bash", "-c", f"samtools reheader {new_header_path} {bam_path} > {tmp_bam_path}"])
    os.replace(tmp_bam_path, bam_path)
    os.remove(new_header_path)
    if os.path.exists(bam_path + ".bai"):
        os.remove(bam_path + ".bai")
    run(["samtools", "index", bam_path])


def main():
    data_dir = sys.argv[1] if len(sys.argv) > 1 else "benchmark/data"

    print("=== windows.bed ===")
    windows_path = os.path.join(data_dir, "windows.bed")
    tmp_path = windows_path + ".tmp"
    with open(windows_path) as src, open(tmp_path, "w") as dst:
        for line in src:
            chrom, start, end = line.rstrip("\n").split("\t")
            dst.write(f"{rename(chrom)}\t{start}\t{end}\n")
    os.replace(tmp_path, windows_path)

    print("=== ref.fa ===")
    migrate_fasta(os.path.join(data_dir, "ref.fa"))

    for name in ("catalog.bed.gz", "confident.bed.gz", "scoring_regions.bed.gz"):
        path = os.path.join(data_dir, name)
        if os.path.exists(path):
            print(f"=== {name} ===")
            migrate_plain_text(path)
            run(["tabix", "-f", "-p", "bed", path])

    for name in ("truth.vcf.gz",):
        path = os.path.join(data_dir, name)
        if os.path.exists(path):
            print(f"=== {name} ===")
            migrate_plain_text(path)
            run(["tabix", "-f", "-p", "vcf", path])

    for bam_path in glob.glob(os.path.join(data_dir, "reads.*.bam")):
        print(f"=== {os.path.basename(bam_path)} ===")
        migrate_bam(bam_path)

    print("Done. Any previously-scored atarva.*.vcf.gz / baseline_history.tsv rows from "
         "before this migration used the old naming and should be re-run.")


if __name__ == "__main__":
    main()
