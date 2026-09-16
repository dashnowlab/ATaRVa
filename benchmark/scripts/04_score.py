#!/usr/bin/env python3
"""Score an ATaRVa VCF against the benchmark subset's truth VCF using Aardvark.

Reuses str-truth-set-v2's approach (docs/truthset_release_check.md in that repo):
Aardvark compares two VCFs by edit distance over haplotype sequences within each
scored region, so two callers that describe the same haplotype with different VCF
records score the same. Only Aardvark's BASEPAIR recall/precision/F1 is reported (not
its GT metric, a stricter hap.py-style per-record match) -- see guarded_metrics()'s
comment for why GT isn't a useful signal here.

Unlike the release check (which validates the truth set itself against an independent
GIAB benchmark), this scores a tool's calls against the truth set directly, so the only
region restriction needed is where the truth set is itself confident (no GIAB benchmark
BED dependency).

Aardvark v1.0.0 only ships a linux-x86_64 binary; this runs it via Docker
(--platform linux/amd64) so it works on Apple Silicon without a Rust toolchain.

Usage:
    python3 04_score.py --tech pacbio --atarva-version 0.7.1+ext0.01
"""
import argparse
import csv
import gzip
import os
import shutil
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import baseline_history

MIN_VARIANT_GAP_BP = 1000
TRUTH_SAMPLE = "syndip"


def run(cmd, **kwargs):
    print(f"> {' '.join(cmd)}", flush=True)
    return subprocess.run(cmd, check=True, **kwargs)


def uppercase_vcf(src_path, dst_gz_path):
    """Uppercase REF/ALT (dipcall carries hg38 soft-masking through into alleles).

    Writes a real BGZF file (via `bgzip`), not plain gzip, so it stays tabix-indexable.
    """
    n_lowercase = 0
    tmp_path = dst_gz_path[:-3] if dst_gz_path.endswith(".gz") else dst_gz_path + ".tmp"
    with gzip.open(src_path, "rt") as src, open(tmp_path, "w") as dst:
        for line in src:
            if line.startswith("#"):
                dst.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            ref, alt = fields[3], fields[4]
            if any(c.islower() for c in ref + alt):
                n_lowercase += 1
            fields[3] = ref.upper()
            fields[4] = alt.upper()
            dst.write("\t".join(fields) + "\n")
    run(["bash", "-c", f"bgzip -f {tmp_path}"])
    return n_lowercase


def decompress(src_gz_path, dst_path):
    with gzip.open(src_gz_path, "rt") as src, open(dst_path, "w") as dst:
        shutil.copyfileobj(src, dst)


def guess_query_sample(vcf_gz_path):
    with gzip.open(vcf_gz_path, "rt") as f:
        for line in f:
            if line.startswith("#CHROM"):
                fields = line.rstrip("\n").split("\t")
                return fields[9]
    raise ValueError(f"No #CHROM header line found in {vcf_gz_path}")


def run_aardvark_docker(aardvark_dir, data_dir, work_dir, reference_path, truth_vcf_path,
                        truth_sample, query_vcf_path, query_sample, regions_path,
                        min_variant_gap, compare_label, output_dir_name):
    # Stage inputs into a plain system temp dir (tempfile.gettempdir(), not work_dir) before
    # mounting into Docker. Bind-mounting straight from a cloud-synced checkout (OneDrive,
    # Dropbox, iCloud Drive) intermittently fails on macOS -- Docker Desktop's file-sharing
    # layer and the sync client's own file-locking step on each other -- surfacing as "Resource
    # deadlock avoided" reading the reference genome. A local, non-synced staging copy sidesteps
    # that regardless of where the repo itself lives.
    final_output_dir = os.path.join(work_dir, output_dir_name)
    if os.path.exists(final_output_dir):
        shutil.rmtree(final_output_dir)

    with tempfile.TemporaryDirectory(prefix="atarva_benchmark_aardvark_") as stage_dir:
        staged = {}
        for name, path in (("reference", reference_path), ("truth", truth_vcf_path),
                           ("query", query_vcf_path), ("regions", regions_path)):
            dest = os.path.join(stage_dir, os.path.basename(path))
            shutil.copy2(path, dest)
            for ext in (".fai", ".tbi"):
                if os.path.exists(path + ext):
                    shutil.copy2(path + ext, dest + ext)
            staged[name] = dest

        stage_output_dir = os.path.join(stage_dir, output_dir_name)
        os.makedirs(stage_output_dir)

        def in_container(path):
            return "/work/" + os.path.relpath(path, stage_dir)

        run([
            "docker", "run", "--rm", "--platform", "linux/amd64",
            "-v", f"{aardvark_dir}:/aardvark",
            "-v", f"{stage_dir}:/work",
            "debian:bookworm-slim",
            "/aardvark/aardvark", "compare",
            "--reference", in_container(staged["reference"]),
            "--truth-vcf", in_container(staged["truth"]), "--truth-sample", truth_sample,
            "--query-vcf", in_container(staged["query"]), "--query-sample", query_sample,
            "--regions", in_container(staged["regions"]),
            "--min-variant-gap", str(min_variant_gap),
            "--compare-label", compare_label,
            "--output-dir", in_container(stage_output_dir),
        ])
        shutil.copytree(stage_output_dir, final_output_dir)

    return os.path.join(final_output_dir, "summary.tsv")


def read_summary_tsv(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def guarded_metrics(summary_rows):
    # Only BASEPAIR (representation-invariant edit distance over the reconstructed haplotype
    # sequence) is reported here, not Aardvark's GT metric (a stricter, hap.py-style per-record
    # match). GT drops whenever ATaRVa's realignment represents the same underlying allele as a
    # differently-shaped VCF record than the truth set does, which is common and largely benign
    # at TR loci specifically -- it isn't a genotyping error, so it isn't a useful signal here.
    metrics = {}
    for row in summary_rows:
        if row["variant_type"] != "ALL" or row["comparison"] != "BASEPAIR":
            continue
        if row["region_label"] != "ALL":
            continue
        for name in ("recall", "precision", "f1"):
            metrics[f"{row['comparison']}/{name}"] = float(row[f"metric_{name}"])
    return metrics


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tech", choices=["pacbio", "ont"], required=True)
    parser.add_argument("--atarva-version", required=True,
                        help="Version tag of the ATaRVa build that produced the VCF being "
                             "scored, e.g. \"0.7.1+ext0.01\" -- from `atarva[-ext] --version`. "
                             "Selects benchmark/data/atarva.<tech>.<version>.vcf.gz as the "
                             "query and tags the appended baseline_history.tsv rows.")
    parser.add_argument("--catalog", default="random_2to6bp_with_extremes",
                        help="Label for the locus catalog in --data-dir -- tags the appended "
                             "baseline_history.tsv rows. Earlier catalog designs (stratified_1000, "
                             "pathogenic_only, random_1000_all_motifs, random_1000_2to6bp) were "
                             "explored and retired; their historical rows remain in "
                             "baseline_history.tsv but their data/scripts were removed.")
    parser.add_argument("--data-dir", default="benchmark/data")
    parser.add_argument("--aardvark-dir",
                        default="benchmark/tools/aardvark-v1.0.0-x86_64-unknown-linux-gnu")
    parser.add_argument("--min-variant-gap", type=int, default=MIN_VARIANT_GAP_BP)
    args = parser.parse_args()

    data_dir = os.path.abspath(args.data_dir)
    aardvark_dir = os.path.abspath(args.aardvark_dir)

    truth_uppercased_path = os.path.join(data_dir, "truth.uppercased.vcf.gz")
    n_lowercase = uppercase_vcf(os.path.join(data_dir, "truth.vcf.gz"), truth_uppercased_path)
    run(["tabix", "-f", "-p", "vcf", truth_uppercased_path])
    print(f"Uppercased truth VCF: {n_lowercase} records had a lowercase base")

    regions_bed_path = os.path.join(data_dir, "scoring_regions.bed")
    decompress(os.path.join(data_dir, "scoring_regions.bed.gz"), regions_bed_path)

    query_vcf_path = os.path.join(data_dir, f"atarva.{args.tech}.{args.atarva_version}.vcf.gz")
    query_sample = guess_query_sample(query_vcf_path)
    print(f"Query sample: {query_sample}")

    summary_path = run_aardvark_docker(
        aardvark_dir=aardvark_dir, data_dir=data_dir, work_dir=data_dir,
        reference_path=os.path.join(data_dir, "ref.fa"),
        truth_vcf_path=truth_uppercased_path, truth_sample=TRUTH_SAMPLE,
        query_vcf_path=query_vcf_path, query_sample=query_sample,
        regions_path=regions_bed_path,
        min_variant_gap=args.min_variant_gap,
        compare_label=f"atarva_benchmark_{args.tech}",
        output_dir_name=f"aardvark_output.{args.tech}.{args.atarva_version}")

    summary_rows = read_summary_tsv(summary_path)
    metrics = guarded_metrics(summary_rows)

    print(f"\n=== ATaRVa {args.atarva_version} vs HG002 truth set "
         f"({args.tech}, {len(summary_rows)} summary rows) ===")
    for key in sorted(metrics):
        print(f"{key:20s} {metrics[key]:.6f}")

    baseline_history.append(args.atarva_version, args.tech, args.catalog, "aardvark_basepair", metrics)


if __name__ == "__main__":
    main()
