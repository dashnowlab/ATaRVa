#!/usr/bin/env python3
"""Score an ATaRVa VCF against the truth set using str-truth-set-v2/str-truth-set's OWN
comparison code, reused as-is, rather than Aardvark.

str-truth-set-v2's actual tool-ranking pipeline (run_tools/run_genotyping_tools.py) does not
use Aardvark at all -- Aardvark is only used by docs/truthset_release_check.md to validate
the truth set itself against an independent GIAB benchmark. The real per-tool comparison
runs a chain of scripts (from broadinstitute/str-truth-set, tool_comparison/scripts/ and
tool_comparison/hail_batch_pipelines/atarva_pipeline.py) that compare ALLELE LENGTH (repeat
copy number), not edit-distance over sequence, and report "% exact match" as the headline
number -- e.g. the tool_comparison_viewer.html rankings.json numbers.

This script reproduces that chain locally on our benchmark subset:
  1. Reverse-map our ATaRVa VCF from pseudo-contig coordinates back to real hg38 coordinates
     (05_reverse_map_vcf.py) -- required because str_analysis.convert_atarva_vcf_to_
     expansion_hunter_json derives each locus's LocusId from CHROM/POS/START/END, which has
     to match the truth genotypes table's real-genome LocusIds for the join to find anything.
  2. str_analysis.convert_atarva_vcf_to_expansion_hunter_json --discard-hom-ref
  3. str_analysis.combine_str_json_to_tsv -> {sample}.variants.tsv.gz
  4. Filter HG002.tandem_repeat_genotypes.tsv.gz (the full genome-wide truth genotype table)
     down to our selected loci
  5. compute_truth_set_tsv_for_comparisons.py on the filtered truth table
  6. add_tool_results_columns.py --tool ATaRVa (left-joins the tool's calls onto every truth
     locus; a locus with no tool call becomes a No Call, not a dropped row)
  7. add_concordance_columns.py --tool ATaRVa --compare-to Truth -> per-variant "Discordant" /
     "OverlappingCIs" / "ExactlyTheSame" label

"% exact match" = ExactlyTheSame / total truth loci -- comparable to rankings.json's "pct".

Requires (installed once into the ATaRVa venv):
    pip install 'git+https://github.com/broadinstitute/str-analysis@3d5e3dc37161d41a1bc6f92afdcf3fc99e81b30a' \
        intervaltree 'pandas<3'
(pandas>=3 defaults to strict string dtypes that these scripts, written pre-pandas-3, don't
tolerate -- see the TypeError this raises if pandas>=3 is active.)

Usage:
    python3 06_score_strts_native.py --tech pacbio --atarva-version 0.7.1+ext0.01
"""
import argparse
import os
import subprocess
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import baseline_history

TRUTH_GENOTYPES_TSV_URL = "https://storage.googleapis.com/str-truth-set-v2/filter_vcf_v2/HG002/HG002.tandem_repeat_genotypes.tsv.gz"
COVERAGE_LABELS = {"pacbio": "30x", "ont": "26x"}

SCRIPTS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "strts_reused")
REVERSE_MAP_SCRIPT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "05_reverse_map_vcf.py")


def run(cmd, **kwargs):
    print(f"> {' '.join(cmd)}", flush=True)
    return subprocess.run(cmd, check=True, **kwargs)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tech", choices=["pacbio", "ont"], required=True)
    parser.add_argument("--atarva-version", required=True,
                        help="Version tag of the ATaRVa build that produced the VCF being "
                             "scored, e.g. \"0.7.1+ext0.01\" -- from `atarva[-ext] --version`. "
                             "Selects benchmark/data/atarva.<tech>.<version>.vcf.gz as input "
                             "and tags the appended baseline_history.tsv rows.")
    parser.add_argument("--catalog", default="random_2to6bp_with_extremes",
                        help="Label for the locus catalog in --data-dir -- tags the appended "
                             "baseline_history.tsv rows. See 04_score.py's --catalog help for "
                             "the convention.")
    parser.add_argument("--data-dir", default="benchmark/data")
    parser.add_argument("--research-dir", default="benchmark/research")
    args = parser.parse_args()

    if pd.__version__.split(".")[0] >= "3":
        sys.exit(f"pandas {pd.__version__} is active; these reused scripts need pandas<3 "
                 f"(strict string-dtype coercion in pandas>=3 breaks them). "
                 f"pip install 'pandas<3' in this venv.")

    data_dir = args.data_dir
    work_dir = os.path.join(data_dir, "strts_comparison")
    os.makedirs(work_dir, exist_ok=True)
    sample = f"reads.{args.tech}.{args.atarva_version}"

    genome_vcf = os.path.join(data_dir, f"atarva.{args.tech}.{args.atarva_version}.genome.vcf")
    run([sys.executable, REVERSE_MAP_SCRIPT,
        os.path.join(data_dir, f"atarva.{args.tech}.{args.atarva_version}.vcf.gz"),
        genome_vcf,
        os.path.join(data_dir, "windows.bed")])
    run(["bash", "-c", f"bgzip -f {genome_vcf}"])
    run(["tabix", "-p", "vcf", genome_vcf + ".gz"])

    run([sys.executable, "-m", "str_analysis.convert_atarva_vcf_to_expansion_hunter_json",
        "--discard-hom-ref", "--sample-id", sample, genome_vcf + ".gz"])

    # the converter names its output by stripping only ".vcf.gz" from the input path
    genome_json = os.path.abspath(genome_vcf[:-len(".vcf")] + ".json")
    run([sys.executable, "-m", "str_analysis.combine_str_json_to_tsv",
        "--output-prefix", sample, genome_json], cwd=work_dir)

    variants_tsv_path = os.path.join(work_dir, f"{sample}.1_json_files.variants.tsv.gz")
    variants_df = pd.read_table(variants_tsv_path)
    variants_df["Coverage"] = COVERAGE_LABELS[args.tech]
    variants_df.to_csv(variants_tsv_path, sep="\t", index=False)

    truth_genotypes_path = os.path.join(args.research_dir, "HG002.tandem_repeat_genotypes.tsv.gz")
    if not os.path.exists(truth_genotypes_path):
        run(["curl", "-fsS", "-o", truth_genotypes_path, TRUTH_GENOTYPES_TSV_URL])

    loci_df = pd.read_csv(os.path.join(data_dir, "loci.bed"), sep="\t", header=None,
                          names=["Chrom", "Start0Based", "End", "Motif", "MotifSize"])
    selected_keys = set(zip(loci_df.Chrom, loci_df.Start0Based, loci_df.End))
    truth_df = pd.read_table(truth_genotypes_path, low_memory=False)
    truth_subset_df = truth_df[[k in selected_keys for k in
                                zip(truth_df.Chrom, truth_df.Start0Based, truth_df.End)]]
    truth_subset_path = os.path.join(work_dir, "truth.tandem_repeat_genotypes.subset.tsv.gz")
    truth_subset_df.to_csv(truth_subset_path, sep="\t", index=False)
    print(f"Truth subset: {len(truth_subset_df)} / {len(selected_keys)} selected loci found in the truth table")

    run([sys.executable, os.path.join(SCRIPTS_DIR, "compute_truth_set_tsv_for_comparisons.py"),
        "--output-dir", work_dir, truth_subset_path])

    for_comparison_path = os.path.join(
        work_dir, os.path.basename(truth_subset_path).replace(".tsv.gz", ".for_comparison.tsv.gz"))
    merged_path = os.path.join(work_dir, f"{sample}.with_ATaRVa_results.tsv.gz")
    run([sys.executable, os.path.join(SCRIPTS_DIR, "add_tool_results_columns.py"),
        "--tool", "ATaRVa", "--output-tsv", merged_path, variants_tsv_path, for_comparison_path])

    concordance_path = os.path.join(work_dir, f"{sample}.with_concordance.tsv")
    run([sys.executable, os.path.join(SCRIPTS_DIR, "add_concordance_columns.py"),
        "--tool", "ATaRVa", "--compare-to", "Truth", "--output-tsv", concordance_path, merged_path])

    concordance_df = pd.read_csv(concordance_path, sep="\t")
    col = "Variant: Concordance: ATaRVa vs Truth"

    def exact_match_pct(df):
        counts = df[col].value_counts()
        total = len(df)
        exact = int(counts.get("ExactlyTheSame", 0))
        return exact, total, 100 * exact / total if total else 0.0

    # A length-only match within a repeat-unit tolerance, alongside the strict (tolerance=0)
    # exact match above. NumRepeats -- not RepeatSize (bp) -- matches what add_concordance_
    # columns.py itself compares. The ATaRVa paper (Sivakumar et al., bioRxiv 2025.05.13.653434)
    # reports both an exact-match figure and a more lenient +/-1bp-tolerant figure as their
    # headline numbers; NaN handling mirrors add_concordance_columns.py's own special case
    # (a missing tool allele only counts as a match when the tool's overall call is hom-ref).
    merged_df = pd.read_table(merged_path)

    def within_tolerance(row, tolerance):
        for allele in (1, 2):
            truth = row[f"NumRepeats: Allele {allele}: Truth"]
            tool = row[f"NumRepeats: Allele {allele}: ATaRVa"]
            if pd.isna(tool):
                if not row.get("IsHomRef: ATaRVa", False):
                    return False
                continue
            if pd.isna(truth) or abs(truth - tool) > tolerance:
                return False
        return True

    def tolerance_pct(df, tolerance):
        total = len(df)
        n_match = int(df.apply(lambda r: within_tolerance(r, tolerance), axis=1).sum())
        return n_match, total, 100 * n_match / total if total else 0.0

    metrics = {}
    print(f"\n=== ATaRVa {args.atarva_version} vs HG002 truth set "
         f"({args.tech}, str-truth-set's own comparison code) ===")

    exact, total, pct = exact_match_pct(concordance_df)
    metrics["ExactMatch/ALL"] = pct
    print(f"ExactMatch/ALL                {pct:.1f}%  ({exact}/{total})")

    bucket_2to6 = concordance_df[(concordance_df.MotifSize >= 2) & (concordance_df.MotifSize <= 6)]
    exact, total, pct = exact_match_pct(bucket_2to6)
    metrics["ExactMatch/2to6bp_motifs"] = pct
    print(f"ExactMatch/2to6bp              {pct:.1f}%  ({exact}/{total})")

    n_match, total, pct = tolerance_pct(merged_df, 1)
    metrics["LengthMatch_within_1_unit/ALL"] = pct
    print(f"LengthMatch (+/-1 unit)/ALL    {pct:.1f}%  ({n_match}/{total})")

    merged_2to6 = merged_df[(merged_df.MotifSize >= 2) & (merged_df.MotifSize <= 6)]
    n_match, total, pct = tolerance_pct(merged_2to6, 1)
    metrics["LengthMatch_within_1_unit/2to6bp_motifs"] = pct
    print(f"LengthMatch (+/-1 unit)/2to6bp {pct:.1f}%  ({n_match}/{total})")

    baseline_history.append(args.atarva_version, args.tech, args.catalog, "strts_native", metrics)


if __name__ == "__main__":
    main()
