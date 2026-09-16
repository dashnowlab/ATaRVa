#!/usr/bin/env python3
"""Select the benchmark's TR locus catalog from the str-truth-set-v2 HG002 catalog.

Design ("random_2to6bp_with_extremes"): a plain (non-stratified) random sample of
TARGET_TOTAL_LOCI loci restricted to 2-6bp motifs (classic STRs -- the motif range most
published tool comparisons report headline numbers for), plus a separate random sample of
N_EXTREMES_PER_SIDE homopolymers (motif_len == 1) and N_EXTREMES_PER_SIDE VNTRs
(motif_len >= VNTR_MOTIF_LEN_MIN) -- the two motif-size extremes the core sample excludes.
All three pools are disjoint (motif-length ranges don't overlap), so no dedup is needed
between them. Every locus is restricted to HG002's dipcall confident regions, since that's
the only region where the truth set's genotypes are considered reliable.

Inputs (research/ dir, already downloaded):
    research/catalog.confident.bed -- str-truth-set-v2 HG002 catalog, confidence-filtered

Output:
    data/loci.bed             -- final locus list, ATaRVa's 5-column catalog format
    data/loci.annotated.tsv   -- same loci with a "source" column, for documentation
"""
import argparse
import os
import random

TARGET_TOTAL_LOCI = 1000
N_EXTREMES_PER_SIDE = 100
VNTR_MOTIF_LEN_MIN = 7
RANDOM_SEED = 42


def load_bed(path):
    rows = []
    with open(path) as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            rows.append(fields)
    return rows


def select_random(pool, n, rng):
    sampled_rows = rng.sample(pool, min(n, len(pool)))
    print(f"Sampled {len(sampled_rows)} (from a pool of {len(pool)})")
    return sampled_rows


def select_random_with_extremes(all_rows, target_total, n_extremes_per_side, seed):
    rng = random.Random(seed)

    core_pool = [r for r in all_rows if 2 <= int(r[4]) <= 6]
    print("2-6bp motifs:", end=" ")
    core_rows = select_random(core_pool, target_total, rng)

    homopolymer_pool = [r for r in all_rows if int(r[4]) == 1]
    print("Homopolymers:", end=" ")
    homopolymer_rows = select_random(homopolymer_pool, n_extremes_per_side, rng)

    vntr_pool = [r for r in all_rows if int(r[4]) >= VNTR_MOTIF_LEN_MIN]
    print("VNTRs:", end=" ")
    vntr_rows = select_random(vntr_pool, n_extremes_per_side, rng)

    return ([(r, "2to6bp") for r in core_rows] +
           [(r, "homopolymer") for r in homopolymer_rows] +
           [(r, "vntr") for r in vntr_rows])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--research-dir", default="benchmark/research")
    parser.add_argument("--data-dir", default="benchmark/data")
    parser.add_argument("--target-total", type=int, default=TARGET_TOTAL_LOCI)
    parser.add_argument("--n-extremes-per-side", type=int, default=N_EXTREMES_PER_SIDE)
    parser.add_argument("--seed", type=int, default=RANDOM_SEED)
    args = parser.parse_args()

    os.makedirs(args.data_dir, exist_ok=True)

    catalog_path = os.path.join(args.research_dir, "catalog.confident.bed")
    all_rows = load_bed(catalog_path)

    final_rows = select_random_with_extremes(all_rows, args.target_total, args.n_extremes_per_side, args.seed)
    final_rows.sort(key=lambda x: (x[0][0], int(x[0][1])))

    loci_bed_path = os.path.join(args.data_dir, "loci.bed")
    annotated_path = os.path.join(args.data_dir, "loci.annotated.tsv")
    with open(loci_bed_path, "w") as bed_out, open(annotated_path, "w") as ann_out:
        ann_out.write("#CHROM\tSTART\tEND\tMOTIF\tMOTIF_LEN\tsource\n")
        for row, source in final_rows:
            bed_out.write("\t".join(row) + "\n")
            ann_out.write("\t".join(row) + f"\t{source}\n")

    print(f"Total loci selected: {len(final_rows)}")
    print(f"Wrote {loci_bed_path} and {annotated_path}")

    source_counts = {}
    for _, source in final_rows:
        source_counts[source] = source_counts.get(source, 0) + 1
    print("Counts by source:", source_counts)

    chrom_counts = {}
    for row, _ in final_rows:
        chrom_counts[row[0]] = chrom_counts.get(row[0], 0) + 1
    print("Per-chromosome counts:", chrom_counts)


if __name__ == "__main__":
    main()
