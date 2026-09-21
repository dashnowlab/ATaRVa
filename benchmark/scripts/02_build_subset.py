#!/usr/bin/env python3
"""Build a small, self-contained benchmark bundle around the selected TR loci.

Each selected locus (data/loci.bed) is padded by FLANK_BP and overlapping/adjacent
padded windows are merged, giving a set of disjoint genomic windows. Each window becomes
its own small "pseudo-contig" in a synthetic reference (named
"<chrom>_<0-based window start>_<window end>"), with coordinates reset to 0. The catalog,
truth VCF, dipcall-confident BED and read alignments are all remapped into this same
coordinate space, so the resulting bundle is a tiny, self-consistent "mini genome" that
any TR genotyper can run against directly -- no whole-genome reference or multi-GB BAM
required.

Reads are fetched per-window with remote `samtools view` region queries (only the bytes
for that window are pulled over HTTP, via the BAM's .bai index) and soft-clipped to the
window boundary, since PacBio HiFi/ONT reads are routinely much longer than a
locus-sized window. The full original SEQ/QUAL is kept (not trimmed) so nothing about
the read itself is altered -- only which portion of it the CIGAR marks as aligned.

Usage:
    python3 02_build_subset.py --data-dir benchmark/data --tech pacbio
    python3 02_build_subset.py --data-dir benchmark/data --tech ont
"""
import argparse
import gzip
import os
import subprocess
import sys

import pysam

REFERENCE_FASTA_URL = "https://storage.googleapis.com/str-truth-set/hg38/ref/hg38.fa"
# str-truth-set-v2's filter step (run_filter_vcf_to_tandem_repeats.py) runs
# str_analysis.filter_vcf_to_tandem_repeats over the genome-wide confident-regions VCF and
# writes this TR-only VCF (one record per TR locus, MOTIF/TR_END INFO fields) as a separate
# output -- this is how that project itself separates TR variation from the surrounding
# background SNPs/indels before any tool comparison, and it's what we want to score
# ATaRVa against too (NOT HG002.high_confidence_regions.vcf.gz, which is genome-wide and
# would count every unrelated flanking small variant as a false negative, since ATaRVa
# only emits one variant per TR locus).
TRUTH_VCF_URL = "https://storage.googleapis.com/str-truth-set-v2/filter_vcf_v2/HG002/HG002.tandem_repeats.vcf.gz"
CATALOG_URL = "https://storage.googleapis.com/str-truth-set-v2/filter_vcf_v2/HG002/HG002.bed.gz"
DIPCALL_CONFIDENT_BED_URL = "https://storage.googleapis.com/str-truth-set-v2/dipcall_pipeline/HG002/HG002.dip.bed.gz"

READ_BAM_URLS = {
    # str-truth-set-v2's own tool rankings use PacBio downsampled to exactly 30x (matches
    # HG002.downsampled_to_30x.bam's own total_depth.txt) and ONT at its native, *not*
    # downsampled, coverage -- HG002.bam's total_depth.txt is 26.36x, which is exactly the
    # "26x" label their rankings.json uses. HG002.downsampled_to_20x.bam is a different,
    # lower-coverage cut they also publish but don't use for the headline ONT numbers.
    "pacbio": "https://storage.googleapis.com/str-truth-set-v2/raw_data/HG002/pacbio/HG002.downsampled_to_30x.bam",
    "ont": "https://storage.googleapis.com/str-truth-set-v2/raw_data/HG002/ONT/HG002.bam",
}

FLANK_BP = 1500

# Aardvark drops any variant not fully contained in a scored region, so like
# str-truth-set-v2's own release check, loci are padded a little before being scored --
# just enough to catch a variant that straddles the locus edge. This is much tighter than
# FLANK_BP (which exists purely to give ATaRVa realignment context): the truth VCF is a
# genome-wide small-variant VCF, not TR-only, so scoring across the full extraction window
# would count every unrelated flanking SNP/indel that ATaRVa never attempts to call as a
# false negative.
SCORING_PADDING_BP = 50

REF_CONSUMING = {0, 2, 3, 7, 8}    # M, D, N, =, X
QUERY_CONSUMING = {0, 1, 4, 7, 8}  # M, I, S, =, X


def run(cmd, **kwargs):
    print(f"> {' '.join(cmd)}", flush=True)
    return subprocess.run(cmd, check=True, **kwargs)


def load_loci(loci_bed_path):
    loci = []
    with open(loci_bed_path) as f:
        for line in f:
            chrom, start, end, motif, motif_len = line.rstrip("\n").split("\t")
            loci.append((chrom, int(start), int(end), motif, int(motif_len)))
    return loci


def build_windows(loci, flank_bp, work_dir):
    """Pad each locus and merge overlapping/adjacent windows. Returns sorted window list."""
    padded_path = os.path.join(work_dir, "padded_loci.bed")
    with open(padded_path, "w") as f:
        for chrom, start, end, _, _ in loci:
            f.write(f"{chrom}\t{max(0, start - flank_bp)}\t{end + flank_bp}\n")

    sorted_path = os.path.join(work_dir, "padded_loci.sorted.bed")
    run(["bash", "-c", f"sort -k1,1 -k2,2n {padded_path} > {sorted_path}"])

    merged_path = os.path.join(work_dir, "windows.bed")
    with open(merged_path, "w") as out:
        run(["bedtools", "merge", "-i", sorted_path], stdout=out)

    windows = []
    with open(merged_path) as f:
        for line in f:
            chrom, start, end = line.rstrip("\n").split("\t")
            windows.append((chrom, int(start), int(end)))

    os.remove(padded_path)
    os.remove(sorted_path)
    return windows


def window_name(chrom, start, end):
    # Letters, not underscores, as separators: ATaRVa v0.7.1's own region-file validation
    # requires the CHROM column to be str.isalnum() (no underscores or other punctuation).
    return f"{chrom}s{start}e{end}"


def build_reference(windows, out_fasta_path):
    """Extract each window's sequence from the remote hg38 FASTA via faidx, one region per call."""
    with open(out_fasta_path, "w") as out:
        for chrom, start, end in windows:
            region = f"{chrom}:{start + 1}-{end}"
            result = subprocess.run(
                ["samtools", "faidx", REFERENCE_FASTA_URL, region],
                check=True, capture_output=True, text=True)
            lines = result.stdout.splitlines()
            out.write(f">{window_name(chrom, start, end)}\n")
            out.write("\n".join(lines[1:]) + "\n")
    run(["samtools", "faidx", out_fasta_path])


def remap_bed_like(src_lines_iter, windows_by_chrom, coord_cols, extra_cols_fn=None):
    """Yield remapped rows for lines whose [start,end) falls inside some window on the same chrom.

    Args:
        src_lines_iter: iterable of split fields, coord_cols[0]/[1] are 0-based start/end (ints already).
        windows_by_chrom: dict chrom -> sorted list of (start, end).
    """
    for fields in src_lines_iter:
        chrom = fields[0]
        start, end = fields[coord_cols[0]], fields[coord_cols[1]]
        for wstart, wend in windows_by_chrom.get(chrom, []):
            if wstart <= start and end <= wend:
                new_fields = list(fields)
                new_fields[0] = window_name(chrom, wstart, wend)
                new_fields[coord_cols[0]] = start - wstart
                new_fields[coord_cols[1]] = end - wstart
                yield new_fields
                break


def build_catalog(loci, windows_by_chrom, out_path):
    rows = ([chrom, start, end, motif, str(motif_len)] for chrom, start, end, motif, motif_len in loci)
    remapped = list(remap_bed_like(rows, windows_by_chrom, (1, 2)))
    tmp_path = out_path[:-3] if out_path.endswith(".gz") else out_path + ".tmp"
    with open(tmp_path, "w") as f:
        for fields in sorted(remapped, key=lambda r: (r[0], r[1])):
            f.write("\t".join(str(x) for x in fields) + "\n")
    run(["bash", "-c", f"bgzip -f {tmp_path}"])
    run(["tabix", "-p", "bed", tmp_path + ".gz"])
    print(f"Catalog: {len(remapped)} / {len(loci)} loci remapped (rest fell outside a window -- unexpected)")


def build_confident_bed(dip_bed_gz_path, windows, windows_by_chrom, out_path):
    """Intersect dipcall-confident regions with our windows, then remap into window coordinates."""
    windows_bed_path = out_path + ".windows.bed"
    with open(windows_bed_path, "w") as f:
        for chrom, start, end in windows:
            f.write(f"{chrom}\t{start}\t{end}\n")

    intersected_path = out_path + ".intersected.bed"
    with open(intersected_path, "w") as out:
        run(["bash", "-c",
             f"gzcat {dip_bed_gz_path} | sort -k1,1 -k2,2n | "
             f"bedtools intersect -a - -b {windows_bed_path} | sort -k1,1 -k2,2n"],
            stdout=out)

    rows = []
    with open(intersected_path) as f:
        for line in f:
            chrom, start, end = line.rstrip("\n").split("\t")
            rows.append([chrom, int(start), int(end)])

    # These are sub-intervals of a window (not necessarily equal to it), so a fields-based
    # remap that requires full containment works directly.
    remapped = list(remap_bed_like(rows, windows_by_chrom, (1, 2)))
    tmp_path = out_path[:-3] if out_path.endswith(".gz") else out_path + ".tmp"
    with open(tmp_path, "w") as f:
        for fields in sorted(remapped, key=lambda r: (r[0], r[1])):
            f.write("\t".join(str(x) for x in fields) + "\n")
    run(["bash", "-c", f"bgzip -f {tmp_path}"])
    run(["tabix", "-p", "bed", tmp_path + ".gz"])
    os.remove(windows_bed_path)
    os.remove(intersected_path)


def build_scoring_regions(loci, windows_by_chrom, dip_bed_gz_path, padding_bp, out_path):
    """Build tight per-locus scoring regions for Aardvark: each locus padded by padding_bp
    (not the full FLANK_BP extraction window), merged, intersected with dipcall-confident
    regions, and remapped into window coordinates. See SCORING_PADDING_BP for why this has
    to be much tighter than the extraction window.
    """
    padded_path = out_path + ".padded_loci.bed"
    with open(padded_path, "w") as f:
        for chrom, start, end, _, _ in loci:
            f.write(f"{chrom}\t{max(0, start - padding_bp)}\t{end + padding_bp}\n")

    sorted_path = out_path + ".padded_loci.sorted.bed"
    run(["bash", "-c", f"sort -k1,1 -k2,2n {padded_path} > {sorted_path}"])

    merged_path = out_path + ".padded_loci.merged.bed"
    with open(merged_path, "w") as out:
        run(["bedtools", "merge", "-i", sorted_path], stdout=out)

    intersected_path = out_path + ".intersected.bed"
    with open(intersected_path, "w") as out:
        run(["bash", "-c",
             f"gzcat {dip_bed_gz_path} | sort -k1,1 -k2,2n | "
             f"bedtools intersect -a {merged_path} -b - | sort -k1,1 -k2,2n"],
            stdout=out)

    rows = []
    with open(intersected_path) as f:
        for line in f:
            chrom, start, end = line.rstrip("\n").split("\t")
            rows.append([chrom, int(start), int(end)])

    remapped = list(remap_bed_like(rows, windows_by_chrom, (1, 2)))
    tmp_path = out_path[:-3] if out_path.endswith(".gz") else out_path + ".tmp"
    with open(tmp_path, "w") as f:
        for fields in sorted(remapped, key=lambda r: (r[0], r[1])):
            f.write("\t".join(str(x) for x in fields) + "\n")
    run(["bash", "-c", f"bgzip -f {tmp_path}"])
    run(["tabix", "-p", "bed", tmp_path + ".gz"])
    for p in (padded_path, sorted_path, merged_path, intersected_path):
        os.remove(p)
    print(f"Scoring regions: {len(remapped)} intervals")


def build_truth_vcf(windows, windows_by_chrom, out_path):
    """Remap the truth VCF's records into window coordinates via remote tabix + pysam."""
    header_written = False
    tmp_vcf_path = out_path[:-3] if out_path.endswith(".gz") else out_path + ".tmp"
    n_written = 0
    with pysam.VariantFile(TRUTH_VCF_URL) as src:
        out_header = src.header.copy()
        out_header.contigs.clear_header()  # rebuild contigs to match our pseudo-contigs
        for chrom, start, end in windows:
            out_header.contigs.add(window_name(chrom, start, end), length=end - start)

        with pysam.VariantFile(tmp_vcf_path, "w", header=out_header) as out:
            for chrom, start, end in windows:
                region = f"{chrom}:{start + 1}-{end}"
                try:
                    fetch_iter = src.fetch(region=region)
                except ValueError:
                    continue
                for rec in fetch_iter:
                    rec_start = rec.start  # 0-based
                    rec_end = rec.stop     # 0-based exclusive
                    if not (start <= rec_start and rec_end <= end):
                        continue  # record straddles the window edge -- drop rather than corrupt it
                    new_rec = out.new_record(
                        contig=window_name(chrom, start, end),
                        start=rec_start - start,
                        stop=rec_end - start,
                        alleles=rec.alleles,
                        id=rec.id,
                        qual=rec.qual,
                        filter=list(rec.filter) if rec.filter else None)
                    for sample_name, sample_data in rec.samples.items():
                        if sample_name not in new_rec.samples:
                            continue
                        for key, value in sample_data.items():
                            try:
                                new_rec.samples[sample_name][key] = value
                            except Exception:
                                pass
                    out.write(new_rec)
                    n_written += 1
    run(["bash", "-c", f"bgzip -f {tmp_vcf_path}"])
    run(["tabix", "-p", "vcf", tmp_vcf_path + ".gz"])
    print(f"Truth VCF: {n_written} records remapped")


def clip_read_to_window(cigartuples, ref_start, wstart, wend):
    """Soft-clip a read's alignment to [wstart, wend). Returns (new_ref_start, new_cigartuples) or None
    if the read does not actually overlap the window. SEQ/QUAL are left untouched by the caller --
    only the CIGAR's boundaries change, so the full original read is preserved under the soft clips.
    """
    ref_pos = ref_start
    out = []
    lead_clip = 0
    trail_clip = 0
    new_ref_start = None
    in_window = False
    finished = False

    for op, length in cigartuples:
        if finished:
            if op in QUERY_CONSUMING:
                trail_clip += length
            continue

        if op not in REF_CONSUMING:
            if not in_window:
                lead_clip += length
            else:
                out.append((op, length))
            continue

        seg_start, seg_end = ref_pos, ref_pos + length
        if seg_end <= wstart:
            if op in QUERY_CONSUMING:
                lead_clip += length
            ref_pos = seg_end
            continue
        if seg_start >= wend:
            if op in QUERY_CONSUMING:
                trail_clip += length
            ref_pos = seg_end
            finished = True
            continue

        kept_start, kept_end = seg_start, seg_end
        if kept_start < wstart:
            before_len = wstart - kept_start
            if op in QUERY_CONSUMING:
                lead_clip += before_len
            kept_start = wstart
        after_len = 0
        if kept_end > wend:
            after_len = kept_end - wend
            kept_end = wend

        kept_len = kept_end - kept_start
        if kept_len > 0:
            if not in_window:
                new_ref_start = kept_start
                in_window = True
            out.append((op, kept_len))
        if after_len > 0:
            if op in QUERY_CONSUMING:
                trail_clip += after_len
            finished = True
        ref_pos = seg_end

    if new_ref_start is None:
        return None

    final_cigar = []
    if lead_clip > 0:
        final_cigar.append((4, lead_clip))
    final_cigar.extend(out)
    if trail_clip > 0:
        final_cigar.append((4, trail_clip))

    return new_ref_start - wstart, final_cigar


def build_bam(windows, bam_url, out_bam_path, ref_fasta_path):
    header_contigs = [{"LN": end - start, "SN": window_name(chrom, start, end)} for chrom, start, end in windows]
    header = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6", "SO": "unsorted"}, "SQ": header_contigs})

    unsorted_path = out_bam_path + ".unsorted.bam"
    n_kept = 0
    n_seen = 0
    seen_query_names_per_window = 0
    with pysam.AlignmentFile(unsorted_path, "wb", header=header) as out:
        with pysam.AlignmentFile(bam_url) as src:
            for chrom, start, end in windows:
                region = f"{chrom}:{start + 1}-{end}"
                tid = out.get_tid(window_name(chrom, start, end))
                try:
                    fetch_iter = src.fetch(region=region)
                except ValueError:
                    continue
                for read in fetch_iter:
                    n_seen += 1
                    if read.is_unmapped or read.is_secondary or read.is_supplementary:
                        continue
                    clipped = clip_read_to_window(read.cigartuples, read.reference_start, start, end)
                    if clipped is None:
                        continue
                    new_ref_start, new_cigar = clipped

                    new_read = pysam.AlignedSegment(header)
                    new_read.query_name = read.query_name
                    new_read.query_sequence = read.query_sequence
                    new_read.flag = read.flag & ~0x1 & ~0x8 & ~0x20  # drop paired/mate-unmapped/mate-reverse bits
                    new_read.reference_id = tid
                    new_read.reference_start = new_ref_start
                    new_read.mapping_quality = read.mapping_quality
                    new_read.cigartuples = new_cigar
                    new_read.next_reference_id = tid
                    new_read.next_reference_start = new_ref_start
                    new_read.template_length = 0
                    new_read.query_qualities = read.query_qualities
                    # MD/NM are only valid against the original reference; drop them and every tag that
                    # encodes original-genome coordinates (SA), regenerating MD/NM via samtools calmd below.
                    new_read.set_tags([(k, v) for k, v in read.get_tags() if k not in ("MD", "NM", "SA")])
                    out.write(new_read)
                    n_kept += 1

    run(["samtools", "sort", "-o", out_bam_path, unsorted_path])
    os.remove(unsorted_path)

    calmd_path = out_bam_path + ".calmd.bam"
    run(["bash", "-c", f"samtools calmd -b {out_bam_path} {ref_fasta_path} > {calmd_path} 2>/dev/null"])
    os.replace(calmd_path, out_bam_path)
    run(["samtools", "index", out_bam_path])
    print(f"BAM: saw {n_seen} region-query hits, kept {n_kept} primary alignments after clipping")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", default="benchmark/data")
    parser.add_argument("--tech", choices=["pacbio", "ont"], required=True)
    parser.add_argument("--flank-bp", type=int, default=FLANK_BP)
    parser.add_argument("--skip-shared", action="store_true",
                        help="Skip building reference/catalog/truth-vcf/confident-bed (already built)")
    args = parser.parse_args()

    loci_bed_path = os.path.join(args.data_dir, "loci.bed")
    loci = load_loci(loci_bed_path)
    windows = build_windows(loci, args.flank_bp, args.data_dir)
    windows_by_chrom = {}
    for chrom, start, end in windows:
        windows_by_chrom.setdefault(chrom, []).append((start, end))
    print(f"{len(loci)} loci -> {len(windows)} merged windows "
          f"(total {sum(e - s for _, s, e in windows):,} bp)")

    ref_fasta_path = os.path.join(args.data_dir, "ref.fa")
    if not args.skip_shared:
        print("=== Building reference ===")
        build_reference(windows, ref_fasta_path)

        print("=== Building catalog ===")
        build_catalog(loci, windows_by_chrom, os.path.join(args.data_dir, "catalog.bed.gz"))

        print("=== Building confident regions BED ===")
        dip_bed_local = os.path.join(args.data_dir, "..", "research", "HG002.dip.bed.gz")
        if not os.path.exists(dip_bed_local):
            run(["curl", "-fsS", "-o", dip_bed_local, DIPCALL_CONFIDENT_BED_URL])
        build_confident_bed(dip_bed_local, windows, windows_by_chrom,
                            os.path.join(args.data_dir, "confident.bed.gz"))

        print("=== Building scoring regions BED ===")
        build_scoring_regions(loci, windows_by_chrom, dip_bed_local, SCORING_PADDING_BP,
                              os.path.join(args.data_dir, "scoring_regions.bed.gz"))

        print("=== Building truth VCF ===")
        build_truth_vcf(windows, windows_by_chrom, os.path.join(args.data_dir, "truth.vcf.gz"))

    print(f"=== Building {args.tech} BAM ===")
    bam_out_path = os.path.join(args.data_dir, f"reads.{args.tech}.bam")
    build_bam(windows, READ_BAM_URLS[args.tech], bam_out_path, ref_fasta_path)

    print("Done.")


if __name__ == "__main__":
    main()
