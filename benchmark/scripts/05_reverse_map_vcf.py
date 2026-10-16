#!/usr/bin/env python3
"""Reverse-map an ATaRVa VCF from our synthetic pseudo-contig coordinates back to real
hg38 genome coordinates.

Needed to reuse str-truth-set's own comparison scripts (str_analysis.convert_atarva_vcf_to_
expansion_hunter_json etc.), which derive each locus's LocusId directly from the VCF's
CHROM/POS and its START/END INFO fields -- those have to match the truth genotypes table's
real-genome LocusIds for the join in add_tool_results_columns.py to find anything.

Plain text-based (not pysam VariantRecord construction): the only edits needed are the
CHROM column, the POS column, and the START/END INFO sub-fields -- everything else in each
line passes through untouched. This sidesteps pysam's reserved-field handling of END/stop,
which fought the record-by-record approach.

Usage: 05_reverse_map_vcf.py <in.vcf.gz> <out.vcf> <windows.bed>
"""
import gzip
import re
import sys


def parse_window_name(name):
    m = re.match(r"^(chr[0-9XYM]+)s(\d+)e(\d+)$", name)
    if not m:
        raise ValueError(f"Unexpected pseudo-contig name: {name}")
    chrom, start, end = m.group(1), int(m.group(2)), int(m.group(3))
    return chrom, start, end


def rewrite_info(info, offset):
    if info in (".", ""):
        return info
    parts = []
    for field in info.split(";"):
        if "=" in field:
            key, value = field.split("=", 1)
            if key in ("START", "END"):
                value = str(int(value) + offset)
            parts.append(f"{key}={value}")
        else:
            parts.append(field)
    return ";".join(parts)


def main():
    in_vcf_gz_path, out_vcf_path, windows_bed_path = sys.argv[1], sys.argv[2], sys.argv[3]

    real_chroms = set()
    with open(windows_bed_path) as f:
        for line in f:
            real_chroms.add(line.split("\t")[0])

    header_lines = []
    data_rows = []  # (real_chrom, real_pos_int, line)
    n_skipped = 0

    with gzip.open(in_vcf_gz_path, "rt") as f:
        for line in f:
            if line.startswith("##contig"):
                continue  # replaced below with real chromosome names
            if line.startswith("#"):
                if line.startswith("#CHROM"):
                    # insert real ##contig lines right before the #CHROM header line
                    header_lines.extend(f"##contig=<ID={c}>\n" for c in sorted(real_chroms))
                header_lines.append(line)
                continue

            fields = line.rstrip("\n").split("\t")
            chrom, pos, info = fields[0], int(fields[1]), fields[7]
            alt = fields[4]
            if alt == ".":
                n_skipped += 1  # no-ALT (hom-ref/no-call) record -- uninformative, drop it
                continue

            real_chrom, wstart, wend = parse_window_name(chrom)
            new_pos = wstart + pos
            fields[0] = real_chrom
            fields[1] = str(new_pos)
            fields[7] = rewrite_info(info, wstart)
            data_rows.append((real_chrom, new_pos, "\t".join(fields) + "\n"))

    data_rows.sort(key=lambda r: (r[0], r[1]))

    with open(out_vcf_path, "w") as out:
        out.writelines(header_lines)
        out.writelines(row[2] for row in data_rows)

    print(f"Wrote {len(data_rows)} reverse-mapped records to {out_vcf_path} ({n_skipped} no-ALT records skipped)")


if __name__ == "__main__":
    main()
