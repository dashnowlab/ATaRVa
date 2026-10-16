#!/usr/bin/env python3
"""
Merge many single-sample tandem-repeat VCFs into one multi-sample VCF,
splitting work into chunks with ~equal numbers of loci.

Strategy
  * inputs: bgzipped + tabix-indexed single-sample VCFs from the SAME locus catalog
  * catalog: catalog file
  * a reference fasta file is required with the catalog is provided 
  * a reference VCF is scanned once to place chunk boundaries every K loci;
    chunks may span contigs (list of (contig, lo, hi) segments, lo <= POS < hi)
  * each worker opens all inputs once, k-way merges its segments with a heap,
    and writes a BGZF-compressed part file
  * parts are concatenated byte-wise (valid BGZF) in chunk order and tabix-indexed
  * locus key = (POS, END, REF, TRID/ID); ALT alleles unioned by sequence, GT remapped
  * missing loci -> ./.

Requirements: pysam.   For thousands of inputs:  ulimit -n 65536
Memory note: each open TabixFile holds its index in RAM -> ~ threads x samples x index size.

Usage:
  python merge_tr_vcfs_chunked.py -l vcf_list.txt -o merged.vcf.gz -t 16
  python merge_tr_vcfs_chunked.py -l list.txt -o out.vcf.gz -t 16 --chunks-per-thread 1 --ref-vcf best.vcf.gz
"""
import argparse
import heapq
import os
import re
import shutil
import sys
import tempfile
from array import array
from multiprocessing import Pool

import pysam

END_RE       = re.compile(r"(?:^|;)END=(\d+)")
TRID_RE      = re.compile(r"(?:^|;)ID=([^;]+)")
GT_SPLIT     = re.compile(r"([/|])")
DEFAULT_KEYS = ['GT', 'AL', 'CN', 'AR', 'SD', 'DP', 'SN', 'SQ', 'MA', 'MR', 'DS', 'MV', 'OL', 'FV', 'PS']
FLUSH_EVERY  = 100

_H = []  # per-process tabix handles, set by init_worker


def merge_parser(subparsers):
    """
    Merge ATaRVa VCF files over specified regions.

    :param subparsers: ArgumentParser subparsers object
    :return: None
    """

    parser = subparsers.add_parser("merge", help="merging multiple ATaRVa VCF files over specified regions", description="Merge ATaRVa VCF files",
                                   formatter_class=argparse.RawDescriptionHelpFormatter)
    parser._action_groups.pop()

    required = parser.add_argument_group('Required arguments')
    required.add_argument('-i', '--vcfs', nargs='+', required=True, metavar="<FILE>", help="ATaRVa output VCF files to be merged.")
    required.add_argument('-r', '--regions',         required=True, metavar='<FILE>', help="bgzipped + tabix-indexed regions file. Prepare: sort with bedtools"
                                                                                           "→ bgzip → tabix index. Include with reference fasta file.")
    required.add_argument('-f', '--fasta',           required=True, metavar='<FILE>', help="input reference fasta file. requed if --regions is provided. "
                                                                                           "Prepare: samtools faidx.")
    required.add_argument('-o', '--out',  type=str,  required=True, metavar='<FILE>', help="name of the output file, output is in vcf format.")

    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument('--contigs',       nargs='+', help='contigs to get merged [chr1 chr12 chr22 ..]. If not mentioned every contigs in the region file will be merged.')
    optional.add_argument('--format-keys',   nargs='+', metavar='<STR>', default=','.join(DEFAULT_KEYS), help=f'FORMAT keys to keep in output. GT is by default part of the output. (default: {','.join(DEFAULT_KEYS)})')
    optional.add_argument('-t', '--threads', type=int,  metavar='<INT>', default=1,    help='number of threads. [default: 1]')

    if (len(sys.argv) == 2) and (sys.argv[1] == 'merge'):
        parser.print_help()
        sys.exit()

    parser.set_defaults(func=run_merge)


def init_worker(paths):
    global _H
    _H = [pysam.TabixFile(p) for p in paths]


def read_header(path):
    tbx     = pysam.TabixFile(path)
    hdr     = list(tbx.header)
    contigs = list(tbx.contigs)
    tbx.close()
    return hdr, contigs


def sample_name(header_lines):
    chrom_line = [h for h in header_lines if h.startswith("#CHROM")][0]
    return chrom_line.rstrip("\n").split("\t")[9]


def segment_iter(tbx, contig, lo, hi):
    """Yield (key, fields) for records with lo <= POS < hi (hi=None -> open end)."""
    try:
        it = tbx.fetch(contig, max(lo - 1, 0), None if hi is None else hi - 1)
    except ValueError:  # contig absent in this file
        return
    for line in it:
        f = line.split("\t", 10)
        pos = int(f[1])
        if pos < lo or (hi is not None and pos >= hi):
            continue  # tabix returns overlapping records; keep only those starting here
        m = END_RE.search(f[7])
        t = TRID_RE.search(f[7])
        locus_id = t.group(1) if t else f[2]
        yield (pos, int(m.group(1)) if m else 0, f[3], locus_id), f


def catalog_iter(tbx, contig, lo, hi, fasta):
    """Yield (key, fields) for records with lo <= POS < hi (hi=None -> open end)."""
    try:
        it = tbx.fetch(contig, max(lo, 0), hi)
    except ValueError:  # contig absent in this file
        return
    for line in it:
        f = line.split("\t")
        pos = int(f[1])
        end = int(f[2])
        motif = f[3]
        locus_id = f'{f[0]}-{pos}-{end}-{motif}'

        if pos < lo or pos > hi:
            continue  # tabix returns overlapping records; keep only those starting here

        ref_seq = fasta.fetch(f[0], pos, end)
        yield (pos+1, int(f[2]), ref_seq, locus_id), f + ["CATALOG"]  # add a sentinel to distinguish catalog records from sample records


def remap_gt(gt, amap):
    parts = GT_SPLIT.split(gt)
    for i in range(0, len(parts), 2):
        a = parts[i]
        if a not in (".", "0", ""):
            parts[i] = str(amap[int(a)])
    return "".join(parts)


def tamatr(args):

    cid, catalog, ref_file, segments, keys, tmpdir = args
    tbx   = pysam.TabixFile(catalog)
    fasta = pysam.FastaFile(ref_file)

    n = len(_H)
    missing_str = ":".join("./." if k == "GT" else "." for k in keys)
    fmt_cache = {}
    out_path = os.path.join(tmpdir, f"chunk_{cid:06d}.part.gz")
    buf = []

    with pysam.libcbgzf.BGZFile(out_path, "wb") as out:
        for segment in segments:
            contig, lo, hi = segment
            lo = lo[0]
            hi = hi[0]
            iters = [catalog_iter(tbx, contig, lo, hi, fasta)] + [segment_iter(h, contig, lo, hi) for h in _H]
            heap = []
            for i, it in enumerate(iters):
                nxt = next(it, None)
                if nxt:
                    heap.append((nxt[0], i, nxt[1]))
            heapq.heapify(heap)

            while heap:
                key = heap[0][0]
                group = []
                locus = None
                while heap and heap[0][0] == key:
                    _, i, f = heapq.heappop(heap)
                    if f[-1] == "CATALOG": locus = f
                    else: group.append((i, f))
                    nxt = next(iters[i], None)
                    if nxt:
                        heapq.heappush(heap, (nxt[0], i, nxt[1]))

                first = None
                alts  = []
                cols  = [missing_str] * n
                if len(group) > 0:
                    first = group[0][1]
                    alts = set()
                    for _, f in group:
                        if f[4] != ".":
                            alts.update(f[4].split(","))
                    alts = sorted(alts, key=lambda s: (len(s), s))
                    alt_idx = {a: j + 1 for j, a in enumerate(alts)}

                    cols = [missing_str] * n
                    for i, f in group:
                        fmt = f[8]
                        sel = fmt_cache.get(fmt)
                        if sel is None:
                            names = fmt.split(":")
                            sel = [names.index(k) if k in names else -1 for k in keys]
                            fmt_cache[fmt] = sel
                        vals = f[9].rstrip("\n").split(":")
                        local = f[4].split(",") if f[4] != "." else []
                        amap = {j + 1: alt_idx[a] for j, a in enumerate(local)}
                        row = []
                        for k, s in zip(keys, sel):
                            v = vals[s] if 0 <= s < len(vals) else "."
                            if k == "GT" and v != ".":
                                v = remap_gt(v, amap)
                            row.append(v)
                        cols[i-1] = ":".join(row)

                elif len(group) == 0:
                    ref_seq = key[-1]
                    first = [locus[0], str(int(locus[1])+1), '.', ref_seq, '0', 'FAIL', '.', ':'.join(keys)]
                    alts = []

                buf.append(
                    "\t".join([first[0], first[1], first[2], first[3], ",".join(alts) or ".",
                               ".", "PASS", first[7], ":".join(keys)])
                    + "\t" + "\t".join(cols) + "\n"
                )
                if len(buf) >= FLUSH_EVERY:
                    out.write("".join(buf).encode())
                    buf.clear()
        if buf:
            out.write("".join(buf).encode())

    tbx.close()
    fasta.close()
    return cid, out_path


def splitfile_threads(regions, contigs, threads):
    """
    split the input region file into per-thread chunks for parallel processing.

    :param tbx:        pysam Tabixfile object
    :param contigs:    list of contigs to process
    :param total_loci: total number of loci across all contigs
    :param threads:    number of parallel threads
    :return:           list of region chunk tuples, one per thread
    """

    tbx = pysam.TabixFile(regions)
    total_loci = 0
    for contig in contigs:
        for _ in tbx.fetch(contig):
            total_loci += 1

    split_size  = total_loci // threads
    if split_size == 0:
        split_size = total_loci
    fetcher     = []
    line_count  = 0
    cur_chunk   = []
    chrom       = None
    start_coord = None

    for contig in contigs:
        first_in_contig = True
        itr = tbx.fetch(contig)

        for row in itr:
            fields = row.split('\t')
            line_count += 1

            if first_in_contig:
                chrom       = fields[0]
                start_coord = (int(fields[1]), int(fields[2]))
                first_in_contig = False

            if len(fetcher) < threads - 1 and line_count >= split_size:
                while True:
                    try: row = next(itr)
                    except StopIteration:
                        end_coord = (int(fields[1]), int(fields[2]))
                        cur_chunk.append([chrom, start_coord, end_coord])
                        fetcher.append(tuple(cur_chunk))
                        cur_chunk       = []
                        start_coord     = None
                        line_count      = 0
                        first_in_contig = True
                        break
                    _fields = row.split('\t')

                    if int(_fields[1]) > int(fields[2]) + 10:
                        end_coord = (int(fields[1]), int(fields[2]))
                        cur_chunk.append([chrom, start_coord, end_coord])
                        fetcher.append(tuple(cur_chunk))
                        cur_chunk       = []
                        start_coord     = (int(_fields[1]), int(_fields[2]))
                        line_count      = 1
                        first_in_contig = True
                        break

        if not first_in_contig:
            end_coord = (int(fields[1]), int(fields[2]))
            cur_chunk.append([chrom, start_coord, end_coord])

    fetcher.append(tuple(cur_chunk)) # pad with empty chunks if needed
    tbx.close()
    return fetcher, total_loci


def run_merge(args):
    
    paths = list(args.vcfs or [])
    if not paths: sys.exit("Input VCFs are not provided!")

    keys = args.format_keys.split(",")
    if keys[0] != "GT":
        keys = ['GT'] + keys
        sys.exit("--format-keys must start with GT")

    catalog = args.regions
    fasta   = args.fasta

    hdr0 = None
    contigs, seen, samples = [], set(), []
    for p in paths:
        h, cs = read_header(p)
        if hdr0 is None: hdr0 = h
        samples.append(sample_name(h))
        for c in cs:
            if c not in seen:
                seen.add(c)
                contigs.append(c)
    if args.contigs:
        contigs = args.contigs.split(",")
    if len(set(samples)) != len(samples):
        sys.exit("Error: Duplicate sample names among inputs")

    avail  = {l.split("ID=")[1].split(",")[0]: l for l in hdr0 if l.startswith("##FORMAT=<ID=")}
    header = [l.rstrip("\n") for l in hdr0
              if l.startswith("##") and not l.startswith("##FORMAT")
              and not l.startswith("##bcftools") and not l.startswith("##source=")]
    header += [avail[k].rstrip("\n") for k in keys if k in avail]
    header.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples))

    fetcher, total_loci = splitfile_threads(catalog, contigs, args.threads)

    print(f"Reference loci: {total_loci}; threads: {args.threads} (~{-(-total_loci // len(fetcher))} loci each)", file=sys.stderr)

    tmpdir = tempfile.mkdtemp(prefix="trmerge_", dir=os.path.dirname(os.path.abspath(args.out)))
    try:
        jobs = [(i, catalog, fasta, segs, keys, tmpdir) for i, segs in enumerate(fetcher)]
        with Pool(min(args.threads, len(jobs)), initializer=init_worker, initargs=(paths,)) as pool:
            results = dict(pool.imap_unordered(tamatr, jobs))

        if not args.out.endswith("gz"):
            args.out += ".gz"

        with pysam.libcbgzf.BGZFile(args.out, "wb") as out:
            out.write(("\n".join(header) + "\n").encode())
        with open(args.out, "ab") as out:
            for i in range(len(fetcher)):  # chunk order == contig/position order
                with open(results[i], "rb") as part:
                    shutil.copyfileobj(part, out)
        pysam.tabix_index(args.out, preset="vcf", force=True)
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)
    print(f"Merged {len(paths)} samples -> {args.out}", file=sys.stderr)
