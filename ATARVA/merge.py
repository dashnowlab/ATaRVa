from math import remainder, sqrt, ceil
import sys, os
import pysam
import timeit as ti

from multiprocessing import Process
from ATARVA.tamatr import chop_tamatar


def merge_parser(subparsers):
    """
    Merge ATaRVa VCF files over specified regions.

    :param subparsers: ArgumentParser subparsers object
    :return: None
    """

    parser = subparsers.add_parser("merge", help="merging multiple ATaRVa VCF files over specified regions", description="Merge ATaRVa VCF files")
    parser._action_groups.pop()

    required = parser.add_argument_group('Required arguments')
    required.add_argument('-i', '--vcfs', nargs='+', required=True, metavar="<FILE>", help="ATaRVa output VCF files to be merged.")
    required.add_argument('-r', '--regions', required=True, metavar="<FILE>", help="sorted, bgzipped and indexed bed file containing regions of interest.")
    required.add_argument('-f', '--fasta', required=True, metavar='<FILE>', help='input reference fasta file. The file should be indexed.')

    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument('--contigs', nargs='+', help='contigs to get merged [chr1 chr12 chr22 ..]. If not mentioned every contigs in the region file will be merged.')
    optional.add_argument('-o', '--output', type=str, metavar='<STR>', default='', help='name of the output file, output is in vcf format.')
    optional.add_argument('-t', '--threads', type=int, metavar='<INT>', default=1,  help='number of threads. [default: 1]')
    optional.add_argument('--mem', type=int, metavar='<INT>', default=4, help='memory in GB for the whole program. [default: 4]')

    if (len(sys.argv) == 2) and (sys.argv[1] == 'merge'):
        parser.print_help()
        sys.exit()

    parser.set_defaults(func=run_merge)


def count_loci(args, tbx):
    """
    Count the number of loci in a bed file.

    :param bedfile: Path to the bed file
    :return: Total number of loci in the bed file
    """

    n = 0
    contigs = []
    if not args.contigs:
        contigs = sorted(tbx.contigs)
        for row in tbx.fetch(): n += 1
    else:
        contigs = sorted(args.contigs)
        for c in sorted(args.contigs):
            for row in tbx.fetch(c): n += 1
    return n, contigs


def shredder(threads):
    """
    Splits the total number of threads into groups for parallel processing.

    :param threads: Total number of threads
    :return: List of thread counts for each group
    """
    ngrps      = ceil(sqrt(threads))
    base, rem  = divmod(threads, ngrps)
    thread_grp = [base + 1] * rem + [base] * (ngrps - rem)
    return len(thread_grp), thread_grp


def run_merge(args):
    """
    Merge ATaRVa VCF files over specified regions.

    :param args: Parsed command-line arguments
    :return: None
    """
    print('atarva', ' '.join(sys.argv[1:]))
    sys.exit()

    start_time = ti.default_timer()

    for arg in vars(args):
        if arg in ['func', 'command']: continue
        print (arg, getattr(args, arg))
    print('\n')

    out = sys.stdout
    if args.output:
        if '.vcf' == args.output[-4:]:
            out = f'{args.output}'[:-4]
        elif args.output[-1]=='/':
            out = args.output + "atarva_merged"
        else:
            out = f'{args.output}'
    else:
        out = "atarva_merged"

    vcf_list = []
    if len(args.vcfs) == 1:
        with open(args.vcfs[0], 'r') as vh:
            for line in vh:
                if line[0] == '#': continue
                line = line.strip()
                vcf_list.append(line)
    else:
       vcf_list = args.vcfs

    tbx  = pysam.Tabixfile(args.regions)
    total_loci, contigs = count_loci(args, tbx)

    threads = args.threads
    threads = threads - 1

    split_point = 5000 if total_loci > 5000 else total_loci // 5
    if split_point == 0:
        split_point = 1
        partition_point = 1
    else:
        partition_point = total_loci//split_point

    split_point_chunks = 0 # to count number of split_point chunks excluding the 'minimum chunks' eg 9920 from 1 contig and 80 from another contig to add up to 10000
    fetcher       = []
    current_split = []
    line_count    = 0

    for contig in contigs:
        init = 0
        for row in tbx.fetch(contig):
            line_count += 1
            row = row.strip().split('\t')
            if init == 0:
                chrom = row[0]
                start_coord = (int(row[1]), int(row[2]))
                init = 1

            if split_point_chunks < partition_point-1:
                if line_count % split_point == 0:
                    end_coord = (int(row[1]), int(row[2]))
                    current_split.append([chrom, start_coord, end_coord])
                    # fetcher.append(tuple(current_split))
                    fetcher.extend(current_split)
                    split_point_chunks += 1
                    line_count = 0
                    current_split = []
                    init = 0

        if init != 0:
            end_coord = (int(row[1]), int(row[2]))
            current_split.append([chrom, start_coord, end_coord])

    fetcher.extend(current_split)
    tbx.close()

    region_file = args.regions
    ref_file    = args.fasta
    if threads > 1:

        nprocs, thread_grps = shredder(threads)
        thread_pool = []
        partition = len(fetcher) // nprocs
        initial = 0
        track = partition
        for tidx in range(nprocs):
            if tidx == nprocs - 1:
                reader_contigs = fetcher[initial : ]
            else:
                reader_contigs = fetcher[initial : track]

            t = Process(target = chop_tamatar, args = (out, region_file, ref_file, vcf_list, reader_contigs,
                                                       tidx, thread_grps[tidx]))
            t.start()
            thread_pool.append(t)

            initial = track
            track += partition

        # joining Threads
        for t in thread_pool: t.join()
        # emptying thread_pool
        thread_pool.clear()
        #sys.exit()
        out = open(f'{out}.vcf', 'a')

        print('Concatenating thread outputs!', file=sys.stderr)
        for tidx in range(nprocs):
            thread_out = f'{out}_thread_{tidx}.vcf'
            print(thread_out)
            with open(thread_out, 'r') as fh:
                # if tidx!=0: next(fh)
                for line in fh:
                    repeat_info = line.strip().split('\t')
                    #print(*repeat_info, file=out, sep='\t')
                    out.write("\t".join(map(str, repeat_info)) + "\n")
            os.remove(thread_out)
        out.close()
        print('Concatenation completed!! ^_^', file=sys.stderr)

    else:
        chop_tamatar(out, region_file, ref_file, vcf_list, fetcher, -1, 0)

    time_now = ti.default_timer()
    sys.stderr.write('CPU time: {} seconds\n'.format(time_now - start_time))
