import sys, os
import pysam
import copy
import threading
import polars as pl
from functools import reduce

INFO_MP_CUTOFF = 0.5

BASE_SCHEMA = {"C": pl.Categorical,
               "S": pl.Int32,
               "E": pl.Int32,
               "R": pl.Categorical,
               "I": pl.Categorical}
FILE_SCHEMA = {"C": pl.Categorical,
               "S": pl.Int32,
               "E": pl.Int32}

COLUMNS   = ["C", # CHROM
             "S", # START
             "E", # END
             "R", # REF
             "I"] # INFO
FILE_COLUMNS = ["C", # CHROM
                "S", # START
                "E"] # END


def extract_names(vcf_files):
    """
    Extract the sample names from the VCF files

    :param vcf_files: List of VCF file paths
    :return: List of sample names
    """

    sample_names = []
    for f in vcf_files:
        with pysam.VariantFile(f) as vcf:
            sample_names.extend(list(vcf.header.samples))
    return sample_names


def joiner(frames, parquet_batch, tidx, outfile):
    """
    Joining multiple dataframes with different samples

    :param frames: List of dataframes to be joined
    :param parquet_batch: Batch number for the output parquet file
    :param tidx: Thread index for multi-threaded processing
    :param outfile: The output file path for the parquet file
    :return: None, writes the joined dataframe to a parquet file
    """

    base = reduce(lambda l, r: l.join(r, on=['C', 'S', 'E'], how='left'), frames)
    df   = base.collect(engine="streaming")
    df.write_parquet(f"{outfile}_reader{tidx}_batch{parquet_batch}.parquet", compression="zstd")


def write_header(out, bam_name, source_vcf_path):
    """
    Writing the VCF header to the output file based on the source VCF and sample names.

    :param out: The output file handle to write the VCF header.
    :param bam_name: List of sample names extracted from the BAM files.
    :param source_vcf_path: Path to the source VCF file to extract contig information and other metadata.
    :return: None, writes the VCF header to the output file
    """

    source_vcf = pysam.VariantFile(source_vcf_path)
    vcf_header = pysam.VariantHeader()

    # command
    vcf_header.add_line(f"##command=Tamatr {' '.join(sys.argv)}")

    for contig, metadata in source_vcf.header.contigs.items():
        vcf_header.contigs.add(contig, length=metadata.length)
    info_mp_cutoff = source_vcf.header.info["MPC"].description
    source_vcf.close()
    
    #sample_name
    for each_sample in bam_name:
        vcf_header.add_sample(each_sample)

    # FILTER
    vcf_header.filters.add('LESS_READS', number=None, type=None, description="Read depth below threshold")

    # INFO
    vcf_header.info.add("AC",    number='A', type="Integer",  description="Number of alternate alleles in called genotypes")
    vcf_header.info.add("AN",    number=1,   type="Integer",  description="Number of alleles in called genotypes")
    vcf_header.info.add("MOTIF", number=1,   type="String",   description="Repeat motif")
    vcf_header.info.add("START", number=1,   type="Integer",  description="Start position of the repeat region in 0-based coordinate system")
    vcf_header.info.add("END",   number=1,   type="Integer",  description="End position of the repeat region")
    vcf_header.info.add("ID",    number=1,   type="String",   description="Locus identifier tag")
    vcf_header.info.add("REFCN", number=1,   type="Integer",  description="Reference allele copy number")
    vcf_header.info.add("CT",    number=1,   type="String",   description="Cluster type")
    vcf_header.info.add("MPC",   number=1,   type="String",   description=f"{INFO_MP_CUTOFF}")
    vcf_header.info.add("AFD",   number=1,   type="String",   description="Allele length frequency for the reads supporting the locus")

    # FORMAT
    vcf_header.formats.add("GT", number=1,   type="String",   description="Genotype")
    vcf_header.formats.add("AL", number=2,   type="Integer",  description="Allele length in base pairs")
    vcf_header.formats.add("CN", number=2,   type="Integer",  description="Motif copy number for each allele")
    vcf_header.formats.add("AR", number='.', type="String",   description="Allele length range")
    vcf_header.formats.add("SD", number='.', type="Integer",  description="Number of reads supporting each haplogroup")
    vcf_header.formats.add("DP", number=1,   type="Integer",  description="Total supporting reads for the repeat locus")
    vcf_header.formats.add("SN", number='.', type="Integer",  description="Number of informative SNPs used for phasing")
    vcf_header.formats.add("SQ", number='.', type="Float",    description="Average Phred-scale base call quality of each informative SNP across the supporting reads")
    vcf_header.formats.add("MA", number='.', type="Float",    description="Average  methylation level for each allele")
    vcf_header.formats.add("MR", number='.', type="Integer",  description="Number of informative reads for methylation scoring for each allele")
    vcf_header.formats.add("DS", number='A', type="String",   description="Motif decomposed sequence for each allele sequence")
    vcf_header.formats.add("MV", number='.', type="String",   description="Base methylation score encoded for visualization for each allele")
    vcf_header.formats.add("PS", number='1', type="String",   description="Phase Set assigned in the phasing process if available in alignment file")

    out.write(str(vcf_header))


def processor(process_df, outfile, tidx, each_thread, nsamples):
    """
    Processing of each locus with information from all samples and writing the output to a VCF file.

    :param process_df: DataFrame containing the loci and sample information to be processed
    :param outfile: The output VCF file path
    :param tidx: Thread index for multi-threaded processing
    :param each_thread: Thread index for the current processing thread
    :param total_samples: Total number of samples being processed
    :return: None, writes the processed information to a VCF file
    """

    out = open(f'{outfile}_reader{tidx}_processor{each_thread}.vcf', 'w')
    for row in process_df.iter_rows(named=True):
        genotyped = 0
        sample_formats = []

        ALTs = []
        ALT_counts = {}

        for vid in range(nsamples):
            sformat = row[f'F{vid:06d}']

            if sformat:
                sformat = sformat.split(':')
            else:
                sample_formats.append('.:.:.:.:.:.:.:.:.')
                continue

            genotyped += 1
            GT = []
            sALTs    = sformat[0].split(',') if sformat[0] != '.' else ""
            seq_lens = [0 if seq=='<DEL>' else len(seq) for seq in sALTs]

            for ALT in sALTs:
                if ALT not in ALTs:
                    ALTs.append(ALT)
                    ALT_counts[ALT] = 1
                else: 
                    ALT_counts[ALT] += 1

            sGT = sformat[1];
            sep = ''
            if '/' in sGT: sGT = sGT.split('/'); sep = '/'
            elif '|' in sGT: sGT = sGT.split('|'); sep = '|'

            for gt in sGT:
                if gt == '0':
                    GT.append('0')
                else:
                    alt_idx = int(gt) - 1
                    GT.append(str(ALTs.index(sALTs[alt_idx]) + 1))

            sformat[1] = sep.join(GT)
            sample_formats.append(':'.join(sformat[1:]))
        if genotyped:
            pass
        else:
            continue
        if ALTs:
            AC = []
            for ALT in ALTs:
                AC.append(str(ALT_counts[ALT]))
            AC = ','.join(AC)
        else:
            AC = '0'
        AN = str(genotyped * 2)
        info = 'AC='+AC+';AN='+AN+';' + row['I']
        ref_seq = row['R']
        start   = row['S']
        chrom   = row['C']
        filter  = '.'
        id      = '.'
        q       = '.'
        alt     = ','.join(ALTs) if ALTs else '.'
        format  = 'GT:AL:CN:LPM:AR:SD:DP:SN:SQ:MA:MR:DS:MV'

        repeat_info = [chrom, start, id, ref_seq, alt, q, filter, info, format, *sample_formats]
        del sample_formats
        tot_tabs = len(repeat_info)
        chunk_size = 100
        for i in range(0, tot_tabs, chunk_size):
            chunk = repeat_info[i:i + chunk_size]
            out.write("\t".join(map(str, chunk)))
            if i<tot_tabs-1:
                out.write("\t")
        out.write("\n")

        del repeat_info
    out.close()


def chop_tamatar(outfile, bedfile, ref_file, vcf_files, contigs, tidx, process_thread):
    """
    The main function that handles the merging of ATaRVa VCF files and writes the merged output to a new VCF file.

    :param outfile: The output VCF file path.
    :param bedfile: The BED file path containing the regions of interest.
    :param ref: The reference genome file path.
    :param vcfs: A list of ATaRVa VCF file paths to be merged.
    :param contigs: A list of contigs to be processed.
    :param tidx: Thread index for multi-threaded processing.
    :param process_thread: Number of threads to be used for processing.

    :return: None, write output the main output VCF or thread specific VCF files.
    """

    n_samples = len(vcf_files)

    tbx = pysam.TabixFile(bedfile)
    ref = pysam.FastaFile(ref_file)

    # read each vcf file as a tabix input
    vcfs = []
    for f in vcf_files:
        vcfs.append(pysam.TabixFile(f))

    if tidx != -1: # multi thread
        if tidx == 0: # first process
            sample_names = extract_names(vcf_files)
            out = open(f'{outfile}.vcf', 'w')
            write_header(out, sample_names, vcf_files[0])
        else:
            out = open(f'{outfile}_thread_{tidx}.vcf', 'w')
    else: # single thread
        sample_names = extract_names(vcf_files)
        out = open(f'{outfile}.vcf', 'w')
        write_header(out, sample_names, vcf_files[0])

    thread_pool = list()
    print('Reader thread = ', tidx)
    print(f'Inside reader{tidx} = length of contig = {len(contigs)}')
    for contig in contigs:

        Chrom, Start, End = contig
        # print("\nReading new block............")
        base_frame = pl.DataFrame().lazy()
        frames = []
        parquet_batch = 0
        fcount = 0

        for vidx, vcf in enumerate(vcfs):

            # Create a deep copy of the DataFrame schema

            if vidx == 0: columns = COLUMNS + [f'F{vidx:06d}']
            else: columns = FILE_COLUMNS + [f'F{vidx:06d}']
            vcf_data = {col: [] for col in columns}

            if vidx == 0:
                schema = copy.deepcopy(BASE_SCHEMA)
                schema[f'F{vidx:06d}'] = pl.Categorical
                for line in tbx.fetch(Chrom, Start[0], End[1]):
                    line  = line.strip().split('\t')
                    chrom = line[0]
                    start = int(line[1])
                    end   = int(line[2])

                    if (start >= Start[0]) and (end <= End[1]):
                        if start == Start[0]:
                            if end == Start[1]: pass
                            else: continue
                        pass
                    elif start < Start[0]: continue
                    elif start >= End[0]:  break

                    motif    = line[3]
                    ref_alen = end - start
                    REFCN    = ref_alen // float(line[4])
                    ID       = line[5] if len(line) > 5 else "."
                    del line

                    ref_string = True
                    has_region = False

                    if Chrom in vcf.contigs:
                        for record in vcf.fetch(chrom, start+1, end):
                            record = record.strip().split('\t')
                            vpos   = int(record[1])
                            if (vpos - 1) != start: # -1 to match with 0-based coord
                                continue

                            info = {a.split('=')[0]: a.split('=')[1] for a in record[7].split(';', 5)}
                            if ((vpos - 1) == start) and (info['END'] == end): # -1 to match with 0-based coord
                                has_region = True
                                vcf_data['C'].append(chrom)
                                vcf_data['S'].append(vpos)
                                vcf_data['E'].append(info['END'])
                                vcf_data['R'].append(record[3])
                                vcf_data['I'].append(f"MOTIF={motif};START={start};END={end};ID={ID};REFCN={REFCN}")
                                if record[9][0] =='.':
                                    vcf_data[f'F{vidx:06d}'].append(None)
                                else:
                                    vcf_data[f'F{vidx:06d}'].append(record[9])
                                del record
                                del info
                            break

                    if not has_region:
                        vcf_data['C'].append(chrom)
                        vcf_data['S'].append(start + 1)
                        vcf_data['E'].append(end)
                        vcf_data['R'].append(ref.fetch(chrom, start, end))
                        vcf_data['I'].append(f"MOTIF={motif};START={start};END={end};ID={ID};REFCN={REFCN}")
                        vcf_data[f'F{vidx:06d}'].append(None)

                fcount += 1

                df = pl.DataFrame(vcf_data, schema_overrides=schema).lazy()
                df = df.unique(subset=['C', 'S', 'E'], keep='first', maintain_order=True)
                base_frame = df.collect().select(['S', 'E']).lazy()

                frames.append(df)
                del df

            else:
                fcount += 1

                schema = copy.deepcopy(FILE_SCHEMA)
                schema[f'F{vidx:06d}']   = pl.Categorical
                vcf_data[f'F{vidx:06d}'] = []

                if fcount >= 200:
                    joiner(frames, parquet_batch, tidx, outfile)
                    parquet_batch += 1
                    del frames
                    frames = [base_frame]
                    fcount = 0

                if Chrom not in vcf.contigs:
                    df = pl.DataFrame(vcf_data, schema=schema).lazy()
                    frames.append(df)
                    continue

                for record in vcf.fetch(Chrom, Start[0], End[1]):
                    record = record.strip().split('\t')

                    if record[9][0] == '.':
                        continue
                    else:
                        vpos = int(record[1])
                        info = {a.split('=')[0]: a.split('=')[1] for a in record[7].split(';', 5)}

                        vcf_data['C'].append(Chrom)
                        vcf_data['S'].append(vpos)
                        vcf_data['E'].append(int(info['END']))
                        vcf_data[f'F{vidx:06d}'].append(record[4] + ':' + record[9])

                df = pl.DataFrame(vcf_data, schema=schema).lazy()
                df = df.unique(subset=['C', 'S', 'E'], keep='first', maintain_order=True)
                frames.append(df)
                del df

        if frames:
            joiner(frames, parquet_batch, tidx, outfile)
            parquet_batch += 1
            del frames

        # print('Done reading & joining!!!!!!!!!')
        if thread_pool:
            # joining previous threads - waiting for previous threads to be over
            for t in thread_pool: t.join()
            thread_pool.clear()

            # print('Concatenating processor files..............')
            for each_thread in range(process_thread):
                thread_out = f'{outfile}_reader{tidx}_processor{each_thread}.vcf'

                with open(thread_out, 'r') as fh:
                    for line in fh:
                        repeat_info = line.strip().split('\t')
                        tot_tabs    = len(repeat_info)
                        chunk_size  = 100
                        for i in range(0, tot_tabs, chunk_size):
                            chunk = repeat_info[i:i + chunk_size]
                            out.write("\t".join(map(str, chunk)))
                            if i<tot_tabs-1:
                                out.write("\t")
                        out.write("\n")
                        del repeat_info

                os.remove(thread_out)

        batch_files = [f"{outfile}_reader{tidx}_batch{batch_val}.parquet" for batch_val in range(parquet_batch)]
        parquet_frames = [pl.read_parquet(f).lazy() for f in batch_files]
        for p_files in batch_files:
            os.remove(p_files)

        if parquet_frames:
            merged = reduce(lambda l, r: l.join(r, on=['C','S','E'], how='left'), parquet_frames)
            whole_df = merged.collect(engine="streaming")

            if process_thread > 0:
                loci_count  = whole_df.shape[0]
                split_count = loci_count // process_thread
                if split_count == 0:
                    split_count = 1
                initial = 0
                track = split_count

                # initializing threads
                for each_thread in range(process_thread):
                    if each_thread+1 == process_thread:
                        process_df = whole_df[initial : ]
                    else:
                        process_df = whole_df[initial : track]

                    t = threading.Thread(target = processor, args = (process_df, outfile, tidx, each_thread, n_samples))
                    t.start()
                    thread_pool.append(t)

                    initial = track
                    track += split_count

            else:
                processor(whole_df, outfile, tidx, 0, n_samples)
                thread_out = f'{outfile}_reader{tidx}_processor{0}.vcf'
                with open(thread_out, 'r') as fh:
                    for line in fh:
                        repeat_info = line.strip().split('\t')
                        out.write("\t".join(map(str, repeat_info)) + "\n")
                os.remove(thread_out)

            del whole_df

    if thread_pool:
        # joining previous threads - waiting for previous threads to be over
        for thread_x in thread_pool:
            # print('waiting for ', thread_x)
            thread_x.join()
        thread_pool.clear()

        # print('Concatenating processor files..............')
        for each_thread in range(process_thread):
            thread_out = f'{outfile}_reader{tidx}_processor{each_thread}.vcf'
            # print('opening ', thread_out)
            with open(thread_out, 'r') as fh:
                for line in fh:
                    repeat_info = line.strip().split('\t')
                    out.write("\t".join(map(str, repeat_info)) + "\n")
            # print('Removing ', thread_out)
            os.remove(thread_out)

    # closing all the opened files
    for vcf in vcfs: vcf.close()
    ref.close()
    tbx.close()
    out.close()
