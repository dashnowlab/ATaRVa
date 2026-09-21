# ATaRVa - a tandem repeat genotyper
![Badge-PyPI](https://img.shields.io/badge/PyPI-v{{VERSION}}-brightgreen)
![Badge-License](https://img.shields.io/badge/License-MIT-blue)

<p>
  <b>NOTE:</b> This is a forked version of ATaRVa. The original repository can be found at <a href="https://github.com/SowpatiLab/ATaRVa">https://github.com/SowpatiLab/ATaRVa</a> There have been the following changes made to the original code:
  <ul>
    <li>The softclipped portions of the read are by default searched for the presence of repeat loci from the catalogs by searching for the flanking sequences.</li>
    <li>Implements a slightly different version of read-haplogrouping based on the informative SNPs. The reads are heirarchically grouped based on the minimising the Qvalue penalty of the informative SNPs.</li>
    <li>Added a new option `--rna` to enable genotyping of tandem repeats from RNA-seq data. This option is designed parse the CIGAR string of the typicaly RNA-seq alignment which include introns and other splicing events encoded as N.</li>
    <li>Added a new option `--instability` which outputs an additional file with allele information from each read at a locus.</li>
</p>

<p align=center>
  <img src="lib/ATaRVa_logo.png" alt="Logo of ATaRVa" width="200"/>
</p>

ATaRVa (pronounced uh-thur-va, IPA: /əθərvə/, Sanskrit: अथर्व) is a technology-agnostic tandem repeat genotyper, specially designed for long read data. The name expands to **A**nalysis of **Ta**ndem **R**epeat **Va**riation, and is derived from the the Sanskrit word _Atharva_ meaning knowledge.


## Motivation
Long-read sequencing propelled comprehensive analysis of tandem repeats (TRs) in genomes. Current long-read TR genotypers are either platform specific or computationally inefficient. ATaRva outperforms existing tools while running an order of magnitude faster. ATaRVa also supports multi-threading, haplotyping, motif decomposition and methylation profiling, making it an invaluable tool for population scale TR analyses.

## Table of contents:

* [Installation](#installation)
  * [PyPI installation](#pypi-installation)
* [Usage](#usage)
  * [`genotype` command](#genotype-command)
    * [Reference genome](#reference-genome)
    * [Alignment file](#alignment-file)
    * [Region file](#region-file)
  * [`merge` command](#merge-command)
* [Changelog](#changelog)
* [Analysis script](#analysis-script)
* [Citation](#citation)
* [Contact](#contact)

## Installation

### Source installation
This version of ATaRVa can be installed from the source code:<br>
It is recommended to install this inside a Python virtual environment.

```bash
# Create a python env
$ python -m venv atarva_env

# Activate the env
$ source atarva_env/bin/activate
$ pip install build

# Download the git repo
$ git clone https://github.com/avvaruakshay/ATaRVa.git

# Install
$ cd ATaRVa
$ python -m build
$ pip install .

# Deactivate the env
$ deactivate
```
Both of the methods add a console command `atarva`, which can be executed from any directory


## Usage
The help message and available subcommands can be accessed using

```bash
$ atarva -h
#  or
$ atarva --help
```
which gives the following output

```
ATaRVa - Analysis of Tandem Repeat Variants
Dashnow lab

Usage:
    atarva [OPTIONS] <COMMAND>

Commands:
  genotype  Tandem Repeat Genotyper
  merge     Merge ATaRVa VCF files

Options:
  -h, --help     Print help
  -v, --version  Print version
```

## `genotype` command
`atarva genotype` accepts read alignments and a set of TR regions of interest and outputs TR genotypes, including the consensus sequence, allele length, and decomposed motifs.

Overview of the ATaRVa worflow:
1. ATaRVa processes the input BAM file read-wise, assuming that most long reads span multiple TR loci.
2. After flank realignment and adjustment of read-wise allele lengths, ATaRVa clusters reads into haplotypes using nearby informative *SNV*s, or applies a *edit-distance* based clustering approach when SNV information is unavailable.
3. It derives consensus allele sequences using partial order alignment, decomposes each TR allele into motif-level representations, and outputs the results in VCF format.

The help message and available options can be accessed using

```bash
$ atarva genotype -h
#  or
$ atarva genotype --help
```
which gives the following output

```
usage: atarva genotype [-h] -f <FILE> -b <FILE> [<FILE> ...] -r <FILE> [-o <FILE>] [--aln-format <STR>] [--rna] [--instability] [--contigs <STR> [<STR> ...]] [--karyotype <STR> [<STR> ...]] [-q <INT>] [--min-reads <INT>]
                           [--max-reads <INT>] [--flank <INT>] [--snp-dist <INT>] [--snp-count <INT>] [--snp-qual <INT>] [--snp-read <FLOAT>] [--phasing-read <FLOAT>] [--haplotag <STR>] [--meth-prob <FLOAT>] [--methviz] [--read-wise]
                           [--locus-wise] [--decompose] [-t <INT>] [-log] [-v]

Tandem Repeat Genotyper

Required arguments:
  -f <FILE>, --fasta <FILE>
                        input reference fasta file
  -b <FILE> [<FILE> ...], --bam <FILE> [<FILE> ...]
                        sample alignment files [SAM | BAM | CRAM]
  -r <FILE>, --regions <FILE>
                        bgzipped + tabix-indexed regions file. Prepare: sort with bedtools → bgzip → tabix index

Optional arguments:
  -o <FILE>, --vcf <FILE>
                        output VCF file [default: stdout]
  --aln-format <STR>    alignment format [cram | bam | sam] [default: bam]
  --skip-softclip       skip softclip processing [default: False]
  --rna                 if the input alignment data is RNA-seq [default: False]
  --instability         generates read level allele information for each locus as TSV [default: False]
  --contigs <STR> [<STR> ...]
                        contigs to genotype e.g. chr1 chr12 [default: all]
  --karyotype <STR> [<STR> ...]
                        sample karyotypes e.g. XY XX
  -q <INT>, --map-qual <INT>
                        minimum mapping quality [default: 5]
  --min-reads <INT>     minimum read coverage at a locus [default: 10]
  --max-reads <INT>     maximum reads per locus [default: 100]
  --flank <INT>         flank length (bp) to search for insertions [default: 10]
  --snp-dist <INT>      max SNP distance from repeat [default: 3000]
  --snp-count <INT>     number of SNPs for phasing [default: 3]
  --snp-qual <INT>      min base quality at SNP position [default: 20]
  --snp-read <FLOAT>    min SNP read fraction [default: 0.2]
  --phasing-read <FLOAT>
                        min phased read cluster fraction [default: 0.4]
  --haplotag <STR>      haplotag for phasing e.g. HP [default: None]
  --meth-prob <FLOAT>   min methylation probability [default: 0.5]
  --methviz             write methylation-encoded sequence to VCF [default: False]
  --read-wise           read-wise genotyping for dense BED regions
  --locus-wise          locus-wise genotyping for sparse BED regions
  --decompose           write motif-decomposed sequence to VCF
  -t <INT>, --threads <INT>
                        number of threads [default: 1]
  -log, --debug_mode    write debug messages to log file
  -v, --version         show program's version number and exit
```

The details of required input files are given below:

### Reference genome
#### `-f or --fasta`
**Expects**: *FILE*<br>
**Default**: *None*<br>
The `-f` or `--fasta` option is used to specify the input FASTA file. The corresponding index file (`.fai`) should be in the same directory. ATaRVa uses [pysam](https://github.com/pysam-developers/pysam)'s `FastaFile` parser to read the input FASTA file.

### Alignment file
#### `-b or --bam`
**Expects**: *FILE*<br>
**Default**: *None*<br>
The `-b` or `--bam` option is used to specify one or more input alignment files in the same format. ATaRVa accepts any of the three alignment formats: SAM, BAM, or CRAM. The alignment file should be sorted by coordinates. The format should be specified using the `--format` option. The corresponding index file (`.bai` or `.csi`) should be located in the same directory. An alignment file can be sorted and indexed using the following commands:

```bash
# to sort the alignment file
$ samtools sort -o sorted_output.bam input.bam

# to generate .bai index file
$ samtools index -b sorted_output.bam
```

An alignment file containing at least one of the following tags is preferred for faster processing: `MD` tag, `CS` tag, or a `CIGAR` string with `=/X` operations.

- The CS tag is generated using the --cs option when aligning reads with the [minimap2](https://github.com/lh3/minimap2) aligner. (`--cs=short` is prefered over `--cs=long`)
- The MD tag can be generated using the --MD option in minimap2.

If the alignment files were generated without any of these tags, you can generate the `MD` tag by running the following command to 

```bash
# input: reference genome fasta file & alignment file
# output: an alignment file with MD tag in it

# for generating MD tag
$ samtools calmd -b aln.bam ref.fa > aln_md.bam
```
### Region file
#### `-r or --regions`
**Expects**: *FILE*<br>
**Default**: *None*<br>
The `-r` or `--regions` option is used to specify the input TR regions file. ATaRVa requires a sorted, bgzipped BED file of TR repeat regions, along with its corresponding tabix-indexed file. The BED file should contain the following columns:

1. Chromosome name where TR is located
2. Start position of the TR
3. End position of the TR
4. Repeat motif
5. Motif length

Below is an example of a repeat region BED file. **NOTE: The BED file should either have no header or a header that starts with `#` symbol. The .gz and .tbi files should be in same directory**

| #CHROM | START | END | MOTIF | MOTIF_LEN |
|--------|-------|-----|-------|-----------|
| chr1   | 10000 | 10467 | TAACCC | 6    |
| chr1   | 10481 | 10497 | GCCC | 4      |
| chr2   | 10005 | 10173 | CCCACACACCACA | 13 |
| chr2   | 10174 | 10604 | ACCCTA | 6    |
| chr17  | 60483 | 60491 | AGA    | 3    |

To sort, bgzip, and index the BED file, use the following commands:

#### Sort
```bash
# input: Unsorted bed file
# output: Sorted bed file

# Sorting the BED file using sort
$ sort -k1,1 -k2,2n input.bed > sorted_output.bed
# or using bedtools
$ bedtools sort -i input.bed > sorted_output.bed
```
#### Bgzip
```bash
# input: Sorted bed file
# output: bgzipped bed file

# To keep the original file unchanged and generate separate gz file
$ bgzip -c sorted_output.bed > sorted_output.bed.gz
# or to compress the original file; converts sorted_output.bed to sorted_output.bed.gz
$ bgzip sorted_output.bed
```
#### Index
```bash
# input: bgzipped bed file
# output: tabix indexed file (.tbi)

# install samtools to use tabix
$ tabix -p bed sorted_output.bed.gz
```
For detailed information on advanced genotyping options, refer to the [Advanced Commands and Usage](./docs/genotype_usage.md) documentation.

## `merge` command
`atarva merge` merges VCF files generated by ATaRVa. The tool is optimized to handle large datasets efficiently by reading and processing multiple files in small chunks, thereby avoiding excessive memory usage and ensuring fast, memory-efficient execution.

The tool requires the following inputs:
- BGZipped and tabix-indexed VCF files
- A BGZipped and tabix-indexed BED file specifying the regions of interest

The help message and available options can be accessed using

```bash
$ atarva merge -h
#  or
$ atarva merge --help
```
which gives the following output

```
usage: atarva merge [-h] -r <FILE> -i <FILE> [<FILE> ...] -f <FILE> [--contigs CONTIGS [CONTIGS ...]] [-o <STR>] [-p <INT>]

Required arguments:
  -r <FILE>, --regions <FILE>
                        input regions file. the regions file should be strictly in bgzipped tabix format. If the regions input file is in bed format. First sort it using bedtools. Compress it using
                        bgzip. Index the bgzipped file with tabix command from samtools package.
  -i <FILE> [<FILE> ...], --vcfs <FILE> [<FILE> ...]
                        text file containing paths to input vcf files to be merged. The text file should list each path on a separate line. The vcf files should be strictly in bgzipped tabix format. If
                        the vcfs input file is in vcf format. First sort it using bcftools. Compress it using bgzip. Index the bgzipped file with tabix command from samtools package.
  -f <FILE>, --fasta <FILE>
                        input reference fasta file. The file should be indexed.

Optional arguments:
  --contigs CONTIGS [CONTIGS ...]
                        contigs to get merged [chr1 chr12 chr22 ..]. If not mentioned every contigs in the region file will be merged.
  -o <STR>, --outname <STR>
                        name of the output file, output is in vcf format.
  -t <INT>, --threads <INT>
                        number of threads. [default: 1]
```
**NOTE: This will merge only those loci that are present in the input BED file** <br>
For detailed information on advanced merging options, refer to the [Tamatr](./docs/merge_usage.md) documentation.

## Changelog

### v0.7.1+ext0.01
* Reports PS (phase set) in the VCF for haplotagged BAM inputs.
* Fixed - Report SQ value as '.' for loci with no informative SNPs.

### v0.7.1+ext
* First release. Extended version of [ATaRVa](https://github.com/SowpatiLab/ATaRVa) with additional features as below.
* Default processing of the softclipped portions of the read to search for the presence of repeat loci from the catalogs by searching for the flanking sequences.
* Implements a slightly different version of read-haplogrouping based on the informative SNPs. The reads are heirarchically grouped based on the minimising the Qvalue penalty of the informative SNPs.
* Added a new option `--rna` to enable genotyping of tandem repeats from RNA-seq data. This option is designed parse the CIGAR string of the typicaly RNA-seq alignment which include introns and other splicing events encoded as N.
* Added a new option `--instability` which outputs an additional file with allele information from each read at a locus

## Contact
For queries or suggestions, please contact:

Akshay Kumar Avvaru - avvaruakshay at gmail dot com
