# CRISP: Comprehensive Read Analysis for Identification of SNVs (and short indels) from Pooled sequencing data

Contact: Vikas Bansal, vibansal@ucsd.edu

---

## Introduction

CRISP is a software program to detect SNPs and short indels from pooled sequencing data. CRISP is designed to detect both rare and common variants by utilizing sequence reads from next-generation sequencing of multiple DNA pools. CRISP uses a cross-pool comparison approach to distinguish sequencing errors from rare variants. Note that the method is not designed for variant detection from a single pool. CRISP has been evaluated on several pooled sequencing datasets (human and bacterial) generated using the Illumina sequencing platform. In principle, it should work for sequence data from other sequencing platforms. The method requires each pool to be sequenced using the same sequencing platform.


## Installation

CRISP is implemented in C and uses **htslib** for reading BAM and CRAM files. Install htslib before compiling:

- **Debian/Ubuntu**: `sudo apt-get install libhts-dev`
- **macOS (Homebrew)**: `brew install htslib`
- **From source**: https://github.com/samtools/htslib

The reference FASTA must be indexed with `samtools faidx` (or equivalent) before running CRISP.

Once htslib is installed, type `make` in the source directory. The compiled executable `CRISP.binary` will be in the `bin/` subfolder. CRISP has been tested on Debian & Ubuntu Linux.


## Running CRISP

The first step is to align the reads from the sequencing experiment to the reference genome (using BWA or similar) and generate a sorted, indexed BAM or CRAM file for each pool. To improve detection of indels, it is recommended to run an indel realigner (such as GATK) to realign reads. Once a single BAM or CRAM file per pool is ready, CRISP can be run to detect variants.

```
./CRISP [options] --bams file_bam_paths --ref reference.fasta --VCF variantcalls.VCF -p poolsize --bed targets.bed > variantcalls.log
```

**Important Notes:**

1. CRISP requires at least four arguments: poolsize (`-p`), reference FASTA (`--ref`), BAM/CRAM file list (`--bams`), and output VCF (`--VCF`).
2. Both BAM and CRAM input files are supported. **CRAM files require `--ref`** to specify the reference FASTA used during alignment.
3. BAM/CRAM files can be specified in two ways:
   - `--bams`: a text file with one BAM/CRAM path per line
   - `--bam`: specify individual files on the command line: `--bam pool1.bam --bam pool2.bam ...`
4. The reference FASTA must be indexed (`samtools faidx reference.fasta`); the `.fai` index must be in the same directory. CRISP validates that contig names and lengths in the BAM/CRAM header match the reference FASTA.
5. For targeted sequencing studies, use `--bed` to restrict variant calling to target regions. Without a bed file, every position in the genome is evaluated (large output).
6. CRISP requires at least two pools to make variant calls; at least 5 pools are recommended. Run CRISP separately for pools sequenced on different instruments.
7. The number of haplotypes per pool (`--poolsize`) is assumed to be the same across pools. For variable pool sizes, specify the pool size per file in the `--bams` input file (see FAQ).
8. CRISP can be parallelized by chromosome/region using `--regions`. BAM/CRAM files must be indexed (`.bam.bai` / `.cram.crai`).


## Command-line arguments

```
--bams                  textfile with list of BAM/CRAM file paths (one per pool)
--bam                   BAM/CRAM file for one pool; repeat for each pool: --bam pool1.bam --bam pool2.bam ...
--ref / --reference     indexed reference FASTA file (required; mandatory for CRAM input)
--bed / --bedfile       BED file with intervals for variant calling
-p / --poolsize         poolsize (number of haploid genomes per pool); for diploid samples: 2 × number of individuals
--VCF                   output VCF file for variant calls
--sm <0/1>              use @RG SM tag from BAM/CRAM header as VCF sample name (default 0: derive from filename)
--mbq <int>             minimum base quality for variant calling, default 10
--mmq <int>             minimum read mapping quality for variant calling, default 20
--regions               region for variant calling, e.g. chr1:654432-763332 (requires indexed BAM/CRAM)
--minc <int>            minimum reads with alternate allele required to call a variant, default 4
--ctpval <float>        contingency table p-value threshold (log10), default -3.5; increase as pool count grows
--filterreads <0/1>     filter reads with excessive mismatches vs reference, default 1 (set 0 to disable)
--EM <0/1>              genotype calling method: 1 = EM algorithm (default), 0 = allele frequency method
--OPE <0/1>             identify overlapping paired-end reads, treat as single read in overlap, default 1
--mincov <int>          minimum average coverage per haploid genome to call variants at a site
--flankingbases <int>   call variants in regions flanking bed targets (use 50–100 for targeted sequencing), default 0
--verbose / --pflag <0/1/2>   output verbosity: 0 = none, 1 = medium (default), 2 = detailed
```


## Output VCF file

CRISP outputs variants to a VCF file. For each variant, allele depths (read counts) for the reference and variant alleles are reported per pool/sample. Multi-allelic variants and indels are included.

The default output (`--EM 1`) reports pool genotypes as allele counts. A script in the `scripts/` subdirectory converts the CRISP VCF to standard pooled genotype VCF format.


## INFO fields

| Field | Description |
|-------|-------------|
| NP | Number of pools with data at this position |
| VP | Number of pools with the variant allele(s) |
| DP | Total read counts (+strand, -strand, bidirectional) across all pools (filtered reads only) |
| CT | Contingency table p-value (log10) for each variant allele, in the order listed in ALT |
| AC | Estimated allele count (across all pools) for each variant allele |
| AF | Estimated allele frequency (EM algorithm) for each variant allele |
| VT | Variant type: SNV, DELETION, or INSERTION |
| VF | EMpass or EMfail; EMfail variants should be discarded |
| EMstats | Log-likelihood statistics from the EM algorithm |
| HWEstats | Hardy-Weinberg equilibrium p-value across all pools |
| MQS | Read counts by mapping quality bin: 0–9, 10–19, 20–39, 40–255 |
| HP | Homopolymer run length (or microsatellite length) at the variant position; indicates ambiguity in indel positioning |
| FLANKSEQ | Flanking sequence context — 10 bases upstream : reference haplotype : 10 bases downstream. Useful for inspecting indels in repetitive sequence. |


## FORMAT fields

| Field | Description |
|-------|-------------|
| GT | Genotype |
| DP | Read depth per pool/sample (all reads) |
| GQ | Genotype quality |
| AC | Maximum-likelihood allele count estimate for ALT allele(s), per pool |
| AF | Variant allele frequency in pool |
| ADf | Forward-strand read counts supporting REF and each ALT allele |
| ADr | Reverse-strand read counts supporting REF and each ALT allele |
| ADb | Overlapping paired-end read counts supporting REF and each ALT allele |


## FILTER values

| Value | Description |
|-------|-------------|
| PASS | Variant passes all filters |
| StrandBias | Strand bias in read distribution between REF and ALT alleles |
| LowDepth | Low average coverage (< 1 filtered read per haplotype across all samples) |
| LowMQ20 | > 20% of reads have mapping quality < 20 |
| LowMQ10 | > 10% of reads have mapping quality < 20 |


## References

- V. Bansal. A statistical method for the detection of variants from next-generation resequencing of DNA pools. *Bioinformatics* 26(12), 2010.
- V. Bansal et al. Efficient and Cost Effective Population Resequencing by Pooling and In-Solution Hybridization. *PLoS ONE*, 2011.
