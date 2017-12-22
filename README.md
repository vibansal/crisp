####################################################################################################

CRISP: Comprehensive Read Analysis for Identification of SNVs (and short indels) from Pooled sequencing data

Contact: Vikas Bansal, vibansal@ucsd.edu

Citation: A statistical method for the detection of variants from next-generation resequencing of DNA pools, V. Bansal. Bioinformatics 26(12), 2010

######################################################################################################

Introduction:
=============

CRISP is a software program to detect SNPs and short indels from pooled sequencing data. CRISP is designed to detect both rare and common variants by utilizing sequence reads from next-generation sequencing of multiple DNA pools. CRISP uses a cross-pool comparison approach to distinguish sequencing errors from rare variants. Note that the method is not designed for variant detection from a single pool. CRISP has been evaluated on several pooled sequencing datasets (human and bacterial) generated using the Illumina sequencing platform. In principle, it should work for sequence data from other sequencing platforms. The method requires each pool to be sequenced using the same sequencing platform. 


Implementation
=============

CRISP is implemented in C and uses the SAMtools API for reading BAM files. 


Running CRISP:
=============


The first step is to align the reads from the sequencing experiment to the reference genome using BWA (or your favorite aligner) and generate sorted BAM file for each pool. To improve detection of indels, it is recommended to run an indel realigner program (such as GATK) to realign reads and generate a realigned BAM file. Once a single bam file for each pool is generated, CRISP can be run on the BAM files to detect variants by comparing against a reference sequence. 

./CRISP [options] --bams file_bam_paths --ref reference.fasta --VCF variantcalls.VCF --poolsize poolsize --bed targets.bed > variantcalls.log


Important Notes:

1. CRISP requires at least four arguments: poolsize, reference fasta file, bamfiles and the output VCF file
2. BAM files can be specified in two ways: 
	(a) --bams: list of all bam files in a single text file
	(b) --bam: specify individual bam files on the command line 
3. The reference sequence file should be indexed using 'samtools faidx' and the index file reference.fasta.fai placed in the same directory as reference.fasta
4. For targeted sequencing studies, it is recommended to use a bed file for variant calling. If a bedfile is not specified, the program will evaluate each base in the genome for variant calling.

4. CRISP requires at least two pools to make variant calls, but at least 5 pools are ideal. CRISP should be run separately on pools sequenced on different sequencing instruments.

5. The number of haplotypes in each pool (--poolsize) is assumed to be the same. For variable poolsizes, the poolsize should be specified inn the input file with the list of bam files (see FAQ for details on the format).

6. Increasing the number of permutations (from 20K to 50K-100K) can slightly improve the accuracy of variant detection at the cost of increasing running time. Choose this parameter based on the number of pools and the total target region sequenced.

7. CRISP can be parallelized by calling variants on specific chromosomes (or regions) using the --regions option. This requires the bam files to be indexed. 


Command-line arguments for CRISP
================================

         --bams         textfile with list of bam file paths (one for each pool)
         --bam          bam file for one pool, specify filename for each pool using --bam pool1.bam --bam pool2.bam .... --bam pooln.bam
         --ref       	Indexed Reference Sequence file (fasta)
	 --bed		bedfile with intervals in which variants will be called
         --poolsize     poolsize (number of haploid genomes in each pool), for diploid genomes: 2 x # individuals
         --VCF       	VCF file to which the variant calls will be output 
         --qvoffset  	quality value offset, 33 for Sanger format, 64 for Illumina 1.3+ format
         --mbq       	minimum base quality to consider a base for variant calling, default 10
         --mmq       	minimum read mapping quality to consider a read for variant calling, default 20
         --regions      region(s) in which variants will be called, e.g chr1:654432-763332. BAM files should be indexed for this option with the index file pooln.bam.bai 
         --minc      	minimum number of reads with alternate allele required for calling a variant, default 4
         --ctpval    	threshold on the contingency table p-value for calling position as variant (specified as log10), default is -3.5, increase this threshold as number of pools increases
         --qvpval    	threshold on the quality values based p-value for calling position as variant (specified as log10), default is -5
         --perms     	maximum number of permutations for calculating contingency table p-value, default 20,000
	 --filterreads  filter reads with excessive number of mismatches (and gaps) compared to the reference sequence, Default is 1. Set to 0 to disable filtering 



Output VCF file:
================

CRISP outputs the variants identified to a VCF file. For each variant, CRISP outputs the allele depth (read counts) for the reference and variant alleles at each poisition in the VCF file. This can be used for calculating allele frequencies or doing case-control association analysis. Multi-allelic variants and indels are also reported.

Update: the latest version of CRISP outputs genotypes for each pool (--EM 1, default) as allele counts. The CRISP VCF can be converted to a standard pooled genotype VCF using a script provided in the scripts sub-directory. 

Description of fields
=====================

NP: Number of Pools With Data
VP: Number of Pools with variant allele(s)
DP: Total number of reads (+strand,-strand) across all pools (filtered reads only)
CT: contingency table p-value for each variant allele in same order as listed in column 5
QVpf (QVpr): quality values based p-value for each variant allele using forward(reverse) strand reads 
MQ: number of reads with mapping qualities <10 | 10-19 |  20-39 |  >=40 
VT: variant type, SNV | DELETION | INSERTION 
HP: length in the ambiguity of positioning of indels (homopolymer length or microsatellite length)
FLANKSEQ: This represents the reference haplotype sequence (length spans the homopolymer or microsatellite tract) with 10 bases either side of the variant position (10bases_upstream:reference_haplotype_sequence:10bases_downstream). This is useful to eyeball indels that occur in long homopolymer or microsatellite tracts.

SB: variant demonstrates strand bias
PASS: variant passes all filters 

AF: frequency of variant alleles(s) in the pool in order listed in column 5
ADf: Number of reads aligned to the forward strand of the genome supporting reference allele and the alternate alleles in the order listed
ADr: Number of reads aligned to the reverse strand of the genome supporting reference allele and the alternate alleles in the order listed

Please note that some of these fields may be removed/updated in future releases. 








