

** CRISP can run in two modes: --EM 0 and --EM 1 option **

--EM 1 corresponds to the new CRISP method that can estimate the 'discrete' pooled genotype for each pool
and uses a likelihood ratio test (in addition to the contingency table analysis) to call variants. This is 
recommended for artificial pooled sequencing studies. 


--EM 0 corresponds to the original CRISP method that uses the contingency table analysis and per pool filters
to identify variant sites and report allele frequencies for each pool. The ctpval and qvpval 
thresholds are only used with EM = 0 option


**running CRISP on targeted sequencing experiments**

It is recommended to specify a bedfile (--bed option) for targeted sequencing experiments. If a bedfile is 
not specified, the program will evaluate each base for variant calling and the output file can be quite
large. If the bedfile is specified, CRISP will only call variants in the regions listed in the bedfile. To call variants in the
region flanking the targeted intervals, --flanking option can be used. With "--flanking 50" option, 
variants will also be called in a 50 bp region flanking every interval. 

***Some requirements before running CRISP***

1. The reference sequence file should be indexed using 'samtools faidx' or a similar 
program and placed in same directory as fasta file with extension .fai
2. The aligned reads for each pool should be in a single bam file that 
is sorted by chromosomal coordinates. If there are multiple bam files for each pool, merge them
using samtools or Picard. 
3. BAM files should be indexed in order to use the --regions option with the 
indexed bam file named as pooln.bam.bai (pooln.bai will not work)
4. Please make sure that the reference sequence file is the same as the one used to align 
the reads in the BAM files


***calling of short insertions/deletions (indels) ***

CRISP can call short indels as well as SNVs. However, it is difficult to call multi-allelic indels (
i.e. indels with multiple variant alleles) from pooled sequence data. For indel analysis, CRISP assumes that indels are left justified, 
--leftalign 1 option can be used to left justify gaps in aligned reads


*** Common sources of Errors while running CRISP ***

If you get a segmentation fault while running CRISP, try to compile CRISP from the source code 
on your system.


