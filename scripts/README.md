

## DESCRIPTION 

Many programs for processing VCF files (such as vcftools) do not handle VCF files for polyploid or pooled DNA samples. Therefore, we have written python programs to process 
VCF files output by CRISP. This directory contains scripts to analyze and post-process the VCF files generated using CRISP or other variant calling programs. 
All these programs take as input a VCF and output summary statistics (and a new VCF file) 

1. convert_pooled_vcf.py convert CRISP VCF (with allele counts) into standard VCF for pooled samples (0/0/0/0/0/0/0/1/1/1) 

2. VCFstats.py: calculate Ti/Tv (transition/transversion) ratio for SNPs and 3n/non-3n indel ratio 

3. VCFstats\_extended.py statistics about allele frequencies and singletons for each pool 

4. SingleVarTest.py single variant association analysis using Fisher's exact test 

5. Generate a VCF file with variants that overlap intervals in a bed file: This can be done using bedtools (http://bedtools.readthedocs.org/en/latest/)


## TO BE ADDED SOON

A. filter variants using LR statistic and chi-square table statistic (user-specified thresholds)

B. INDELS: calculate fraction of indels in homopolymer runs and identify regions with clusters of variants in low-complexity sequence 

C. identify somatic mutations in tumor-normal sequencing using list of tumor-normal pairs 

D. split tri-allelic variants into multiple variants for ease of annotation and rare variant association testing

