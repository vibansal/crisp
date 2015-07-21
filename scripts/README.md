This directory contains scripts to analyze and post-process the VCF files generated using CRISP or other variant 
calling programs. All these programs take as input a VCF and output summary statistics of a new VCF


VCFstats.py: calculate Ti/Tv (transition/transversion) ratio for SNPs and 3n/non-3n indel ratio 

convert_pooled_vcf.py convert CRISP VCF (with allele counts) into standard VCF for pooled samples (0/0/0/0/0/0/0/1/1/1) 




## TO BE ADDED SOON ### 

0. filter variants using LR statistic and chi-square table statistic (determine thresholds)
2. calculate fraction of indels in homopolymer runs...
3. for pooled sequencing: # of singletons in each pool 
4. merge multi-base SNPs (MNPs) and complex variants into a single variant (using proximity and genotype analysis)
5. split tri-allelic variants into multiple variants for ease of annotation and rare variant association testing 
7. identify regions with clusters of variants in low-complexity sequence 
8. generate a VCF file with variants that overlap intervals in a bed file (for exome analysis) 
9. identify somatic mutations in tumor-normal sequencing: VCF file + list of tumor-normal pairs 

