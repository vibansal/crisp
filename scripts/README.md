This directory contains scripts to analyze and post-process the VCF files generated using CRISP or other variant 
calling programs. All these programs take as input a VCF and output summary statistics of a new VCF

1. calculate Ti/Tv (transition/transversion) ratio for SNPs and 3n/non-3n indel ratio 
2. calculate fraction of indels in homopolymer runs...
3. for pooled sequencing: # of singletons in each pool 
4. merge multi-base SNPs (MNPs) and complex variants into a single variant (using proximity and genotype analysis)
5. identify regions with clusters of variants in low-complexity sequence 
6. convert CRISP VCF (with allele counts) into standard VCF for pooled samples (0/0/0/0/0/0/0/1/1/1) 
