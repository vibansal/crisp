

## DESCRIPTION 

Many programs for processing VCF files (such as vcftools) do not handle VCF files for polyploid or pooled DNA samples. Therefore, we have written python programs to process 
VCF files output by CRISP. This directory contains scripts to analyze and post-process the VCF files generated using CRISP or other variant calling programs. 
All these programs take as input a VCF and output summary statistics (and a new VCF file). Note that the VCFstats\_extended.py and SingleVarTest.py should be run on the standard VCF file generated using convert\_pooled\_vcf.py

1. *convert_pooled_vcf.py* convert CRISP VCF (with allele counts) into standard VCF for pooled samples (0/0/0/0/0/0/0/1/1/1) 

2. *VCFstats.py*: calculate Ti/Tv (transition/transversion) ratio for SNPs and 3n/non-3n indel ratio 

3. *VCFstats\_extended.py*: statistics about allele frequencies and singletons for each pool 

4. *SingleVarTest.py*: single variant association analysis using Fisher's exact test 

5. Generate a VCF file with variants that overlap intervals in a bed file: This can be done using bedtools (http://bedtools.readthedocs.org/en/latest/)


### These scripts work with older versions of python (tested with v2.7) and will not work with new (3+) version of python.
