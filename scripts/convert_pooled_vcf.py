#! /usr/bin/env python
import sys, os, glob, string, subprocess,time, math, re, compiler

## coded july 15 2015, code to convert CRISP pooled allele count to genotype 0/0/0/0/0/1 format so that it is compatible with standard tools

## if the sum of allele counts is more than poolsize, program ignores that variant. 

if len(sys.argv) < 3: 
	print >>sys.stderr, "\n\nPROGRAM to convert CRISP VCF for pooled samples to a standard VCF where genotypes are in the form 0/0/0/0/0/1"; 
	print >>sys.stderr, "\nHOW TO RUN: python convert_pooled_VCF.py VCF_file_from_CRISP.vcf bam_file_list poolsize(optional) > NEW_VCF_file.vcf "; 
	print >>sys.stderr, "\nNOTE: poolsize should be specified as the third argument when all pools have the same pool size\n\n";
	sys.exit();
if len(sys.argv) > 3: 
	poolsize = int(sys.argv[3]); 
else: poolsize = -1; 
	
bam_file_list = sys.argv[2]; 
File = open(bam_file_list);
pool_sizes = []; 
for line in File: 
	pool = line.strip().split();
	if len(pool) >=2 and 'PS=' in pool[1]: pool_sizes.append(int(pool[1].split('=')[1])); 
	elif poolsize < 0: print >>sys.stderr, "pool size information is not available, it should be provided on the command line (same pool size for all pools) or using the bam_file_list file"; sys.exit(); 
File.close(); 
print >>sys.stderr, "pool sizes:", pool_sizes; #sys.exit();

variants=0; 
File = open(sys.argv[1]); # pooled VCF file
 
for line in File:

	if line[0] == '#': 
		print line,
		continue;

	var = line.strip().split('\t'); chrom = var[0]; position = int(var[1]); alleles = var[4].split(',');

	valid_variant = 1;
	GENOTYPES = []; 
	for i in xrange(9,len(var)):
		genotype = var[i].split(':'); counts = genotype[0].split(','); 

		if poolsize > 0: 
			refsum = poolsize; GV = [];

			if genotype[0] == '.': 
				for c in xrange(poolsize): GV.append('.');
			else: 
				for c in xrange(len(counts)): refsum -= int(counts[c]); 
				if refsum < 0: 
					print >>sys.stderr, "poolsize is smaller than allele counts, please provide correct pool size, terminating program";
					print >>sys.stderr, refsum,var[i],var[0:5]
					valid_variant = 0; break; 
					sys.exit();
				for c in xrange(refsum): GV.append('0'); 
				for a in xrange(len(counts)): 
					for c in xrange(int(counts[a])):  GV.append(`(a+1)`); 
			#print 'genotype','/'.join(GV);					
		else: ## variable pool size so genotype is different 
			GV = [];
			if genotype[0] == '.': 
				for c in xrange(pool_sizes[i-9]): GV.append('.');
				#print >>sys.stderr,'ERROR',genotype
				#pass;
				## what should we do, need input bam file list 
			else: 
				for a in xrange(len(counts)): 
					for c in xrange(int(counts[a])):  GV.append(`(a)`); 
				if len(GV) > pool_sizes[i-9]: 
					valid_variant = 0; break; 
					print >>sys.stderr, "poolsize is smaller than allele counts, please provide correct pool size, terminating program"; 
					sys.exit();
			
		#sys.stdout.write('\t' + '/'.join(GV) + ':' +  ':'.join(genotype[1:]));
		GENOTYPES.append('/'.join(GV) + ':' +  ':'.join(genotype[1:]));
	
	if valid_variant ==0: print >>sys.stderr, "filtering variant due to alleles",var[0:5]; continue;
	col8 = var[8].split(':')[1:]; new_col8 = 'GT:' + ':'.join(col8); newVec = var[0:8]; newVec.append(new_col8); 

	print "\t".join(newVec + GENOTYPES).strip() 
	#sys.stdout.write(var[8]);
	#print >>sys.stderr, newVec + GENOTYPES
	#print "\t".join(GENOTYPES).strip() 
	#sys.stdout.write('\n')

File.close();

