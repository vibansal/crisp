#! /usr/bin/env python
# AUTHOR VIKAS BANSAL 
### code for doing single variant association tests for pooled sequencing data ##

#python xx.py CRISP.genes/variants.SRR.vcf ../Phenotypes-pools/192pools.phenotypes 40

## code modified to handle variable pool sizes for phase3 data, june 16 2015 

import sys, os, glob, string, subprocess,time, math, re, compiler, random
#from scipy.stats import chi2  # chi square test p-value 
from combinatorial import fet

PRINTGENOTYPES =0;
USE_MEAN_AC = 0;
epsilon = 1e-6

## read the phenotype file with poolID (G1P54), 2,3 etc
def read_phenotypefile(statusfile):
	## read the case control status file  ### 
	poolDX = {}; casepools = 0; controlpools = 0; totalpools =0; earlyonset=0; DC=0;
	File = open(statusfile);
	for line in File: 
		## allow multiple phenotype codes for each pool separated by commalist 
		pool = line.strip().split(); phenotype = pool[1].split(',');  
		poolDX[pool[0]] = [int(phenotype[0]),int(pool[2])]; ## 2nd var is poolsize

		if int(phenotype[0]) == 0: controlpools +=1; 
		if int(phenotype[0]) >= 1: casepools +=1; 
		if int(phenotype[0]) == 2: earlyonset +=1; 
		if len(phenotype) >=2 and int(phenotype[1]) == 3: DC +=1; 
		totalpools +=1; 
	File.close();
	print >>sys.stderr, "totalpools",totalpools,"CASEpools",casepools,"controlpools",controlpools,'earlyonset',earlyonset,'T2D-complication',DC;
	return poolDX; 


def determine_poolphenotype(vcffile,poolstatustable):
	poolDX = [];  # phenotype status for each pool
	File = open(vcffile);
	for line in File:
		if line[0] == '#' and line[1] == '#': continue;
		variant = line.strip().split('\t');
		if variant[0] == '#CHROM': 
			for i in xrange(9,len(variant)):
				sampleid = variant[i].split('/')[-1].split('.')[0];
				#print sampleid;
				try: status = poolstatustable[sampleid]; 
				except KeyError: status= [-1]; 
				poolDX.append(status);
				print >>sys.stderr, status,
			#samples = len(variant)-9; print 'samples',samples,len(variant);
			print >>sys.stderr;
			continue;
	File.close();
	return poolDX; 



############################################################################################
# this uses the VCF format where we have allele counts (ML) instead of readcounts
# ignore multi-allelic variants for now
def calculate_FETpvalues_allelecounts(vcffile,poolstatustable):
	poolDX = determine_poolphenotype(vcffile,poolstatustable);  # phenotype status for each pool
	variants =0; trivariants =0;

	File = open(vcffile);
	for line in File:
		if line[0] == '#': continue;
		variant = line.strip().split('\t');
		samples = len(variant)-9;
		chrom = variant[0]; position = int(variant[1]); refallele = variant[3]; varalleles = variant[4].split(','); 

		if len(varalleles) ==2: triallelic=1; trivariants +=1; 
		else: triallelic=0; variants +=1; 

		if len(varalleles) >= 3 or len(varalleles) ==4: 
			#print >>sys.stderr, '##triallelic',chrom,position,refallele,varalleles,variant[2];
			continue; # ignore multi-allelic variants for now

		## healthy (controls), D (disease), E (early onset), DC (diabetic complications)  
		H0 = 0.0; H1 = 0.0; H2 = 0.0; D0 = 0.0; D1 = 0.0; D2 = 0.0; E0 = 0.0; E1= 0.0; E2 = 0; 
		casepools=0; controlpools=0;
		for i in xrange(samples):
			if poolDX[i][0] == -1 or variant[i+9].split(':')[0] == '.': continue;  # ignore pool for case control analysis 
			try: 
				poolsize = poolDX[i][1]
				genotypes = variant[i+9].split(':');
				if triallelic==0: 
					MLAC = int(genotypes[0].split(',')[1]); #meanAC = float(genotypes[2]); QAC = float(genotypes[1]); 
					QAC = float(genotypes[1]); 
					MLAC2 = 0;
				else: 
					MLAC = int(genotypes[0].split(',')[1]); #meanAC = float(genotypes[2]); QAC = float(genotypes[1]); 
					QAC = float(genotypes[1]); 
					MLAC2 = int(genotypes[0].split(',')[2]); #meanAC = float(genotypes[2]); QAC = float(genotypes[1]); 

				#else: varAF = float(genotypes[3]);
 
				if poolDX[i][0] == 0: H0 += poolsize; H1 += MLAC; H2 += MLAC2;
				elif poolDX[i][0] >= 1:  D0 +=poolsize; D1 += MLAC; D2 += MLAC2;
				if poolDX[i][0] == 2: E0 += poolsize; E1 += MLAC; E2 += MLAC2;

			except IndexError: print 'Exception',i,samples,genotypes; 


		if (H1 + D1 >= 4): pvalue = fet(int(round(H1,0)),int(round(H0,0)),int(round(D1,0)),int(round(D0,0)));
		else: pvalue = [0,0,0];

		## calculate p-values between early onset and controls
		if (H1+E1 >=4): pvalue1 = fet(int(round(H1,0)),int(round(H0,0)),int(round(E1,0)),int(round(E0,0)));
		else: pvalue1 = [0,0,0]; 

		## calculate p-value between early onset and late-onset
		if (D1+E1 >=4): pvalue2 = fet(int(round(D1,0)),int(round(D0,0)),int(round(E1,0)),int(round(E0,0)));
		else: pvalue2 = [0,0,0]; 
		

		print '%0.2f %.1f/%.1f %.1f/%.1f %.1f/%.1f %.1f/%.1f'  %(pvalue[0],H1,H0,D1,D0,E1,E0,DC1,DC0),#'Control',float(H1)/(H0+0.001),'Case',float(D1)/(D0+0.001),
		print '%0.2f %0.2f' %(pvalue1[0],pvalue2[0]),
                print '%0.4f %0.4f %0.4f %0.4f' %(H1/(H0+epsilon),D1/(D0+epsilon),E1/(E0+epsilon),DC1/(DC0+epsilon)),
		#sys.exit()

		if triallelic ==1: 
			if H2 + D2 >=4: pvalue_tri = fet(int(round(H2,0)),int(round(H0,0)),int(round(D2,0)),int(round(D0,0)));
			else: pvalue_tri = [0,0,0];
			print 'TRIALLELIC:%0.2f:%.1f:%.1f'  %(pvalue_tri[0],H2,D2),
		else: print '-',

		print '%3s %9s %s %10s %10s %8s %10s' %(variant[0],variant[1],variant[2],variant[3],variant[4],variant[5],variant[6]),
		print variant[7], 
		print;
	print >>sys.stderr, "variants evaluated",variants,"triallelic or more",trivariants;

############################################################################################
if len(sys.argv) < 3: print >>sys.stderr, "python pooled-association.py VCFfile phenotypefile"; sys.exit();
## poolsize is inferred from VCF file itself, works on new VCF format 0/0/0...../1

vcffile=sys.argv[1];  statusfile = sys.argv[2]; 
poolDX = read_phenotypefile(statusfile); 
calculate_FETpvalues_allelecounts(vcffile,poolDX);



############################################################################################

	
	



