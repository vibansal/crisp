## analyze VCF file to calculate statistics on Ti/Tv and indel lengths and indels in homopolymers

#! /usr/bin/env python
import sys, os, glob, string, gzip

"""
Ti/Tv and indel statistics for all samples together 
Ti/Tv and indel statistics per pool | singletons, missing genotypes per pool 
allele frequency spectrum across all samples 
separate count for multi-allelic indels in homopolymer runs or microsatellite sequences, they are not counted with other indels 

"""
epsilon = 1e-8;
	
def print_INDEL_lengths(INDEL_lengths):
	ilist = [];
	for indel in INDEL_lengths.iterkeys(): ilist.append([indel,INDEL_lengths[indel]]); 
	ilist.sort();
	print 'deletions:',
	for i in xrange(len(ilist)-1,0,-1): 
		if ilist[i][0] < 0: print `ilist[i][0]` + ':' + `ilist[i][1]`, 
	print '\ninsertions:',
	for i in xrange(len(ilist)): 
		if ilist[i][0] > 0: print `ilist[i][0]` + ':' + `ilist[i][1]`, 
	print '\n';
	
def print_sample_stats(TsTvcounts,genotype_stats,pool_sizes):
	for i in xrange(len(TsTvcounts)): 
		sampleid = TsTvcounts[i][7]; Ts = float(TsTvcounts[i][0]); Tv = float(TsTvcounts[i][1]); 
		del_3n = float(TsTvcounts[i][3]); del_rest = float(TsTvcounts[i][4]); 
		ins_3n = float(TsTvcounts[i][5]); ins_rest = float(TsTvcounts[i][6]); 
		
		print 'SAMPLE',sampleid,'SNVs:',Ts+Tv,round(float(TsTvcounts[i][0])/(TsTvcounts[i][1]+epsilon),3),
		if pool_sizes[i] > 2: print 'poolsize',pool_sizes[i],'missing',round(float(genotype_stats[i][1])/genotype_stats[i][0],3),genotype_stats[i][2],'SING:',genotype_stats[i][3],'DOUB:',genotype_stats[i][5],
		else: print 'missing',round(float(genotype_stats[i][1])/genotype_stats[i][0],3),genotype_stats[i][2],genotype_stats[i][3],genotype_stats[i][4],

		if del_3n + del_rest + ins_3n + ins_rest > 0: print 'INDELS',del_3n+del_rest+ins_3n+ins_rest,'DEL:',del_3n,del_3n+del_rest,round(del_3n/(del_3n+del_rest+epsilon),3),'INS:',ins_3n,ins_3n+ins_rest,round(ins_3n/(ins_3n+ins_rest+epsilon),3);
		else: print '\n',

##################################################################################################################################

def variantstats(vfile1):
        print >>sys.stderr, "\ncalculating variant statistics for VCF file",vfile1;

	## global variant statistics 
        SNPs=0; indels=0; indels_3n=0; deletions=0; insertions=0; MNPs=0; varsites=0; Ts=0; Tv=0; AG =0; CT =0; AC =0; AT = 0; CG =0; GT = 0;
	multi_allelic_indels = 0; multi_allelic_variants = [0,0,0]; tri_allelic_SNVs =0;
	INDEL_lengths = {}; 

	if vfile1.endswith('.gz'): File1 = gzip.open(vfile1, 'r')
        elif vfile1 != '-' and vfile1 != 'sys.stdin': File1 = open(vfile1,'r');
        else: File1 = sys.stdin;

	# statistics for each pool/sample
	TsTvcounts = []; ## indexed by sample-ID/col
	genotype_stats = []; # singeltons, missing,etc
	pool_sizes = [];
	total_allele_count = 0; AFS = [0]; AF_hist = [0.0];
	
	############################################################################################################

	CHROM_line_found = 0; columns = -1; 

	st= [0,0,0,0];

        for line in File1:
                var = line.strip().split('\t');  
		if var[0] == '#CHROM': 
			CHROM_line_found = 1; 
			for i in xrange(9,len(var)): 
				TsTvcounts.append([0,0,0,0,0,0,0,var[i]]); ## Ts, Tv, indels, indels_3n, insertions,deletions per sample
				genotype_stats.append([0,0,0,0,0,0]); ## total_genotypes, missing, reference_homozygous, singletons, homozygotes
			columns = len(var); 
				
                if line[0] == '#' or len(var) < columns: continue;
		#if varsites > 20000: break;
 
		if CHROM_line_found ==0: 
			print >>sys.stderr, "VCF file does not have CHROM header line, please provide valid VCF file \n";
			sys.exit();
		
                varsites +=1;
                if varsites > 0 and varsites%10000==0: print >>sys.stderr, 'processed',varsites,'variant records';

                refallele = var[3]; altalleles =var[4].split(',');
		GENOTYPES = []; 
		for i in xrange(9,len(var)): 
			counts = [0,0];  ## missing, reference
			for j in xrange(len(altalleles)): counts.append(0); 
			genotype = var[i].split(':')[0].split('/'); # list = ['0,'0',,,,,,,'1'.....,'2'] 
			if varsites ==1: 
				pool_sizes.append(len(genotype)); ## calculate poolsize using length of genotype vector 
				total_allele_count += len(genotype);  
				for t in xrange(len(genotype)): AFS.append(0); AF_hist.append(0.0); 

			for g in genotype: 
				if g == '.': counts[0] +=1; 
				elif g == '0': counts[1] +=1; 
				else: counts[int(g)+1] +=1; 
			GENOTYPES.append(counts); 

			genotype_stats[i-9][0] +=1;
			if counts[0] > 0: genotype_stats[i-9][1] +=1;  
			elif counts[1] == len(genotype): genotype_stats[i-9][2] +=1; 
			elif counts[2] == len(genotype): genotype_stats[i-9][4] +=1; 
			elif counts[2] ==1 or (len(altalleles) ==2 and counts[2] ==0 and counts[3] ==1): genotype_stats[i-9][3] +=1; 
			elif counts[2] ==2: genotype_stats[i-9][5] +=1; 

		# update allele frequency spectrum
		if len(altalleles) ==1:
			ref = 0; alt = 0; missing =0; 
			for i in xrange(9,len(var)): missing += GENOTYPES[i-9][0];  ref += GENOTYPES[i-9][1]; alt += GENOTYPES[i-9][2]; 
			if missing ==0: AFS[alt] +=1; AFS[0] += 1; 
			if missing*10 < total_allele_count: allele_freq = float(alt)/(alt+ref); AF_hist[int(allele_freq*total_allele_count)] +=1; 
		

		## detect multi-allelic indels TAA -> TA,TAAA,T	
		if len(altalleles) >= 2 and len(refallele) != len(altalleles[0]) and len(refallele) != len(altalleles[1]): 
			multi_allelic_indels +=1; 
			continue;  # if we remove this statement, then multi-allelic indels also get counted among other 

		if len(altalleles) == 2: 
			if len(refallele) == len(altalleles[0]) and len(refallele) == len(altalleles[1]): tri_allelic_SNVs +=1; 
			multi_allelic_variants[0] +=1;  
		elif len(altalleles) == 3: 
			if len(refallele) == len(altalleles[0]) and len(refallele) == len(altalleles[1]) and len(refallele) == len(altalleles[2]): tri_allelic_SNVs +=1; 
			multi_allelic_variants[1] +=1;  
		elif len(altalleles) == 4: multi_allelic_variants[2] +=1;  

		## else
		############################################################################################################
		## vartype: 0 = Ts(SNP), 1 = Tv, 2 = MNP, 3 = DEL_3, 4 = DEL, 5 = INS_3, 6 = INS
		st[3] += len(altalleles);
                for i in xrange(len(altalleles)):
                        rl = len(refallele); al = len(altalleles[i]); vartype = -1; 
                        if rl == al and len(refallele) ==1: vartype = 0; st[2] +=1; 
			elif rl == al: vartype = 2; st[2] +=1; ## MNP
                        elif rl < al:
				st[0] +=1; 
                                if (al-rl)%3 ==0: vartype = 5;
				else: vartype = 6;
				try: INDEL_lengths[al-rl] +=1;
				except KeyError: INDEL_lengths[al-rl] =1;
                        elif rl > al:
				st[1] +=1; 
                                if (rl-al)%3 ==0: vartype = 3;
				else: vartype = 4;
				try: INDEL_lengths[al-rl] +=1;
				except KeyError: INDEL_lengths[al-rl] =1;

                        if rl == al: 
				if (refallele == 'A' and altalleles[i] == 'G') or (refallele == 'G' and altalleles[i] == 'A'): Ts +=1; AG +=1; vartype = 0;
				elif (refallele == 'C' and altalleles[i] == 'T') or (refallele == 'T' and altalleles[i] == 'C'): Ts +=1; CT +=1; vartype = 0

				elif (refallele == 'A' and altalleles[i] == 'C') or (refallele == 'C' and altalleles[i] == 'A'): Tv +=1; AC +=1; vartype = 1
				elif (refallele == 'A' and altalleles[i] == 'T') or (refallele == 'T' and altalleles[i] == 'A'): Tv +=1; AT +=1; vartype = 1
				elif (refallele == 'C' and altalleles[i] == 'G') or (refallele == 'G' and altalleles[i] == 'C'): Tv +=1; CG +=1; vartype = 1
				elif (refallele == 'G' and altalleles[i] == 'T') or (refallele == 'T' and altalleles[i] == 'G'): Tv +=1; GT +=1; vartype = 1
				else: Tv +=1; vartype = 1;
	
			if vartype ==0 or vartype ==1: SNPs +=1; 
			elif vartype ==2: MNPs +=1; 
			elif vartype == 3 or vartype ==4: indels +=1; deletions +=1; 
			elif vartype == 5 or vartype ==6: indels +=1; insertions +=1; 
			if vartype == 3 or vartype ==5 : indels_3n +=1; 
			
				
			for j in xrange(9,len(var)): 
				counts = GENOTYPES[j-9]; 
				if counts[i+2] > 0 and vartype >=0: TsTvcounts[j-9][vartype] +=1; 
				
		############################################################################################################
			

        if vfile1 != '-' and vfile1 != 'sys.stdin': File1.close();

        print 'VARSITES:%d SNPs:%d INDELS:%d insertions:%d deletions:%d 3n-indels:%d:%0.3f MNPs:%d' %(varsites,SNPs,indels,insertions,deletions,indels_3n,float(indels_3n)/(indels+0.001),MNPs);
	print 'Multi-allelic variants',multi_allelic_variants[0],multi_allelic_variants[1],multi_allelic_variants[2],'SNVs',tri_allelic_SNVs,'Indels:',multi_allelic_indels;
        if Ts+Tv >0: print 'Ts/Tv: %0.3f %d/%d Transitions: A-G:%d C-T:%d Transversions: A-C:%d A-T:%d C-G:%d G-T:%d\n' %(float(Ts)/Tv,Ts,Tv,AG,CT,AC,AT,CG,GT);
	else: print;

	print >>sys.stderr, "stats",st;

	print_INDEL_lengths(INDEL_lengths);

	print 'AFS',AFS[0],
	for i in xrange(1,min(100,total_allele_count),1): 
		f = float(AFS[i])/AFS[0]
		print `i` + ':' + `AFS[i]` + ':' + `round(f,3)`,
	#for i in xrange(1,min(100,total_allele_count),1): print `i` + ':' + `AFS[i]` + ':' + `AF_hist[i]`,
	print '\n';

	## print table of statistics per sample
	print "Variant Statistics for each sample/pool in VCF file\n";
	print_sample_stats(TsTvcounts,genotype_stats,pool_sizes);

variantstats(sys.argv[1]);

