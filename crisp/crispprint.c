
#include<time.h>
/* this is now the code for original CRISP printing, no allele count estimation */
/* functions to print CRISP program options, VCF header and variant in VCF format */

// also calculate allele frequency of variant alleles across all pools combined (assumption is that pool size is constant)

// output variant calls and allele counts to VCF file (this is done once for each variant, even if it is multi-allelic) 
int print_pooledvariant(struct VARIANT* variant,FILE* vfile,REFLIST* reflist,int current,ALLELE* maxvec)
{
	variant->previousbase = 'N';
        if (current >=0 && variant->position >=1)  variant->previousbase = toupper(reflist->sequences[current][variant->position-1]);
        if (variant->previousbase != 'A' && variant->previousbase != 'C'  && variant->previousbase != 'G' && variant->previousbase != 'T') variant->previousbase = 'N';
        // one base prior to current base, used for representing indels in VCF file

	int i=0,j=0,k=0; int depthf=0,depthr = 0; int delallele = 0, insallele = 0, subs =0, longestdel=0,coverage=0;
	int printhomseq=0;
	// need to account for both deletions and insertions at same locus +A, -A  which should be output as AA A,AAA
	for (i=0;i<variant->varalleles;i++) 
	{ 
		if (variant->itb[variant->alleles[i]][0] == '-') delallele++; 
		else if (variant->itb[variant->alleles[i]][0] == '+') insallele++; 
		else subs++;
		if (strlen(variant->itb[variant->alleles[i]]) > longestdel && variant->itb[variant->alleles[i]][0] == '-') 
		{ 
			k = i; longestdel = strlen(variant->itb[variant->alleles[i]]); 
		} 
	} 
	if (insallele +delallele ==0)
	{
		fprintf(vfile,"%s\t%d\t.\t%c\t%s",variant->chrom,variant->position+1,variant->refb,variant->itb[variant->alleles[0]]);
		strcpy(variant->type,"SNV");
		for (i=1;i<variant->varalleles;i++) fprintf(vfile,",%s",variant->itb[variant->alleles[i]]);
	}
	else if (delallele > 0 && variant->previousbase != '0') 
	{
		if (insallele ==0 && subs==0) strcpy(variant->type,"DELETION");
		else if (subs > 0 && insallele ==0) strcpy(variant->type,"SNV,DELETION");
		else if (insallele > 0 && subs==0) strcpy(variant->type,"DELETION,INSERTION");
		else strcpy(variant->type,"SNV,DELETION,INSERTION");

		fprintf(vfile,"%s\t%d\t.\t%c",variant->chrom,variant->position,variant->previousbase); 
		for (j=1;j<longestdel;j++) fprintf(vfile,"%c",variant->itb[variant->alleles[k]][j]); fprintf(vfile,"\t");

		// print alt allele 
		for (i=0;i<variant->varalleles;i++)
		{
			fprintf(vfile,"%c",variant->previousbase);
			if (variant->itb[variant->alleles[i]][0] == '+')// insertion 
			{
				for (j=1;j<strlen(variant->itb[variant->alleles[i]]);j++) fprintf(vfile,"%c",variant->itb[variant->alleles[i]][j]);
				for (j=1;j<longestdel;j++) fprintf(vfile,"%c",variant->itb[variant->alleles[k]][j]);
			}
			else if (variant->itb[variant->alleles[i]][0] == '-')  // deletion 
			{
				// bug fixed, the longest allele needs to be used to print the variant allele for deletions March 23 2012
				for (j=strlen(variant->itb[variant->alleles[i]]);j<longestdel;j++) fprintf(vfile,"%c",variant->itb[variant->alleles[k]][j]); 
			}
			else // substitution
			{
				// don't print the reference base
				fprintf(vfile,"%c",variant->itb[variant->alleles[i]][0]);
				for (j=2;j<longestdel;j++) fprintf(vfile,"%c",variant->itb[variant->alleles[k]][j]);
			}
			if (i < variant->varalleles-1) fprintf(vfile,",");
		}
	}
	else if (variant->previousbase != '0') // flag for no reference sequence 
	{
		fprintf(vfile,"%s\t%d\t.\t%c",variant->chrom,variant->position,variant->previousbase); 
		if (subs > 0) { fprintf(vfile,"%c\t",variant->refb); strcpy(variant->type,"SNV,INSERTION"); } 
		else { fprintf(vfile,"\t");  strcpy(variant->type,"INSERTION"); } 
		for (i=0;i<variant->varalleles;i++)
		{
			fprintf(vfile,"%c",variant->previousbase);
			if (variant->itb[variant->alleles[i]][0] == '+') 
			{
				for (j=1;j<strlen(variant->itb[variant->alleles[i]]);j++) fprintf(vfile,"%c",variant->itb[variant->alleles[i]][j]);
				if (subs > 0) fprintf(vfile,"%c",variant->refb); 
			}
			else fprintf(vfile,"%c",variant->itb[variant->alleles[i]][0]);
			if (i < variant->varalleles-1) fprintf(vfile,",");
		}
	}
	//	for (i=0;i<variant->varalleles;i++) fprintf(vfile,"%2.1f,%2.1f,%2.1f,VP=%d,HP=%d;",variant->ctpval[i],variant->qvpvaluef[i],variant->qvpvaluer[i],variant->nvariants[i],variant->HPlength[variant->alleles[i]]); 

	int K=0;
	for (i=0;i<variant->samples;i++) 
	{
		K += variant->ploidy[i];
		depthf += variant->indcounts[i][variant->refbase]; depthr += variant->indcounts[i][variant->refbase+maxalleles];
		for (j=0;j<variant->varalleles;j++) depthf += variant->indcounts[i][variant->alleles[j]];
		for (j=0;j<variant->varalleles;j++) depthr += variant->indcounts[i][variant->alleles[j]+maxalleles];
	}
	
	int lowcoverage = 0;	if (depthf + depthr <= 2*K) lowcoverage = 1; 

	// quality value of variant 
	//if (variant->varalleles ==1) 
#if MLMETHOD ==1
	fprintf(vfile,"\t"); 
	for (i=0;i<variant->varalleles;i++)
	{
		fprintf(vfile,"%d",(int)(variant->crispvar[i].delta*10));
		if (i<variant->varalleles-1)fprintf(vfile,",");
	}
	fprintf(vfile,"\t"); 
#endif

	fprintf(vfile,"\t%d\t",(int)(-variant->ctpval[0]*10.0)); 

#if MLMETHOD ==1
	double reduction =0,reductionf=0,reductionr=0;
	for (i=0;i<variant->varalleles;i++)
	{
		reduction=variant->crispvar[i].delta - variant->crispvar[i].deltaf-variant->crispvar[i].deltar;
		reductionf=variant->crispvar[i].delta - variant->crispvar[i].deltaf;
		reductionr=variant->crispvar[i].delta - variant->crispvar[i].deltar;

		if (reduction >= 0 && variant->crispvar[i].delta >= 2) fprintf(vfile,"EMpass,");
		else if (reduction >= -0.5 && variant->crispvar[i].delta >= 2 && variant->crispvar[i].allele >=4) fprintf(vfile,"EMpass,");
		else if (reduction >= -1 && variant->crispvar[i].delta >= 3) fprintf(vfile,"EMpass1,");
		else if (reduction >= -2 && variant->crispvar[i].delta >= 5) fprintf(vfile,"EMpass2,");
		// special consideration for indels 
		else if (reduction >= -3 && variant->crispvar[i].delta >= 10 && variant->crispvar[i].allele >=4) fprintf(vfile,"EMpass3,");
		else if (reductionf >= 3 && reductionr >=3) fprintf(vfile,"EMpass-low,");

		// bug fixed aug 7 6pm 2012
		//else if (variant->crispvar[variant->alleles[i]].deltaf >= 0 && variant->crispvar[variant->alleles[i]].deltar >= 0 && variant->crispvar[variant->alleles[i]].delta >= 1.5) fprintf(vfile,"EMpass,");
		//else if (reductionf <= 2 && reductionr <= 2 && variant->crispvar[variant->alleles[i]].delta >= 4) fprintf(vfile,"EMpass,");
		//else if (reductionf <= 0.5 && reductionr <= 0.5 && variant->crispvar[variant->alleles[i]].delta >= 3) fprintf(vfile,"EMpass,");
		else fprintf(vfile,"EMfail,"); 
	}
#endif

	// low coverage lable, single strand label, strand-bias label....
	int MQflag = 0; double totalcount = variant->MQcounts[0]+variant->MQcounts[1]+variant->MQcounts[2]+variant->MQcounts[3];
	if (variant->MQcounts[0] + variant->MQcounts[1] > 0.2*totalcount) MQflag = 20;
	else if (variant->MQcounts[0] + variant->MQcounts[1] > 0.1*totalcount) MQflag = 10;

	int filter = 0;  double DOF = variant->samples-1; double deltac =0;
	if (MQflag >0) 
	{
		fprintf(vfile,"LowMQ%d",MQflag); 
		if (variant->strandpass >0)  fprintf(vfile,";StrandBias");
		if (lowcoverage ==1) fprintf(vfile,";LowDepth");
		filter = 1; 
	}
	else if (variant->strandpass >0)  
	{
		fprintf(vfile,"StrandBias");
		if (lowcoverage ==1) fprintf(vfile,";LowDepth");
		filter = 1; 
	}
	else if (lowcoverage ==1) { fprintf(vfile,"LowDepth"); filter = 1; } 
	for (i=0;i<variant->varalleles;i++)
	{
		deltac = variant->chi2stats[i][3]-DOF - (variant->chi2stats[i][1]-DOF) - (variant->chi2stats[i][2]-DOF); 
		if ((variant->chi2stats[i][3]+5 < variant->chi2stats[i][1]) || (variant->chi2stats[i][3]+5 < variant->chi2stats[i][2])) 
		//if (deltac <-5)
		{
			//if (filter ==1) fprintf(stdout,";"); 
			//if (filter > 1) fprintf(stdout,",CHI2_SB"); else fprintf(stdout,"CHI2_SB");
			//filter++;
		}
	}
	if (filter ==0) fprintf(vfile,"PASS");

	fprintf(vfile,"\tNP=%d;DP=%d,%d;VT=%s;CT=",variant->samples,depthf,depthr,variant->type);
	/*
	   for (i=0;i<variant->varalleles;i++) 
	   { 
	   if (i < variant->varalleles-1) fprintf(vfile,","); else fprintf(vfile,";AC="); 
	   } */
	for (i=0;i<variant->varalleles;i++) 
	{ 
	   	fprintf(vfile,"%.1f",variant->ctpval[i]); 
		//fprintf(vfile,"%d",variant->crispvar[i].allelecount); 
		if (i < variant->varalleles-1) fprintf(vfile,","); else fprintf(vfile,";VP=");
	}
	for (i=0;i<variant->varalleles;i++) 
	{ 
		fprintf(vfile,"%d",variant->crispvar[i].variantpools); if (i < variant->varalleles-1) fprintf(vfile,","); else fprintf(vfile,";QVpf="); 
	}
	for (i=0;i<variant->varalleles;i++) { fprintf(vfile,"%.1f",variant->qvpvaluef[i]); if (i < variant->varalleles-1) fprintf(vfile,","); else fprintf(vfile,";QVpr="); } 
	for (i=0;i<variant->varalleles;i++) { fprintf(vfile,"%.1f",variant->qvpvaluer[i]); if (i < variant->varalleles-1) fprintf(vfile,","); } 

	fprintf(vfile,";MQS=%d,%d,%d,%d;",variant->MQcounts[0],variant->MQcounts[1],variant->MQcounts[2],variant->MQcounts[3]);

	if (insallele + delallele > 0)
	{
		printhomseq=0;
		fprintf(vfile,"HP=");  
		for (i=0;i<variant->varalleles;i++) 
		{ 
			fprintf(vfile,"%d",variant->HPlength[variant->alleles[i]]); if (i < variant->varalleles-1) fprintf(vfile,",");  
			if (variant->alleles[i] >=4 && variant->itb[variant->alleles[i]][0] == '+' && variant->HPlength[variant->alleles[i]] > printhomseq) printhomseq= variant->HPlength[variant->alleles[i]];
			else if (variant->alleles[i] >=4 && variant->itb[variant->alleles[i]][0] == '-' && variant->HPlength[variant->alleles[i]] + strlen(variant->itb[variant->alleles[i]])-1 > printhomseq) printhomseq= variant->HPlength[variant->alleles[i]] + strlen(variant->itb[variant->alleles[i]])-1;
		} 
		if (printhomseq >=0 && variant->position-11  >=0 && variant->position+printhomseq+10 < reflist->lengths[current]) 	
		{
			fprintf(vfile,";FLANKSEQ=");  
			for (j=-11;j<-1;j++) fprintf(vfile,"%c",tolower(reflist->sequences[current][variant->position+j]));
			fprintf(vfile,":");
			for (j=-1;j<printhomseq;j++) fprintf(vfile,"%c",toupper(reflist->sequences[current][variant->position+j])); 
			fprintf(vfile,":");
			for (j=printhomseq;j<printhomseq+10;j++) fprintf(vfile,"%c",tolower(reflist->sequences[current][variant->position+j]));
		} 
	}
	else //print flanking sequence even for SNPs
	{
		fprintf(vfile,"FLANKSEQ=");  
		j=-10; if (variant->position +j <0) j = 0;
		for (j=-10;j<0;j++) fprintf(vfile,"%c",tolower(reflist->sequences[current][variant->position+j]));
		fprintf(vfile,":");
		for (j=0;j<1;j++) fprintf(vfile,"%c",toupper(reflist->sequences[current][variant->position+j])); 
		fprintf(vfile,":");
		for (j=1;j<11 && j+variant->position<reflist->lengths[current];j++) fprintf(vfile,"%c",tolower(reflist->sequences[current][variant->position+j]));
	}

#if MLMETHOD ==1
	fprintf(vfile,";EMstats=");
	for (i=0;i<variant->varalleles;i++)
	{
		fprintf(vfile,"%0.6f:%0.2f:%0.2f:%0.2f",variant->crispvar[i].AF_ML,variant->crispvar[i].delta,variant->crispvar[i].deltaf,variant->crispvar[i].deltar);	
		if (i<variant->varalleles-1)fprintf(vfile,",");
	}
	fprintf(vfile,";BFGSstats=");
	for (i=0;i<variant->varalleles;i++)
	{
		fprintf(vfile,"%0.2f:%0.5f:%0.5f:%0.5f:%0.5f",variant->crispvar[i].deltaBFGS,variant->crispvar[i].E[0],variant->crispvar[i].E[1],variant->crispvar[i].E[2],variant->crispvar[i].E[3]);	
		if (i<variant->varalleles-1)fprintf(vfile,",");
	}
#endif
	// transition/transversion label

	// change GT to allele frequency may 7 2012
	//for (i=1;i<variant->varalleles;i++) fprintf(vfile,":A%d",i+1);
	//fprintf(vfile,":R0"); for (i=0;i<variant->varalleles;i++) fprintf(vfile,":R%d",i+1);

	int b0=0;
	fprintf(vfile,"\tGT:AF:ADf:ADr"); 
	//fprintf(vfile,"\tAC:QV:mAC:ADf:ADr:VP\t"); 
	for (i=0;i<variant->samples;i++)
	{
		coverage = variant->indcounts[i][variant->refbase] + variant->indcounts[i][variant->refbase+maxalleles];
		for (j=0;j<variant->varalleles;j++) coverage += variant->indcounts[i][variant->alleles[j]] + variant->indcounts[i][variant->alleles[j]+maxalleles];
		j=0;
		fprintf(vfile,"\t0/0:%0.3f",(float)(variant->indcounts[i][variant->alleles[j]] + variant->indcounts[i][variant->alleles[j]+maxalleles])/(coverage+0.001));
		for (j=1;j<variant->varalleles;j++) fprintf(vfile,",%0.3f",(float)(variant->indcounts[i][variant->alleles[j]] + variant->indcounts[i][variant->alleles[j]+maxalleles])/(coverage+0.001));

		// pool presence probability 
		//fprintf(vfile,"\t0/0:%0.3f",variant->poolpv[i][0]); for (j=1;j<variant->varalleles;j++) fprintf(vfile,",%0.3f",variant->poolpv[i][j]);
		// distribute the bidirectional reads between the two strands dec 5 2012
		b0 = variant->indcounts[i][variant->refbase+2*maxalleles]; 
		fprintf(vfile,":%d",variant->indcounts[i][variant->refbase]+b0/2); 
		for (j=0;j<variant->varalleles;j++) 
		{
			b0 = variant->indcounts[i][variant->alleles[j]+2*maxalleles]; 
			fprintf(vfile,",%d",variant->indcounts[i][variant->alleles[j]]+b0/2);	
		}

		b0 = variant->indcounts[i][variant->refbase+2*maxalleles]-variant->indcounts[i][variant->refbase+2*maxalleles]/2;
		fprintf(vfile,":%d",variant->indcounts[i][variant->refbase+maxalleles]+b0);
		for (j=0;j<variant->varalleles;j++) 
		{
			b0 = variant->indcounts[i][variant->alleles[j]+2*maxalleles]-variant->indcounts[i][variant->alleles[j]+2*maxalleles]/2;
			fprintf(vfile,",%d",variant->indcounts[i][variant->alleles[j]+maxalleles]+b0);
		}
	}
	fprintf(vfile,"\n");
	return 1;
}

void print_crispheader(struct OPTIONS* options)
{
	FILE* vfile = options->vfile; int i=0;
	// also print program options, poolsize, qvoffset, etc 
	time_t now; time(&now); 
	fprintf(options->vfile,"##fileformat=VCFv4.0\n##fileDate=%s",ctime(&now));
	fprintf(vfile,"##source=CRISP_V0.1\n");
	fprintf(vfile,"##reference=%s\n",options->fastafile);
	fprintf(vfile,"##options: poolsize=%d,bamfiles=%d,qvoffset=%d,bedfile=%s,min_base_quality=%d\n",options->POOLSIZE,options->bamfiles,QVoffset,options->bedfile,MINQ);
	fprintf(vfile,"##INFO=<ID=NP,Number=1,Type=Integer,Description=\"Number of Pools With Data\">\n");
	fprintf(vfile,"##INFO=<ID=DP,Number=2,Type=Integer,Description=\"Total number of reads (+strand,-strand) across all pools (filtered reads only)\">\n");
	fprintf(vfile,"##INFO=<ID=CT,Number=.,Type=Float,Description=\"contingency table p-value for each variant allele in same order as listed in column 5\">\n");
	//fprintf(vfile,"##INFO=<ID=QVpf,Number=.,Type=Float,Description=\"quality values based p-value for each variant allele using forward strand reads\">\n");
	//fprintf(vfile,"##INFO=<ID=QVpr,Number=.,Type=Float,Description=\"quality values based p-value for each variant allele using reverse strand reads\">\n");
	fprintf(vfile,"##INFO=<ID=VP,Number=.,Type=Integer,Description=\"Number of Pools with variant allele(s)\">\n");
	fprintf(vfile,"##INFO=<ID=HP,Number=.,Type=Integer,Description=\"Ambiguity in positioning of indel (homopolymer run length or microsatellite length)\">\n");
	fprintf(vfile,"##INFO=<ID=MQ,Number=.,Type=Integer,Description=\"# of reads with mapping qualities 0-9,10-19,20-39,40-255\">\n");
	fprintf(vfile,"##FILTER=<ID=StrandBias,Description=\"strand bias in distribution of reads between reference and variant alleles on two strands\">\n");
	fprintf(vfile,"##FILTER=<ID=LowDepth,Description=\"low average coverage: less than 1 (filtered) read per haplotype across all samples \">\n");
	fprintf(vfile,"##FILTER=<ID=LowMQ20,Description=\" >20 percent of reads have mapping quality score less than 20\">\n");
	fprintf(vfile,"##FILTER=<ID=LowMQ10,Description=\" >10 percent of reads have mapping quality score less than 20\">\n");
	fprintf(vfile,"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
	//fprintf(vfile,"##FORMAT=<ID=PP,Number=.,Type=Float,Description=\"Probability that the pool contains the variant allele(s), one value for each variant allele listed in column 5\">\n");
	fprintf(vfile,"##FORMAT=<ID=AF,Number=.,Type=Float,Description=\"variant allele frequency in pool, one value for each variant allele listed in column 5\">\n");
	fprintf(vfile,"##FORMAT=<ID=ADf,Number=.,Type=Integer,Description=\"Number of reads aligned to the forward strand of the genome supporting reference allele and the alternate alleles in the order listed\">\n");
	fprintf(vfile,"##FORMAT=<ID=ADr,Number=.,Type=Integer,Description=\"Number of reads aligned to the reverse strand of the genome supporting reference allele and the alternate alleles in the order listed\">\n");
	fprintf(options->vfile,"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t");
	for (i=0;i<options->bamfiles-1;i++) fprintf(options->vfile,"%s\t",options->bamfilelist[i]);
	fprintf(options->vfile,"%s\n",options->bamfilelist[options->bamfiles-1]);

}

void print_crispoptions()
{
	fprintf(stderr,"\nCRISP: statistical method to identify SNVs and indels from pooled DNA sequencing data (requires multiple samples (pools), ideally >=5 samples)\n\n");
	fprintf(stderr,"./CRISP [options] --bams file_bam_paths --ref reference.fasta --VCF variantcalls.VCF -p poolsize > variantcalls.log\n\n");
	fprintf(stderr,"Options:\n");
	fprintf(stderr,"         --bams         textfile with list of bam file paths (one bam file for each pool)\n");
	fprintf(stderr,"         --bam          bam file for one pool, specify file for each pool using --bam pool1.bam --bam pool2.bam .... --bam pooln.bam\n");
	fprintf(stderr,"         --ref       	Indexed Reference Sequence file (fasta)\n");
	fprintf(stderr,"         -p/--poolsize <int>    poolsize (number of haploid genomes in each pool), for diploid genomes: 2 x # individuals\n");
	fprintf(stderr,"         --VCF       	VCF file to which the variant calls will be output \n");
	fprintf(stderr,"         --qvoffset  	quality value offset, 33 for Sanger format, 64 for Illumina 1.3+ format\n");
	// change this to --phred33 or phred 64 
	fprintf(stderr,"         --mbq       	minimum base quality to consider a base for variant calling, default 13\n");
	fprintf(stderr,"         --mmq       	minimum read mapping quality to consider a read for variant calling, default 20\n");
	fprintf(stderr,"         --regions      region(s) in which variants will be called, e.g chr1:654432-763332. BAM files should be indexed for using this option.\n");
	fprintf(stderr,"         --minc      	minimum number of reads with alternate allele required for calling a variant, default 4\n");
	fprintf(stderr,"         --ctpval    	threshold on the contingency table p-value for calling position as variant (specified as log10), default is -3.5\n");
	fprintf(stderr,"         --qvpval    	threshold on the quality values based p-value for calling position as variant (specified as log10), default is -5\n");
	fprintf(stderr,"         --perms  <int> maximum number of permutations for calculating contingency table p-value, default 20000\n");
	fprintf(stderr,"         --filterreads 0/1 filter reads with excessive number of mismatches (and gaps) compared to the reference sequence, Default is 1. Set to 0 to disable filtering\n");
	fprintf(stderr,"         --verbose  <int> amount of information to output to log file, 0: no output, 1: medium (default), 2: detailed\n");
	fprintf(stderr,"         --OPE  <int> identify overlapping paired-end reads and treat as single read, default is 1, set to 0 to disable this (can slow the program for high-coverage datasets)\n");
	//fprintf(stderr,"         --seqerate  	estimate of average sequencing error rate. Default value is -1 which uses individual base qualities to estimate p-values\n");
	//fprintf(stderr,"         --maxm      	maximum number of mismatches allowed for read to be considered for snp calling. Default value is 3\n");

	fprintf(stderr,"\nNotes:\n\n");
	fprintf(stderr," 1. CRISP requires poolsize and reference fasta file for making variant calls\n");
	fprintf(stderr," 2. CRISP requires at least two pools to make variant calls, but 5 or more pools are ideal\n");
	fprintf(stderr," 2.1. The aligned reads for each pool should be a single bam file that is sorted by chromosomal coordinates\n");
	fprintf(stderr," 3. The reference sequence file should be indexed using 'samtools faidx' or a similar program and placed in same directory as fasta file with extension .fai\n");
	fprintf(stderr," 4. The ploidy of each pool is assumed to be the same, this will be changed in future releases\n");
	fprintf(stderr," 5. For human resequencing, if a bedfile is not specified, the program will evaluate each base for variant calling\n");
	fprintf(stderr," 6. Please set the quality value offset correctly: 33 for quality values encoded in Sanger format and 64 for Illumina 1.3+ format. Setting this value incorrectly can result in significant overcalling/undercalling of variants. Default is 33\n");
	fprintf(stderr," 7. Please make sure that the reference sequence file is the same as the one used to align the reads in the BAM files and that the BAM files are coordinate sorted, the program will not check this\n");
	fprintf(stderr," 8. For indel analysis, CRISP assumes that indels are left justified\n");
	fprintf(stderr," 9. the program uses the samtools API for reading bam files \n");
	fprintf(stderr," 10. BAM files should be indexed for using the --regions option with the indexed bam file present in the same directory with name pooln.bam.bai\n\n");
	//fprintf(stderr," 11. The bed file can be used to re-call variants at specified positions in the genome\n\n");
	//fprintf(stderr," 7. The value of the parameter --maxm should be choosen based on average length of reads (3 for < 75-bp reads, 4 for < 100 bp reads...)\n\n");

}


