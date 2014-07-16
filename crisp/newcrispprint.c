
/* functions to print variant and genotypes or allele counts to VCF file */

/* TODO 
1. genotype quality for tri-allelic variants is not correct for now
2. low coverage pools -> genotype = ./. or .
3. for single-end reads, don't print bidirectional read counts...

*/
#include<time.h>
#include "crisp_header_options.c" // functions to print CRISP program options, VCF header

// as currently implemented, coverage/depth is only calculated using variant alleles (not seq errors or third/fourth alleles)... check this
// allele counts are printed for reads aligned to forward, reverse and bidirectional (covered from both ends) (3 counts for each allele)

// also calculate allele frequency of variant alleles across all pools combined (assumption is that pool size is constant)

int print_VCF_first_5cols(struct VARIANT* variant,FILE* vfile);
void print_VCF_allelefreqs(struct VARIANT* variant,FILE* vfile);
void print_VCF_genotypes_diploid(struct VARIANT* variant, FILE* vfile);
void print_VCF_genotypes_pooled(struct VARIANT* variant,FILE* vfile);
void print_flanking_sequence(struct VARIANT* variant, FILE* vfile,REFLIST* reflist,int current, int is_indel_variant);
int print_pooledvariant(struct VARIANT* variant,FILE* vfile,REFLIST* reflist,int current,ALLELE* maxvec,int EMflag);


int print_VCF_first_5cols(struct VARIANT* variant,FILE* vfile)
{
	int i=0,j=0,k=0; int delallele = 0, insallele = 0, subs =0, longestdel=0;
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
	return insallele + delallele;
}

// CODE for printing CRISP allele frequencies, no genotypes, old version of CRISP
void print_VCF_allelefreqs(struct VARIANT* variant,FILE* vfile)
{
	int i=0,j=0,b0=0,coverage=0;
	fprintf(vfile,"\tGT:AF:DP:ADf:ADr"); 
	for (i=0;i<variant->samples;i++)
        {
                coverage = variant->indcounts[i][variant->refbase] + variant->indcounts[i][variant->refbase+maxalleles] + variant->indcounts[i][variant->refbase+2*maxalleles];
                for (j=0;j<variant->varalleles;j++) coverage += variant->indcounts[i][variant->alleles[j]] + variant->indcounts[i][variant->alleles[j]+maxalleles] + variant->indcounts[i][variant->alleles[j]+2*maxalleles];

                j=0; fprintf(vfile,"\t0/0:%0.3f",(float)(variant->indcounts[i][variant->alleles[j]] + variant->indcounts[i][variant->alleles[j]+maxalleles] + variant->indcounts[i][variant->alleles[j]+2*maxalleles])/(coverage+0.001));
                for (j=1;j<variant->varalleles;j++) fprintf(vfile,",%0.3f",(float)(variant->indcounts[i][variant->alleles[j]] + variant->indcounts[i][variant->alleles[j]+maxalleles]+variant->indcounts[i][variant->alleles[j]+2*maxalleles])/(coverage+0.001));
		fprintf(vfile,":%d",coverage);
		fprintf(vfile,":%d",variant->readdepths[i]);

		// add 1/2 of bidirection reads to the two strands and print coverage for each allele on each strand
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
}
// for diploid genotypes, also output PL: genotype likelihood for all three or 6 possible genotypes 
void print_VCF_genotypes_diploid(struct VARIANT* variant, FILE* vfile)
{
	int i=0,j=0,k=0,coverage=0,bidir=0;

	if (variant->counts[variant->refbase+2*maxalleles] + variant->counts[variant->alleles[0]+2*maxalleles] >= 2) bidir = 1;
	if (bidir == 1) fprintf(vfile,"\tGT:GQ:DP:ADf:ADr:ADb\t");  else fprintf(vfile,"\tGT:GQ:DP:ADf:ADr\t");
	

	for (i=0;i<variant->samples;i++)
	{
		if ((double)variant->readdepths[i]/variant->ploidy[i] <= 1) fprintf(vfile,"./.:."); 
		else if (variant->varalleles ==1) // print genotypes for bi-allelic variants
		{
			if (variant->crispvar[0].AC[i] ==0) fprintf(vfile,"0/0:%d",(int)(0.5+variant->crispvar[0].QV[i]));
			else if (variant->crispvar[0].AC[i] ==1) fprintf(vfile,"0/1:%d",(int)(0.5+variant->crispvar[0].QV[i]));
			else if (variant->crispvar[0].AC[i] ==2) fprintf(vfile,"1/1:%d",(int)(0.5+variant->crispvar[0].QV[i]));
		}
		else if (variant->varalleles ==2) // print genotypes for multi-allelic variants
		{
			if (variant->crispvar[0].AC[i] ==0 && variant->crispvar[1].AC[i] ==0) fprintf(vfile,"0/0:%d",(int)(0.5+variant->crispvar[0].QV[i]));
			else if (variant->crispvar[0].AC[i] ==1 && variant->crispvar[1].AC[i] ==0) fprintf(vfile,"0/1:%d",(int)(0.5+variant->crispvar[0].QV[i]));
			else if (variant->crispvar[0].AC[i] ==2 && variant->crispvar[1].AC[i] ==0) fprintf(vfile,"1/1:%d",(int)(0.5+variant->crispvar[0].QV[i]));
			else if (variant->crispvar[0].AC[i] ==1 && variant->crispvar[1].AC[i] ==1) fprintf(vfile,"1/2:%d",(int)(0.5+variant->crispvar[0].QV[i]));
			else if (variant->crispvar[0].AC[i] ==0 && variant->crispvar[1].AC[i] ==1) fprintf(vfile,"0/2:%d",(int)(0.5+variant->crispvar[1].QV[i]));
			else if (variant->crispvar[0].AC[i] ==0 && variant->crispvar[1].AC[i] ==2) fprintf(vfile,"2/2:%d",(int)(0.5+variant->crispvar[1].QV[i]));
			else fprintf(vfile,"./.:.");
		}
		else fprintf(vfile,"./.:.");
		//fprintf(vfile,":%d",(int)(variant->crispvar[0].QV[i]+0.5)); 
		//for (j=1;j<variant->varalleles;j++) fprintf(vfile,",%.1f",variant->crispvar[j].QV[i]); 

		//fprintf(vfile,":%d",coverage);
		fprintf(vfile,":%d",variant->readdepths[i]);
		fprintf(vfile,":%d",variant->indcounts[i][variant->refbase]); 
		for (j=0;j<variant->varalleles;j++) fprintf(vfile,",%d",variant->indcounts[i][variant->alleles[j]]);
		fprintf(vfile,":%d",variant->indcounts[i][variant->refbase+maxalleles]);
		for (j=0;j<variant->varalleles;j++) fprintf(vfile,",%d",variant->indcounts[i][variant->alleles[j]+maxalleles]);
		if (bidir ==1) 
		{
			fprintf(vfile,":%d",variant->indcounts[i][variant->refbase+2*maxalleles]);
			for (j=0;j<variant->varalleles;j++) fprintf(vfile,",%d",variant->indcounts[i][variant->alleles[j]+2*maxalleles]);
		}
		fprintf(vfile,"\t");
	}
	fprintf(vfile,"\n");
	
}


void print_VCF_genotypes_pooled_full(struct VARIANT* variant,FILE* vfile)
{
	int i=0,j=0,k=0,coverage=0,bidir=0;
	if (variant->counts[variant->refbase+2*maxalleles] + variant->counts[variant->alleles[0]+2*maxalleles] >= 2) bidir = 1;
	if (bidir == 1) fprintf(vfile,"\tGT:GQ:AC:DP:ADf:ADr:ADb\t");  else fprintf(vfile,"\tGT:GQ:AC:DP:ADf:ADr\t");

	int AC1 = 0, AC2 = 0; int printflag = 0;
	for (i=0;i<variant->samples;i++)
	{
		AC1 += variant->crispvar[0].AC[i]; 
		if (variant->varalleles >= 2) AC2 += variant->crispvar[1].AC[i]; 
		if (variant->crispvar[0].AC[i] > 0 && variant->varalleles >= 2 && variant->crispvar[1].AC[i] > 0) printflag++;
	}

	for (i=0;i<variant->samples;i++)
	{
		if ((double)variant->readdepths[i]/variant->ploidy[i] <= 1) 
		{
			for (j=0;j<variant->ploidy[i];j++) fprintf(vfile,"./"); fprintf(vfile,".:0:.");
		}
		else if (variant->ploidy[0] > 2 && AC2 ==0)  // print reference allele count as well.
		{
			if (variant->ploidy[i]-variant->crispvar[0].AC[i] > 0) fprintf(vfile,"0");
			for (j=0;j<variant->ploidy[i]-variant->crispvar[0].AC[i];j++) fprintf(vfile,"/0"); 
			fprintf(vfile,":%d:0",(int)(variant->crispvar[0].QV[i]+0.5)); 
		}
		else if (variant->ploidy[0] > 2)  // pool size is same
		{
			fprintf(vfile,"%d",variant->crispvar[0].AC[i]); 
			for (j=1;j<variant->varalleles;j++) fprintf(vfile,",%d",variant->crispvar[j].AC[i]); 
			fprintf(vfile,":%d",(int)(variant->crispvar[0].QV[i]+0.5)); 
		}
		else 
		{
			for (j=0;j<variant->ploidy[i];j++) fprintf(vfile,"./"); fprintf(vfile,".:0:.");
		}

		//for (j=1;j<variant->varalleles;j++) fprintf(vfile,",%.1f",variant->crispvar[j].QV[i]); 
		//fprintf(vfile,":"); fprintf(vfile,"%0.1f",variant->crispvar[0].meanAC[i]); 	
		//for (j=1;j<variant->varalleles;j++) fprintf(vfile,",%.1f",variant->crispvar[j].meanAC[i]); 
		coverage = variant->indcounts[i][variant->refbase] + variant->indcounts[i][variant->refbase+maxalleles] + variant->indcounts[i][variant->refbase+2*maxalleles];
		for (j=0;j<variant->varalleles;j++) coverage += variant->indcounts[i][variant->alleles[j]] + variant->indcounts[i][variant->alleles[j]+maxalleles] + variant->indcounts[i][variant->alleles[j]+2*maxalleles];	

		//fprintf(vfile,":%d",coverage);
		fprintf(vfile,":%d",variant->readdepths[i]);
		fprintf(vfile,":%d",variant->indcounts[i][variant->refbase]); 
		for (j=0;j<variant->varalleles;j++) fprintf(vfile,",%d",variant->indcounts[i][variant->alleles[j]]);
		fprintf(vfile,":%d",variant->indcounts[i][variant->refbase+maxalleles]);
		for (j=0;j<variant->varalleles;j++) fprintf(vfile,",%d",variant->indcounts[i][variant->alleles[j]+maxalleles]);
		if (bidir ==1) 
		{
			fprintf(vfile,":%d",variant->indcounts[i][variant->refbase+2*maxalleles]);
			for (j=0;j<variant->varalleles;j++) fprintf(vfile,",%d",variant->indcounts[i][variant->alleles[j]+2*maxalleles]);
		}
		//if (variant->crispvar[0].AC[i] > 0) fprintf(vfile,":="); else fprintf(vfile,":."); 
		fprintf(vfile,"\t");
	}
	fprintf(vfile,"\n");
}

/* CODE for printing columns 9,10..... | GENOTYPE COLUMSN */
void print_VCF_genotypes_pooled(struct VARIANT* variant,FILE* vfile)
{
	int i=0,j=0,k=0,coverage=0,bidir=0;
	if (variant->counts[variant->refbase+2*maxalleles] + variant->counts[variant->alleles[0]+2*maxalleles] >= 2) bidir = 1;
	if (bidir == 1) fprintf(vfile,"\tMLAC:GQ:DP:ADf:ADr:ADb\t");  else fprintf(vfile,"\tMLAC:GQ:DP:ADf:ADr\t");
	for (i=0;i<variant->samples;i++)
	{
		if ((double)variant->readdepths[i]/variant->ploidy[i] <= 1) fprintf(vfile,"."); 
		else if (variant->ploidy[0] > 2 && variant->varpoolsize > 0)  // print reference allele count as well.
		{
			k = variant->ploidy[i]; for (j=0;j<variant->varalleles;j++) k -= variant->crispvar[j].AC[i]; 
			if (k < 0) // total allele counts do not add up, can happen for tri/quatra-allelic variants 
			{
				fprintf(vfile,".");
			}
			else
			{
				fprintf(vfile,"%d,",k); fprintf(vfile,"%d",variant->crispvar[0].AC[i]); 
				for (j=1;j<variant->varalleles;j++) fprintf(vfile,",%d",variant->crispvar[j].AC[i]); 
			}
		}
		else if (variant->ploidy[0] > 2) 
		{
			fprintf(vfile,"%d",variant->crispvar[0].AC[i]); 
			for (j=1;j<variant->varalleles;j++) fprintf(vfile,",%d",variant->crispvar[j].AC[i]); 
		}
		else fprintf(vfile,".");

		fprintf(vfile,":%d",(int)(variant->crispvar[0].QV[i]+0.5)); 
		//for (j=1;j<variant->varalleles;j++) fprintf(vfile,",%.1f",variant->crispvar[j].QV[i]); 
		//fprintf(vfile,":"); fprintf(vfile,"%0.1f",variant->crispvar[0].meanAC[i]); 	
		//for (j=1;j<variant->varalleles;j++) fprintf(vfile,",%.1f",variant->crispvar[j].meanAC[i]); 
		coverage = variant->indcounts[i][variant->refbase] + variant->indcounts[i][variant->refbase+maxalleles] + variant->indcounts[i][variant->refbase+2*maxalleles];
		for (j=0;j<variant->varalleles;j++) coverage += variant->indcounts[i][variant->alleles[j]] + variant->indcounts[i][variant->alleles[j]+maxalleles] + variant->indcounts[i][variant->alleles[j]+2*maxalleles];	

		//fprintf(vfile,":%d",coverage);
		fprintf(vfile,":%d",variant->readdepths[i]);
		fprintf(vfile,":%d",variant->indcounts[i][variant->refbase]); 
		for (j=0;j<variant->varalleles;j++) fprintf(vfile,",%d",variant->indcounts[i][variant->alleles[j]]);
		fprintf(vfile,":%d",variant->indcounts[i][variant->refbase+maxalleles]);
		for (j=0;j<variant->varalleles;j++) fprintf(vfile,",%d",variant->indcounts[i][variant->alleles[j]+maxalleles]);
		if (bidir ==1) 
		{
			fprintf(vfile,":%d",variant->indcounts[i][variant->refbase+2*maxalleles]);
			for (j=0;j<variant->varalleles;j++) fprintf(vfile,",%d",variant->indcounts[i][variant->alleles[j]+2*maxalleles]);
		}
		//if (variant->crispvar[0].AC[i] > 0) fprintf(vfile,":="); else fprintf(vfile,":."); 
		fprintf(vfile,"\t");
	}
	fprintf(vfile,"\n");
}

// print flanking bases for each SNP/indel, use to evaluate sequence context of variant, homopolymer length etc
void print_flanking_sequence(struct VARIANT* variant, FILE* vfile,REFLIST* reflist,int current, int is_indel_variant)
{
	int i=0,j=0,printhomseq=0;
	if (is_indel_variant > 0)
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
}

// what if some pool has very little data | reads compared to others -> mark as missing genotype, NP=190 -> reduce...

// output variant calls and allele counts to VCF file (this is done once for each variant, even if it is multi-allelic) 
int print_pooledvariant(struct VARIANT* variant,FILE* vfile,REFLIST* reflist,int current,ALLELE* maxvec,int EMflag)
{
	variant->previousbase = 'N';
        if (current >=0 && variant->position >=1)  variant->previousbase = toupper(reflist->sequences[current][variant->position-1]);
	if (variant->previousbase != 'A' && variant->previousbase != 'C'  && variant->previousbase != 'G' && variant->previousbase != 'T') variant->previousbase = 'N';
        // one base prior to current base, used for representing indels in VCF file

	int K = 0; int i=0,j=0; int depth[3] = {0,0,0}; int lowcoverage=0; int MQflag =0;  
	int is_indel_variant = print_VCF_first_5cols(variant,vfile); 

	double varscores[variant->varalleles]; double bestscore=0; // variant quality to print in column 6
	for (i=0;i<variant->varalleles;i++) 
	{
		if (EMflag) varscores[i] = variant->crispvar[i].delta*10; else varscores[i] = variant->ctpval[0]*(-10);
		if (varscores[i] > bestscore) bestscore = varscores[i];
	}
	fprintf(vfile,"\t%d\t",(int)bestscore); 

	for (i=0;i<variant->samples;i++) 
	{
		K += variant->ploidy[i];
		depth[0] += variant->indcounts[i][variant->refbase]; depth[1] += variant->indcounts[i][variant->refbase+maxalleles];
		depth[2] += variant->indcounts[i][variant->refbase+2*maxalleles];
		for (j=0;j<variant->varalleles;j++) 
		{
			depth[0] += variant->indcounts[i][variant->alleles[j]]; 
			depth[1] += variant->indcounts[i][variant->alleles[j]+maxalleles];
			depth[2] += variant->indcounts[i][variant->alleles[j]+2*maxalleles];
		}
	}
	if (depth[0] + depth[1] + depth[2] <= 2*K) lowcoverage = 1;  
	// lowdepth is defined as less than 2x reads per haplotype on average // below this depth the contingency table p-value filter is not powerful
	
	double totalcount = variant->MQcounts[0]+variant->MQcounts[1]+variant->MQcounts[2]+variant->MQcounts[3];
	if (variant->MQcounts[0] + variant->MQcounts[1] > 0.2*totalcount) MQflag = 20;
	else if (variant->MQcounts[0] + variant->MQcounts[1] > 0.1*totalcount) MQflag = 10;
	
	// FILTER : PASS if this position has passed all filters, i.e. a call is made at this position. Otherwise, if the site has not passed all filters, a semicolon-separated list of codes for filters that fail. e.g. “q10;s50” might indicate that at this site the quality is below 10 and the number of samples with data is below 50% of the total number of samples.
	// low coverage lable, single strand label, strand-bias label....
	char FILTER[1024]; int filters = 0; //strcpy(FILTER,"\0");
	if (MQflag >0) filters++; if (variant->strandpass > 0 && EMflag ==0) filters++; else if (lowcoverage ==1) filters++;

	if (MQflag >0) 
	{
		fprintf(vfile,"LowMQ%d",MQflag); 
		if (variant->strandpass >0 && EMflag ==0)  fprintf(vfile,";StrandBias");
		if (lowcoverage ==1) fprintf(vfile,";LowDepth");
	}
	else if (variant->strandpass >0 && EMflag ==0)  
	{
		fprintf(vfile,"StrandBias");
		if (lowcoverage ==1) fprintf(vfile,";LowDepth");
	}
	else if (lowcoverage ==1) fprintf(vfile,"LowDepth");
		
	char* filter_strings[] =  { "EMpass", "EMpass1", "EMpass2", "EMpass_indel", "EMfail" }; 
	int varfilters[variant->varalleles]; int PASS_variant=0;

	// use chisquare for deltaf thresholds, p=0.01 = 6.64/2 | p= 0.001, 10.83/2
	for (i=0;i<variant->varalleles;i++)
	{
		if (variant->crispvar[i].deltaf < 3.32 && variant->crispvar[i].delta >= 2) varfilters[i] = 0;
		else if (variant->crispvar[i].deltaf < 5.42 && variant->crispvar[i].delta >= 5) varfilters[i] = 1;
		else if (variant->crispvar[i].deltaf < 5.42 && variant->crispvar[i].delta >= 10) varfilters[i] = 2;
		else if (variant->crispvar[i].deltaf < 5.42 && variant->crispvar[i].delta >= 5 && variant->crispvar[i].allele >=4) varfilters[i] = 3;
		else varfilters[i] = 4;
		if (varfilters[i] != 4) PASS_variant++; 
	}
	if (filters ==0 && PASS_variant > 0) fprintf(vfile,"PASS"); else if (filters ==0) fprintf(vfile,".");

	fprintf(vfile,"\tNP=%d;DP=%d,%d,%d;VT=%s;CT=",variant->samples,depth[0],depth[1],depth[2],variant->type);
	for (i=0;i<variant->varalleles;i++) 
	{ 
		fprintf(vfile,"%.1f",variant->ctpval[i]); 
		if (i < variant->varalleles-1) fprintf(vfile,","); else fprintf(vfile,";VP="); 
	} 
	for (i=0;i<variant->varalleles;i++) 
	{
		if (EMflag ==1) fprintf(vfile,"%d",variant->varpools[i]);  // varpools is number of pools with variant using read count filters... 
		else fprintf(vfile,"%d",variant->crispvar[i].variantpools); 
		if (i < variant->varalleles-1) fprintf(vfile,","); 
	}
	if (EMflag ==1) // allele count is estimated only by EM algorithm
	{
		fprintf(vfile,";VF=");
		for (i=0;i<variant->varalleles;i++)
		{
			fprintf(vfile,filter_strings[varfilters[i]]); if (i < variant->varalleles-1) fprintf(vfile,","); 
		}
		fprintf(vfile,";AC=");
		for (i=0;i<variant->varalleles;i++) 
		{ 
			fprintf(vfile,"%d",variant->crispvar[i].allelecount); if (i < variant->varalleles-1) fprintf(vfile,","); 
		}
		//AF : allele frequency for each ALT allele in the same order as listed: use this when estimated from primary data, not called genotypes
		fprintf(vfile,";AF=");
		for (i=0;i<variant->varalleles;i++) 
		{ 
			fprintf(vfile,"%0.5f",variant->crispvar[i].AF_ML); if (i < variant->varalleles-1) fprintf(vfile,","); 
		}
		
		fprintf(vfile,";EMstats=");
		for (i=0;i<variant->varalleles;i++)
		{
			// print individual quality scores in INFO field separately.
			fprintf(vfile,"%0.2f:%0.2f",variant->crispvar[i].delta,variant->crispvar[i].deltaf);
			if (i<variant->varalleles-1)fprintf(vfile,",");
		}
		fprintf(vfile,";HWEstats=");
		for (i=0;i<variant->varalleles;i++)
		{
			fprintf(vfile,"%0.1f",-10*variant->crispvar[i].deltar);	
			if (i<variant->varalleles-1)fprintf(vfile,",");
		}
	}

	// change to MQS (Mapping quality statistics)	
	fprintf(vfile,";MQS=%d,%d,%d,%d;",variant->MQcounts[0],variant->MQcounts[1],variant->MQcounts[2],variant->MQcounts[3]);

	print_flanking_sequence(variant,vfile,reflist,current,is_indel_variant);

	if (variant->IUPAC ==1) fprintf(vfile,";IUPAC=%c",reflist->sequences[current][variant->position]); // 07/10/13 added as flag if reference base is ambig.

	// print genotypes/allele frequencies
	if (EMflag ==0) print_VCF_allelefreqs(variant,vfile); 
	else if (variant->ploidy[0] ==2 || variant->ploidy[0] ==1) print_VCF_genotypes_diploid(variant,vfile);
	else print_VCF_genotypes_pooled(variant,vfile);

	return 1;
}

/*
fprintf(vfile,";BFGSstats=");
for (i=0;i<variant->varalleles;i++)
{
	fprintf(vfile,"%0.2f:%0.5f:%0.5f:%0.5f:%0.5f",variant->crispvar[i].deltaBFGS,variant->crispvar[i].E[0],variant->crispvar[i].E[1],variant->crispvar[i].E[2],variant->crispvar[i].E[3]);	
	if (i<variant->varalleles-1)fprintf(vfile,",");
}*/
// transition/transversion label

// change GT to allele frequency may 7 2012
//for (i=1;i<variant->varalleles;i++) fprintf(vfile,":A%d",i+1);
//fprintf(vfile,":R0"); for (i=0;i<variant->varalleles;i++) fprintf(vfile,":R%d",i+1);

	//fprintf(vfile,"\tGT:AF:ADf:ADr");  fprintf(vfile,"\tAC:QV:mAC:ADf:ADr:ADb\t"); 
