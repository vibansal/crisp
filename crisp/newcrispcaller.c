/* functions that have CODE for new CRISP caller using genotype likelihood approach (EM based) */

#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include "crispcaller.h"
#include "../FET/contables.h"
#include "../FET/chisquare.h"
#include "allelecounts.h"

//#define MLMETHOD 1 // change this to 0 to disable EM method and optimizations.

//#include "../BFGSmethod/bfgs.h" // for using BMC paper code 
//#include "../Lbfgsb.3.0/lbfgsb.h"  //uncomment this if want to use own BFGS method function adaptation
#include "../FET/chisquare-stratified.c" // functions for calculating chi-square based p-value (stratified counts)
#include "crispEM.c" // functions for finding maximum likelihood pooled genotypes

int USE_BASE_QVS =1;
int MAXITER = 20000; int MINPOS = 0, MAXPOS =0;
double thresh1 = -3.5, MINQpv = -5, thresh3 = -2;  //three pvalue thresholds 
double SER = -1; double RELAX = 0.75; int MIN_READS = 4; double  Qmin = 23; // error rate of 0.005
int FAST_FILTER =1;

#include "newcrispprint.c" // print CRISP program options and CRISP VCF record
#include "pooledFET.c" // calculate FET p-value using variant->indcounts (allele counts)
#include "poolstats.c"
#include "oldcrispmethod.c"

int OUTPUTGENOTYPES = 0; 
int EMflag =1;

int evaluate_variant(READQUEUE* bq, struct VARIANT* variant,int allele1,int allele2,int allele3,double paf[])
{
	//betanorm = log10(alpha_1+beta_1-1) + ncr(alpha_1+beta_1-2,alpha_1-1); betanorm *= -1;
	// x ^ (alpha-1). (1-x)^ (beta-1)
	int i=0,k=0;  int permflag =0;
	double p[2]={0.0,0.0};
	double p1=0,p2=0,ep; double C1=0,C2=0,C3=0,N=0;
	double pv1=0,pv2=0,chi2stats[6] = {0,0,0,0,0,0};
	int K=0;

	variant->crispvar[variant->varalleles].AF_ML=0; variant->crispvar[variant->varalleles].delta=0;
	variant->crispvar[variant->varalleles].deltaf=0; variant->crispvar[variant->varalleles].deltar=0;
		
	// calculate initial estimate of population allele frequency
	for (i=0;i<variant->samples;i++)
	{
		K += variant->ploidy[i];
		C1=0,C2=0,C3=0; 
		for (k=0;k<3;k++) C1 += variant->indcounts[i][allele1+k*maxalleles]; 
		for (k=0;k<3;k++) C2 += variant->indcounts[i][allele2+k*maxalleles]; 
		for (k=0;k<3;k++) C3 += variant->indcounts[i][allele3+k*maxalleles]; 
		p1 = C2/(C1+C2+C3+0.001); p1 *= variant->ploidy[i]; // p1 is estimate of number of variant alleles in the pool (or diploid sample)

		if (p1 >= 0.2 && C2 >=4 ) p[0] += p1; // extra condition added that could work as fast filter for potential variants 
		else if (p1 >= 0.4 && C2 >=3 ) p[0] += p1;
		else if (variant->ploidy[i] ==2 && p1 >= 0.25 && C2 >=1) p[0] += p1;  // for diploid genomes loose filter 

		if (C3 >= 2) 
		{
			p2 = C3/(C1+C2+C3+0.001); p1 *= variant->ploidy[i];
			if (p2 >= 0.25 && C3 >=4 ) p[1] += p2; 
			else if (p2 >= 0.4 && C3 >=3 ) p[1] += p2;
			else if (variant->ploidy[i] ==2 && p2 >= 0.25 && C3 >= 1) p[1] += p2;  // for diploid genomes loose filter 
		}
	}
	//fprintf(stdout,"%d %d:%d %d:%d %f %d\n",i,variant->indcounts[i][allele2],variant->indcounts[i][allele2+maxalleles],variant->indcounts[i][allele1],variant->indcounts[i][allele1+maxalleles],p1,variant->ploidy[i]);
	double e01f_e10r=0,e01r_e10f=0;
	// FILTERS FOR REMOVING NON-VARIANT SITES || p[0] <= 0.5 filter loses variant sites in exome data
	if ( p[0] < 0.25 || (variant->ploidy[0] ==2 && p[0] <= 0.4)) return 0;
	else //if (p[0] <= variant->samples ) 
	{
		// if p is small, calculate the chi-square p-value to get another filtering step for rare variants in pooled data
		populate_contable(variant,allele1,allele2,variant->refbase,allele3);
		//chi2pvalue_stratified(variant,allele1,allele2,allele3,&pv1,chi2stats);

		// for low-coverage diploid data, the exact permutation test is called 
		// add additional check where we call permutation test if average coverage per allele is less than 2-3x 06/20/13
		permflag =0;
		// these two function calls do not account for bidirectional reads !! BUG 06/20/13
		if (variant->ploidy[0] > 2 && CHISQ_PERMUTATION ==1) 
		{
			permflag = pvalue_contable_iter_stranded(variant->ctable_stranded,variant->ctable_weighted,variant->strata,variant->samples,MAXITER,variant->newtable,&pv1,&pv2,chi2stats);
			if (permflag ==0) variant->ctpval[variant->varalleles] = pv2; else variant->ctpval[variant->varalleles] = pv1;
			variant->chi2pval[variant->varalleles] = pv2;
		}
		else if (variant->ploidy[0] ==2 && CHISQ_PERMUTATION ==1) 
		{
			struct acounts* atable = calloc(sizeof(struct acounts),variant->samples); 
			pvalue_lowcoverage(variant->ctable_stranded,variant->samples,MAXITER,atable,&pv1);
			variant->ctpval[variant->varalleles] = pv1; variant->chi2pval[variant->varalleles] = 0;
			free(atable); 
		}
		else 
		{
			if (PIVOTSAMPLE ==0) chi2pvalue_stratified(variant,allele1,allele2,allele3,&pv1,chi2stats);
			else chi2pvalue_stratified_pivotsample(variant,allele1,allele2,allele3,&pv1,chi2stats);
			
			variant->ctpval[variant->varalleles] = pv1; variant->chi2pval[variant->varalleles] = pv1;
		}

		e01r_e10f = (double)variant->filteredreads[5]/(variant->filteredreads[5]+variant->filteredreads[6]+variant->counts[allele2+2*maxalleles]+variant->counts[allele1+2*maxalleles]+0.1);
		e01f_e10r = (double)variant->filteredreads[6]/(variant->filteredreads[5]+variant->filteredreads[6]+variant->counts[allele2+2*maxalleles]+variant->counts[allele1+2*maxalleles]+0.1);
		fprintf(stdout,"FL=%d,%d,E01f_E10r=%0.4f,E01r_E10f=%0.4f SPV=%0.2f:%0.2f\t",variant->filteredreads[5],variant->filteredreads[6],e01f_e10r,e01r_e10f,chi2stats[3],chi2stats[4]);

		if (allele2 >=4) 
		{
			fprintf(stdout,"indel %d AL:%s:%f %s:%f FET:%0.2f:%0.2f\n",allele2,variant->itb[allele2],p[0],variant->itb[allele3],p[1],pv1,pv2);
			int hplength = variant->HPlength[allele2]/(strlen(variant->itb[allele2])-1);
			if (FAST_FILTER ==1 && variant->ctpval[variant->varalleles] >= -3 && variant->ploidy[0] > 2 && hplength >= 4 && (p[0]/K) <= 0.05) return 0;
			else if (FAST_FILTER ==1 & variant->ctpval[variant->varalleles] >= -3 && variant->ploidy[0] > 2 && hplength >= 8 && (p[0]/K) <= 0.1) return 0;
			// special filter for indels in homopolymer runs added 10/23/13
		}
		else fprintf(stdout,"AL:%s:%f %s:%f FET:%0.2f:%0.2f\n",variant->itb[allele2],p[0],variant->itb[allele3],p[1],pv1,pv2);

		// filter low-frequency variants that are not significant using chi-square distribution test
		if (FAST_FILTER ==1)
		{
			if (variant->ctpval[variant->varalleles] >= -3 && variant->ploidy[0] > 2 && p[0] <= variant->samples/4) return 0;
			else if (variant->ctpval[variant->varalleles] >= -2 && variant->ploidy[0] > 2 && p[0] <= variant->samples/2) return 0;
			else if (pv1  >= -2 && pv2 >= -2 && variant->ploidy[0] ==2 && p[0] <= 5 && p[0] <= variant->samples/4) return 0; // specialfilter for diploid
		}
	}
	paf[0] = p[0]; paf[1] = p[1]; 
	return 1;
	//else fprintf(stdout,"allele count estimate %s:%0.4f %s:%0.4f\n",variant->itb[allele2],p[0],variant->itb[allele3],p[1]);
}	


// code added 03/26/2015 to output indel allele counts to text file, --EM 2 option has to be used 
// also output SNP alleles, useful to see if some regions need haplotype analysis to examine combinations...
void output_allele_counts(REFLIST* reflist,int current,struct VARIANT* variant,ALLELE* maxvec,char* vtype)
{
	FILE* vfile = stdout;
	int allele=0,allele1 = maxvec[0].al, allele2 = maxvec[1].al, allele3 = maxvec[2].al, allele4=maxvec[3].al;
	int i=0,j=0; int printhomseq =0;
	for (j=1;j<4;j++)
	{
		allele = maxvec[j].al;
		if (allele >=4 &&  variant->itb[allele][0] == '+' && variant->HPlength[allele] > printhomseq) printhomseq= variant->HPlength[allele];
		else if (allele >=4 &&  variant->itb[allele][0] == '-' && (variant->HPlength[allele] + strlen(variant->itb[allele])-1 > printhomseq)) printhomseq= variant->HPlength[allele]+strlen(variant->itb[allele])-1;
	}

	// also print flanking sequence for indel allele 
	fprintf(vfile,"\ncandidate %s %s %d ",vtype,variant->chrom,variant->position); 
	for (j=0;j<4;j++)
	{
		allele = maxvec[j].al;
		if (maxvec[j].ct < 3 && j > 0 ) continue; // ignore alleles with less than 3 reads 
		fprintf(vfile,"%s ",variant->itb[allele]);
	}
	fprintf(vfile,"FLANKSEQ=");
	if (variant->position-21 >=0 && variant->position + printhomseq+20 <= reflist->lengths[current])
	{
		for (j=-21;j<-1;j++) fprintf(vfile,"%c",tolower(reflist->sequences[current][variant->position+j]));
		fprintf(vfile,":");
		for (j=-1;j<printhomseq;j++) fprintf(vfile,"%c",toupper(reflist->sequences[current][variant->position+j]));
		fprintf(vfile,":");
		for (j=printhomseq;j<printhomseq+20;j++) fprintf(vfile,"%c",tolower(reflist->sequences[current][variant->position+j]));
		fprintf(vfile,"\n");
	}
	else fprintf(vfile,".");

	for (j=0;j<4;j++)
	{
		allele = maxvec[j].al;
		if (maxvec[j].ct < 3  && j > 0) continue; // ignore alleles with less than 3 reads but always print reference 
		fprintf(vfile,"candidate "); //allele %s ",variant->itb[allele]);
		for (i=0;i<variant->samples;i++)
		{
			fprintf(vfile,"%d,%d,%d ",variant->indcounts[i][allele],variant->indcounts[i][allele+maxalleles],variant->indcounts[i][allele+2*maxalleles]);
		}
		fprintf(vfile,"\n");
	}
}

// aug 23 2012, return 0 -> no variant, 1 for SNP, 2 for indel, 3 for multi-allelic // maintain Ti/Tv ratio
// function called in main for pooled variant calling
int newCRISPcaller(REFLIST* reflist,int current,int position,READQUEUE* bq,struct BAMFILE_data* bd,struct VARIANT* variant,ALLELE* maxvec,FILE* vfile)
{
	int i=0,j=0,al=0,potentialvar=0,is_variant =0,variantpools=0;
	int allele=0, allele1 = maxvec[0].al, allele2 = maxvec[1].al, allele3 = maxvec[2].al, allele4=-1;
	double paf[2] = {0.0,0.0};
	int maxambiguity=0; int indelflag = 0;

	for (al=1;al<=5;al++) 
	{
		if (maxvec[al].ct < 3 ) continue; // ignore alleles with less than 3 reads, for low coverage data -> 2 reads could be good enough 07/03/13 !!
		allele = maxvec[al].al; variant->HPlength[allele] =0;
		// if there are no indels then the allele bases may not be even set
		if (current >=0 && allele >=4 && (variant->itb[allele][0] == '+' || variant->itb[allele][0] == '-')) calculate_indel_ambiguity(variant,reflist,current,allele);
		if (variant->HPlength[allele] > maxambiguity) maxambiguity = variant->HPlength[allele];
	}
	// we are removing ambigous bases even for SNP calling | BUG june 13 2013 
	//if (maxambiguity >0 ) remove_ambiguous_bases(bq,variant,maxambiguity,1); // this function is called after allele counts are determined, OPE issue

	int variant_output = 0;
	variant->strandpass =0; variant->varalleles=0;
	//need fast filtering function to evaluate potential alleles for being true variant 
	for (al=1;al<=5;al++)
	{
		if (maxvec[al].ct < 3 ) continue; // ignore alleles with less than 3 reads 
		allele1=maxvec[0].al; allele2 = maxvec[al].al;  allele3 = maxvec[al ==1 ? al+1: 1].al;
		variant->ctpval[variant->varalleles] = 0; potentialvar = 0; indelflag =0; is_variant =0;

		// remove ambigous reads from allele counts for indels prior to evaluating variant 
		if ((allele2 >= 4 && variant->HPlength[allele2] > 0) || (allele3 >=4 && variant->HPlength[allele3] >= 100)) 
		{
			indelflag =variant->HPlength[allele2]; 	remove_ambiguous_bases(bq,variant,indelflag,1);
		}

		potentialvar = evaluate_variant(bq,variant,allele1,allele2,allele3,paf); 
		if (potentialvar ==1) 
		{
			if (EMflag >=1) 
			{
				is_variant = compute_GLL_pooled(bq,variant,allele1,allele2,allele3,paf);
				if (is_variant ==1) 
				{
					variantpools = calculate_pool_statistics(bq,bd,variant,allele1,allele2); // function that uses non-genotype based filters to estimate # of variant pools 
					// if the number of pools passing filters =0 -> label variant as low confidence 
					variant->varpools[variant->varalleles] = variantpools; // new variable to store this info
					variant->alleles[variant->varalleles] = allele2;  
					variant->nvariants[variant->varalleles] = variantpools; // use it to store # of pools that pass filter
					variant->varalleles++;
				}
				if (EMflag >= 2 && allele2 >=4 && variant_output ==0) 
				{
					//fprintf(stdout,"indel candidate \n"); // this is an indel, we just print it to VCF  for realignment here 
					output_allele_counts(reflist,current,variant,maxvec,"INDEL"); // added 03/26/2015
					variant_output = 1; // only do it once per site
				}
			}
			else if (EMflag ==0) // use old pooled variant calling method
			{
				is_variant =0; variantpools =0; 
				is_variant = pooled_variantcaller(bq,bd,variant,allele1,allele2); 
			}
		}
		if (indelflag > 0) remove_ambiguous_bases(bq,variant,indelflag,-1); 
	}
	if (PFLAG >=1) fprintf(stdout,"\n");

	if (variant->varalleles >= 1)  // print VCF variant for pooled sequencing
	{
		if (EMflag ==3 && variant_output ==0) output_allele_counts(reflist,current,variant,maxvec,"SNP");  // also do it for SNPs

		// BUG here the allele counts will have been changed back to original count 07/02/13, run remove_ambiguous_bases again with max HP length
		indelflag =0;
		for (i=0;i<variant->varalleles;i++) 
		{
			if (variant->alleles[i] >= 4 && variant->HPlength[variant->alleles[i]] > indelflag) indelflag = variant->HPlength[variant->alleles[i]];
		}
		if (indelflag > 0) 
		{
			remove_ambiguous_bases(bq,variant,indelflag,1); fprintf(stdout,"\n");
		}
		//if (variant->alleles[0] >= 0) print_indel_haplotypes(bq,variant,variant->alleles[0]); 
			
		print_pooledvariant(variant,vfile,reflist,current,maxvec,EMflag);
		return 1;
	}
	return 0;
}


