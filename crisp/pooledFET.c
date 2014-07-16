
/*
   this file has the functions for calculating p-value for pooled contingency table directly from allele counts variant->indcounts 
   1. for indels, there is no distribution based on base quality values (maybe use average QV of neighboring bases ?)
   2. we can calculate p-value for multiple alleles simultaneously, conditional p-values keeping allele counts for other alleles fixed...
 */

#include<stdio.h>
#include<math.h>
#include<stdlib.h>

// calculate permutation based p-value for variant allele counts
// general permutation that works for multiple bins and multiple alleles...
// conditional test where we permute the third allele keeping second allele constant and vice versa....
// bins should be based on strand (+/-), first/2nd read (1/2) and quality value (30+ and remaining)
// for indels, we shouldn't do quality values based bins but mapping based bins
// to speed up permutation based p-value:  stop as soon as clranew exceeds clra
int FETpvalue_variant(struct VARIANT* variant,int maxiter,double* permutationpvalue,int allele1,int allele2,int allele3,int allele4)
{
	int i=0,j=0,k=0,r=0,b=0,pvalue=0,offset=0,iter=0,R=0;
//	fprintf(stdout,"calling new permutation test alleles %d:%d:%d:%d\n",allele1,allele2,allele3,allele4);
	double clra =0, clranew=0;
	int** table_index = (int**)malloc(sizeof(int*)*12);
	// we have 12 bins (2 for each strand x (3 for QV + 3 for R12_Single)) within which we permute 
	int** newtable = calloc(variant->samples,sizeof(int*)); 
	static int allelecounts[12][4]; // for total, allele2, allele3, allele4
	double* TLL = calloc(maxiter,sizeof(double)); // random table likelihoods
	for (i=0;i<12;i++) 
	{
		for (j=0;j<4;j++) allelecounts[i][j] = 0;
	}

	for (i=0;i<variant->samples;i++) 
	{
		for (j=0;j<6;j++) 
		{
			allelecounts[j][0] += variant->indcounts_binned[i][allele1][j] + variant->indcounts_binned[i][allele2][j] + variant->indcounts_binned[i][allele3][j];
			allelecounts[j][1] += variant->indcounts_binned[i][allele2][j];
			allelecounts[j][2] += variant->indcounts_binned[i][allele3][j];
			allelecounts[j+6][0] += variant->indcounts_binned[i][allele1+maxalleles][j] + variant->indcounts_binned[i][allele2+maxalleles][j] + variant->indcounts_binned[i][allele3+maxalleles][j];
			allelecounts[j+6][1] += variant->indcounts_binned[i][allele2+maxalleles][j];
			allelecounts[j+6][2] += variant->indcounts_binned[i][allele3+maxalleles][j];
		}
		R = variant->indcounts[i][allele1]+variant->indcounts[i][allele1+maxalleles] + variant->indcounts[i][allele2]+variant->indcounts[i][allele2+maxalleles] + variant->indcounts[i][allele3]+variant->indcounts[i][allele3+maxalleles]; 
		newtable[i] = calloc(4,sizeof(int)); newtable[i][0] = R; 
		//clra += multinomial(R,variant->indcounts[i][allele2]+variant->indcounts[i][allele2+maxalleles],variant->indcounts[i][allele3]+variant->indcounts[i][allele3+maxalleles]);
		clra += ncr(R,variant->indcounts[i][allele2]+variant->indcounts[i][allele2+maxalleles]);
	}

	for (j=0;j<12;j++) table_index[j] = (int*)malloc(sizeof(int)*(allelecounts[j][0]+1));
	for (j=0;j<12;j++)
	{
		offset =0; 
		for (i=0;i<variant->samples;i++) 
		{
			if (j < 6) R = variant->indcounts_binned[i][allele1][j] + variant->indcounts_binned[i][allele2][j] + variant->indcounts_binned[i][allele3][j]; 
			else R = variant->indcounts_binned[i][allele1+maxalleles][j-6] + variant->indcounts_binned[i][allele2+maxalleles][j-6] + variant->indcounts_binned[i][allele3+maxalleles][j-6];
			for (k=0;k<R;k++) table_index[j][offset+k] =i; offset += R;
		}

	}
	for (iter=1;iter<=maxiter;iter++)
	{
		for (i=0;i<variant->samples;i++) { newtable[i][1] = 0; newtable[i][2] = 0; newtable[i][3] = 0; }
		clranew=0;
		// do the bins for quality values
		for (j=0;j<3;j++)
		{
			for (i=0;i<allelecounts[j][1];i++) 
			{
				// freedom to permute among all alleles rather than just allele1/2 // check this 
				r = (int)(drand48()*(allelecounts[j][0]-i))+i; 
				b = table_index[j][(int)r]; table_index[j][(int)r] = table_index[j][i]; table_index[j][i]=b;
				//fprintf(stdout,"clra %f %d:%d\n",clranew,newtable[i][0],newtable[b][1]);
				clranew += log10(newtable[b][0]-newtable[b][1])- log10(newtable[b][1]+1); newtable[b][1]++;
			}
			for (i=0;i<allelecounts[j+6][1];i++) 
			{
				r = (int)(drand48()*(allelecounts[j+6][0]-i))+i; 
				b = table_index[j+6][(int)r]; table_index[j+6][(int)r] = table_index[j+6][i]; table_index[j+6][i]=b;
				clranew += log10(newtable[b][0]-newtable[b][1])- log10(newtable[b][1]+1); newtable[b][1]++;
			}
		}
		if (clranew <= clra+epsilon) 	pvalue +=1;  
		TLL[iter-1] = clranew;
		if (pvalue >= 10 || (iter <= 1000 && pvalue >= 5)) break; 
	}

	//qsort(TLL,iter-1,sizeof(double),doublecomp);
	//fprintf(stdout,"\nTLL %d iter %d clra %f | ",variant->position,iter,clra);
	//for (i=0;i<iter-1;i++) fprintf(stdout,"%f ",TLL[i]); fprintf(stdout,"\n");
	// sort the TLL array, find peak of PDF, take left half of it, it is approximated as N(0,1) with mean at maxima PDF

	free(TLL);	
	for (i=0;i<12;i++) free(table_index[i]); free(table_index);
	for (i=0;i<variant->samples;i++) free(newtable[i]); free(newtable);
	//fprintf(stdout,"clra %f mean %f sq %f %d\n",clra,clrasum,clrasumsq,k+1);
	fprintf(stdout,"table %d NEWTEST:%d/%d/%f\t",variant->position,pvalue,iter,clra);
	*permutationpvalue = log10(pvalue+1)-log10(iter+1); 
	return 1;
}

// this works for only two strands (indels)
int FETpvalue_indels(struct VARIANT* variant,int maxiter,double* permutationpvalue,int allele1,int allele2,int allele3)
{
	int i=0,j=0,k=0,r=0,b=0,pvalue=0,offset=0,iter=0,R=0;
	double clra =0, clranew=0;
	int** table_index = calloc(2,sizeof(int*));
	int** newtable = calloc(variant->samples,sizeof(int*)); 
	static int allelecounts[2][4]; // for total, allele2, allele3, allele4
	for (i=0;i<2;i++) 
	{
		for (j=0;j<4;j++) allelecounts[i][j] = 0;
	}

	for (i=0;i<variant->samples;i++) 
	{
		allelecounts[0][0] += variant->indcounts[i][allele1] + variant->indcounts[i][allele2] + variant->indcounts[i][allele3];
		allelecounts[0][1] += variant->indcounts[i][allele2];
		allelecounts[0][2] += variant->indcounts[i][allele3];
		allelecounts[1][0] += variant->indcounts[i][allele1+maxalleles] + variant->indcounts[i][allele2+maxalleles] + variant->indcounts[i][allele3+maxalleles];
		allelecounts[1][1] += variant->indcounts[i][allele2+maxalleles];
		allelecounts[1][2] += variant->indcounts[i][allele3+maxalleles];
		R = variant->indcounts[i][allele1]+variant->indcounts[i][allele1+maxalleles] + variant->indcounts[i][allele2]+variant->indcounts[i][allele2+maxalleles] + variant->indcounts[i][allele3]+variant->indcounts[i][allele3+maxalleles]; 
		newtable[i] = calloc(4,sizeof(int)); newtable[i][0] = R; 
		//clra += multinomial(R,variant->indcounts[i][allele2]+variant->indcounts[i][allele2+maxalleles],variant->indcounts[i][allele3]+variant->indcounts[i][allele3+maxalleles]);
		clra += ncr(R,variant->indcounts[i][allele2]+variant->indcounts[i][allele2+maxalleles]);
	}

	for (j=0;j<2;j++)
	{
		table_index[j] = (int*)malloc(sizeof(int)*(allelecounts[j][0]+1));
		offset =0; 
		for (i=0;i<variant->samples;i++) 
		{
			if (j ==0) R = variant->indcounts[i][allele1] + variant->indcounts[i][allele2] + variant->indcounts_binned[i][allele3]; 
			else R = variant->indcounts[i][allele1+maxalleles] + variant->indcounts[i][allele2+maxalleles] + variant->indcounts[i][allele3+maxalleles];
			for (k=0;k<R;k++) table_index[j][offset+k] =i; offset += R;
		}

	}
	for (iter=1;iter<=maxiter;iter++)
	{
		for (i=0;i<variant->samples;i++) { newtable[i][1] = 0; newtable[i][2] = 0; newtable[i][3] = 0; }
		clranew=0;
		for (j=0;j<2;j++)
		{
			for (i=0;i<allelecounts[j][1];i++) 
			{
				// freedom to permute among all alleles rather than just allele1/2 // check this 
				r = (int)(drand48()*(allelecounts[j][0]-i))+i; 
				b = table_index[j][(int)r]; table_index[j][(int)r] = table_index[j][i]; table_index[j][i]=b;
				//fprintf(stdout,"clra %f %d:%d\n",clranew,newtable[i][0],newtable[b][1]);
				clranew += log10(newtable[b][0]-newtable[b][1])- log10(newtable[b][1]+1); newtable[b][1]++;
			}
		}
		if (clranew <= clra+epsilon) 	pvalue +=1;  
		if (pvalue >= 10 || (iter <= 1000 && pvalue >= 5)) break; 
	}

	for (i=0;i<2;i++) free(table_index[i]); free(table_index);
	for (i=0;i<variant->samples;i++) free(newtable[i]); free(newtable);
	//fprintf(stdout,"clra %f mean %f sq %f %d\n",clra,clrasum,clrasumsq,k+1);
	fprintf(stdout,"table %d NEWTEST:%d/%d/%f\t",variant->position,pvalue,iter,clra);
	*permutationpvalue = log10(pvalue+1)-log10(iter+1); 
	return 1;
}


/*
for chi-square we need to collapse bins when we have low count in some bin 
first check strand bin (if small collapse onto other strand) 
then check within strand bins: Q10-20 onto Q20-A30 and so on
*/

int chi2pvalue_variant(struct VARIANT* variant,int maxiter,double* chi2pvalue,int allele1,int allele2,int allele3,int allele4)
{
	int i=0,j=0;
	static int allelecounts[12][4]; // for total, allele2, allele3, allele4
	static int parentbin[12]; 
	static double E[12];
	double correction = 0.5;
	int MINREADS = 50; //int collapsef =0,collapser=0; 
	for (i=0;i<12;i++) 
	{
		for (j=0;j<4;j++) allelecounts[i][j] = 0;
		parentbin[i] = i; // if it is collapsed, then we change this parent
	}
	for (i=0;i<variant->samples;i++) 
	{
		for (j=0;j<6;j++) 
		{
			allelecounts[j][0] += variant->indcounts_binned[i][allele1][j] + variant->indcounts_binned[i][allele2][j];// + variant->indcounts_binned[i][allele3][j];
			allelecounts[j][1] += variant->indcounts_binned[i][allele2][j];
			allelecounts[j][2] += variant->indcounts_binned[i][allele3][j];
			allelecounts[j+6][0] += variant->indcounts_binned[i][allele1+maxalleles][j] + variant->indcounts_binned[i][allele2+maxalleles][j];// + variant->indcounts_binned[i][allele3+maxalleles][j];
			allelecounts[j+6][1] += variant->indcounts_binned[i][allele2+maxalleles][j];
			allelecounts[j+6][2] += variant->indcounts_binned[i][allele3+maxalleles][j];
		}
	}

	if (variant->counts[allele1]+variant->counts[allele2] + variant->counts[allele3] < MINREADS) 
	{
		// combine reads from first 6 bins with last 6 bins 
		for (j=0;j<6;j++) 
		{ 
			for (i=0;i<4;i++) { allelecounts[j+6][i] += allelecounts[j][i];  allelecounts[j][i] = 0; }
			parentbin[j] = parentbin[j+6]; 
		} 
	}
	else if (variant->counts[allele1+maxalleles]+variant->counts[allele2+maxalleles] + variant->counts[allele3+maxalleles] < MINREADS) 
	{
		for (j=0;j<6;j++) 
		{ 
			for (i=0;i<4;i++) { allelecounts[j][i] += allelecounts[j+6][i];  allelecounts[j+6][i] = 0;} 
			parentbin[j+6] = parentbin[j]; 
		} 
	} 
	for (j=0;j<2;j++)
	{
		if (allelecounts[j][0] < MINREADS) 
		{
			for (i=0;i<4;i++) { allelecounts[j+1][i] += allelecounts[j][i]; allelecounts[j][i] = 0; }
			parentbin[j] = parentbin[j+1]; 
		}
		if (allelecounts[j+6][0] < MINREADS) 
		{
			for (i=0;i<4;i++) { allelecounts[j+7][i] += allelecounts[j+6][i]; allelecounts[j+6][i] = 0;}
			parentbin[j+6] = parentbin[j+7]; 
		}
	}
	for (j=0;j<12;j++) 
	{
		if (parentbin[j]==j) E[j] = ((double)allelecounts[j][1]+correction)/(allelecounts[j][0]+correction);
	}
	double delta=0,tf=0,chi2stat=0,chi2pval=0; 
	for (i=0;i<variant->samples;i++)
	{
		delta = 0; tf = 0; 
		for (j=0;j<3;j++) 
		{
			delta += variant->indcounts_binned[i][allele2][j]- E[parentbin[j]]*(variant->indcounts_binned[i][allele1][j] + variant->indcounts_binned[i][allele2][j]);
			tf += E[parentbin[j]]*(1.0-E[parentbin[j]])*(variant->indcounts_binned[i][allele1][j] + variant->indcounts_binned[i][allele2][j]);
			delta += variant->indcounts_binned[i][allele2+maxalleles][j]- E[parentbin[j+6]]*(variant->indcounts_binned[i][allele1+maxalleles][j] + variant->indcounts_binned[i][allele2+maxalleles][j]);
			tf += E[parentbin[j+6]]*(1.0-E[parentbin[j+6]])*(variant->indcounts_binned[i][allele1+maxalleles][j] + variant->indcounts_binned[i][allele2+maxalleles][j]);

		}
		//if(i==0)fprintf(stdout,"chi %f delta %f tf %f\n",chi2stat,delta,tf);
		chi2stat += delta*delta/tf; 
	}
	for (j=0;j<3;j++) 
	{
	//	fprintf(stdout,"\ncounts:%d %d %d/%d %f |",j,parentbin[j],allelecounts[j][1],allelecounts[j][0],E[parentbin[j]]);
	//	fprintf(stdout,"%d %d %d/%d %f\n",j+6,parentbin[j+6],allelecounts[j+6][1],allelecounts[j+6][0],E[parentbin[j+6]]);
	}
	if (chi2stat > 0) chi2pval = kf_gammaq((double)(variant->samples-1)/2,(double)chi2stat/2);
	fprintf(stdout,"statistic %f %f\n",chi2stat,chi2pval);
	
}


/* general function that can handle multiple alleles (3-4) instead of just 2 as well as multiple bins */
struct ALLELECOUNTS
{
	int pools; int alleles; int bins; // bins correspond to strand, quality values 
	// data for N pools 
	int*** counts; // corresponding to counts for each allele and each bin and N samples 
	int*** randcounts; // counts for doing permutation test
	int* totalcounts; // totalcounts across all alleles and all bins for each pool
	int* bincounts; // counts of reads in each bin across all pools
};
int doublecomp(const void* elem1, const void* elem2)
{
    if(*(const double*)elem1 < *(const double*)elem2)
        return -1;
    return *(const double*)elem1 > *(const double*)elem2;
}
