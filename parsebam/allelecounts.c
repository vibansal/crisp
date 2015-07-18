#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include "variant.h"
#include "allelecounts.h"
#include "calculatelikelihoods.c"
//#include "contables.h"
//int USE_QV = 1; // added dec 21 2011 to allow for weighted table use in chi-square calculation int QVset = 0;

// compare two basecall (used for determining unique start site reads for variants)
int bcall_compare(const void* a1, const void* a2)
{
	struct alignedread **a = (struct alignedread **)a1; struct alignedread **b = (struct alignedread **)a2;
	if ( (*a)->strand == (*b)->strand ) return (*a)->l1 + (*a)->delta - (*b)->l1 - (*b)->delta;
	else return (*a)->strand - (*b)->strand;
}

// print variant reads for SOLID data july 30 2012 
void print_variantreads(READQUEUE* bq,struct BAMFILE_data* bamfiles_data,struct VARIANT* variant,int allele1,int allele2,int sample)
{
	int i=0,j=0;
	struct alignedread* bcall = bamfiles_data[sample].first;
	for (bcall = bamfiles_data[sample].first; bcall != NULL && bcall->position <= variant->position;bcall = bcall->nextread)
	{
		if (bcall->lastpos <= variant->position || bcall->filter == '1') continue; 
		if (bcall->quality[bcall->l1+bcall->delta] >= MINQ+QVoffset && bcall->type ==0 && allele2 <4)
		{
			if (toupper(bcall->sequence[bcall->l1+bcall->delta]) == variant->itb[allele2][0]) 
			{
				fprintf(stdout,"%s:%d:%d:%d:%d:%d:%d:",bcall->readid,bcall->strand,bcall->position,bcall->CM,bcall->IS,bcall->quality[bcall->l1+bcall->delta]-QVoffset,bcall->mquality);
				for (i=0;i<bcall->cigs;i++) fprintf(stdout,"%d%c",bcall->cigarlist[i]>>4,INT_CIGAROP[bcall->cigarlist[i]&0xf]);
				fprintf(stdout,"\t");
			}
		}
	}
	fprintf(stdout,"\n");
}

// for each sample calculate # of unique start site reads and middle reads for indels
void calculate_uniquereads(READQUEUE* bq,struct BAMFILE_data* bamfiles_data,struct VARIANT* variant,int allele1,int allele2,int sample,int* unique,int* middle)
{
	int UFLAG=1;
	int i=0,j=0,reads=0; *unique =1; *middle = 0;
	// we know the number of reads for allele2 from variant->indcounts so we can use that to allocate a static array 
	reads= variant->indcounts[sample][allele2] + variant->indcounts[sample][allele2+maxalleles] + variant->indcounts[sample][allele2+2*maxalleles];
	struct alignedread** bpointers = (struct alignedread**)malloc(sizeof(struct alignedread*)*(reads+1)); 
	struct alignedread* bcall = bamfiles_data[sample].first;
	for (bcall = bamfiles_data[sample].first; bcall != NULL && bcall->position <= variant->position;bcall = bcall->nextread)
	{
		if (bcall->lastpos <= variant->position || bcall->filter == '1') continue; 
		if (bcall->quality[bcall->l1+bcall->delta] >= MINQ+QVoffset && bcall->type ==0 && allele2 <4 && i < reads)
		{
			//for (j=0;j<bcall->fcigs;j+=2) fprintf(stdout,"%d%c",bcall->fullcigarlist[j],bcall->fullcigarlist[j+1]);
			//fprintf(stdout,"%s:%c:%d:%c:%d ",bcall->readid,bcall->strand,bcall->position,bcall->filter,bcall->l1+bcall->delta);
			if (toupper(bcall->sequence[bcall->l1+bcall->delta]) == variant->itb[allele2][0]) bpointers[i++] = bcall;
		}
		else if ((bcall->type > 0) && allele2 >=4 && i < reads) 
		{
			if (variant->itb[allele2][0] == '+' && strlen(variant->itb[allele2])-1 == bcall->type) bpointers[i++] = bcall;
		}
		else if ((bcall->type < 0) && allele2 >=4 && i < reads) 
		{
			if (variant->itb[allele2][0] == '-' && strlen(variant->itb[allele2])-1 == -1*bcall->type) bpointers[i++] = bcall;
		}
	}
	reads = i;
	if (reads > 2) qsort(bpointers,reads,sizeof(struct alignedread*),bcall_compare); 
	for (i=0;i<reads-1;i++) 
	{
		if (bpointers[i]->strand != bpointers[i+1]->strand || (bpointers[i]->l1+bpointers[i]->delta) != (bpointers[i+1]->l1+bpointers[i+1]->delta)) (*unique)++; 
	}
	for (i=0;i<reads;i++)
	{
		if ( bpointers[i]->l1+bpointers[i]->delta >= 5 && bpointers[i]->alignedbases-bpointers[i]->l1-bpointers[i]->delta >= 5) (*middle)++;
	}
	//if (allele2 >=4 && UFLAG) fprintf(stdout,"UAR:%d %d,%d:%d,%d ",sample,variant->indcounts[sample][allele1],variant->indcounts[sample][allele1+maxalleles],variant->indcounts[sample][allele2],variant->indcounts[sample][allele2+maxalleles]);
	//if(allele2 >= 4 && UFLAG && sample == variant->samples-1)fprintf(stdout,"\n");
	free(bpointers);
}

void sort_allelefreqs(struct VARIANT* variant,ALLELE* maxvec)
{
	double PAF = 0, WS = 0,pi,wi,TAC,p01; int i=0,j=0,k=0;
	for (i=0;i<maxalleles;i++)
	{
		PAF = 0; WS = 0; maxvec[i].vars = 0.0;
		for (j=0;j<variant->samples;j++)
		{
			pi = variant->indcounts[j][i]+ variant->indcounts[j][i+maxalleles]+variant->indcounts[j][i+2*maxalleles]; 
			TAC =0; for (k=0;k<maxalleles;k++) TAC += variant->indcounts[j][k] + variant->indcounts[j][k+maxalleles] + variant->indcounts[j][k+2*maxalleles]; 
			if (TAC ==0) continue; 
			p01= 0;//fake 
			//p01 = minreads_pvalue(TAC,pi,0.5/(double)poolsize);
			if (pi >= 2) maxvec[i].vars += pow(10,p01);

			pi /= TAC; wi = TAC*2/(TAC+1); PAF += pi*wi; WS += wi; 
		}
		PAF /= WS; 
		if (PAF <= 0.00001) PAF = 0.00001; // PAF is always greater than zero
		variant->allelefreqs[i] = PAF; 
		//if (PAF >= 0.0001) fprintf(stdout," PAF %d:%f %f\t",i,PAF,maxvec[i].vars);
	}

}

// new function for picall allele sorting, feb 12 2012, refbase is always first allele in sorted list of alleles 
int sort_allelecounts_diploid(struct VARIANT* variant, ALLELE* maxvec,int refbase)
{
	// alleles should be considered in order of population allele frequency rather than absolute read counts, dec 30 2011 
	sort_allelefreqs(variant,maxvec);
	int i=0,j=0,temp,r; double t1;
	for (i=0;i<maxalleles;i++) 
	{ 
		maxvec[i].ct = variant->counts[i]+variant->counts[i+maxalleles]+variant->counts[i+2*maxalleles]; 
		maxvec[i].al = i; 
	}
	// vars changed to population allele frequency... for sorting variant alleles feb 23 2012
	for (i=0;i<maxalleles;i++) maxvec[i].vars = variant->allelefreqs[i]; 

	if (refbase >= 1 && variant->IUPAC ==0) 
	{
		t1 =maxvec[0].vars; maxvec[0].vars = maxvec[refbase].vars; maxvec[refbase].vars = t1;
		temp =maxvec[0].ct; maxvec[0].ct = maxvec[refbase].ct; maxvec[refbase].ct = temp;
		temp =maxvec[0].al; maxvec[0].al = maxvec[refbase].al; maxvec[refbase].al = temp;
	}

	// sort the remaining alleles by counts
	r = 1; if (variant->IUPAC ==1) r = 0;
	for (i=r;i<maxalleles-1;i++)
	{
		for (j=i+1;j<maxalleles;j++)
		{
			if (maxvec[i].vars < maxvec[j].vars)  
				// changed jan 4 2012, order of sorting alleles, still a heuristic...
				// updated feb 10 2012 so that refbase always ends up as first in sorted list
			{
				t1 =maxvec[i].vars; maxvec[i].vars = maxvec[j].vars; maxvec[j].vars = t1;
				temp =maxvec[i].ct; maxvec[i].ct = maxvec[j].ct; maxvec[j].ct = temp;
				temp =maxvec[i].al; maxvec[i].al = maxvec[j].al; maxvec[j].al = temp;
			}
		}
	}
	if (variant->IUPAC ==1 && maxvec[0].al < 4) 
	{ 
		variant->refbase = maxvec[0].al; variant->refb = variant->itb[maxvec[0].al][0]; 
	} 
	else if (variant->IUPAC ==1) return 0;
	//	fprintf(stdout," position %d ALLELES %d %d %d %f %f\n",variant->position,refbase,maxvec[0].al,maxvec[1].al,maxvec[1].vars,maxvec[2].vars);
	return 1;
}

// this is for CRISP only 
int sort_allelecounts(struct VARIANT* variant, ALLELE* maxvec,int refbase)
{
	// alleles should be considered in order of population allele frequency rather than absolute read counts, dec 30 2011 
	sort_allelefreqs(variant,maxvec);
	int i=0,j=0,temp,r=0; double t1;
	//int potalleles[8] = {0,0,0,0,0,0,0,0}; 
	for (i=0;i<maxalleles;i++) 
	{ 
		maxvec[i].ct = variant->counts[i]+variant->counts[i+maxalleles]+variant->counts[i+2*maxalleles];
		maxvec[i].al = i; 
		if (i == refbase) r =i; 
	} 

	// changing code july 10 2012, reference is always first allele in maxvec 
	if (r >= 1 && variant->IUPAC ==0)
	{
		t1 =maxvec[0].vars; maxvec[0].vars = maxvec[r].vars; maxvec[r].vars = t1; // bug fixed jan 5 2012 
		temp = maxvec[0].ct; maxvec[0].ct = maxvec[r].ct; maxvec[r].ct = temp;
		temp = maxvec[0].al; maxvec[0].al = maxvec[r].al; maxvec[r].al = temp;
	}
	r = 1; if (variant->IUPAC ==1) r = 0;
	for (i=r;i<maxalleles-1;i++)
	{
		for (j=i+1;j<maxalleles;j++)
		{
			if (maxvec[i].ct < maxvec[j].ct )  // changed jan 4 2012, order of sorting alleles, still a heuristic...
			{
				t1 =maxvec[i].vars; maxvec[i].vars = maxvec[j].vars; maxvec[j].vars = t1;
				temp =maxvec[i].ct; maxvec[i].ct = maxvec[j].ct; maxvec[j].ct = temp;
				temp =maxvec[i].al; maxvec[i].al = maxvec[j].al; maxvec[j].al = temp;
			}
		}
	}
	if (variant->IUPAC ==1 && maxvec[0].al < 4) 
	{ 
		variant->refbase = maxvec[0].al; variant->refb = variant->itb[maxvec[0].al][0]; 
	} 
	else if (variant->IUPAC ==1) return 0;
	//fprintf(stdout,"IUPAC %d %d \n",variant->IUPAC,variant->refbase);
	// ensure that the reference base is always either the first or second allele in sorted order in maxvec, WHY ?? jan 4 2012 
	/*
	r=0; for (i=0;i<maxalleles;i++) { if (maxvec[i].al == refbase) r=i; }
	if (r >=2  && r <0) // old code
	{
		t1 =maxvec[1].vars; maxvec[1].vars = maxvec[r].vars; maxvec[r].vars = t1; // bug fixed jan 5 2012 
		temp = maxvec[1].ct; maxvec[1].ct = maxvec[r].ct; maxvec[r].ct = temp;
		temp = maxvec[1].al; maxvec[1].al = maxvec[r].al; maxvec[r].al = temp;
	}*/
	return 1;
}

// this function is no longer correct since bidirectional bases are not accounted for nov 30 2012
// allele3 is also added to counts, is this correct 06/28/2013 
int populate_contable(struct VARIANT* variant,int allele1,int allele2,int refbase,int allele3)
{
	int i=0,asymptotic_table=1,c0=0,c1=0;
	for (i=0;i<variant->samples;i++)
	{
		variant->ctable_stranded[i][4] = 0; variant->ctable_stranded[i][5]=0;
		c0 = variant->indcounts[i][allele1] + variant->indcounts[i][maxalleles+allele1]+variant->indcounts[i][allele1+2*maxalleles];
		if (allele3 != -1) c0 += variant->indcounts[i][allele3] + variant->indcounts[i][maxalleles+allele3] + variant->indcounts[i][allele3+2*maxalleles];
		c1 = variant->indcounts[i][allele2] + variant->indcounts[i][maxalleles+allele2]+variant->indcounts[i][allele2+2*maxalleles];
		if (c0 < 5 || c1 < 5) asymptotic_table = 0;
		if (allele1 == refbase || (allele1 != refbase && allele2 != refbase))  //  c1 is always alternate base  
		{
			variant->ctable[i][0] = c0+c1; variant->ctable[i][1] = c1;
			variant->ctable_stranded[i][0] = variant->indcounts[i][allele1] + variant->indcounts[i][allele2];
			// adding allele3 counts to allele1 
			if (allele3 != -1) variant->ctable_stranded[i][0] += variant->indcounts[i][allele3];
			variant->ctable_stranded[i][2] = variant->indcounts[i][allele1+maxalleles] + variant->indcounts[i][allele2+maxalleles];
			if (allele3 != -1) variant->ctable_stranded[i][2] += variant->indcounts[i][allele3+maxalleles];
			variant->ctable_stranded[i][1] = variant->indcounts[i][allele2];
			variant->ctable_stranded[i][3] = variant->indcounts[i][allele2+maxalleles];
			variant->ctable_stranded[i][4] = variant->indcounts[i][allele1+2*maxalleles] + variant->indcounts[i][allele2+2*maxalleles];
			variant->ctable_stranded[i][5] = variant->indcounts[i][allele2+2*maxalleles]; // added for bidirectional allele counts..

			variant->ctable_weighted[i][0] = variant->Qhighcounts[i][allele1] + variant->Qhighcounts[i][allele2];
			if (allele3 != -1) variant->ctable_weighted[i][0] += variant->Qhighcounts[i][allele3];
			variant->ctable_weighted[i][2] = variant->Qhighcounts[i][allele1+maxalleles] + variant->Qhighcounts[i][allele2+maxalleles];
			if (allele3 != -1) variant->ctable_weighted[i][2] += variant->Qhighcounts[i][allele3+maxalleles];
			variant->ctable_weighted[i][1] = variant->Qhighcounts[i][allele2];
			variant->ctable_weighted[i][3] = variant->Qhighcounts[i][allele2+maxalleles];
		}
		else
		{
			//fprintf(stderr,"ref allele is 2nd \n");
			variant->ctable[i][0] = c0+c1; variant->ctable[i][1] = c0;
			variant->ctable_stranded[i][0] = variant->indcounts[i][allele1] + variant->indcounts[i][allele2];
			if (allele3 != -1) variant->ctable_stranded[i][0] += variant->indcounts[i][allele3];
			variant->ctable_stranded[i][1] = variant->indcounts[i][allele1];
			variant->ctable_stranded[i][2] = variant->indcounts[i][allele1+maxalleles] + variant->indcounts[i][allele2+maxalleles];
			if (allele3 != -1) variant->ctable_stranded[i][2] += variant->indcounts[i][allele3+maxalleles];
			variant->ctable_stranded[i][3] = variant->indcounts[i][allele1+maxalleles];
			variant->ctable_stranded[i][4] = variant->indcounts[i][allele1+2*maxalleles] + variant->indcounts[i][allele2+2*maxalleles];
			variant->ctable_stranded[i][5] = variant->indcounts[i][allele1+2*maxalleles]; // added for bidirectional allele counts..

			variant->ctable_weighted[i][0] = variant->Qhighcounts[i][allele1] + variant->Qhighcounts[i][allele2];
			if (allele3 != -1) variant->ctable_weighted[i][0] += variant->Qhighcounts[i][allele3];
			variant->ctable_weighted[i][1] = variant->Qhighcounts[i][allele1];
			variant->ctable_weighted[i][2] = variant->Qhighcounts[i][allele1+maxalleles] + variant->Qhighcounts[i][allele2+maxalleles];
			if (allele3 != -1) variant->ctable_weighted[i][2] += variant->Qhighcounts[i][allele3+maxalleles];
			variant->ctable_weighted[i][3] = variant->Qhighcounts[i][allele1+maxalleles];

		}
	}
	return asymptotic_table;
}

