#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>

/* functions to initialize the VARIANT data structure (once for allocation and once for every potential variant site) */

void init_variant(struct VARIANT* variant,int bamfiles,int samples)
{
	int i=0,j=0;
	variant->samples = samples; variant->bamfiles= bamfiles;
	variant->itb = (char**)malloc(sizeof(char*)*maxalleles);
	//variant->haplotypes = (HAPLOTYPE*)malloc(sizeof(HAPLOTYPE)*maxalleles); variant->haps = 0;
	for (i=0;i<maxalleles;i++) variant->itb[i] = (char*)malloc(1024);
	strcpy(variant->itb[0],"A"); strcpy(variant->itb[1],"C"); strcpy(variant->itb[2],"G"); strcpy(variant->itb[3],"T");  
	//strcpy(variant->itb[4],"D");

	variant->NCarray = calloc(variant->samples,sizeof(double*)); 
	for (i=0;i<variant->samples;i++) 
	{
		variant->NCarray[i] = calloc(variant->ploidy[i]+1,sizeof(double)); 
		variant->NCarray[i][0] = 0.0; variant->NCarray[i][variant->ploidy[i]] = 0.0;
		for (j=1;j<variant->ploidy[i];j++) variant->NCarray[i][j] = variant->NCarray[i][j-1] + log10(variant->ploidy[i]-j+1) - log10(j);
	}

	for (i=0;i<8;i++)
	{
		variant->crispvar[i].AC = calloc(sizeof(int),variant->samples);
		variant->crispvar[i].QV = calloc(sizeof(double),variant->samples);
		variant->crispvar[i].meanAC = calloc(sizeof(double),variant->samples);
		variant->crispvar[i].varAC = calloc(sizeof(double),variant->samples);
	}	

	variant->chrom = (char*)malloc(1024);
	variant->strata = (int*)malloc(sizeof(int*)*variant->samples);
	variant->readdepths = (int*)malloc(sizeof(int*)*variant->samples);
	variant->indcounts = (int**)malloc(sizeof(int*)*variant->samples);
	variant->indcounts_binned = (int***)malloc(sizeof(int**)*variant->samples);
	variant->Qhighcounts = (double**)malloc(sizeof(double*)*variant->samples);
	variant->stats = (double**)malloc(sizeof(double*)*variant->samples);
	variant->qpvalue = (double**)malloc(sizeof(double*)*variant->samples);
	variant->ctable = (int**)malloc(sizeof(int*)*variant->samples);		
	variant->newtable = (int**)malloc(sizeof(int*)*variant->samples);		
	variant->ctable_stranded = (int**)malloc(sizeof(int*)*variant->samples);
	variant->ctable_weighted = (double**)malloc(sizeof(double*)*variant->samples);
	variant->genotypes = (int*)malloc(sizeof(int)*variant->samples);
	variant->poolpv = (double**)malloc(sizeof(double*)*variant->samples);
	variant->genotypeQV = (double*)malloc(sizeof(double)*variant->samples);
	variant->GENLL = (double**)malloc(sizeof(double*)*variant->samples);
	variant->GENLLf = (double**)malloc(sizeof(double*)*variant->samples);
	variant->GENLLr = (double**)malloc(sizeof(double*)*variant->samples);
	variant->GENLLb = (double**)malloc(sizeof(double*)*variant->samples);
	variant->gpriors = (double**)malloc(sizeof(double*)*variant->samples);
	for (i=0;i<variant->samples;i++) 
	{
		variant->strata[i] = 0;
		variant->genotypes[i] =0;
		variant->readdepths[i] = 0; // used for new caller
		variant->ctable[i] = (int*)malloc(sizeof(int)*6);
		variant->newtable[i] = (int*)malloc(sizeof(int)*4);
		variant->ctable_stranded[i] = (int*)malloc(sizeof(int)*6);
		variant->ctable_weighted[i] = (double*)malloc(sizeof(double)*12); // 2 alleles, 2 strands and 3 quality value bins
		variant->qpvalue[i] = (double*)malloc(sizeof(double)*4);
		variant->indcounts[i] = (int*)malloc(sizeof(int)*4*maxalleles); // for two strands, bidirectional counts in two bins
		variant->indcounts_binned[i] = (int**)malloc(sizeof(int*)*2*maxalleles); 
		for (j=0;j<2*maxalleles;j++) 
		{
			variant->indcounts_binned[i][j] = (int*)malloc(sizeof(int)*6); // 6 bins for each strand and individual
		}

		variant->Qhighcounts[i] = (double*)malloc(sizeof(double)*6*maxalleles); // 3 bins for quality values
		variant->stats[i] = (double*)malloc(sizeof(double)*2*maxalleles);
		variant->poolpv[i] = (double*)malloc(sizeof(double)*maxalleles);
		if (variant->ploidy[i] ==2) 
		{
			variant->GENLL[i] = (double*)malloc(sizeof(double)*10); // maximum of 4 alleles 
			variant->GENLLf[i] = (double*)malloc(sizeof(double)*10); // maximum of 4 alleles 
			variant->GENLLr[i] = (double*)malloc(sizeof(double)*10); // maximum of 4 alleles 
			variant->GENLLb[i] = (double*)malloc(sizeof(double)*10); // maximum of 4 alleles 
			variant->gpriors[i] = (double*)malloc(sizeof(double)*10); // maximum of 4 alleles 
		}
		else // for pooled variant calling
		{
			variant->GENLL[i] = (double*)malloc(sizeof(double)*(variant->ploidy[i]+1)*5); 
			variant->GENLLf[i] = (double*)malloc(sizeof(double)*(variant->ploidy[i]+1)*5); 
			variant->GENLLr[i] = (double*)malloc(sizeof(double)*(variant->ploidy[i]+1)*5); 
			variant->GENLLb[i] = (double*)malloc(sizeof(double)*(variant->ploidy[i]+1)*5); // maximum of 4 alleles 
			variant->gpriors[i] = (double*)malloc(sizeof(double)*(variant->ploidy[i]+1)*5); 
		}
		//variant->genotypeQV[i] = (int*)malloc(sizeof(int)*maxalleles);
		for (j=0;j<4*maxalleles;j++) variant->indcounts[i][j] = 0;
		for (j=0;j<6*maxalleles;j++) variant->Qhighcounts[i][j] = 0.0;
		for (j=0;j<2*maxalleles;j++) variant->stats[i][j] = 0;

	}
	variant->counts = (int*)malloc(sizeof(int)*4*maxalleles);
	variant->tstats = (double*)malloc(sizeof(double)*2*maxalleles);
	variant->allelefreqs = (double*)malloc(sizeof(double)*maxalleles); 

}

// initialize the counts and determine depth of coverage for each sample
int initialize_variant(REFLIST* reflist,int current, int position,struct VARIANT* variant)
{
	int i=0,j=0,p=0,q=0; 

	variant->refb = toupper(reflist->sequences[current][position]); // set the refbase 
	if (variant->refb != 'A' && variant->refb != 'C' && variant->refb != 'G' && variant->refb != 'T' && variant->refb != 'N') variant->IUPAC = 1;
	else variant->IUPAC = 0;
	variant->refbase = BTI[variant->refb]; // what if this is IUPAC ambig. base 
	strcpy(variant->chrom,reflist->names[current]); variant->position = position;  // add 1 to account for 0 offset
	for (i=0;i<10;i++) { variant->HPlength[i] = 0; }
	variant->varalleles =0;
	for (i=0;i<maxalleles;i++) variant->allelefreqs[i] = 0.0;
	for (q=0;q<4*maxalleles;q++) variant->counts[q] = 0;  
	for (q=0;q<2*maxalleles;q++) variant->tstats[q] = 0.0;  
	// if reference allele is 'N', no need to call variant at this position may 3 2012
	for (i=0;i<4;i++) variant->MQcounts[i] =0; 
	for (i=0;i<8;i++) variant->filteredreads[i]=0; 
	for (i=0;i<variant->samples;i++)
	{
		variant->readdepths[i] =0;
		for (q=0;q<4*maxalleles;q++) variant->indcounts[i][q] = 0; 
		for (q=0;q<2*maxalleles;q++) 
		{
			for (p=0;p<6;p++) variant->indcounts_binned[i][q][p] = 0; 
		}
		for (q=0;q<2*maxalleles;q++) variant->stats[i][q] = 0.0;
		for (q=0;q<6*maxalleles;q++) variant->Qhighcounts[i][q] = 0.0;
	}
}


void init_poolsizes(struct VARIANT* variant, struct OPTIONS* options,int VC_METHOD)
{
        int i=0;
        if ((VC_METHOD ==0 || VC_METHOD ==3) && options->varpoolsize ==0)
        {
                for (i=0;i<options->bamfiles;i++) variant->ploidy[i] = options->POOLSIZE; variant->varpoolsize = 0;
                fprintf(stderr,"poolsize for each sample is %d \n",variant->ploidy[0]);
        }
        else if ((VC_METHOD ==0 || VC_METHOD ==3) && options->varpoolsize ==1)
        {
                for (i=0;i<options->bamfiles;i++) variant->ploidy[i] = options->ploidy[i];
                variant->varpoolsize =0;
                for (i=0;i<options->bamfiles-1;i++)
                {
                        if (variant->ploidy[i] != variant->ploidy[i+1]) variant->varpoolsize = 1;
                }
                fprintf(stderr,"poolsizes are variable and read from file: %d...%d\n",options->ploidy[0],options->ploidy[options->bamfiles-1]);
        }
        else
        {
                for (i=0;i<options->bamfiles;i++) variant->ploidy[i] = 2;  variant->varpoolsize =0;
        }
        //for (i=0;i<options->bamfiles;i++) fprintf(stdout,"i %d %d \n",i,variant.ploidy[i]);
        //        //if (options->samples > 0) init_variant(&variant,options->bamfiles,options->samples);  // init should be done with number of samples not bamfiles 
}
