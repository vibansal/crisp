
#ifndef INC_allelecounts_H
#define INC_allelecounts_H
#include <stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include "readfasta.h"
#include "variant.h"

//extern int USE_QV, QVset; // added dec 21 2011 to allow for weighted table use in chi-square calculation

//int compute_counts(struct pileup* pile,int* counts,float* Qhighcounts,int* indcounts,float* stats,float* tstats,int** poscounts,int posfilter,int sample,struct VARIANT* variant);
//int indelanalysis(struct pileup** pilelist,int noofsamples,int* counts,int** indcounts,float** Qhighcounts,int*** poscounts,char** indel);

//int compute_GLLs(struct VARIANT* variant,int posfilter,int allele1, int allele2,int allele3);
void compute_GLLs(READQUEUE* bq,struct BAMFILE_data* bd,struct VARIANT* variant,int allele1,int allele2,int allele3);

//void compute_GLL_pooled(READQUEUE* bq, struct VARIANT* variant,int allele1,int allele2,int allele3);

void print_variantreads(READQUEUE* bq,struct BAMFILE_data* bamfiles_data,struct VARIANT* variant,int allele1,int allele2,int sample);

void calculate_uniquereads(READQUEUE* bq,struct BAMFILE_data* bd, struct VARIANT* variant,int allele1,int allele2,int sample,int* unique,int* middle);

int correct_counts_indels(struct VARIANT* variant,REFLIST* reflist,int allele1,int allele2);

int sort_allelecounts(struct VARIANT* variant, ALLELE* maxvec,int refbase);
int sort_allelecounts_diploid(struct VARIANT* variant, ALLELE* maxvec,int refbase);

int populate_contable(struct VARIANT* variant,int allele1,int allele2,int refbase,int allele3);


#endif
