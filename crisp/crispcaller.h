/***** CODE FOR CRISP pooled variant caller *****/

#ifndef INC_crispcaller_H
#define INC_crispcaller_H
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include "../variant.h"
#include "../bamsreader.h"
#include "../allelecounts.h"

extern double thresh1,MINQpv,thresh3,SER,RELAX;  //three pvalue thresholds 
extern int MAXITER; // max number of permutations
extern int MINPOS, MAXPOS, MIN_READS;
extern int USE_BASE_QVS; // default set to 1, if set to 0, use empirical error rates even for SNP calling
extern int OUTPUTGENOTYPES;
//extern int CHISQ_PERMUTATION;

int newCRISPcaller(REFLIST* reflist,int current,int position,READQUEUE* bq,struct BAMFILE_data* bd,struct VARIANT* variant,ALLELE* maxvec,FILE* vfile);

int populate_contable(struct VARIANT* variant,int allele1,int allele2,int refbase,int allele3);

int CRISPcaller(REFLIST* reflist,int current,int position,READQUEUE* bq,struct BAMFILE_data* bamfiles_data,struct VARIANT* variant,ALLELE* maxvec,FILE* vfile);

int call_genotypes(READQUEUE* bq,struct BAMFILE_data* bamfiles_data,struct VARIANT* variant,int allele1,int allele2);

// output variant calls and allele counts to VCF file (this is done once for each variant, even if it is multi-allelic) 
//int print_pooledvariant(struct VARIANT* variant,FILE* vfile);

int pooled_variantcaller(READQUEUE* bq,struct BAMFILE_data* bamfiles_data,struct VARIANT* variant, int allele1, int allele2);

void print_crispheader(struct OPTIONS* options);

void print_crispoptions();

void print_crispoptions_additional();


#endif
