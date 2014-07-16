#ifndef INC_variant_H
#define INC_variant_H
#include <stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include "bamsreader.h"
#include "common.h"

struct VCF_ALLELE
{
        char* bases; uint16_t length;
        char type;  // S/I/D/Complex
        uint32_t** allelecounts; // counts for this allele on +/-/b strand
        uint16_t HPlength;
        uint32_t score; double chi2pval;
        uint32_t AC; uint32_t VP; double AF; double delta; double strand_score; double HWE_pvalue;
        uint16_t varfilter; // EMpass, EMfail...
	char discard; 

};

struct VCF_RECORD // 96 bytes 
{
        uint32_t samples; uint32_t chrom; uint32_t position;
        uint32_t* readdepths; uint16_t varalleles;
        double QUAL;  // best variant score
        char* FILTER;  // filter string 
        char bidirectional; // 0/1 
        char* flankingseq;
        char lowcoverage; char IUPAC;
        uint32_t MQcounts[4];  uint16_t delalleles, insalleles; uint32_t total_depths[3];
	struct VCF_ALLELE* allelelist; // list of alleles, one for reference base
};


struct CRISPVAR // to store info about crisp variantallele
{
	int allele;
	double FETpvalue; double AF_ML; 
	double delta,deltaf,deltar;
	double deltaBFGS,E[6]; // BFGS likelihoods and error matrix
	int allelecount; int variantpools; // additional info for each variant 
	int* AC; double* QV; double* meanAC; double* varAC; // information about CRISP variant 

	// store info about error rates (both strands) and allele counts that can be used in VCF 
};

struct VARIANT // encapsulates all information about a single variant 
{
	struct OPTIONS* options; // added on aug 7 2012
	int samples; int bamfiles;
	char* chrom; int position; char refb; int refbase; 
	char** itb; // bases corresponding to each allele, A is 0, C is 1, G is 2, T is 3, then indel alleles 
	//HAPLOTYPE* haplotypes; int haps; // no of haplotypes
	char previousbase;
	int* counts; // total counts
	int* readdepths; // added may 25 2012
	int** indcounts; // individual allele counts, store total allele counts as well as for +strand,-strand,bidirectional,filtered
	int*** indcounts_binned; // individual allele counts binned based on strand, read1_2, quality values...
	double** Qhighcounts; 
	double** stats; double* tstats; double** qpvalue;
	int** ctable; int** ctable_stranded; double** ctable_weighted;
	int** newtable;
	int* strata; // 0/1/2/3 corresponding to heterogeneity in sample sequencing 
	char type[256]; int totaldepth; int variants;
	double paf[10]; // population allele frequency
	double ctpval[10]; double qvpvaluef[10]; double qvpvaluer[10]; double chi2pval[10]; double chi2stats[10][4];
	int nvariants[10]; int alleles[10]; int varpools[10];
	int varalleles; // how many variant alleles at this chrom apart from the reference base 
	double** GENLL; //genotype likelihoods for 10 possible SNP genotypes
	double** GENLLf; double** GENLLr; double** GENLLb; // forward & reverse strand genotype likelihoods
	//double** GENLLpooled; // pooled genotype likelihood for two alleles (0,1,2,......p+1) possible genotypes
	double** gpriors; // priors for genotype likelihoods based on allele frequencies 
	int* genotypes; // genotypes for diploid samples are of the form: allele1 allele2 
	// if only two alleles genotype = 0,1,2,...poolsize 
	// for triallelic: genotype = 0*poolsize + poolsize*allele2 + allele3 
	double** poolpv; // probability that the variant allele is present in a given pool, added dec 16 2011 
	double* genotypeQV; // difference between best posterior genotype and 2nd best genotype
	double* allelefreqs; // added dec 30, stores allele frequency for each allele calculated using Neilsen formula 
	// some information about homopolymer run length for each allele AAAAAA / TTC^n
	int HPlength[10];
	int strandpass,hetpass,u2pass,u3pass; 
	int MQcounts[4]; // may 7 2012, counts of reads with mapping quality of 0,<20,<30,all, for both reference and alternate allele 
	int filteredreads[8]; // due to high # of mismatches
	int insidetarget; // whether variant is inside an interval in bed file or not
	int* ploidy; // poolsize 
	// 0-> reference base  1:allele1  2:allele2 3: allele3 
	int allele1,allele2,allele3; // current alleles for variant evaluation 
	struct CRISPVAR crispvar[8];
	int varpoolsize; // variable is 1 if the pool size is variable and is read from bamlist file 06/14/2013
	double** NCarray; double** MCarray; 
	short IUPAC; // 0 if reference is A|C|T|G|N and 1 otherwise
};

struct INDEL { char* bases; int sample; char strand; int length; int reads; }; // l is length

struct GLL // genotype likelihood prior, likelihood using integration on forward and reverse strands
{
	double LLtotal, LLf, LLr, prior, ef,er; // error rates for reference allele 
	double LLsum; // sum over multiple genotypes
	double ILLe_f,ILLe_r;
};




void init_variant(struct VARIANT* variant,int bamfiles,int samples);

int initialize_variant(REFLIST* reflist,int current, int position,struct VARIANT* variant);

void calculate_allelecounts(REFLIST* reflist,int current,READQUEUE* bq,struct BAMFILE_data* bd, struct VARIANT* variant);

//void compute_PAF(struct VARIANT* variant);

////////////////////////////////////// INDEL ANALYSIS CODE ///////////////////////////////////

int indel_cmp(const void *a,const void *b);
int remove_ambiguous_bases(READQUEUE* bq,struct VARIANT* variant,int HPlength,int delta); // delta =+1 bases removed | delta = -1 bases added back
int calculate_indel_ambiguity(struct VARIANT* variant,REFLIST* reflist,int current,int allele);

void init_poolsizes(struct VARIANT* variant, struct OPTIONS* options,int VC_METHOD);
#endif
