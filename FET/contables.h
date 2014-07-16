
#ifndef INC_contables_H
#define INC_contables_H
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#define epsilon 1e-6 // comparison of two doubleing point numbers is not working correctly for some reason Dec 15 2010
extern int MAXN;
extern int CHISQ_PERMUTATION;

struct acounts { int N,R,N0,N1,N2,R0,R1,R2; int r,r0,r1,r2; int sumL,sumR;  }; // r0 and r1 correspond to random allele counts 
// also stores partial sums sumL (sum from 1 to this column) and sumR (sum from this column to end)

typedef struct COUNT
{ 
	int n,n0,n1,n2,r0,r1,r2,r; // strand information  
} COUNT; 

int ctable_cmp(const void *a,const void *b);
int ctable_cmp_full(const void *a,const void *b);
int ctable_cmp_descending(const void *a,const void *b);

int int_cmp(const void *a, const void *b);

int COUNT_cmp(const void *a, const void *b);

double factorial(int n);

double ncr(int n,int r);

double multinomial(int n,int r1,int r2);

double fetsided(int c1,int r1,int c2,int r2,double* pv1,double* pv2);

double fet(int c1,int r1,int c2,int r2);// c1 <= r1 and c2 <= r2 // function simplified on July 22 2011 

// cumulative probability of observing A or less reads in a sample of R reads with P(A) = e 
double minreads_pvalue(int R,int A,double e);
double minreads_pvalue_dirichlet(int R,int A,int poolsize,double cov); // standard_deviation/mean 

double binomial_test(int R, int A, double e,double* pv0,double* pv1);

double binomial_pvalue(int R, int A, double e);

void compute_NCRtable(double** NCRtable);

double correctionfactor(struct acounts* ctable,int size);

double compute_pvalue_permutation(int** ctable,int size,double cl,int iter,int** newtable,int fast,struct acounts* atable,int st);

// recursive algorithm for calculating the p-value of a contingency table 2 x k
// size is number of columns, R is row sum for two rows, A is row sum, prob is probability of partial table, st is index into ctable (for strandednedness) 
double compute_pvalue_exact(int** ctable,int size,int offset,int R,int A,double prob,double** NCRtable,int st);

double pvalue_contable_iter(int** ctable,int size,int iter,double** NCRtable,int** newtable,int fast, struct acounts* atable,int st);

///////////////////////////////// joint contingency table ppvalue implemented july 15 2011, last modified dec 16 2011

int pvalue_contable_iter_stranded(int** ctable,double** ctable_weighted,int strata[],int size,int iter,int** newtable,double* permutationpvalue,double* chisqpvalue,double chi2stats[]);

void chernoff_bound(int* counts,double* stats,int refbase,int altbase,double SER, double* P0,double* P1,int maxalleles);


// specific method for low coverage variant p-value calculation, dec 20 2011
//int pvalue_lowcoverage(int** ctable,int size,int iter,struct acounts* atable,double* permutationpvalue);

//int pvalue_lowcoverage_slow(int** ctable,int size,int iter,struct acounts* atable,double* permutationpvalue);

#endif
