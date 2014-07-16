#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include<string.h>
#include<ctype.h>
#include "chisquare.h"
#include "contables.h"
#include "../common.h"
int MAXN = 2001; // maxsize of NCRtable, this will be used for exact p-value computation for small k x 2 tables
int CHISQ_PERMUTATION = 0;

int ctable_cmp_full(const void *a,const void *b)
{
	const struct acounts *ia = (const struct acounts*)a;
	const struct acounts *ib = (const struct acounts*)b;
	if (ia->N == ib->N) 
	{
		if (ia->N0 == ib->N0) 
		{
			if (ia->r0 == ib->r0) 
			{
				if (ia->r1 == ib->r1) return (ia->r2 -ib->r2); else return (ia->r1 - ib->r1);
			}
			else return (ia->r0-ib->r0);
		}
		else return ia->N0-ib->N0;
	}
	else return (ia->N -ib->N);
}

int ctable_cmp(const void *a,const void *b)
{
	const struct acounts *ia = (const struct acounts*)a;
	const struct acounts *ib = (const struct acounts*)b;
	if (ia->N == ib->N) 
	{
		return ( (ia->r0+ia->r1+ia->r2) - (ib->r0+ib->r1+ib->r2)); // consistent between the two sorting methods
	}
	else return (ia->N -ib->N);
}

int ctable_cmp_descending(const void *a,const void *b)
{
	const struct acounts *ia = (const struct acounts*)a; const struct acounts *ib = (const struct acounts*)b;
	return (ib->N -ia->N);
}

int int_cmp(const void *a, const void *b) 
{ 
	const int *ia = (const int *)a; const int *ib = (const int *)b; 
	return (abs(*ia)-abs(*ib));  // change to absolute value 
}

int COUNT_cmp(const void *a, const void *b)
{
	const COUNT *ia = (const COUNT *)a; const COUNT *ib = (const COUNT *)b;
	if (ia->n == ib->n)
	{
		return ( (ia->r0+ia->r1) - (ib->r0+ib->r1)); // consistent between the two sorting methods
		/*
		if (ia->r0 == ib->r0) 
		{
			return (ia->r1 - ib->r1);
		}
		else return (ia->r0-ib->r0);
		*/
	}
	else return ( ia->n-ib->n );
}


double factorial(int n)
{
	double f=0;
	while (n > 0) f += log10(n--); 
	return f;
}

double ncr(int n,int r)
{
	if (r ==0 || n ==r) return 0;
	double ll =0; int i=0; 
	if (2*r > n) r = n-r;
	for (i=0;i<r;i++)  ll += log10((double)(n-i)/(i+1));
	return ll;
}

double multinomial(int n,int r1,int r2)
{
	double ll =0; int i=0; 
	for (i=0;i<r1;i++)  ll += log10((double)(n-i)/(i+1));
	for (i=0;i<r2;i++) ll += log10((double)(n-i-r1)/(i+1));
	return ll;
}

// samtools has FET test: https://github.com/lh3/samtools/blob/master/bcftools/fet.c 

double fetsided(int c1,int r1,int c2,int r2,double* PVleft,double* PVright) // c1 <= r1 and c2 <= r2 // function simplified on July 22 2011 
{
	double pvalue0 = -10000000; double  pvalue1 =-10000000; double pvalue = -10000000;
	double ptable = ncr(r1,c1) + ncr(r2,c2); double DEN= ncr(r1+r2,c1+c2);
	int minc =0, maxc = c1+c2;
	if (r2 < c1+c2) minc = c1+c2-r2;
	if (r1 < c1+c2) maxc = r1;
	int c11 = minc; int c21 = c1+c2-c11;
	double  p1 = ncr(r1,c11); double  p2 =  ncr(r2,c21); 

	while (c11 <= maxc)
	{
		if (c11 <= c1)
		{
			if (p1 + p2 <= pvalue0) pvalue0 += log10(1+pow(10,p1+p2-pvalue0));
			else  pvalue0 = p1+p2 + log10(1.0+pow(10,pvalue0-p1-p2));
		}

		if (c11 >= c1)
		{
			if (p1 + p2 <= pvalue1) pvalue1 += log10(1+pow(10,p1+p2-pvalue1));
			else pvalue1 = p1+p2 + log10(1.0+pow(10,pvalue1-p1-p2));
		}

		if (r1-c11 > 0) p1 += log10(r1-c11) - log10(c11+1);
		if (c21 > 0) p2 += log10(c21) - log10(r2-c21+1);
		c11 +=1; c21 -= 1;
	}
	pvalue0 -= DEN; pvalue1  -= DEN; *PVleft =pvalue0; *PVright =pvalue1;
	return pvalue; 
}

double fet(int c1,int r1,int c2,int r2) // c1 <= r1 and c2 <= r2 // function simplified on July 22 2011 
{
	int minc =0, maxc = c1+c2;
	if (r2 < c1+c2) minc = c1+c2-r2;
	if (r1 < c1+c2) maxc = r1;
	double ptable = ncr(r1,c1) + ncr(r2,c2); 	
	int c11 = minc; int c21 = c1+c2-minc;

	double  p1 = ncr(r1,c11); double  p2 =  ncr(r2,c21); 
	double pvalue0 = -10000000; double  pvalue1 =-10000000; double pvalue = -10000000;
	int flag =0;

	while (c11 <= maxc)
	{
		if (p1+p2 <= ptable + epsilon) 
		{
			if (p1 + p2 <= pvalue) pvalue += log10(1.0+pow(10,p1+p2-pvalue));
			else pvalue = p1 + p2 + log10(1.0+pow(10,pvalue-p1-p2)); 
			flag = 1; 
		}
		if (c11 <= c1)
		{
			if (p1 + p2 <= pvalue0) pvalue0 += log10(1+pow(10,p1+p2-pvalue0));
			else  pvalue0 = p1+p2 + log10(1.0+pow(10,pvalue0-p1-p2));
		}

		if (c11 >= c1)
		{
			if (p1 + p2 <= pvalue1) pvalue1 += log10(1+pow(10,p1+p2-pvalue1));
			else pvalue1 = p1+p2 + log10(1.0+pow(10,pvalue1-p1-p2));
		}

		if (r1-c11 > 0) p1 += log10(r1-c11) - log10(c11+1);
		if (c21 > 0) p2 += log10(c21) - log10(r2-c21+1);
		c11 +=1; c21 -= 1;
		//printf(" p1+p2 %f ptable %f %d %f \n",p1+p2,ptable-p1-p2,c11,sum);
	}
	if (flag ==0) pvalue = ptable; // pvalue was not updated so we set it to probability of observed table, sept 20 2012
	pvalue -= ncr(r1+r2,c1+c2);  return pvalue; 
	// we can also return pvalue0 pvalue1, one sided tests....
}


// cumulative probability of observing A or less reads in a sample of R reads with P(A) = e where e is from a beta distribution with mean = 1/poolsize and variance = ???
// when we pool DNA from multiple individuals, alpha = 
double minreads_pvalue_dirichlet(int R,int A,int poolsize,double cov) // standard_deviation/mean 
{
	// find parameters of beta distribution 
	double alpha = cov*cov*(poolsize-1) -1; alpha /= poolsize; double beta = (poolsize-1)*alpha; 
	alpha= 10; beta = 10*poolsize-10; 
	int r=0; 
	double betaI = ncr((int)alpha+(int)beta-2,(int)alpha-1) +log10(alpha+beta-1); 
	double pvlog = betaI - ncr(R+(int)alpha+(int)beta-2,r+(int)alpha-1) - log10(R+(int)alpha+(int)beta-1);
	double sum = pvlog; 
	for (r=1;r<=A;r++)
	{
		pvlog += log10(R-r+1) - log10(r) - (log10(R-r+(int)beta) - log10(r+(int)alpha-1)); 
		if (sum > pvlog) sum += log10(1.0+pow(10,pvlog-sum)); 
		else sum = pvlog + log10(1.0+pow(10,sum-pvlog)); 
	}
//	printf("sum %f %f\n",sum,sum);
	return sum;
}


// cumulative probability of observing A or less reads in a sample of R reads with P(A) = e 
double minreads_pvalue(int R,int A,double e)
{
	double e1 = log10(e); double e2 = log10(1.0-e); 
	double ll = R*e2; double pvlog = ll; double sum = ll;  int r=0;
	for (r=1;r<=A;r++)
	{ 
		pvlog += log10(R-r+1) - log10(r) + e1-e2;
		if (pvlog < sum) sum += log10(1.0+pow(10,pvlog-sum));
		else sum = pvlog + log10(1.0+pow(10,sum-pvlog));
	}
//	printf("sum %f \n",sum);
	return sum;
}

double binomial_test(int R, int A, double e,double* pv0,double* pv1)
{
	double e1 = log10(e); double e2 = log10(1.0-e); int r=0;
	double ll = ncr(R,A) + A*e1 + (R-A)*e2; double pvlog = ll; *pv1 = ll;
	for (r=A+1;r<R+1;r++)
	{
		pvlog += log10(R-r+1) - log10(r) + e1-e2;
		if (pvlog > *pv1) *pv1 = pvlog + log10(1.0+pow(10,*pv1-pvlog));
		else *pv1 += log10(1.0+pow(10,pvlog-*pv1));
	}
	pvlog = R*e2; *pv0 = pvlog; 
	for (r=1;r< A+1;r++)
	{
		pvlog += log10(R-r+1) - log10(r) + e1-e2; 
		if (pvlog > *pv0) *pv0 = pvlog + log10(1.0+pow(10,*pv0-pvlog));
		else *pv0 += log10(1.0+pow(10,pvlog-*pv0));
	}
	return 1;
}


double binomial_pvalue(int R, int A, double e)
{
	double e1 = log10(e); double e2 = log10(1.0-e); int r=0;
	double ll = ncr(R,A) + A*e1 + (R-A)*e2; double pvlog = ll; double sum = ll;
	for (r=A+1;r<R+1;r++)
	{
		pvlog += log10(R-r+1) - log10(r) + e1-e2;
		sum += log10(1.0 + pow(10,pvlog-sum));
	}
	return sum;
}

void compute_NCRtable(double** NCRtable)
{
	// NCRtable has to be already memory allocated 
	int i=0, j=0; double value = 0; 
	for (i=1;i<MAXN;i++)
	{
		for (j=0;j<MAXN;j++) NCRtable[i][j] = 0;
		NCRtable[i][0] = 0; value =0;
		for (j=1;j<i;j++)
		{
			value += log10(i-j+1) - log10(j); NCRtable[i][j] = value; 
		}
	}
}


///////////////////////////////// joint contingency table pvalue implemented july 15 2011, last modified dec 16 2011
int pvalue_contable_iter_stranded(int** ctable,double** ctable_weighted,int* strata, int size,int maxiter,int** newtable,double* permutationpvalue,double* chisqpvalue,double chi2stats[])
{
	int i=0,j=0,offset=0; 
	int  k=0,onesf=0,readsf=0,r=0,b=0; int onesr =0,readsr=0; int pvalue=0;
	double clra =0, clranew=0; double deltac=0; double DOF = size-1;
	//double chisqpvalue1=0;
	//fprintf(stdout,"\n");
	//for (i=0;i<size;i++)
	//{
		//fprintf(stdout,"counts %d %d %d %d \t",ctable[i][0],ctable[i][1],ctable[i][2],ctable[i][3]);
		//fprintf(stdout,"countsW %f %f %f %f \n",ctable_weighted[i][0],ctable_weighted[i][1],ctable_weighted[i][2],ctable_weighted[i][3]);
		//for (j=0;j<4;j++) ctable_weighted[i][j] = ctable[i][j]; 
	//}
	int computechi2; //double chistat[4];
	
	computechi2 = chi2pvalue(ctable_weighted,strata,size,chisqpvalue,chi2stats); 
	//if (PFLAG >=2) fprintf(stdout,"chi:%.1f:%.2f",chistat[0],*chisqpvalue);
	if (PFLAG >=2) 
	{
		fprintf(stdout,"chi:%.1f,%.1f,%.1f,%.1f:%.2f",chi2stats[0],chi2stats[1],chi2stats[2],chi2stats[3],*chisqpvalue);
                deltac =chi2stats[3]-DOF - (chi2stats[1]-DOF) - (chi2stats[2]-DOF);

		if ((chi2stats[3]+3 < chi2stats[1]) || (chi2stats[3]+3 < chi2stats[2])) fprintf(stdout,",SB-CHI2");
	}

	if ( (*chisqpvalue <= 0) && computechi2 ==1 && CHISQ_PERMUTATION ==0) 
	//if ( (*chisqpvalue <= -7) && computechi2 ==1 && CHISQ_PERMUTATION ==0) 
	{
		// if statistic is really high, then we skip the part below and speed up everything...
		*permutationpvalue = 1; // no permutations done
		if (PFLAG >=1) fprintf(stdout,"\t");
		return 0;
	}
//	else { fprintf(stdout,"doing_permutations "); }

	for (i=0;i<size;i++) {  onesf += ctable[i][1]; readsf += ctable[i][0]; onesr += ctable[i][3]; readsr += ctable[i][2]; }
	for (i=0;i<size;i++) { newtable[i][0]= ctable[i][0]; newtable[i][2] = ctable[i][2]; newtable[i][1] = 0; newtable[i][3] =0; }
	for (i=0;i<size;i++) clra += ncr(ctable[i][0]+ctable[i][2],ctable[i][1]+ctable[i][3]);
	//for (i=0;i<size;i++) clra += ncr(ctable[i][0],ctable[i][1]) + ncr(ctable[i][2],ctable[i][3]);

	int* table_index;
//	fprintf(stdout,"reads %d %d \n",readsf,readsr); for (i=0;i<size;i++) fprintf(stdout,"%d %d %d %d %d \n",i,ctable[i][0],ctable[i][1],ctable[i][2],ctable[i][3]);
	table_index = (int*)malloc(sizeof(int)*(readsf+readsr+1)); 
	for (i=0;i<size;i++) {  for (j=0;j<ctable[i][0];j++) table_index[offset+j] = i;  offset += ctable[i][0]; }
	for (i=0;i<size;i++) { 	for (j=0;j<ctable[i][2];j++) table_index[offset+j] = i;  offset += ctable[i][2]; }

	//double clrasum =0; double clrasumsq =0;
	for (k=1;k<=maxiter;k++)
	{
		//continue; fprintf(stdout,"test \n");
		clranew = 0; for (i=0;i<size;i++) { newtable[i][1] = 0; newtable[i][3] = 0; } 
		for (i=0;i<onesf;i++)
		{
			r = (int)(drand48()*(readsf-i))+i; 
			b = table_index[(int)r]; table_index[(int)r] =  table_index[i]; table_index[i] =b; 
			clranew += log10( (double)(newtable[b][0]+newtable[b][2]-newtable[b][1])/(newtable[b][1]+1)); newtable[b][1]++; 
		}
		for (i=0;i<onesr;i++)
		{
			r = (int)(drand48()*(readsr-i))+i + readsf; 
			b = table_index[(int)r]; table_index[(int)r] =  table_index[i+readsf]; table_index[i+readsf] =b; 
			clranew += log10( (double)(newtable[b][2]+newtable[b][0]-newtable[b][3]-newtable[b][1])/(newtable[b][1]+newtable[b][3]+1)); newtable[b][3]++;
		}
		//printf("table %f %f %d\n",clranew,clra,k);
		if (clranew <= clra+epsilon) 	pvalue +=1; 
		//clrasum += clranew; clrasumsq += clranew*clranew;
		if (pvalue >= 10) break;
		if (k <= 1000 && pvalue >= 5) break;
		//if (k >= 100 && pvalue >= 5) break;
		//fprintf(stdout," %f ",clranew);
	}
	//clrasum /= (k+1); clrasumsq /= (k+1); clrasumsq -= clrasum*clrasum; clrasumsq = sqrt(clrasumsq);
	//fprintf(stdout,"clra %f mean %f sq %f %d\n",clra,clrasum,clrasumsq,k+1);
	if (PFLAG >=2) fprintf(stdout,":%d/%d\t",pvalue,k);
	free(table_index);
	*permutationpvalue = log10(pvalue+1)-log10(k+1); 
	return 1;
	//	return (double)(pvalue+1)/(k+2); // apply Yates correction
}


void chernoff_bound(int* counts,double* stats,int refbase,int altbase,double SER, double* P0,double* P1,int maxalleles)
{
	double mu0,mu1,delta0,delta1,p0,p1;
	if (SER == -1)
	{
		mu0 = stats[altbase]+stats[refbase]; mu1 = stats[altbase+maxalleles]+stats[refbase+maxalleles];
		if (counts[altbase] ==0 || mu0 ==0) delta0 = 1;
		else delta0 = counts[altbase]/mu0-1;
		if (counts[altbase+maxalleles] ==0 || mu1 ==0) delta1 = 1;
		else delta1 = counts[altbase+maxalleles]/mu1-1;

		if (delta0 <= 0) p0 = 0;
		else p0 = mu0*(delta0*log10(2.718281)-(1+delta0)*log10(1+delta0));
		if (delta1 <= 0) p1 = 0;
		else    p1 = mu1*(delta1*log10(2.718281)-(1+delta1)*log10(1+delta1));
	}
	else
	{
		p0 = binomial_pvalue(counts[refbase]+counts[altbase],counts[altbase],SER);
		p1 = binomial_pvalue(counts[refbase+maxalleles]+counts[altbase+maxalleles],counts[altbase+maxalleles],SER);
	}
	*P0 = p0; *P1 = p1;
}

