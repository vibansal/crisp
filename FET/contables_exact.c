#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include<string.h>
#include<ctype.h>

/* exact and permutation tests for calculating p-value of contingency table, not used anymore 06/20/2013 */

double compute_pvalue_permutation(int** ctable,int size,double cl,int iter,int** newtable,int fast,struct acounts* atable,int st)
{
	double pvalue = 0, clra =0;
	int  i=0,j=0,k=0,offset=0,ones=0,reads=0,r=0,b=0;
	//	int** newtable = (int**)malloc(size*sizeof(int*));
	//	for (i=0;i<size;i++) { newtable[i] = (int*)malloc(2*sizeof(int)); newtable[i][0]= ctable[i][0]; newtable[i][1] = 0; }
	for (i=0;i<size;i++) { newtable[i][0]= ctable[i][0+st]; newtable[i][1] = 0; }

	for (i=0;i<size;i++) {  ones += ctable[i][1+st]; reads += ctable[i][0+st]; }
	int* table_index = (int*)malloc(sizeof(int)*reads); 
	for (i=0;i<size;i++)
	{
		for (j=0;j<ctable[i][0+st];j++) table_index[offset+j] = i; 
		offset += ctable[i][0+st];
	}

	//	struct acounts* atable = (struct acounts*)malloc(sizeof(struct acounts)*size);
	for(i=0;i<size;i++) { atable[i].N = ctable[i][0+st]; atable[i].R = ctable[i][1+st]; }
	//	for (i=0;i<size-1;i++) { if (atable[i].n != atable[i+1].n) { printf("%d %d | ",atable[i].n,icells); icells = 1; } else icells++; } printf("\n");
	double factor =0, ifactor =0; // factn = factorial(size);
	//ifactor = factn + correctionfactor(atable,size); 
	for (k=1;k<iter+1;k++)
	{
		for (i=0;i<size;i++) newtable[i][1] = 0; 
		clra =0;
		for (i=0;i<ones;i++)
		{
			r = (int)(drand48()*(reads-i))+i; 
			b = table_index[(int)r]; 	table_index[(int)r] =  table_index[i]; table_index[i] =b; 
			clra += log10(newtable[b][0]-newtable[b][1]) - log10(newtable[b][1]+1); newtable[b][1] +=1; 
		}

		for(i=0;i<size;i++) { atable[i].N = newtable[i][0]; atable[i].R = newtable[i][1]; }
		//factor = factn + correctionfactor(atable,size);
		//  	printf("table %f %f %d\n",cl,clra,k);
		if (clra + factor <= cl+ifactor+epsilon) 	pvalue +=1; 
		//	if (clra  <= cl) 	pvalue +=1; 
		if (fast ==1)
		{
			if (k <= 20 && pvalue >= 2) break;
			if (k <= 100 && pvalue >= 3) break;
			//if (k <= 1000 && pvalue >= 3) break;
		}
	}
	free(table_index); //free(atable);
	//	for (i=0;i<size;i++) free(newtable[i]); free(newtable);
	return (double)(pvalue+1)/(k+1); // apply Yates correction
	//else  return (double)(pvalue)/k;
}

// recursive algorithm for calculating the p-value of a contingency table 2 x k
// size is number of columns, R is row sum for two rows, A is row sum, prob is probability of partial table, st is index into ctable (for strandednedness) 
double compute_pvalue_exact(int** ctable,int size,int offset,int R,int A,double prob,double** NCRtable,int st)
{
	int i=0; double a=0, cp=0,sum=0;

	if (A > R) return -1;
	if (size ==1)
	{
		//printf("size = 1 A> R %d %d \n",A,R);
		if (ctable[offset][0+st] > 0) cp = NCRtable[ctable[offset][0+st]][A];
		else cp = 0; 
		if (cp > prob)  return -1;
		else return cp;
	}
	else
	{
		sum =0;
		for (i=0;i<=ctable[offset][0+st] && i <= A;i++)  
		{
			if (ctable[offset][0+st] > 0) cp = NCRtable[ctable[offset][0+st]][i];
			else cp = 0; 
			a = compute_pvalue_exact(ctable,size-1,offset+1,R-ctable[offset][0+st],A-i,prob-cp,NCRtable,st);
			if (a > -1) sum += pow(10,a + cp-prob);
		}

		if (sum > 0) return log10(sum) + prob;
		return -1; 
	}
}

// table with lowest probability can be constructed by sorting columns by column sums and packing the 
double pvalue_contable_iter(int** ctable,int size,int iter,double** NCRtable,int** newtable,int fast, struct acounts* atable,int st)
{
	int ones = 0; int reads = 0; int  i=0; double p=0;
	double cl =0;
	for (i=0;i<size;i++) cl += ncr(ctable[i][0+st],ctable[i][1+st]);
	int maxcol = 0;
	for (i=0;i<size;i++) 
	{  
		ones += ctable[i][1+st]; reads += ctable[i][0+st]; 
		if (ctable[i][0+st] > maxcol) maxcol = ctable[i][0+st];
	}

	if (ones < 70 && maxcol < MAXN && size >= 2 && size < 8)
	{
		//	printf("computing p-value using exact method %d %d \n",reads,ones);
		p = compute_pvalue_exact(ctable,size,0,reads,ones,cl,NCRtable,st);
		//printf("%f %d %d %f %f \n",cl,reads,ones,p,ncr(reads,ones));  
		return pow(10,p-ncr(reads,ones));
	}
	else
	{
		//	printf("computing p-value using permutation method %d %d \n",reads,ones);
		p= compute_pvalue_permutation(ctable,size,cl,iter,newtable,fast,atable,st);
		return p;
	}
}

