#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include<string.h>

/* the code below considers all tables that have identical bin counts as equal, i.e. (2,0) (2,1) (2,0) (3,1)  ===  (2,1) (2,0) (2,0) (3,1) == (2,0) (2,0) (2,1) (3,1)  | 3 tables equivalent 
therefore when we generate a random table, we have to determine the number of possible tables that are in same class by permutation 
ifactor = multiplicative factor that accounts for this... 
 why do we need emptyfactor | if bin is empty -> C(n,0) = 1 | 06/20/2013
emptyfactor = (k1)! (k2)! where there are k1 bins with 2 reads, k2 bins with 3 reads in total.... 

correction factor for bins of size 'n'  =  k1!/ k1_0! k1_1! k1_2! ... k1_n!
where there are k1 bins of size 'n', k1_0 = bins with no alt.reads, k1_1 = bins with 1 alt.reads 
k1! is same for all permutations... 
k1_1! k1_2! how do we get k1_0! since in nebins this are not accounted for

correction factor should be multiplied by table probability...
*/

// methods updated to handle bidirectional reads 06/28/2013
// specific method for low coverage variant p-value calculation, dec 20 2011
int pvalue_lowcoverage(int** ctable,int size,int iter,struct acounts* atable,double* permutationpvalue)
{
	int i=0,j=0,offset=0; int den=0; int num=0;
	int  k=0,onesf=0,readsf=0,r=0,b=0; int onesr =0,readsr=0,onesb=0,readsb=0;
	double clra =0, clranew=0; unsigned int pvalue =0, pvalue_1 = 0;
	double factor =0, ifactor =0; double clrasum =0; double clrasumsq =0;
	int maxcoverage =0;
	int nebins =0;
	COUNT* NEBcounts; int* ballindex; int* table_index; int* bincounts;

	for (i=0;i<size;i++) 
	{ 
		atable[i].r0 = ctable[i][1]; atable[i].r1 = ctable[i][3]; atable[i].r2 = ctable[i][5];
		atable[i].N0 = ctable[i][0]; atable[i].N1 = ctable[i][2]; atable[i].N2 = ctable[i][4];
		atable[i].N = ctable[i][0]+ctable[i][2]+ctable[i][4]; atable[i].R = ctable[i][1]+ctable[i][3] + ctable[i][5];
		onesf += ctable[i][1]; readsf += ctable[i][0]; onesr += ctable[i][3]; readsr += ctable[i][2]; onesb += ctable[i][5]; readsb += ctable[i][4]; 
		if (atable[i].R > 0) clra += ncr(atable[i].N,atable[i].R);
		if (atable[i].N > maxcoverage) maxcoverage= atable[i].N;
	}
	if (onesf + onesr + onesb < 2) { *permutationpvalue = 0; return 0; } 

	//fprintf(stdout,"maxcov %d \n",maxcoverage);
	// bincounts is number of bins with given number of reads...
	bincounts = (int*)malloc(sizeof(int)*(maxcoverage+1));	 // important to add plus 1 since 'MAXCOV' is also included, atable[i].N can be equal to maxcoverage
	for (i=0;i<maxcoverage+1;i++) bincounts[i] = 0; for (i=0;i<size;i++) bincounts[atable[i].N]++; 

	/*
	qsort(atable,size,sizeof(struct acounts),ctable_cmp); i=1; den =2;
	while (i < size)
	{
		if (atable[i-1].N != atable[i].N) ifactor += factorial(bincounts[atable[i-1].N]); 
		//fprintf(stdout,"%d:%d ",atable[i-1].N,atable[i-1].r0+atable[i-1].r1);
		if (atable[i-1].N == atable[i].N && (atable[i-1].r) == (atable[i].r) ) ifactor -=  log10(den++); 
		else den =2;
		i++;
	}
	ifactor += factorial(bincounts[atable[size-1].N]);   // if last bin is same as previous bin then we need this boundary factor, else this should be 1 
	*/

	table_index = (int*)malloc(sizeof(int)*(readsf+readsr+readsb+1)); 
	
	offset = 0;
	for (i=0;i<size;i++) {  for (j=0;j<atable[i].N0;j++) table_index[offset+j] = i;  offset += atable[i].N0; }
	for (i=0;i<size;i++) { 	for (j=0;j<atable[i].N1;j++) table_index[offset+j] = i;  offset += atable[i].N1; }
	for (i=0;i<size;i++) { 	for (j=0;j<atable[i].N2;j++) table_index[offset+j] = i;  offset += atable[i].N2; }
	ballindex = (int*)malloc(sizeof(int)*(onesf+onesr+onesb)); // bin where the alternate allele (ball) is placed 
	NEBcounts = (COUNT*)malloc(sizeof(COUNT)*(onesf+onesr+onesb)); // non-empty bin counts  

	// need to reduce the loop complexity to less than O(n) -> O(ones) 
	// code updated dec 27 2011, check for end conditions to make sure it is working propoerly
	for (k=1;k<=iter;k++)
	{
		for (i=0;i<onesf;i++)
		{
			r = (int)(drand48()*(readsf-i))+i; 
			b = table_index[r]; table_index[r] =  table_index[i]; table_index[i] =b; ballindex[i] = b; 
		}
		for (i=0;i<onesr;i++)
		{
			r = (int)(drand48()*(readsr-i))+i + readsf; 
			b = table_index[r]; table_index[r] =  table_index[i+readsf]; table_index[i+readsf] =b; ballindex[i+onesf] = b; 
		}
		for (i=0;i<onesb;i++)
		{
			r = (int)(drand48()*(readsb-i))+i + readsf+readsr; 
			b = table_index[r]; table_index[r] =  table_index[i+readsf+readsr]; table_index[i+readsf+readsr] =b; ballindex[i+onesf+onesr] = b; 
		}

		qsort(ballindex,onesf+onesr+onesb,sizeof(int),int_cmp);  // sort the list of bins that in which '1's were assigned to generate a set of unique counts
		nebins =0;
		for (i=0;i<onesf+onesr+onesb;i++)
		{
			if (i > 0 && ballindex[i] == ballindex[i-1]) NEBcounts[nebins-1].r++; 
			else 
			{ 
				NEBcounts[nebins].r = 1; NEBcounts[nebins].n = atable[ballindex[i]].N;
				nebins++;
			}
		} 
		clranew=0; for (i=0;i<nebins;i++) clranew += ncr(NEBcounts[i].n,NEBcounts[i].r); // basic table probability 
		
		/*
		if (nebins > 2) qsort(NEBcounts,nebins,sizeof(COUNT),COUNT_cmp); 
		i=1; den = 2; factor=0; num= bincounts[NEBcounts[0].n];
		while (i < nebins)
		{
			//fprintf(stdout,"%d:%d:%d,%d |",bincounts[NEBcounts[i].n],NEBcounts[i].n,NEBcounts[i].r0,NEBcounts[i].r1);
			if (NEBcounts[i].n == NEBcounts[i-1].n  && (NEBcounts[i].r) ==(NEBcounts[i-1].r)) factor -= log10(den++);
			else den = 2;
			if (NEBcounts[i-1].n == NEBcounts[i].n) num--;
			else 
			{
				factor += factorial(bincounts[NEBcounts[i-1].n]) - factorial(num-1);  // 06/20/13 change
				num = bincounts[NEBcounts[i].n];
			}
			i++;
		}
		factor += factorial(bincounts[NEBcounts[i-1].n]) - factorial(num-1);  // 06/20/13 change
		if (clranew + factor <= clra+ifactor+epsilon) 
		{
			fprintf(stdout," exceed %f:%f %f:%f ",clranew,factor,clra,ifactor);  fprintf(stdout,"nebins %d %f delta %f\n",nebins,factor,clra+ifactor-clranew-factor);
		}
		else
		{
			fprintf(stdout," ..... %f:%f %f:%f ",clranew,factor,clra,ifactor);  fprintf(stdout,"nebins %d %f\n",nebins,factor);
			
		}*/
		factor =0; ifactor = 0;

		if (clranew + factor <= clra+ifactor) 	pvalue +=1; 
		if (clranew + factor >= clra+ifactor) 	pvalue_1 +=1; 
		clrasum += clranew+factor; clrasumsq += (clranew+factor)*(clranew+factor);
		if ((pvalue >= 10 && pvalue_1 >=10) || (k <= 100 && pvalue >= 5 && pvalue_1 >=5)) break; // increased these values on may 4 2012
		//fprintf(stdout,"%f ",clranew+factor);
	}
	clrasum /= (k); clrasumsq /= (k); clrasumsq -= clrasum*clrasum; clrasumsq = sqrt(clrasumsq);
	if (pvalue_1 < pvalue) 
	{
		*permutationpvalue = log10((double)pvalue_1+1)-log10(k+1);
		if (PFLAG >=1 && *permutationpvalue <= -3) fprintf(stdout,"PV1 ");
	}
	else *permutationpvalue = log10((double)pvalue+1)-log10(k+1);
	if (PFLAG >=1) fprintf(stdout,"LL %.2f %.2f/%.2f %d/%d %d\t",clra+ifactor,clrasum,clrasumsq,pvalue,k,pvalue_1);
	//ultra important: free memory in same order as it was allocated
	free(bincounts); free(table_index); free(ballindex); free(NEBcounts);
	return 1;
}

// extend it to account for multiple populations or sub-groups .... ctable is actually a list of tables...
// non-weighted test without correction factor should be done and if it fails -> do more sensitive test with correction factor...
// specific method for low coverage variant p-value calculation, dec 20 2011
// this method is exact | does not use optimizations 
int pvalue_lowcoverage_stranded(int** ctable,int size,int iter,struct acounts* atable,double* permutationpvalue)
{
	int i=0,j=0,offset=0; int den=0;
	int  k=0,onesf=0,readsf=0,r=0,b=0; int onesr =0,readsr=0,onesb=0,readsb=0;
	double clra =0, pvalue = 0, clranew=0; double CT = 0;
	double factor =0, ifactor =0; double clrasum =0; double clrasumsq =0;
	double sstat =0; 

	for (i=0;i<size;i++) 
	{ 
		// ctable[i][1] += ctable[i][3]; ctable[i][3] = 0; ctable[i][0] += ctable[i][2]; ctable[i][2] = 0; 
		atable[i].r0 = ctable[i][1]; atable[i].r1 = ctable[i][3]; atable[i].r2 = ctable[i][5];
		atable[i].N0 = ctable[i][0]; atable[i].N1 = ctable[i][2]; atable[i].N2 = ctable[i][4];
		atable[i].N = ctable[i][0]+ctable[i][2]+ctable[i][4]; atable[i].R = ctable[i][1]+ctable[i][3] + ctable[i][5];
		onesf += ctable[i][1]; readsf += ctable[i][0]; onesr += ctable[i][3]; readsr += ctable[i][2]; onesb += ctable[i][5]; readsb += ctable[i][4]; 
		CT = ncr(ctable[i][0]+ctable[i][2]+ctable[i][4],ctable[i][1]+ctable[i][3]+ctable[i][5]);	clra += CT;
	}
	if (onesf + onesr + onesb < 2) { *permutationpvalue = 0; return 0; } 
	qsort(atable,size,sizeof(struct acounts),ctable_cmp_full); i=1; den =2;
	while (i++ < size)
	{
		//fprintf(stdout,"%d:%d ",atable[i-1].N,atable[i-1].r0+atable[i-1].r1);
		if (atable[i-1].N == atable[i].N && atable[i-1].N0 == atable[i].N0 && atable[i].r0 == atable[i-1].r0 && atable[i].r1 ==atable[i-1].r1 && atable[i].r2 == atable[i-1].r2) ifactor -= log10(den++); 
		else den =2;
	}

	int* table_index = (int*)malloc(sizeof(int)*(readsf+readsr+readsb)); 
	for (i=0;i<size;i++) {  for (j=0;j<ctable[i][0];j++) table_index[offset+j] = i;  offset += ctable[i][0]; }
	for (i=0;i<size;i++) { 	for (j=0;j<ctable[i][2];j++) table_index[offset+j] = i;  offset += ctable[i][2]; }
	for (i=0;i<size;i++) { 	for (j=0;j<ctable[i][4];j++) table_index[offset+j] = i;  offset += ctable[i][4]; }

	// permutation takes O(nlog(n) + ones) where 'n' is # of samples 
	for (k=1;k<=iter;k++)
	{
		clranew = 0; 
		for (i=0;i<size;i++) 
		{ 
			atable[i].r0 = 0; atable[i].r1 = 0; atable[i].r2 =0; atable[i].N = ctable[i][0]+ctable[i][2]+ctable[i][4]; 
			atable[i].N0 = ctable[i][0]; 	atable[i].R = ctable[i][1]+ctable[i][3]+ctable[i][5]; 
		}
		for (i=0;i<onesf;i++)
		{
			r = (int)(drand48()*(readsf-i))+i; 
			b = table_index[(int)r]; table_index[(int)r] =  table_index[i]; table_index[i] =b; 
			clranew += log10(atable[b].N-atable[b].r0) - log10(atable[b].r0+1); 
			atable[b].r0++;
		}
		for (i=0;i<onesr;i++)
		{
			r = (int)(drand48()*(readsr-i))+i + readsf; 
			b = table_index[(int)r]; table_index[(int)r] =  table_index[i+readsf]; table_index[i+readsf] =b; 
			clranew += log10(atable[b].N-atable[b].r1-atable[b].r0) - log10(atable[b].r0+atable[b].r1+1); 
			atable[b].r1++;
		}
		for (i=0;i<onesb;i++)
		{
			r = (int)(drand48()*(readsb-i))+i + readsf+readsr; 
			b = table_index[(int)r]; table_index[(int)r] =  table_index[i+readsf+readsr]; table_index[i+readsf+readsr] =b; 
			clranew += log10(atable[b].N-atable[b].r1-atable[b].r0-atable[b].r2) - log10(atable[b].r0+atable[b].r1+atable[b].r2+1); 
			atable[b].r2++;
		}

		qsort(atable,size,sizeof(struct acounts),ctable_cmp_full); i=1; den =2; factor =0;
		while (i++ < size)
		{
			if (atable[i-1].N == atable[i].N && atable[i-1].N0 == atable[i].N0 && atable[i].r0 == atable[i-1].r0 && atable[i].r1 ==atable[i-1].r1 && atable[i-1].r2 == atable[i].r2) factor -= log10(den++); 
			else den =2;
		}
		if (clranew + factor <= clra+ifactor+epsilon) 	pvalue +=1; 
		clrasum += clranew+factor; clrasumsq += (clranew+factor)*(clranew+factor);
		if (pvalue >= 10 || (k <= 100 && pvalue >= 5)) break;
	}
	clrasum /= (k); clrasumsq /= (k); clrasumsq -= clrasum*clrasum; clrasumsq = sqrt(clrasumsq);
	*permutationpvalue = log10(pvalue+1)-log10(k+1);
	fprintf(stdout,"LL %.2f %.2f/%.2f %d/%d\t",clra+ifactor,clrasum,clrasumsq,(int)pvalue,k);
	//	fprintf(stdout,"\t");
	free(table_index); //free(atable);
	return 1;

	//	return (double)(pvalue+1)/(k+2); // apply Yates correction
}

// this function is not correct, does not account for bidirectional read counts
int pvalue_lowcoverage_stranded1(int** ctable,int size,int iter,struct acounts* atable,double* permutationpvalue)
{
	int i=0,j=0,offset=0; int den=0;
	int  k=0,onesf=0,readsf=0,r=0,b=0; int onesr =0,readsr=0;
	double clra =0, pvalue = 0, clranew=0; double CT = 0;
	double factor =0, ifactor =0; double clrasum =0; double clrasumsq =0;
	double sstat[5] ={0,0,0,0,0}; double s=0;

	for (i=0;i<size;i++) 
	{ 
		// ctable[i][1] += ctable[i][3]; ctable[i][3] = 0; ctable[i][0] += ctable[i][2]; ctable[i][2] = 0; 
		
		onesf += ctable[i][1]; readsf += ctable[i][0]; onesr += ctable[i][3]; readsr += ctable[i][2]; 
		atable[i].r0 = ctable[i][1]; atable[i].r1 = ctable[i][3]; atable[i].r2 = 0;
		atable[i].N0 = ctable[i][0]; atable[i].N1 = ctable[i][2]; atable[i].N2 = 0;
		atable[i].N = ctable[i][0]+ctable[i][2]; atable[i].R = ctable[i][1]+ctable[i][3]; 
		CT = ncr(ctable[i][0]+ctable[i][2],ctable[i][1]+ctable[i][3]);	
		clra += CT;
	}
	for (i=0;i<size;i++)
	{
		if (atable[i].N0 == 0 || atable[i].N1 == 0) continue;
		s = (double)atable[i].r0/(double)atable[i].N0 -  (double)atable[i].r1/atable[i].N1; s *= s;
		s /= 1.0/atable[i].N0 + 1.0/atable[i].N1; 
		sstat[0] +=s; 
	}

	qsort(atable,size,sizeof(struct acounts),ctable_cmp_full); i=1; den =2;
	while (i++ < size)
	{
		//fprintf(stdout,"%d:%d ",atable[i-1].N,atable[i-1].r0+atable[i-1].r1);
		//if (atable[i-1].N == atable[i].N && (atable[i].r0+atable[i].r1 == atable[i-1].r0+atable[i-1].r1)) ifactor -= log10(den++); 
		if (atable[i-1].N == atable[i].N && atable[i-1].N0 == atable[i].N0 && atable[i].r0 == atable[i-1].r0 && atable[i].r1 ==atable[i-1].r1) ifactor -= log10(den++); 
		else den =2;
	}

	int* table_index = (int*)malloc(sizeof(int)*(readsf+readsr)); 
	for (i=0;i<size;i++) {  for (j=0;j<ctable[i][0];j++) table_index[offset+j] = i;  offset += ctable[i][0]; }
	for (i=0;i<size;i++) { 	for (j=0;j<ctable[i][2];j++) table_index[offset+j] = i;  offset += ctable[i][2]; }

	for (k=1;k<=iter;k++)
	{
		clranew = 0; 
		for (i=0;i<size;i++) 
		{ 
			atable[i].r0 = 0; atable[i].r1 = 0; atable[i].N = ctable[i][0]+ctable[i][2]; atable[i].N0 = ctable[i][0];
			atable[i].R = ctable[i][1]+ctable[i][3]; 
		}
		for (i=0;i<onesf;i++)
		{
			r = (int)(drand48()*(readsf-i))+i; 
			b = table_index[(int)r]; table_index[(int)r] =  table_index[i]; table_index[i] =b; 
			clranew += log10(atable[b].N-atable[b].r0) - log10(atable[b].r0+1); 
			atable[b].r0++;
		}
		for (i=0;i<onesr;i++)
		{
			r = (int)(drand48()*(readsr-i))+i + readsf; 
			b = table_index[(int)r]; table_index[(int)r] =  table_index[i+readsf]; table_index[i+readsf] =b; 
			clranew += log10(atable[b].N-atable[b].r1-atable[b].r0) - log10(atable[b].r0+atable[b].r1+1); 
			atable[b].r1++;
		}

		// calculate_sstat
		sstat[1] = 0; sstat[2] = 0;
		for (i=0;i<size;i++)
		{
			if (atable[i].N0 == 0 || atable[i].N1 == 0) continue;
			s = (double)atable[i].r0/(double)atable[i].N0 -  (double)atable[i].r1/atable[i].N1; s *= s; s /= 1.0/atable[i].N0 + 1.0/atable[i].N1; 
			//s = (double)atable[i].r0/(double)atable[i].N0 -  (double)ctable[i][3]/atable[i].N1; s *= s; s /= 1.0/atable[i].N0 + 1.0/atable[i].N1; 
			sstat[1] +=s; 
			s = (double)ctable[i][1]/(double)atable[i].N0 -  (double)atable[i].r1/atable[i].N1; s *= s; s /= 1.0/atable[i].N0 + 1.0/atable[i].N1; 
			sstat[2] +=s; 
		}
		if (sstat[1] < sstat[0]) sstat[3] +=1;	if (sstat[2] < sstat[0]) sstat[4] +=1;
		//fprintf(stdout,"%f %f %f \n",sstat[0],sstat[1],sstat[2]);

		qsort(atable,size,sizeof(struct acounts),ctable_cmp_full); i=1; den =2; factor =0;
		while (i++ < size)
		{
			if (atable[i-1].N == atable[i].N && atable[i-1].N0 == atable[i].N0 && atable[i].r0 == atable[i-1].r0 && atable[i].r1 ==atable[i-1].r1) factor -= log10(den++); 
			//if (atable[i-1].N == atable[i].N && (atable[i].r0+atable[i].r1 == atable[i-1].r0+atable[i-1].r1)) factor -= log10(den++); 
			else den =2;
		}
		if (clranew + factor <= clra+ifactor+epsilon) 	pvalue +=1; 
		clrasum += clranew+factor; clrasumsq += (clranew+factor)*(clranew+factor);
		if (pvalue >= 10 || (k <= 100 && pvalue >= 5)) break;
	}
	clrasum /= (k); clrasumsq /= (k); clrasumsq -= clrasum*clrasum; clrasumsq = sqrt(clrasumsq);
	*permutationpvalue = log10(pvalue+1)-log10(k+1);
	sstat[3] /= k; sstat[4] /=k;
	fprintf(stdout,"LL %.2f %.2f/%.2f %d/%d SPV %f:%f:%f\t",clra+ifactor,clrasum,clrasumsq,(int)pvalue,k,sstat[3],sstat[4],sstat[0]);
	//	fprintf(stdout,"\t");
	free(table_index); //free(atable);
	return 1;

	//	return (double)(pvalue+1)/(k+2); // apply Yates correction
}

//int main(char** argc,int argv) { return 1; }
