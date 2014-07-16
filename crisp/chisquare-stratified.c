/* code to calculate chi-square heterogeneity p-value for multiple pools */
#include<stdio.h>
#include<math.h>
#include<stdlib.h>

// this code analyzes bidirectional allele counts as separate bin in chi-square calculation
// strata has information about samples sequenced on same run/platform etc (0,1,2,3), added nov 16 2012 
// we should stratify by read 1/2 for Illumina data in addition to strand 
// dynamic stratification also by position in read if there is sufficient data... 
int chi2pvalue_stratified(struct VARIANT* variant,int allele1,int allele2,int allele3,double* chisqpvalue,double chistat[])
{
	int i=0; int index=0; int a=0,b=0,c=0;
	double deltaf=0,deltar=0,deltab=0; //double d=0;
	double pvalue[3] = {0,0,0}; double SF,RUGF;
	double statisticjoint[2] = {0,0}; double statisticf[2]={0,0},statisticr[2]={0,0};
	double statnew[2]={0,0},pairwisestat[2]={0,0};

	double  Ef[3]={0,0,0}; double Er[3] = {0,0,0};  double Eb[3] = {0,0,0};
	int computechi2 =1;//computechi2f=1,computechi2r=1;
	*chisqpvalue =1;
	double MINREADS = 25; double yatescorr =0.25;
	double counts[6] = {0,0,0,0,0,0};	double counts1[6] = {0,0,0,0,0,0};
	double w1=0,w2=0,w3=0,nfactor;
	int DOF = variant->samples-1; // degress of freedom, if some pool has no reads, we should reduce the degrees of freedom
	
	for (i=0;i<variant->samples;i++) 
	{  
		//fprintf(stdout,"%d %.0f:%.0f %.0f:%.0f \n",i,ctable[i][0],ctable[i][1],ctable[i][2],ctable[i][3]);
		a= variant->indcounts[i][allele1] + variant->indcounts[i][allele2]+variant->indcounts[i][allele3];
		b= variant->indcounts[i][allele1+maxalleles] + variant->indcounts[i][allele2+maxalleles]+variant->indcounts[i][allele3+maxalleles];
		c= variant->indcounts[i][allele1+2*maxalleles] + variant->indcounts[i][allele2+2*maxalleles]+variant->indcounts[i][allele3+2*maxalleles];
		counts[0] += a; counts[1] += variant->indcounts[i][allele2];
		counts[2] += b; counts[3] += variant->indcounts[i][allele2+maxalleles];
		counts[4] += c; counts[5] += variant->indcounts[i][allele2+2*maxalleles];
	}
	// problem occurs if Ef/Er = 1.0
	Ef[0] = (counts[1]+yatescorr)/(counts[0]+yatescorr); Er[0] = (counts[3]+yatescorr)/(counts[2]+yatescorr); Eb[0] = (counts[5]+yatescorr)/(counts[4]+yatescorr);
	if (Ef[0] >=0.999) Ef[0] = 0.999;  if (Er[0] >=0.999) Er[0] = 0.999; if (Eb[0] >=0.999) Eb[0] = 0.999;  

	//fprintf(stdout,"Ef %0.5f %0.5f Er %0.5f %0.5f %0.0f/%0.0f %0.0f/%0.0f\n",Ef[0],Ef[1],Er[0],Er[1],counts1[1],counts1[0],counts1[3],counts1[2]);

	for (i=0;i<variant->samples;i++) 
	{
		index =0;
		a= variant->indcounts[i][allele1] + variant->indcounts[i][allele2]+variant->indcounts[i][allele3];
		b= variant->indcounts[i][allele1+maxalleles] + variant->indcounts[i][allele2+maxalleles]+variant->indcounts[i][allele3+maxalleles];
		c= variant->indcounts[i][allele1+2*maxalleles] + variant->indcounts[i][allele2+2*maxalleles]+variant->indcounts[i][allele3+2*maxalleles];

		if (computechi2 ==1 && a+b+c >= 1 && counts[0]+counts[2]+counts1[0]+counts1[2]+counts[4]+counts1[4] >= MINREADS) // at least two reads
		{
			deltaf =0; deltar=0; deltab=0; nfactor = sqrt(a*a+b*b+c*c);
			w1 = a/nfactor; w2 = b/nfactor; w3 = c/nfactor;
			if (a>0) deltaf = (variant->indcounts[i][allele2]-Ef[index]*a)/sqrt((1.0-Ef[index])*Ef[index]*a);
			if (b > 0) deltar = (variant->indcounts[i][allele2+maxalleles]-Er[index]*b)/sqrt((1.0-Er[index])*Er[index]*b); 
			if (c > 0) deltab = (variant->indcounts[i][allele2+2*maxalleles]-Eb[index]*c)/sqrt((1.0-Eb[index])*Er[index]*c); 
			statisticjoint[index] += pow(w1*deltaf+w2*deltar+w3*deltab,2);  // weighted sum of normal random variables ^2 = chi-square 
			statisticf[index] += pow(deltaf,2); statisticr[index] += pow(deltar,2);
			pairwisestat[index] += 2*w1*w2*deltaf*deltar;
		}
		else
		{
			// return that statistic was not computed
		}
	}
	//if (computechi2 ==1 && statisticjoint > 0) *chisqpvalue = gammq( (double)(DOF-1)/2,statisticjoint/2); else *chisqpvalue = 0;
	SF = statisticjoint[0]; 
	RUGF = kf_gammaq( (double)(DOF-1)/2,(double)SF/2);
	if (SF > 0) pvalue[2] = RUGF/log(10); else pvalue[2] = 0; *chisqpvalue = pvalue[2];
	chistat[0] = statisticjoint[0];  
	chistat[1] = statisticf[0]; chistat[2] = statisticr[0];  chistat[3] = statnew[0]; chistat[4] = pairwisestat[0];

	return computechi2;
}


// PIVOTSAMPLE is sample that divides the list of samples into two groups sequenced on two different runs...
// first implemented for PF 50 pools: 40 on one run and 10 on another
int chi2pvalue_stratified_pivotsample(struct VARIANT* variant,int allele1,int allele2,int allele3,double* chisqpvalue,double chistat[])
{
	int i=0; int index=0; int size =variant->samples;
	int a=0,b=0,c=0;
	double delta=0,deltaf=0,deltar=0,deltab=0,tf; //double d=0;
	double pvalue[3] = {0,0,0}; double SF,RUGF;
	double statisticjoint[2] = {0,0}; double statistic[2]={0,0}; double statisticf[2]={0,0},statisticr[2]={0,0};
	double statnew[2]={0,0},pairwisestat[2]={0,0};

	double  Ef[3]={0,0,0}; double Er[3] = {0,0,0};  double Eb[3] = {0,0,0};
	int computechi2 =1;//computechi2f=1,computechi2r=1;
	*chisqpvalue =1;
	double MINREADS = 25;
	double yatescorr =0.25;
	//int partitions =1;  for (i=0;i<size;i++) { if (strata[i] > partitions) partitions = strata[i]+1; } 

	double counts[6] = {0,0,0,0,0,0};	double counts1[6] = {0,0,0,0,0,0};
	double w1=0,w2=0,w3=0,nfactor;
	int DOF = size-2; // degress of freedom, if some pool has no reads, we should reduce the degrees of freedom
	if (PIVOTSAMPLE ==0) DOF++; 
	int DOF1 = PIVOTSAMPLE-1; int DOF2 = size-PIVOTSAMPLE-1; 
	
	for (i=0;i<size;i++) 
	{  
		//fprintf(stdout,"%d %.0f:%.0f %.0f:%.0f \n",i,ctable[i][0],ctable[i][1],ctable[i][2],ctable[i][3]);
		a= variant->indcounts[i][allele1] + variant->indcounts[i][allele2]+variant->indcounts[i][allele3];
		b= variant->indcounts[i][allele1+maxalleles] + variant->indcounts[i][allele2+maxalleles]+variant->indcounts[i][allele3+maxalleles];
		c= variant->indcounts[i][allele1+2*maxalleles] + variant->indcounts[i][allele2+2*maxalleles]+variant->indcounts[i][allele3+2*maxalleles];
		if (i < PIVOTSAMPLE) 
		{ 
			counts[0] += a; counts[1] += variant->indcounts[i][allele2];
			counts[2] += b; counts[3] += variant->indcounts[i][allele2+maxalleles];
			counts[4] += c; counts[5] += variant->indcounts[i][allele2+2*maxalleles];
		} 
		else 
		{
			counts1[0] += a; counts1[1] += variant->indcounts[i][allele2];
			counts1[2] += b; counts1[3] += variant->indcounts[i][allele2+maxalleles];
			counts1[4] += c; counts1[5] += variant->indcounts[i][allele2+2*maxalleles];
		}
	}
	// problem occurs if Ef/Er = 1.0
	Ef[0] = (counts[1]+yatescorr)/(counts[0]+yatescorr); Er[0] = (counts[3]+yatescorr)/(counts[2]+yatescorr); Eb[0] = (counts[5]+yatescorr)/(counts[4]+yatescorr);
	Ef[1] = (counts1[1]+yatescorr)/(counts1[0]+yatescorr); Er[1] = (counts1[3]+yatescorr)/(counts1[2] + yatescorr); Eb[1] = (counts1[5]+yatescorr)/(counts1[4]+yatescorr); 

	if (Ef[0] >=0.999) Ef[0] = 0.999; if (Ef[1] >=0.999) Ef[1] = 0.999; 
	if (Er[0] >=0.999) Er[0] = 0.999; if (Er[1] >=0.999) Er[1] = 0.999; 
	if (Eb[0] >=0.999) Eb[0] = 0.999; if (Eb[1] >=0.999) Eb[1] = 0.999; 

	if (PIVOTSAMPLE > 0 && PFLAG >=2) fprintf(stdout,"newf Ef %0.5f %0.5f Er %0.5f %0.5f %0.0f/%0.0f %0.0f/%0.0f %0.0f/%0.0f\n",Ef[0],Ef[1],Er[0],Er[1],counts[1],counts[0],counts[3],counts[2],counts[5],counts[4]);
	//fprintf(stdout,"Ef %0.5f %0.5f Er %0.5f %0.5f %0.0f/%0.0f %0.0f/%0.0f\n",Ef[0],Ef[1],Er[0],Er[1],counts1[1],counts1[0],counts1[3],counts1[2]);

	for (i=0;i<size;i++) 
	{
		if (i < PIVOTSAMPLE) index = 0; else index = 1; 
		a= variant->indcounts[i][allele1] + variant->indcounts[i][allele2]+variant->indcounts[i][allele3];
		b= variant->indcounts[i][allele1+maxalleles] + variant->indcounts[i][allele2+maxalleles]+variant->indcounts[i][allele3+maxalleles];
		c= variant->indcounts[i][allele1+2*maxalleles] + variant->indcounts[i][allele2+2*maxalleles]+variant->indcounts[i][allele3+2*maxalleles];

		if (computechi2 ==1 && a+b+c >= 1 && counts[0]+counts[2]+counts1[0]+counts1[2]+counts[4]+counts1[4] >= MINREADS) // at least two reads
		{
			deltaf =0; deltar=0; deltab=0; nfactor = sqrt(a*a+b*b+c*c);
			w1 = a/nfactor; w2 = b/nfactor; w3 = c/nfactor;
			if (a>0) deltaf = (variant->indcounts[i][allele2]-Ef[index]*a)/sqrt((1.0-Ef[index])*Ef[index]*a);
			if (b > 0) deltar = (variant->indcounts[i][allele2+maxalleles]-Er[index]*b)/sqrt((1.0-Er[index])*Er[index]*b); 
			if (c > 0) deltab = (variant->indcounts[i][allele2+2*maxalleles]-Eb[index]*c)/sqrt((1.0-Eb[index])*Er[index]*c); 
			//statnew[index] += pow(w1*deltaf+w2*deltar,2); 
			statisticjoint[index] += pow(w1*deltaf+w2*deltar+w3*deltab,2); 
			statisticf[index] += pow(deltaf,2); statisticr[index] += pow(deltar,2);
			pairwisestat[index] += 2*w1*w2*deltaf*deltar;
			if (i==size-1 && PIVOTSAMPLE ==0 && PFLAG >=20) fprintf(stdout,"pool %d %f %f %f %f stats %f:%f:%f\n",i,w1*deltaf+w2*deltar,deltaf,deltar,pairwisestat[index],statnew[index],statisticf[index],statisticr[index]);

			if ((counts1[0]+counts[0]) >= MINREADS || (counts1[2]+counts[2]) >= MINREADS)
			{
				// this is correct as it generates a chi-square statistic for 2xk tables 
				// but it is unweighted 
				delta = variant->indcounts[i][allele2]-Ef[index]*a + variant->indcounts[i][allele2+maxalleles]-Er[index]*b; 
				tf = Ef[index]*(1.0-Ef[index])*a + Er[index]*(1.0-Er[index])*b; 
				delta = w1*(variant->indcounts[i][allele2]-Ef[index]*a) + w2*(variant->indcounts[i][allele2+maxalleles]-Er[index]*b);
				tf = w1*w1*Ef[index]*(1.0-Ef[index])*a + w2*w2*Er[index]*(1.0-Er[index])*b; 
				//statisticjoint[index] += delta*delta/tf; 
				statnew[index] += delta*delta/tf; 
			//fprintf(stdout,"tf %f Ef %f Er %f ct %d %d delta %f\n",tf,Ef,Er,ctable[i][0],ctable[i][2],delta);
			}
		}
		//else DOF--;
	}
	//if (computechi2 ==1 && statisticjoint > 0) *chisqpvalue = gammq( (double)(DOF-1)/2,statisticjoint/2); else *chisqpvalue = 0;
	SF = statisticjoint[0] + statisticjoint[1]; 
	RUGF = kf_gammaq( (double)(DOF-1)/2,(double)SF/2);
	if (SF > 0) pvalue[2] = RUGF/log(10); else pvalue[2] = 0; *chisqpvalue = pvalue[2];
	chistat[0] = statisticjoint[0]+statisticjoint[1]; 
	chistat[1] = statisticf[0]+statisticf[1]; chistat[2] = statisticr[0]+statisticr[1]; 
	chistat[3] = statnew[0]+statnew[1]; chistat[4] = pairwisestat[0]+pairwisestat[1];

	if (PIVOTSAMPLE > 0)
	{
		fprintf(stdout,"chi2 %0.2f:%.2f ",SF,pvalue[2]);//,statnew,chi2pvr);
		RUGF = kf_gammaq( (double)(DOF1-1)/2,(double)statisticjoint[0]/2);
		if (statisticjoint[0] > DOF1-1) pvalue[0] = RUGF/log(10); else pvalue[0] = 0;
		RUGF = kf_gammaq( (double)(DOF2-1)/2,(double)statisticjoint[1]/2);
		if (statisticjoint[1] > DOF2-1) pvalue[1] = RUGF/log(10); else pvalue[1] = 0;
		fprintf(stdout,"%d:%0.2f:%.2f %d:%0.2f:%.2f ",DOF1,statisticjoint[0],pvalue[0],DOF2,statisticjoint[1],pvalue[1]);
		if (pvalue[0] < pvalue[2] && pvalue[0] < pvalue[1]) 
		{
			*chisqpvalue = pvalue[0]; 
			chistat[0] = statisticjoint[0]; chistat[1] = statisticf[0]; chistat[2] = statisticr[0]; 
			chistat[3] = statnew[0]; chistat[4] = pairwisestat[0];
		}
		else if (pvalue[1] < pvalue[2]) 
		{
			*chisqpvalue = pvalue[1];
			chistat[0] = statisticjoint[1]; chistat[1] = statisticf[1]; chistat[2] = statisticr[1]; 
			chistat[3] = statnew[1]; chistat[4] = pairwisestat[1];
		}
		fprintf(stdout,"\n");
	}

	return computechi2;
}

