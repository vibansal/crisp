#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include "chisquare.h"
#include "kfunc.c"
#include "../common.h"
//int PIVOTSAMPLE = 0; 

// changed on aug 22 2012 to use samtools kfunc.c for calculating chi-square pvalue using upper-incomplete gamma function
// samtools kfunc.c modified to return the logarithm directly instead of doing exponent
//gcc -I../samtools -lm chisquare.c

// need to calculate stranded-chi square statistics for very low-frequency variant detection

// strata has information about samples sequenced on same run/platform etc (0,1,2,3), added nov 16 2012 
int chi2pvalue(double** ctable,int* strata, int size,double* chisqpvalue,double chistat[])
{
	int i=0; int index=0;
	double delta=0,deltaf=0,deltar=0,tf; //double d=0;
	double pvalue[3] = {0,0,0}; double SF,RUGF;
	double statisticjoint[2] = {0,0}; double statistic[2]={0,0}; double statisticf[2]={0,0},statisticr[2]={0,0};
	double statnew[2]={0,0},pairwisestat[2]={0,0};

	double  Ef[3]={0,0,0}; double Er[3] = {0,0,0}; 
	int computechi2 =1;//computechi2f=1,computechi2r=1;
	*chisqpvalue =1;
	double MINREADS = 25;
	double yatescorr =0.25;
	//int partitions =1;  for (i=0;i<size;i++) { if (strata[i] > partitions) partitions = strata[i]+1; } 

	double counts[4] = {0,0,0,0};	double counts1[4] = {0,0,0,0};
	double w1=0,w2=0;
	int DOF = size-2; // degress of freedom, if some pool has no reads, we should reduce the degrees of freedom
	if (PIVOTSAMPLE ==0) DOF++; 
	int DOF1 = PIVOTSAMPLE-1; int DOF2 = size-PIVOTSAMPLE-1; 
	
	for (i=0;i<size;i++) 
	{  
		//fprintf(stdout,"%d %.0f:%.0f %.0f:%.0f \n",i,ctable[i][0],ctable[i][1],ctable[i][2],ctable[i][3]);
		if (i < PIVOTSAMPLE) 
		{ 
			counts[0] += ctable[i][0]; counts[1] += ctable[i][1]; counts[2] += ctable[i][2]; counts[3] += ctable[i][3]; 
		} 
		else 
		{
			counts1[0] += ctable[i][0]; counts1[1] += ctable[i][1]; counts1[2] += ctable[i][2]; counts1[3] += ctable[i][3]; 
		}
	}

	//double w1 = sqrt(readsf)/sqrt(readsf+readsr); double w2 = sqrt(readsr)/sqrt(readsf+readsr); 
	//if (onesf+onesr >= 1 && readsf+readsr-onesf-onesr >= 1 && readsf+readsr >= MINREADS) computechi2 =1; else computechi2 = 0;
	//E = ((double)(onesf+onesr)+yatescorr)/((double)readsf+(double)readsr+yatescorr); 
	// problem occurs if Ef/Er = 1.0
	Ef[0] = (counts[1]+yatescorr)/(counts[0]+yatescorr); Er[0] = (counts[3]+yatescorr)/(counts[2]+yatescorr);
	Ef[1] = (counts1[1]+yatescorr)/(counts1[0]+yatescorr); Er[1] = (counts1[3]+yatescorr)/(counts1[2] + yatescorr);
	if (Ef[0] >=0.999) Ef[0] = 0.999; if (Ef[1] >=0.999) Ef[1] = 0.999; 
	if (Er[0] >=0.999) Er[0] = 0.999; if (Er[1] >=0.999) Er[1] = 0.999; 

	if (PIVOTSAMPLE > 0) fprintf(stdout,"Ef %0.5f %0.5f Er %0.5f %0.5f %0.0f/%0.0f %0.0f/%0.0f\n",Ef[0],Ef[1],Er[0],Er[1],counts[1],counts[0],counts[3],counts[2]);
	//fprintf(stdout,"Ef %0.5f %0.5f Er %0.5f %0.5f %0.0f/%0.0f %0.0f/%0.0f\n",Ef[0],Ef[1],Er[0],Er[1],counts1[1],counts1[0],counts1[3],counts1[2]);

	for (i=0;i<size;i++) 
	{
		if (i < PIVOTSAMPLE) index = 0; else index = 1; 

		//printf("%f %f \n",ctable[i][0],ctable[i][1]);
		if (computechi2 ==1 && ctable[i][0]+ctable[i][2] >= 1 && counts[0]+counts[2]+counts1[0]+counts1[2] >= MINREADS) // at least two reads
		{
			//R = ctable[i][1] + ctable[i][3]; N = ctable[i][0]+ctable[i][2]; statistic += (R-E*N)*(R-E*N)/(N*E*(1.0-E));
			deltaf =0; deltar=0;
			w1 = ctable[i][0]/sqrt(ctable[i][0]*ctable[i][0] + ctable[i][2]*ctable[i][2]); 
			w2 = sqrt(1.0-w1*w1);	
			//w1 = sqrt(ctable[i][0])/sqrt(ctable[i][0] + ctable[i][2]); w2 = sqrt(1.0-w1*w1);	
			if (ctable[i][0] > 0) deltaf = (ctable[i][1]-Ef[index]*ctable[i][0])/sqrt((1.0-Ef[index])*Ef[index]*ctable[i][0]);
			if (ctable[i][2] > 0) deltar = (ctable[i][3]-Er[index]*ctable[i][2])/sqrt((1.0-Er[index])*Er[index]*ctable[i][2]); 
			//statnew[index] += pow(w1*deltaf+w2*deltar,2); 
			statisticjoint[index] += pow(w1*deltaf+w2*deltar,2); 
			statisticf[index] += pow(deltaf,2); statisticr[index] += pow(deltar,2);
			pairwisestat[index] += 2*w1*w2*deltaf*deltar;
			if (i==size-1 && PIVOTSAMPLE ==0) fprintf(stdout,"pool %d %f %f %f %f stats %f:%f:%f\n",i,w1*deltaf+w2*deltar,deltaf,deltar,pairwisestat[index],statnew[index],statisticf[index],statisticr[index]);

			if ((counts1[0]+counts[0]) >= MINREADS || (counts1[2]+counts[2]) >= MINREADS)
			{
				// this is correct as it generates a chi-square statistic for 2xk tables 
				// but it is unweighted 
				delta = ctable[i][1]-Ef[index]*ctable[i][0] + ctable[i][3]-Er[index]*ctable[i][2]; 
				tf = Ef[index]*(1.0-Ef[index])*ctable[i][0] + Er[index]*(1.0-Er[index])*ctable[i][2]; 
				delta = w1*(ctable[i][1]-Ef[index]*ctable[i][0]) + w2*(ctable[i][3]-Er[index]*ctable[i][2]);
				tf = w1*w1*Ef[index]*(1.0-Ef[index])*ctable[i][0] + w2*w2*Er[index]*(1.0-Er[index])*ctable[i][2]; 
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

	//if (statistic > 0) chi2pvf = gammq( (double)(DOF-1)/2,statistic/2); else chi2pvf = 0;
	//if (statnew > 0) chi2pvr = gammq( (double)(DOF-1)/2,statnew/2); else chi2pvr = 0;
//	fprintf(stdout,"statistics %0.1f %0.1f %0.1f\n",statisticf,statisticr,statisticjoint); 
	//if (computechi2f ==1 && statisticf > 0) chi2pvf = gammq( (double)(size-1)/2,statisticf/2); else chi2pvf  =0; 
	//if (computechi2r ==1 && statisticr > 0) chi2pvr = gammq( (double)(size-1)/2,statisticr/2); else chi2pvr = 0;
	//fprintf(stdout,"chi:%.1f:%.2f:%.1f:%.1f",statisticjoint,*chisqpvalue,chi2pvf,chi2pvr);//,statnew,chi2pvr);
	//if (PFLAG >=1) fprintf(stdout,"chi:%.1f:%.2f",statisticjoint,*chisqpvalue);
	return computechi2;
}

//#define CHI_MAIN 1

#ifdef CHI_MAIN
int main(char** argc, int argv) 
{  
	double x = 600, y= 95; double pv = 0;
	int i=0;
	//for (i=1;i<=1000000;i++) pv = gammq(y,x); 
	//for (i=1;i<=1000000;i++) pv = log10(kf_gammaq(y,x)+TINY);
	//printf("chiswuare pvalue %lg \n",gammq(y,x)); 
	printf("chiswuare pvalue %lg \n",kf_gammaq(y,x)/log10); 
	return 1; 
}
#endif 
			
			/** full chi-square calculation
			delta =0;
			w1 = ctable[i][0]/sqrt(ctable[i][0]*ctable[i][0] + ctable[i][2]*ctable[i][2]); w2 = sqrt(1.0-w1*w1);	
			if (ctable[i][0] > 0) delta += w1*(ctable[i][1]-Ef*ctable[i][0])/sqrt(Ef*ctable[i][0]);
			if (ctable[i][2] > 0) delta += w2*(ctable[i][3]-Er*ctable[i][2])/sqrt(Er*ctable[i][2]); 
			statnew += delta*delta;
			
			delta =0;
			if (ctable[i][0] > 0) delta += w1*(ctable[i][0]-ctable[i][1]-(1.0-Ef)*ctable[i][0])/sqrt((1.0-Ef)*ctable[i][0]);
			if (ctable[i][2] > 0) delta += w2*(ctable[i][2]-ctable[i][3]-(1.0-Er)*ctable[i][2])/sqrt((1.0-Er)*ctable[i][2]); 
			statnew += delta*delta;
			**/
