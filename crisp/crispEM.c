#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#define MINP 1e-6
#define LOWLL -100000
// need to extend indel method for multi-allelic indels in homopolymer runs... +CA, -CA 
// incorporate heterozygosity bias for indels, reference reads are aligned better
// pooled variance where the frequency of rare allele can be below 1/poolsize -> 1/2*poolsize 
/* functions for finding best population allele frequency and genotypes for pooled sequence data */

struct VARIANT* BFGSvariant; // global variable for use by BFGS function optimizer 
double THETA=0.001;  // mutation rate per site
double ALPHA =0.001;
double THRESH = 2;
int STRAND_CALC=2;

//int alpha_1 =1; int beta_1 =25;
//double betanorm =0;
//int USEPRIORS =1;

//#include "BFGScode.c"
#include "../FET/lowcovFET.c"

double AGILENT_REF_BIAS = 0.50;  // reference bias = fraction of reads for reference at heterozygous sites...

int calculate_pooled_likelihoods_indel(struct VARIANT* variant,int allele1, int allele2, int allele3,double E01[], double E10[]);
void calculate_pooled_likelihoods_snp(READQUEUE* bq,struct VARIANT* variant,int allele1, int allele2, int allele3,double E01[],double E10[]);
void calculate_genotypes(struct VARIANT* variant,int allele1, int allele2, int allele3,double p,double* HWEpvalue);


#include "LRstatistic.c"

// calculate bias in favor of mapping reference allele for long indels | based on the length of the indel
double indel_bias(struct VARIANT* variant,int allele) 
{
	int indellength = strlen(variant->itb[allele])-1; double het = 0;
	double READLENGTH = 100; double MINFLANK =5;
        if (variant->itb[allele][0] == '+' && indellength >= 1) 
        {
                het += (MINFLANK+indellength/2)/(double)(READLENGTH-MINFLANK-indellength/2);
        }
        else if (indellength >= 1 && variant->itb[allele][0] == '-')
        {
                het += (MINFLANK)/(double)(READLENGTH-MINFLANK);
        }
	return het;
}

// now model e01, e10 => could be different 
// for third allele, add E02 extra parameter, E12 = E10 = E20 
// s1 and s2 are range of samples for which likelihood is calculated [s1,s2) -> added sept 4 2013
int calculate_pooled_likelihoods_indel(struct VARIANT* variant,int allele1, int allele2, int allele3,double E01[], double E10[])
{
	int i=0,sample=0; double ef,er,eb,f;
	int indellength = strlen(variant->itb[allele2])-1; double het = 0.5;
	MINFLANK += variant->HPlength[allele2]; 
	if (variant->itb[allele2][0] == '+') het += 0.5*(MINFLANK+indellength/2)/(double)(READLENGTH-MINFLANK-indellength/2);
	else if (variant->itb[allele2][0] == '-') het += 0.5*(MINFLANK)/(double)(READLENGTH-MINFLANK);
	//fprintf(stdout,"het %f HP %d %d %f:%f:%f\n",het,variant->HPlength[allele2],indellength,epf,epr,epb);
	MINFLANK -= variant->HPlength[allele2];
	for (i=0;i<3;i++) 
	{
		if (E01[i] < MINP) E01[i]= MINP; else if (E01[i] > 1.0-MINP) E01[i] = 1.0-MINP; 
		if (E10[i] < MINP) E10[i]= MINP; else if (E10[i] > 1.0-MINP) E10[i] = 1.0-MINP; 
	}
	// product of error probabilities for bidirectional reads instead of independent estimate
	E01[2] = E01[0]*E01[1];  E10[2] = E10[0]*E10[1]; 
	double het_bias = 0;
	het_bias = indel_bias(variant,allele2); //fprintf(stdout,"het bias %f %d \n",het_bias,strlen(variant->itb[allele2])-1);

	for (sample=0;sample<variant->samples;sample++)
	{
		for (i=0;i<=variant->ploidy[sample];i++)
		{
			f = (double)i/variant->ploidy[sample]; f /= (1+het_bias); // reduce f to account for indel bias 
			ef = (1.0-E01[0])*(1.0-f) + E10[0]*f; er = (1.0-E01[1])*(1.0-f) + E10[1]*f; eb = (1.0-E01[2])*(1.0-f) + E10[2]*f;
			// add factor for reference bias for indels..
			variant->GENLLf[sample][i] = variant->indcounts[sample][allele1]*log10(ef) + variant->indcounts[sample][allele2]*log10(1.0-ef);
			variant->GENLLr[sample][i] = variant->indcounts[sample][allele1+maxalleles]*log10(er) + variant->indcounts[sample][allele2+maxalleles]*log10(1.0-er);
			variant->GENLLb[sample][i] = variant->indcounts[sample][allele1+2*maxalleles]*log10(eb) + variant->indcounts[sample][allele2+2*maxalleles]*log10(1.0-eb);
			variant->GENLL[sample][i] = variant->GENLLf[sample][i] + variant->GENLLr[sample][i] + variant->GENLLb[sample][i]; 
		}
	}
	return 1;
}


//this function can be modified to use error rates ef/er for allele2 instead of base quality values
void calculate_pooled_likelihoods_snp(READQUEUE* bq,struct VARIANT* variant,int allele1, int allele2, int allele3,double E01[],double E10[])
{
	int i=0,sample=0; int base; double f,e,ep,le; 
	struct alignedread* bcall; //fprintf(stdout,"bcall %s %d\n",bcall->readid,bcall->position);
	for (sample=0;sample<variant->samples;sample++)
	{
		for (i=0;i<=variant->ploidy[sample];i++) 
		{
			variant->GENLL[sample][i] = 0.0; variant->GENLLf[sample][i]= 0.0; variant->GENLLr[sample][i] = 0.0; variant->GENLLb[sample][i]=0.0;
		}
	}

	//  table for log10(e) and log10(1-e) for all e that correspond to quality value 0-60 and f = 0 to 1 (grid) will speed up things
	for (bcall = bq->first; bcall != NULL && bcall->position <= variant->position; bcall = bcall->next)
	{
		if (bcall->lastpos < variant->position) continue;
		if (bcall->filter == '0' && bcall->type ==0)
		{
			base = BTI[bcall->sequence[bcall->l1+bcall->delta]]; 
			ep = pow(0.1,((double)bcall->quality[bcall->l1+bcall->delta]-QVoffset)/10);
			//if (ep < 0.001) ep = 0.001;
			for (i=0;i<=variant->ploidy[bcall->sampleid];i++)
			{
				//for i=1, f actually is a beta distribution (alpha,(variant->ploidy-1)*alpha), so we need to integrate over that... 
				// do this for pooled samples with for which the allele frequency is > 0 and < 2/variant->ploidy, discrete 100 bin integration...
				//f = (double)i/variant->ploidy[bcall->sampleid];  // assumes equal proportion of DNA from each pooled sample, 
				f = (double)i*(1.0-AGILENT_REF_BIAS); f /= f + (double)(variant->ploidy[bcall->sampleid]-i)*AGILENT_REF_BIAS; 
				e = (1.0-ep)*(1.0-f) + ep*f;
				if (base == allele1) 
				{
					le = log10(e); variant->GENLL[bcall->sampleid][i] += le;
					if (bcall->bidir ==1) variant->GENLLb[bcall->sampleid][i] += le;
					else if (bcall->strand ==0) variant->GENLLf[bcall->sampleid][i] += le;
					else variant->GENLLr[bcall->sampleid][i] += le;
				}
				else if (base == allele2) 
				{
					le = log10(1.0-e); variant->GENLL[bcall->sampleid][i] += le;
					if (bcall->bidir ==1) variant->GENLLb[bcall->sampleid][i] += le;
					else if (bcall->strand ==0) variant->GENLLf[bcall->sampleid][i] += le;
					else variant->GENLLr[bcall->sampleid][i] += le; 
					//if (bcall->sampleid ==19 && i ==1 && bcall->strand ==1) fprintf(stdout,"%f %f base %d %f \n",variant->GENLLr[bcall->sampleid][0],variant->GENLLr[bcall->sampleid][1],base,ep);
				}
			}
		}
	}
	//for (sample=0;sample<variant->samples;sample++) fprintf(stdout,"LL %f %d\n",variant->GENLL[sample][0],variant->ploidy[sample]);
}

// once we have estimated 'p', assign the ML genotype for each pool, also calculate HWE statistic 
void calculate_genotypes(struct VARIANT* variant,int allele1, int allele2, int allele3,double p,double* HWEpvalue)
{
	int i=0,j=0;
	double HWEstat = 0; double ecount =0; double DOF =0;
	int maxploidy =0; for (i=0;i<variant->samples;i++) { if (variant->ploidy[i] > maxploidy) maxploidy = variant->ploidy[i]; } 
        double* genotype_counts = calloc(sizeof(double),(maxploidy+1)); for (j=0;j<=maxploidy;j++) genotype_counts[j] = 0;

	int bestgenotype =0,sbgenotype =0; 
	double temp=0,Lb=0,bgLL=-100000,sbgLL=-100000;  
	//double meanAF=0,varianceAF=0;//mean allele frequency of pool and its variance 
	// calculate HWE statistic 

	variant->crispvar[variant->varalleles].variantpools = 0; variant->crispvar[variant->varalleles].allelecount = 0; 
	for (i=0;i<variant->samples;i++)
	{
		Lb=-100000; bestgenotype =0; sbgenotype=0; bgLL=-100000; sbgLL =-100000;
		for (j=0;j<=variant->ploidy[i];j++) 
		{
			temp = variant->GENLL[i][j] + variant->NCarray[i][j] + log10(p)*j + log10(1.0-p)*(variant->ploidy[i]-j); 
			if (temp > bgLL) 
			{
				sbgenotype = bestgenotype; sbgLL = bgLL; bestgenotype = j; bgLL = temp;
			}
			else if (temp > sbgLL)
			{
				sbgenotype = j; sbgLL = temp;
			}
			if (temp > Lb || j ==0) Lb = temp + log10(1.0+pow(10,Lb-temp));
			else Lb += log10(1.0+pow(10,temp-Lb));
		}
		for (j=0;j<=variant->ploidy[i];j++) 
		{
			temp = variant->GENLL[i][j] + variant->NCarray[i][j] + log10(p)*j + log10(1.0-p)*(variant->ploidy[i]-j); 
                        genotype_counts[j] += pow(10,temp-Lb);
		}
		//meanAF -= Lb; varianceAF -=Lb; varianceAF += log10(1.0-pow(10,2*meanAF-varianceAF));
		//LLtotalf += Lf; LLtotalr += Lr;
		variant->crispvar[variant->varalleles].AC[i] = bestgenotype; 
		variant->crispvar[variant->varalleles].QV[i] =10*(bgLL-sbgLL); 
		//variant->crispvar[variant->varalleles].meanAC[i] = pow(10,meanAF); 
		//variant->crispvar[variant->varalleles].varAC[i] = pow(10,varianceAF); 
		if (variant->crispvar[variant->varalleles].AC[i] > 0) variant->crispvar[variant->varalleles].variantpools++;
		variant->crispvar[variant->varalleles].allelecount += variant->crispvar[variant->varalleles].AC[i];
		//if (OUTPUTGENOTYPES==1 && vardelta >= THRESH) fprintf(stdout,"%d:%0.1f:%0.4f:%0.2f ",bestgenotype,10*(bgLL-sbgLL),pow(10,meanAF),varianceAF);
	}

	if (variant->varpoolsize ==0)  // calculate HWE only if all pool sizes are equal 
	{
		i=0; // dummy index	
		for (j=0;j<=variant->ploidy[i];j++)
		{
			ecount = pow(10,variant->NCarray[i][j] + j*log10(p) + (variant->ploidy[i]-j)*log10(1.0-p) + log10(variant->samples));
			if (ecount >= 2)
			{
				HWEstat += pow(genotype_counts[j]-ecount,2)/ecount; DOF +=1;
				fprintf(stdout,"%d:%.1f:%.1f ",j,genotype_counts[j],ecount);
			}
		}
		if (HWEstat > 0 && DOF >=2) *HWEpvalue = kf_gammaq( (double)(DOF-1)/2,(double)HWEstat/2)/log(10);  else *HWEpvalue =0; 
		fprintf(stdout,"HWEcounts %0.2f %0.0f %0.2f\n",HWEstat,DOF-1,*HWEpvalue);       
	}
	//free(genotype_counts); // check free memory here ??
}

// calculate genotype likelihood for all pools for this allele frequency 'p' (subtract likelihood for G=0) 
// Pr (Data | p) = product (Pr (data_i |p)) = \prod [ \sum_k Pr(data_i | k var allele) Pr (k var alleles | p)  ]
double EMmethod(struct VARIANT* variant,int allele1, int allele2, int allele3,double p)
{
	// for two variant alleles, need p1,p2
	int i=0,j=0,k=0,K=0,iter=0; 
	double E01[3] = {0.001,0.001,0.001}; double E10[3] = {0.001,0.001,0.001}; 
	double T00[3] = {0,0,0}; double T01[3] = {0,0,0}; double T10[3] = {0,0,0}; double T11[3] = {0,0,0};
	double f=0;
	double pnew=0,pnewf=0,pnewr=0,pnewb=0; double deltaf=0,deltar=0; double Lf=0,Lr=0,Lj=0,Lb=0;
	double LLtotal =0,LLnullf=0,LLnullr=0,LLnullb=0;
	double LLtotalf=0,LLtotalr=0,LLtotalb=0,LLtotalprev=0,LLtotalprevf=0,LLtotalprevr=0,LLtotalprevb=0,LLnullprev=0;
	double LLrefsub=0,LLrefadd=0;
	double temp=0,tempf=0,tempr=0,tempb=0; double LLnull=0;
	double yatescorr = 1.0; // yates correction 
	double LLnull_noprior[3]; 
	double pf=p,pr=p,pb=p; double SDELTA=0; double vardelta =0; double HWEpvalue =0;

	ALPHA = 0.0; for (i=0;i<variant->samples;i++) K += variant->ploidy[i]; for (j=1;j<=K;j++) ALPHA += 1.0/j; 
	ALPHA *= THETA; // prior for non-reference genotypes 

	int maxiter = 1000; double convergence_delta = 0.01; int exitloop =0;
	double alpha = 1.5,beta = 50;	 // priors for EM pseudocounts, change them for indels in homopolymer runs
        double log_p[8];

	//what is the likelihood function that we are maximizing (with prior or without)

	// EM main loop
	while (iter++ < maxiter && exitloop ==0)
	{
		LLtotalprev = LLtotal; 	LLtotalprevf = LLtotalf; LLtotalprevr= LLtotalr; LLtotalprevb = LLtotalb; LLnullprev = LLnull;
		pnew = 0; pnewf = 0; pnewr= 0;
		LLtotal = 0; LLtotalf = 0; LLtotalr= 0; LLtotalb=0; LLnull = 0; LLnullf=0; LLnullr=0;
		for (i=0;i<3;i++) { T00[i]= T01[i] = T10[i] = T11[i]= 0; } 
		LLnull_noprior[0] = 0.0;
		log_p[0] = log10(p); log_p[1]= log10(1.0-p); log_p[2] = log10(pf); log_p[3]= log10(1.0-pf);
                log_p[4] = log10(pr); log_p[5]= log10(1.0-pr); log_p[6] = log10(pb); log_p[7]= log10(1.0-pb);
		for (i=0;i<variant->samples;i++)
		{
			Lb = -100000;Lf = -100000;Lr = -100000; Lj = -100000;
			LLnull_noprior[0] += variant->GENLL[i][0];
			// for low frequency variants, we don't need to go from 0 to variant->ploidy...we can truncate
			for (j=0;j<=variant->ploidy[i];j++) 
			{
				temp = variant->GENLL[i][j] + variant->NCarray[i][j] + log_p[0]*j + log_p[1]*(variant->ploidy[i]-j);
                                tempf = variant->GENLLf[i][j] + variant->NCarray[i][j] + log_p[2]*j + log_p[3]*(variant->ploidy[i]-j);
                                tempr = variant->GENLLr[i][j] + variant->NCarray[i][j] + log_p[4]*j + log_p[5]*(variant->ploidy[i]-j);
                                tempb = variant->GENLLr[i][j]+ variant->GENLLf[i][j] + variant->NCarray[i][j] + log_p[6]*j + log_p[7]*(variant->ploidy[i]-j);

				if (j==0) LLnull += temp; if (j==0) LLnullf += tempf; if (j==0) LLnullr += tempr; 
				if (j==0) Lj = temp; else if (temp > Lj) Lj = temp + log10(1.0+pow(10,Lj-temp)); else Lj += log10(1.0+pow(10,temp-Lj));
				if (j==0) Lf = tempf; else if (tempf > Lf) Lf = tempf + log10(1.0+pow(10,Lf-tempf)); else Lf += log10(1.0+pow(10,tempf-Lf));
				if (j==0) Lr= tempr; else if (tempr > Lr) Lr = tempr + log10(1.0+pow(10,Lr-tempr)); else Lr += log10(1.0+pow(10,tempr-Lr));
				if (j==0) Lb= tempb; else if (tempb > Lb) Lb = tempb + log10(1.0+pow(10,Lb-tempb)); else Lb += log10(1.0+pow(10,tempb-Lb));
			}
			LLtotal += Lj; LLtotalf += Lf; LLtotalr += Lr; LLtotalb += Lb;
			LLtotalr += variant->GENLLf[i][0]+variant->GENLLb[i][0]; 
			LLtotalf += variant->GENLLr[i][0]+ variant->GENLLb[i][0]; 
			LLtotalb += variant->GENLLb[i][0];
			for (j=0;j<=variant->ploidy[i];j++) 
			{
				f = (double)j/variant->ploidy[i]; 
				temp = variant->GENLL[i][j] + variant->NCarray[i][j] + log_p[0]*j + log_p[1]*(variant->ploidy[i]-j);
                                tempf = variant->GENLLf[i][j] + variant->NCarray[i][j] + log_p[2]*j + log_p[3]*(variant->ploidy[i]-j);
                                tempr = variant->GENLLr[i][j] + variant->NCarray[i][j] + log_p[4]*j + log_p[5]*(variant->ploidy[i]-j);
                                tempb = variant->GENLLr[i][j]+ variant->GENLLf[i][j] + variant->NCarray[i][j] + log_p[6]*j + log_p[7]*(variant->ploidy[i]-j);
				pnew += (double)j*pow(10,temp-Lj); 
				pnewf += (double)j*pow(10,tempf-Lf); pnewr += (double)j*pow(10,tempr-Lr); pnewb += (double)j*pow(10,tempb-Lb);

				for (k=0;k<3;k++) T00[k] += (double)variant->indcounts[i][allele1+k*maxalleles]*pow(10,temp-Lj)*(1.0-f)*(1.0-E01[k])/( (1.0-f)*(1.0-E01[k]) + f*E10[k]);
                                for (k=0;k<3;k++) T10[k] += (double)variant->indcounts[i][allele1+k*maxalleles]*pow(10,temp-Lj)*f*E10[k]/( (1.0-f)*(1.0-E01[k]) + f*E10[k]);

                                for (k=0;k<3;k++) T11[k] += (double)variant->indcounts[i][allele2+k*maxalleles]*pow(10,temp-Lj)*f*(1.0-E10[k])/( (1.0-f)*E01[k] + f*(1.0-E10[k]));
                                for (k=0;k<3;k++) T01[k] += (double)variant->indcounts[i][allele2+k*maxalleles]*pow(10,temp-Lj)*(1.0-f)*E01[k]/ ( (1.0-f)*E01[k] + f*(1.0-E10[k]));
			}
		}
		pnew /= K; pnewf /= K; pnewr /= K; pnewb /= K;

		for (i=0;i<3;i++)
                {
			//LLtotal += log10(1-E01[i])*50 + log10(E01[i])*0.5; LLtotal += log10(1-E10[i])*50 + log10(E10[i])*0.5; 
                        //if (T01[i]+T00[i] < 5) { T01[i] += 0.002; T00[i] += 1; } 
                        //if (T11[i]+T10[i] < 5) { T10[i] += 0.002; T11[i] += 1; } 
			yatescorr =0;
			if (T01[i]+T00[i] < 25) E01[i] = 0.005; else E01[i] = (T01[i]+yatescorr)/(T01[i]+T00[i]+yatescorr); 
			if (T11[i]+T10[i] < 25) E10[i] = 0.005; else E10[i] = (T10[i]+yatescorr)/(T10[i] + T11[i]+yatescorr);
		
			E01[i] = (T01[i]+alpha-1)/(T01[i]+T00[i]+alpha+beta-1);  E10[i] = (T10[i]+alpha-1)/(T10[i] + T11[i]+beta+alpha-1);
			if (E10[i] > 0.05) E10[i] = 0.05; 
			//if (E01[i] > 0.05) E01[i] = 0.05;  
                }

		// for SNPs, only variant allele probability should be used, i.e. Pr(observing variant | true = reference) due to errors 
		if (allele2 >=4 || USE_BASE_QVS==0) calculate_pooled_likelihoods_indel(variant,allele1,allele2,allele3,E01,E10);
		if (pnew < MINP) pnew= MINP; if (1.0-pnew <=MINP) pnew = 1.0-MINP;
		if (pnewf < MINP) pnewf= MINP; if (1.0-pnewf <=MINP) pnewf = 1.0-MINP;
		if (pnewr < MINP) pnewr= MINP; if (1.0-pnewr <=MINP) pnewr = 1.0-MINP;
		if (pnewb < MINP) pnewb= MINP; if (1.0-pnewb <=MINP) pnewb = 1.0-MINP;
		//if (LLtotal < LLtotalprev && LLtotalprev < -1) fprintf(stdout,"LLdecrease! %f \n",LLtotalprev-LLtotal);
		if (iter >=2 && fabsf(LLtotal-LLtotalprev) <= 0.0001 && fabsf(LLtotalf-LLtotalprevf) <= 0.01 && fabsf(LLtotalr-LLtotalprevr) <=0.01) exitloop = 1;
		if (PFLAG >=1 && (iter < 10 || iter%20 ==0 || exitloop ==1) ) 
		{
			fprintf(stdout,"AF %0.6f %0.6f LLtotal %0.6f LLnull %0.6f %0.2f:%0.2f:%0.2f E01 %0.5f:%0.5f:%0.5f E10 %0.5f:%0.5f:%0.5f \n",p,pnew,LLtotal,LLnull,LLtotalf,LLtotalr,LLtotalb,E01[0],E01[1],E01[2],E10[0],E10[1],E10[2]);
			//fprintf(stdout,"AF %0.6f %0.6f LLtotal %0.3f LLnull %0.2f %0.2f:%0.2f:%0.2f E01 %0.5f:%0.5f:%0.5f E10 %0.5f:%0.5f:%0.5f \n",p,pnew,LLtotal,LLnull,LLtotalf,LLtotalr,LLtotalb,E01[0],E01[1],E01[2],E10[0],E10[1],E10[2]);
			if (exitloop ==1) fprintf(stdout,"AF %0.6f %0.6f LLtotal %0.3f LLnull %0.2f %0.2f:%0.2f:%0.2f\n",p,pnew,LLtotalprev,LLnullprev,LLtotalprevf,LLtotalprevr,LLtotalprevb);
			//fprintf(stdout,"T vals %f %f %f %f \n",T01[0],T00[0],T10[0],T11[0]);
		}
		p = pnew; pf = pnewf; pr = pnewr; pb = pnewb; 
		//if (iter >=2 && fabsf(LLtotal-LLtotalprev) <= convergence_delta) break; 
	}
	if (iter < maxiter) fprintf(stdout,"convergence at iter %d\n",iter); else fprintf(stdout,"EM didnot converge after %d iters\n",maxiter);

        SDELTA = (LLtotalf> LLtotalr) ? LLtotalf-LLtotal-log10(2.0) : LLtotalr-LLtotal-log10(2.0); deltaf = SDELTA;

	if (STRAND_CALC ==20) 
	{
		SDELTA = (LLtotalf> LLtotalr) ? LLtotalf-LLtotal : LLtotalr-LLtotal;
		// check if bidirectional genotype likelihoods are much better for reference genotypes 
		if (variant->counts[allele1+2*maxalleles] +variant->counts[allele2+2*maxalleles] >= 5 && (LLtotalb - LLtotal > SDELTA) && allele2 <4 && SOLID ==0) SDELTA = LLtotalb-LLtotal;  // only for SNPs for now since for indels, OPE reads are not properly handled
		deltaf = SDELTA;
	}
	if (PFLAG >=1) fprintf(stdout,"%0.6f:%0.6f:%0.6f LL %0.2f:%0.2f:%0.2f:%0.2f SDELTA %0.1f \n",p,pf,pr,LLtotal,LLtotalf,LLtotalr,LLtotalb,SDELTA);

	// subtract the likelihood of G0 (reference genotypes) from the total sum and add sum with different prior 
	//LLrefsub = LLnull_noprior[0] + log10(1.0-p)*variant->ploidy*variant->samples;  
	LLrefadd = LLnull_noprior[0] + log10(1.0-ALPHA); 
	LLtotal += log10(ALPHA);// - log10(1.0-pow(1.0-p,variant->ploidy*variant->samples)); 
	if (LLrefadd < LLtotal) LLtotal += log10(1.0+pow(10,LLrefadd-LLtotal));
	else LLtotal = LLrefadd + log10(1.0+pow(10,LLtotal-LLrefadd));
	vardelta = LLtotal - LLrefadd;  // -log10(posterior probability of reference genotype configuration) 
	//fprintf(stdout,"sub %f add %f %f\n",LLrefsub,LLrefadd,LLtotal);

	if (PFLAG >=1) fprintf(stdout,"delta %0.2f p %0.4f LL %0.2f %0.2f ",vardelta,p,LLtotal,LLnull_noprior[0]);
	if (vardelta >= 2 && PFLAG >=1) fprintf(stdout,"pvariant\n");
	if (vardelta >=2) calculate_genotypes(variant,allele1,allele2,allele3,p,&HWEpvalue);
	// if poolsizes are variable, HWE statistic is not calculated 

	//store information about allele freq and deltas for crisp vcf 
	variant->crispvar[variant->varalleles].AF_ML=p; variant->crispvar[variant->varalleles].delta=vardelta;
	variant->crispvar[variant->varalleles].deltaf=deltaf; variant->crispvar[variant->varalleles].deltar=HWEpvalue;
	variant->crispvar[variant->varalleles].allele =allele2;
	return p;
}

// main function called from crispcaller.c for detection of variants 
// calculate likelihood for pooled variant by finding the allele frequency 'p' that maximizes Pr(Data over all pools |p) using EM algorithm 
// initial estimate of 'p' is provided using read count estimate
int compute_GLL_pooled(READQUEUE* bq, struct VARIANT* variant,int allele1,int allele2,int allele3,double p[])
{
	int i=0;
	int K=0; for (i=0;i<variant->samples;i++) K += variant->ploidy[i]; 
	p[0] /= K; if (1.0-p[0] < MINP) p[0]=1.0-MINP; p[1] /= K; if (1.0-p[1] < MINP) p[1]=1.0-MINP;
	double yatescorr = 1.0;
	double E01[3] = {0.001,0.001,0.001}; double E10[3] = {0.001,0.001,0.001}; 

	if (allele2 >=4 || USE_BASE_QVS==0) // likelihood calculation for indels
	{
		E01[0]  = (double)(variant->counts[allele2]+yatescorr)/(variant->counts[allele1]+variant->counts[allele2]+yatescorr);
		E01[1]  = (double)(variant->counts[allele2+maxalleles]+yatescorr)/(variant->counts[allele1+maxalleles]+variant->counts[allele2+maxalleles]+yatescorr);
		E01[2]  = (double)(variant->counts[allele2+2*maxalleles]+yatescorr)/(variant->counts[allele1+2*maxalleles]+variant->counts[allele2+2*maxalleles]+yatescorr);
		// check if ef,er are too high that may lead to EM algorithm converging to local optima
		if (E01[0] >=0.005 || E01[1] >=0.005) 
		{
			E01[0]=0.005; E01[1]=0.005; E01[2] = 0.005; 
		}
		E10[0] = E01[0]; E10[1]= E01[1]; E10[2] = E01[2]; 
		calculate_pooled_likelihoods_indel(variant,allele1,allele2,allele3,E01,E10);
	}
	else if (allele2 < 4) calculate_pooled_likelihoods_snp(bq,variant,allele1,allele2,allele3,E01,E10);

	// call EM algorithm for variant calling //p[0]  = 0.0001; 
	double p_EM= EMmethod(variant,allele1,allele2,allele3,p[0]); 

	// code for calculating case-control likelihood ratio statistic, added sept 4 2013
	if (variant->options->association ==1 && p_EM >= 0.0001 && variant->crispvar[variant->varalleles].delta >=3)
	{
		double pvalue = calculate_LRstatistic(variant,allele1,allele2,allele3,p_EM,E01,E10); 
	}
	// prior for genotyping each pool should be based on case/control status -> otherwise undercall rare variants that are only present in cases (or controls)

	p[0] = p_EM;
	if (variant->crispvar[variant->varalleles].delta >=2 ) return 1; else return 0;
}
