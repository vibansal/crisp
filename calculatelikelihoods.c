#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>

// macro function for updating genotype likelihoods for each allele
void update_glls(double** GENLL,int i,int allele,double ep,double ep1,int allele1,int allele2,int allele3)
{
	if (allele == allele1) 
	{
		GENLL[i][0] += ep1; GENLL[i][1] += log10(0.5); GENLL[i][3] += log10(0.5); 
		GENLL[i][2] += ep;  GENLL[i][4] += ep;  GENLL[i][5] += ep;
	}
	else if (allele == allele2) 
	{
		GENLL[i][2] += ep1; GENLL[i][1] += log10(0.5); GENLL[i][4] += log10(0.5); 
		GENLL[i][0] += ep;  GENLL[i][3] += ep;  GENLL[i][5] += ep;
	}
	else if (allele == allele3) 
	{
		GENLL[i][5] += ep1; GENLL[i][3] += log10(0.5); GENLL[i][4] += log10(0.5); 
		GENLL[i][0] += ep;  GENLL[i][1] += ep;  GENLL[i][2] += ep;
	}
}

// calculate genotype likelihoods using three alleles
void compute_GLLs(READQUEUE* bq,struct BAMFILE_data* bamfiles_data,struct VARIANT* variant,int allele1,int allele2,int allele3)
{
	// 0 -> A1,A1   1 = A1,A2  2 = A2,A2  3 = A1,A3  4 = A2,A3  5 = A3,A3 
	int b=0,i=0,j=0,c=0; double qv;
	int allele; double ep,ep1;
	struct alignedread* bcall;
	for (i=0;i<variant->samples;i++) 
	{
		for (j=0;j<10;j++) { variant->GENLL[i][j] = 0.0;  variant->GENLLf[i][j] = 0.0; }
	}

	for (b=0;b<variant->samples;b++) 
	{ 
		//if (variant->samples == variant->bamfiles) i = b; else i = variant->options->BAM_TO_SAMPLE[b];
		i = b; 
		for (bcall=bamfiles_data[b].first; bcall != NULL && bcall->position <= variant->position;bcall=bcall->nextread)
		{
			if (bcall->lastpos <= variant->position || bcall->filter == '1' || bcall->type !=0) continue;
			allele = BTI[bcall->sequence[bcall->l1+bcall->delta]]; 
			qv = (double)bcall->quality[bcall->l1+bcall->delta]-QVoffset; 
			if (bcall->mquality < (int)qv) qv = bcall->mquality; // take min of base quality and mapping quality 
			if (qv > 30 && SOLID ==1) qv = 30; // cap the maximum base quality to 30 for SOLID DATA
		
			//if (qv < 1) 	fprintf(stderr,"ERROR qv %d %d %s %s %d %f\n",bcall->l1+bcall->delta,bcall->readlength,bcall->quality,bcall->readid,bcall->sampleid,qv);
			qv /=10; ep1 = log10(1.0-pow(0.1,qv)); ep = -1*qv;
			
			if (bcall->bidir ==1) update_glls(variant->GENLLb,i,allele,ep,ep1,allele1,allele2,allele3);
			else if (bcall->strand == 0) update_glls(variant->GENLLf,i,allele,ep,ep1,allele1,allele2,allele3);
			else update_glls(variant->GENLLr,i,allele,ep,ep1,allele1,allele2,allele3);
		}
		ep = 0.01;
		if (allele2 >=4) // indel allele so we use the counts of indelallele to update genotype-likelihood 
		{
			c = variant->indcounts[i][allele2+maxalleles];
			variant->GENLLr[i][2] += c*log10(1.0-ep);
			variant->GENLLr[i][1] += c*log10(0.5); variant->GENLL[i][4] += c*log10(0.5); 
			variant->GENLLr[i][0] += c*log10(ep);  variant->GENLL[i][3] += c*log10(ep);  variant->GENLL[i][5] += c*log10(ep);
			c = variant->indcounts[i][allele2];
			variant->GENLLf[i][2] += c*log10(1.0-ep);
			variant->GENLLf[i][1] += c*log10(0.5); variant->GENLLf[i][4] += c*log10(0.5); 
			variant->GENLLf[i][0] += c*log10(ep);  variant->GENLLf[i][3] += c*log10(ep);  variant->GENLLf[i][5] += c*log10(ep);
			c = variant->indcounts[i][allele2+2*maxalleles];
			variant->GENLLb[i][2] += c*log10(1.0-ep);
			variant->GENLLb[i][1] += c*log10(0.5); variant->GENLLf[i][4] += c*log10(0.5); 
			variant->GENLLb[i][0] += c*log10(ep);  variant->GENLLf[i][3] += c*log10(ep);  variant->GENLLf[i][5] += c*log10(ep);
		}
		if (allele3 >=4)
		{
			c = variant->indcounts[i][allele3+maxalleles];
			variant->GENLLr[i][5] += c*log10(1.0-ep);
			variant->GENLLr[i][3] += c*log10(0.5); variant->GENLL[i][4] += c*log10(0.5);
			variant->GENLLr[i][0] += c*log10(ep);  variant->GENLL[i][1] += c*log10(ep);  variant->GENLL[i][2] += c*log10(ep);
			c = variant->indcounts[i][allele3];
			variant->GENLLf[i][5] += c*log10(1.0-ep);
			variant->GENLLf[i][3] += c*log10(0.5); variant->GENLLf[i][4] += c*log10(0.5);
			variant->GENLLf[i][0] += c*log10(ep);  variant->GENLLf[i][1] += c*log10(ep);  variant->GENLLf[i][2] += c*log10(ep);
			c = variant->indcounts[i][allele3+2*maxalleles];
			variant->GENLLb[i][5] += c*log10(1.0-ep);
			variant->GENLLb[i][3] += c*log10(0.5); variant->GENLLf[i][4] += c*log10(0.5);
			variant->GENLLb[i][0] += c*log10(ep);  variant->GENLLf[i][1] += c*log10(ep);  variant->GENLLf[i][2] += c*log10(ep);
		}
	}

	for (i=0;i<variant->samples;i++) 
	{
		for (j=0;j<10;j++) { variant->GENLL[i][j] = variant->GENLLf[i][j]+ variant->GENLLr[i][j]+variant->GENLLb[i][j]; }
	}
}

