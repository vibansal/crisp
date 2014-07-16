/* modified on October 26 2012 to be used for low-frequency variant calling so that complexity of algorithm is not highly dependent on poolsize */

#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include "crispcaller.h"
#include "../FET/contables.h"
#include "../FET/chisquare.h"
#include "../allelecounts.h"

int USE_BASE_QVS =1;

#include "crispprint.c" // print CRISP program options and CRISP VCF record
#include "../FET/chisquare-stratified.c" // calculate FET p-value using variant->indcounts (allele counts)

int MAXITER = 20000; int MINPOS = 0, MAXPOS =0;
double thresh1 = -3.5, MINQpv = -5, thresh3 = -2;  //three pvalue thresholds 
double SER = -1; double RELAX = 0.75; int MIN_READS = 4; double  Qmin = 23; // error rate of 0.005
//int CHISQ_PERMUTATION = 0; // if = 1, permutation test is used instead of the asymptotic chi-square test
int OUTPUTGENOTYPES = 0; 

// aug 23 2012, return 0 -> no variant, 1 for SNP, 2 for indel, 3 for multi-allelic // maintain Ti/Tv ratio
// function called in main for pooled variant calling
int CRISPcaller(REFLIST* reflist,int current,int position,READQUEUE* bq,struct BAMFILE_data* bd,struct VARIANT* variant,ALLELE* maxvec,FILE* vfile)
{
	int i=0,j=0,al=0,pv=0,sample=0,potentialvar=0;
	int qpvars =0; double qpsum[2];
	double pvalue_all[2]; double ser,tablepvaluef; double pvaluelog=0,pvalchi2log=0,chi2stats[6] = {0,0,0,0,0,0};
	int asymptotic_table;  
	int permflag = 0; // if permutations are done for calculating p-value then this flag is set to 1 
	int allele=0, allele1 = maxvec[0].al, allele2 = maxvec[1].al, allele3 = maxvec[2].al;
	// indel ambiguity bug here corrected july26 2012
	int maxambiguity=0;
	for (i=1;i<=5;i++) 
	{
		if (maxvec[i].ct < 3 ) continue; // ignore alleles with less than 3 reads 
		allele = maxvec[i].al; variant->HPlength[allele] =0;
		// if there are no indels then the allele bases may not be even set
		if (current >=0 && allele >=4 && (variant->itb[allele][0] == '+' || variant->itb[allele][0] == '-')) calculate_indel_ambiguity(variant,reflist,current,allele);
		if (variant->HPlength[allele] > maxambiguity) maxambiguity = variant->HPlength[allele];
	}
	if (maxambiguity >0 )remove_ambiguous_bases(bq,variant,maxambiguity,1);

	//need fast filtering function to evaluate potential alleles for being true variant 
	for (al=1;al<=5;al++)
	{
		if (maxvec[al].ct < 3 ) continue; // ignore alleles with less than 3 reads 
		allele1=maxvec[0].al; allele2 = maxvec[al].al;   allele3 = maxvec[al ==1 ? al+1: 1].al;
		potentialvar=0;
		for (sample=0;sample<variant->samples;sample++)
        	{
	                if (variant->indcounts[sample][allele2] + variant->indcounts[sample][allele2+maxalleles]+variant->indcounts[sample][allele2+2*maxalleles] >=3) potentialvar++;
	                else if (variant->indcounts[sample][allele2] + variant->indcounts[sample][allele2+maxalleles]+variant->indcounts[sample][allele2+2*maxalleles] >=2 && variant->ploidy[sample] ==2) potentialvar++;
		}
		if (potentialvar==0) continue;

		// compute contingency table for alleles 1 and 2 
		pvaluelog =0;
		if (al > 1 && allele2 >= 4 && maxvec[1].al < 4 && maxvec[0].al < 4) asymptotic_table = populate_contable(variant,allele1,allele2,variant->refbase,maxvec[1].al); // bug fixed jan 5 2012
		else asymptotic_table = populate_contable(variant,allele1,allele2,variant->refbase,-1);
                if (PIVOTSAMPLE ==0) chi2pvalue_stratified(variant,allele1,allele2,allele3,&pvaluelog,chi2stats);
	        else chi2pvalue_stratified_pivotsample(variant,allele1,allele2,allele3,&pvaluelog,chi2stats);
		// chi2pvalue_stratified(variant,allele1,allele2,allele3,&pvaluelog,chi2stats);
		
		//permflag = 0; pvaluelog = 0;
		//permflag = pvalue_contable_iter_stranded(variant->ctable_stranded,variant->ctable_weighted,variant->strata,variant->samples,MAXITER,variant->newtable,&tablepvaluef,&pvalchi2log,chi2stats);
		//if (permflag ==0) pvaluelog = pvalchi2log; else pvaluelog = tablepvaluef;
		//if (PFLAG >=1)fprintf(stdout,"%s:%d:%0.2f ",variant->itb[allele2],permflag,pvaluelog);
		//fet(variant->counts[allele2],variant->counts[allele2]+variant->counts[allele1],variant->counts[allele2+maxalleles],variant->counts[allele2+maxalleles]+variant->counts[allele1+maxalleles]);
		variant->ctpval[variant->varalleles] = pvaluelog; variant->chi2pval[variant->varalleles] = pvaluelog;
		for (i=0;i<5;i++) variant->chi2stats[variant->varalleles][i] = chi2stats[i]; 
		
		// chernoff bound should be computed for alternate allele, FIX ERROR HERE dec14/2010
		if (allele1 >= 4 || allele2 >=4) ser = 0.05; else ser = SER;
		qpvars =0; qpsum[0] = qpsum[1] =0; pvalue_all[0] = 0; pvalue_all[1] = 0;
		if (allele1 < 4 && allele2 < 4)
		{
			for (i=0;i<variant->samples;i++)
			{
				chernoff_bound(variant->indcounts[i],variant->stats[i],allele1,allele2,ser,&variant->qpvalue[i][0],&variant->qpvalue[i][1],maxalleles);
				if (variant->qpvalue[i][0] <= MINQpv && variant->qpvalue[i][1] <= MINQpv)
				{
					//fprintf(stdout,"qpstats %d %d;%d %d:%d %f %f \n",i,variant->indcounts[i][allele1],variant->indcounts[i][allele1+maxalleles],variant->indcounts[i][allele2],variant->indcounts[i][allele2+maxalleles],variant->qpvalue[i][0],variant->qpvalue[i][1]);
					qpvars +=1; 
					qpsum[0] += variant->qpvalue[i][0]; qpsum[1] += variant->qpvalue[i][1];
				}
			}
			chernoff_bound(variant->counts,variant->tstats,allele1,allele2,pow(0.1,0.1*Qmin),&pvalue_all[0],&pvalue_all[1],maxalleles);
		}
		variant->qvpvaluef[variant->varalleles] = pvalue_all[0]; variant->qvpvaluer[variant->varalleles] = pvalue_all[1];
		/* chernoff bound calculation */		
		//fprintf(stderr,"stats %f %f \n",pvalue_all[0],pvalue_all[1]);

		pv =0;  pv= pooled_variantcaller(bq,bd,variant,allele1,allele2,qpvars,pvalue_all,qpsum,al,asymptotic_table);
		variant->paf[variant->varalleles] = 0;
	}
	if (PFLAG >=1) fprintf(stdout,"\n");
	// print VCF variant for pooled sequencing
	if (variant->varalleles >= 1) 
	{
		print_pooledvariant(variant,vfile,reflist,current,maxvec);
		return 1;
	}
	return 0;
}

// evaluate each pool for the variant, do binomial test, stranded test and unique reads test 
// call_genotypes is actually a misnomer for this function, change it at some point 
int call_genotypes(READQUEUE* bq,struct  BAMFILE_data* bd,struct VARIANT* variant,int allele1,int allele2)
{
	//if (variant->ploidy >2) MIN_READS = 4;
	double ctpval = variant->ctpval[variant->varalleles];
	double ratio = 0; double p=0;
	int variants =0; int a=0; int r=0,i=0;
	int uniquess =0; int middle = 0;
	double lls;
	//double lls[variant->ploidy+1]; 
	double best=-10000,secondbest=-10000; int bestcount =0; int sbcount =0; int pass;
	double sum =0;
	double pvallow =0,pvallowf=0,pvallowr=0,pvallowf1,pvallowr1,minpvlog,minpvlog1;
	//double pvbeta=0; double minfraction = 1.0/variant->ploidy; 
	double Ef =0.001, Er=0.001;
	variant->strandpass = 0;
	double pvstrand[variant->samples]; int spools =0;
	int C1,N1,C2,N2,C3,N3;
	//variant->hetpass = hetpass; variant->u2pass = u2pass; variant->u3pass = u3pass;

	for (a=0;a<variant->samples;a++)
	{
		variant->genotypes[a] = 0; pvstrand[a] = 0;
		if (variant->ploidy[a] > 2) variant->poolpv[a][variant->varalleles] = 0;// this variable is only malloc in pooled CRISP code not picall code, dec 20, 2011
		// use binary search or ML estimation to find allele count estimate, 07/02/2013
		best = secondbest = -10000; bestcount = sbcount =0; sum=0;
		for(i=0;i<=variant->ploidy[a] && i <=200;i++) // added i < 200 filter for cancer genomes analysis where ploidy = 10,000
		{
			ratio = (double)i/variant->ploidy[a]; p = ratio*(1-Ef) + (1-ratio)*Ef; lls =0;
			lls += log10(1-p)*(variant->ctable_stranded[a][0]-variant->ctable_stranded[a][1]) + log10(p)*variant->ctable_stranded[a][1];
			lls += log10(1-p)*(variant->ctable_stranded[a][2]-variant->ctable_stranded[a][3]) + log10(p)*variant->ctable_stranded[a][3];
			lls += log10(1-p)*(variant->ctable_stranded[a][4]-variant->ctable_stranded[a][5]) + log10(p)*variant->ctable_stranded[a][5];
			if (lls > best && best == -10000) { best = lls; bestcount = i; }
			else if (lls > best && best != -10000) { secondbest = best; sbcount = bestcount; best = lls; bestcount = i; }
			else if (lls > secondbest) { secondbest = lls; sbcount = i; }
		}
		//for(i=0;i<variant->ploidy+1;i++) sum += pow(10,lls[i]-best);
		//sum = log10(sum) + best; best -= sum; secondbest -= sum;
		variant->ctable[a][2] = bestcount; variant->ctable[a][3] = (int)(10*best-10*secondbest); variant->ctable[a][4] = sbcount;
		//for(i=0;i<variant->ploidy+1;i++) lls[i] -= sum;
		uniquess = 0; middle = 0; 
		if (variant->indcounts[a][allele2] + variant->indcounts[a][allele2+maxalleles] + variant->indcounts[a][allele2+2*maxalleles] < 2) continue;
		calculate_uniquereads(bq,bd,variant,allele1,allele2,a,&uniquess,&middle);

		//if (variant->ploidy ==1) RELAX = 0.95;  // added nov 30 2011 for ION TORRENT haploid data, default is 0.8 instead of 0.75
		//pvbeta = minreads_pvalue_dirichlet(variant->ctable[a][0],variant->ctable[a][1],variant->ploidy,0.1);
		pvallow = minreads_pvalue(variant->ctable[a][0],variant->ctable[a][1],1.0/variant->ploidy[a]);
		
		C1= variant->indcounts[a][allele2]; N1 = C1 + variant->indcounts[a][allele1];
		C2= variant->indcounts[a][allele2+maxalleles]; 	N2 = C2 + variant->indcounts[a][allele1+maxalleles];
		C3 = variant->indcounts[a][allele2+2*maxalleles];  N3 = C3 + variant->indcounts[a][allele1+2*maxalleles];
		pvstrand[a] = fet(C1,N1,C2,N2);
		//C1= variant->indcounts[a][allele2]; N1 = C1 + variant->indcounts[a][allele1];
		//C2= variant->indcounts[a][allele2+maxalleles]; N2 = C2 + variant->indcounts[a][allele1+maxalleles];

		//pvstrand[a] =  fet(variant->ctable_stranded[a][1],variant->ctable_stranded[a][0],variant->ctable_stranded[a][3],variant->ctable_stranded[a][2]); 
		// we want the allele counts on each strand to be not too different // we also want the allele counts to be at least 1/variant->ploidy[a] 
		pvallowf = minreads_pvalue(variant->ctable_stranded[a][0],variant->ctable_stranded[a][1],1.0/variant->ploidy[a]);
		pvallowr = minreads_pvalue(variant->ctable_stranded[a][2],variant->ctable_stranded[a][3],1.0/variant->ploidy[a]);
		minpvlog = pvallowf + pvallowr;

		pvallowf1 = minreads_pvalue(variant->ctable_stranded[a][0],variant->ctable_stranded[a][1],RELAX/variant->ploidy[a]);
		pvallowr1 = minreads_pvalue(variant->ctable_stranded[a][2],variant->ctable_stranded[a][3],RELAX/variant->ploidy[a]);
		minpvlog1 = pvallowf1 + pvallowr1;

		// HACK2 for indels, pooled vs individual should be different
		if (variant->ploidy[a] > 2) 
		{
			pass =0;
			if ( minpvlog >= thresh3 || (minpvlog1 >= thresh3 && ctpval-minpvlog <= thresh1 && variant->ploidy[a] >= 10) || (allele2 >=4 && minpvlog1 >= thresh3 )) pass++;
			else if (minpvlog < thresh3 && (pvallow + pvstrand[a] >= -3)) pass++;
			if  ( (variant->ctable_stranded[a][1] > 0 && variant->ctable_stranded[a][3]> 0 && uniquess >= 3) || uniquess >= 4 ) pass++;
			if (variant->ctable[a][1] >= MIN_READS) pass++;
			if (middle >0) pass++;
			if (pass ==4) { variants +=1; variant->genotypes[a] = 1; }  // added sept30 2011
			variant->poolpv[a][variant->varalleles] = pow(10,minpvlog); 
			//variant->poolsv[a][variant->varalleles] = pvstrand; 
		}
		else
		{
			if (minpvlog >= thresh3 && variant->ctable_stranded[a][1] > 0 && variant->ctable_stranded[a][3] > 0 && uniquess >= 2) variants +=1;
			else if (minpvlog >= thresh3 && variant->ctable[a][1] >= MIN_READS && uniquess >= 2) variants +=1;
			else if (minpvlog1 >= thresh3 && allele2 >= 4 && variant->ctable[a][1] >= MIN_READS && uniquess >= 2) variants +=1;
			else if (minpvlog1 >= thresh3 && allele2 < 4 && variant->ctable[a][1] >= 5 && uniquess >= 3) variants +=1;		
		}
		if (minpvlog1 >= thresh3) 	
		{
			if (spools ==0 && PFLAG >=1) fprintf(stdout,"\npoolstats %s:%d:%s:%s |",variant->chrom,variant->position,variant->itb[allele1],variant->itb[allele2]);
			//printf("S:%d (%d,%d) (%d,%d) (%d,%d) | ",a,variant->ctable_stranded[a][0],variant->ctable_stranded[a][1],variant->ctable_stranded[a][2],variant->ctable_stranded[a][3],uniquess,middle);
			if (PFLAG >=1) fprintf(stdout," S:%d %d:%d %d:%d %d:%d (%d,%d) %1.2f %1.2f %1.2f ",a,N1,C1,N2,C2,N3,C3,uniquess,middle,pvallow,pvstrand[a],minpvlog);
			if (PFLAG >=1) fprintf(stdout, "|");
			spools++; // variant pools for strand pvalue threshold
		}
		else pvstrand[a] = 0; 
	}
	for (a=0;a<variant->samples;a++)
	{
		if (pvstrand[a] < (log10(0.01) -log10(1.0+spools)) ) 
		{ 
			variant->strandpass++; if (PFLAG >=1) fprintf(stdout,"strandbias:%d:%0.2f:%d ",a,pvstrand[a],spools); 
		}
	}
	//if (variants ==1 && ctpval >= -2) variants =0;  // HACK1 if variants = 1 and ctpval is high, ignore 
	//printf("PAF: %0.3f ",variant->paf[variant->varalleles]);
	//free(lls);

	return variants;
}

void print_bincounts(struct VARIANT* variant,int allele1,int allele2)
{
	int i=0,j=0;
	for (i=0;i<variant->samples;i++)
	{
		if (variant->indcounts[i][allele2]+variant->indcounts[i][allele2+maxalleles]+variant->indcounts[i][allele2+2*maxalleles]<3) continue;
		fprintf(stdout," pool %d ",i);
		for (j=0;j<6;j++) fprintf(stdout,"%d,%d | ",variant->indcounts_binned[i][allele1][j],variant->indcounts_binned[i][allele2][j]);

		fprintf(stdout,"- ");
		for (j=0;j<6;j++) fprintf(stdout,"%d,%d | ",variant->indcounts_binned[i][allele1+maxalleles][j],variant->indcounts_binned[i][allele2+maxalleles][j]);
		fprintf(stdout,"\n");
	}
}

// impose condition that the pool with qpvars ==1 should also pass other filters, nov 17 2012 !!

// heuristic method to call variant using chi-square p-value and chernoff bound p-values for each pool and overall 
int pooled_variantcaller(READQUEUE* bq,struct BAMFILE_data* bd,struct VARIANT* variant, int allele1, int allele2,int qpvars,double* pvalue_all,double* qpsum,int al,int asymptotic_table)
{
	// strict criteria for strandedness, middle of read filter .... at least one read...
	// if #unique start sites too low as compared to non-ref reads: red flag..... close to indel, PCR duplicate..

	int low_coverage_strand = 1;
	int variants=0; double pvaluelog = variant->ctpval[variant->varalleles];  
	int i=0,j=0;
	int C1= variant->counts[allele2]; int N1 = C1 + variant->counts[allele1];
	int C2= variant->counts[allele2+maxalleles]; int N2 = C2 + variant->counts[allele1+maxalleles];
	int C3= variant->counts[allele2+2*maxalleles]; int N3 = C3 + variant->counts[allele1+2*maxalleles];
	double freq[3]; freq[0] = (double)C1/(N1+1); freq[1] = (double)C2/(N2+1); freq[2] =(double)C3/(N3+1);
	double strandpvallog = fet(C1,N1,C2,N2);
	//fprintf(stdout,"C3 N3 %d %d %f\n",C3,N3,freq[2]);

	if (variant->counts[allele2] + variant->counts[allele1] < 0.2*(variant->counts[allele2+maxalleles]+variant->counts[allele1+maxalleles])) low_coverage_strand = 1;
	else if (variant->counts[allele2+maxalleles] + variant->counts[allele1+maxalleles] < 0.2*(variant->counts[allele2]+variant->counts[allele1])) low_coverage_strand = 1;
	//fprintf(stdout,"lc %d qpvars %d pvaluelog %f \n",low_coverage_strand,qpvars,pvaluelog);
	// nov 17 2012, filter below is highly heurisic, qpvars >=2
	int variant_type =0;
	if (variant->ctpval[variant->varalleles] <= thresh1) variant_type =1; 
	else if (variant->ctpval[variant->varalleles] <= -3 && pvalue_all[0] <= -5 && pvalue_all[1] <= -5 ) variant_type =2; 
	else if (variant->ctpval[variant->varalleles] <= -2 && pvalue_all[0] + pvalue_all[1] <= -10 && strandpvallog >= -1) variant_type =3; 
	else if (low_coverage_strand ==1 && strandpvallog >= -1 && pvalue_all[0] + pvalue_all[1] <= -20 && variant->counts[allele2] >= 5 && variant->counts[allele2+maxalleles] >=5) variant_type =4;
	else if (qpvars >= 2) variant_type=5;

	if (variant_type > 0)
	{
		//printf("QVp %+2.1f %+2.1f %2.1f %2.1f strandP %+2.1f | ",qpsum[0],qpsum[1],pvalue_all[0],pvalue_all[1],strandpvallog);
		if (PFLAG >=1) fprintf(stdout,"QV %.1f:%.1f:%.1f:%.1f SD %+0.1f:%0.4f:%0.4f:%0.4f ",qpsum[0],qpsum[1],pvalue_all[0],pvalue_all[1],strandpvallog,freq[0],freq[1],freq[2]);
		//fprintf(stdout,"SPV:%.3f:%.3f ",SPVf,SPVr);
		// this is the only place where the poolsize affects complexity of method 
		variants = call_genotypes(bq,bd,variant,allele1,allele2);
		if (PFLAG >=1) 
		{
			if (allele1 >= 4 || (allele2 >=4 && allele1 == variant->refbase))  fprintf(stdout," indel variants-%d %d\t",variant_type,variants);  
			else fprintf(stdout," snp variants-%d %d\n",variant_type,variants); 
		}
		if (variants > 0)
		{
			if (PFLAG>=20) print_bincounts(variant,allele1,allele2);
			if (allele1 == variant->refbase || (allele1 != variant->refbase && allele2 !=variant->refbase)) variant->alleles[variant->varalleles] = allele2;  
			else variant->alleles[variant->varalleles] = allele1; variant->nvariants[variant->varalleles] = variants;
			// added on nov 13 2012 for new crisp output
			variant->crispvar[variant->varalleles].variantpools = variants; variant->crispvar[variant->varalleles].allelecount = 0; 
			variant->varalleles++; 
			return 1;
		}
	}
	else if (PFLAG >=2)
	{
		fprintf(stdout,"QV %.1f:%.1f:%.1f:%.1f\t",qpsum[0],qpsum[1],pvalue_all[0],pvalue_all[1]);
	}
	return 0;
}



/*
for (i=0;i<variant->samples && i==0;i++)
{
	fprintf(stdout,"+ ");
	for (j=0;j<6;j++) fprintf(stdout,"%d,%d | ",variant->indcounts_binned[i][allele1][j],variant->indcounts_binned[i][allele2][j]);

	fprintf(stdout,"- ");
	for (j=0;j<6;j++) fprintf(stdout,"%d,%d | ",variant->indcounts_binned[i][allele1+maxalleles][j],variant->indcounts_binned[i][allele2+maxalleles][j]);
	fprintf(stdout," pool %d\n",i);
}
*/
// chr12_122445542 non 0 chr12 122445542 T T A 45 1819 24 88 G 0 | CT -1.2 chisq-1 -1.166328 QV -47.7 -51.4 -16.9 -34.0 SD -15.8 | S:0 (345,9) (17,6) (3,0) -0.01 -4.79 S:1 (385,12) (15,3) (2,0) -0.02 -1.83 S:2 (381,10) (19,6) (4,0) -0.01 -4.53 S:3 (287,8) (22,8) (4,0) -0.00 -5.81 S:4 (421,6) (15,1) (2,0) -0.83 -0.66 PAF: 0.000 snp variants 2
//if (pvaluelog  <= thresh1 || (qpvars >= 2) || (pvalue_all[0] <= -5 && pvalue_all[1] <= -5 && pvaluelog <= -3) || (pvaluelog <= -2 && pvalue_all[0] + pvalue_all[1] <= -10 && strandpvallog >= -1 ) || (low_coverage_strand ==1 && strandpvallog >= -1 && pvalue_all[0] + pvalue_all[1] <= -20 && variant->counts[allele2] >= 5 && variant->counts[allele2+maxalleles] >=5) )
