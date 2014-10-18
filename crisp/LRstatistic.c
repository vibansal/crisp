
/* code for calculating likelihood ratio statistic for association using pooled data, code added october 22 2013 */

// only calculate likelihood for a subset of pools (case/control), E01, E10 are already determined for indels...
// much cleaner code than full likelihood EMmethod 
double calculate_likelihood_EM(struct VARIANT* variant,int allele1, int allele2, int allele3,double p,double* dataLikelihood,double E01[],double E10[],int pool_flag)
{
	int i=0,j=0,k=0,K=0,iter=0; 
	double pnew=0;  double LLtotal =0,LLtotalprev=0; double temp=0; double Lj=0,lp0=0,lp1=0;
	int maxiter = 1000; double convergence_delta = 0.01; int exitloop =0;

	for (i=0;i<variant->samples;i++) 
	{
		if ( (variant->options->phenotypes[i] & pool_flag) ==0 ) continue;
		K += variant->ploidy[i];
	}
	//fprintf(stdout,"K %d %d\n",K,pool_flag);

	while (iter++ < maxiter && exitloop ==0)
	{
		pnew =0; LLtotal = 0; lp0 = log10(p); lp1 = log10(1.0-p); 
		for (i=0;i<variant->samples;i++)
		{
			if ( (variant->options->phenotypes[i] & pool_flag) ==0 ) continue;
			for (j=0;j<=variant->ploidy[i];j++) 
			{
				temp = variant->GENLL[i][j] + variant->NCarray[i][j] + lp0*j + lp1*(variant->ploidy[i]-j); 
				if (j ==0) Lj =temp; else if (temp > Lj) Lj = temp + log10(1.0+pow(10,Lj-temp)); else Lj += log10(1.0+pow(10,temp-Lj));
			}
			LLtotal += Lj; 
			for (j=0;j<=variant->ploidy[i];j++) 
			{
				temp = variant->GENLL[i][j] + variant->NCarray[i][j] + lp0*j + lp1*(variant->ploidy[i]-j); 
				pnew += (double)j*pow(10,temp-Lj); 
			}
		}
		pnew /= K; 
		if (allele2 >=4 || USE_BASE_QVS==0) calculate_pooled_likelihoods_indel(variant,allele1,allele2,allele3,E01,E10); // for indels or if we don't use base quality values
		if (pnew < MINP) pnew= MINP; if (1.0-pnew <=MINP) pnew = 1.0-MINP;
		if (iter >=2 && fabsf(LLtotal-LLtotalprev) <= convergence_delta) exitloop = 1;
		p = pnew; LLtotalprev = LLtotal; 	
	}
	*dataLikelihood = LLtotal; 	
	return p;
}

// encoding is 1:control, 2:case, 4:earlyonset 0:missing
double calculate_LRstatistic(struct VARIANT* variant,int allele1, int allele2, int allele3,double p,double E01[],double E10[])
{
        double dataLL,dataLL_cases,dataLL_controls,dataLL_cases_early,dataLL1;
	double pcombined = calculate_likelihood_EM(variant,allele1,allele2,allele3,p,&dataLL,E01,E10,7);
	double pcontrols = calculate_likelihood_EM(variant,allele1,allele2,allele3,p,&dataLL_controls,E01,E10,1);
	double pcases = calculate_likelihood_EM(variant,allele1,allele2,allele3,p,&dataLL_cases,E01,E10,6);
	double OR = pcases*(1.0-pcontrols)/((1.0-pcases)*pcontrols);
	double LLR = (dataLL_cases+dataLL_controls-dataLL);
	double pvalue = kf_gammaq(0.5,LLR*log(10))/log(10);

	// only do this when the number of pools with phenotype = 4 is sufficiently large
	double pcases1 = calculate_likelihood_EM(variant,allele1,allele2,allele3,p,&dataLL_cases_early,E01,E10,4);
	double pcombined1 = calculate_likelihood_EM(variant,allele1,allele2,allele3,p,&dataLL1,E01,E10,5);
	double OR1 = pcases1*(1.0-pcontrols)/((1.0-pcases1)*pcontrols);
	double LLR1 = (dataLL_cases_early+dataLL_controls-dataLL1);
	double pvalue1 = kf_gammaq(0.5,LLR1*log(10))/log(10);

	fprintf(stdout,"ASSOC %s %d %s %s | ",variant->chrom,variant->position+1,variant->itb[allele1],variant->itb[allele2]);
	fprintf(stdout,"%0.4f %0.4f %0.4f %0.4f LL %0.2f %0.2f %0.2f LLR %0.2f:%0.3f:%0.2f \t",pcombined,pcontrols,pcases,pcases1,dataLL,dataLL_controls,dataLL_cases,LLR,OR,pvalue);
	fprintf(stdout,"LL %0.2f %0.2f %0.2f LLR %0.2f:%0.3f:%0.2f \n",dataLL1,dataLL_controls,dataLL_cases_early,LLR1,OR1,pvalue1);
	return pvalue;
}

