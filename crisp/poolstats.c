
// function that uses non-genotype based statistics (p-value for strand bias, p-value for min allele fraction, unique start site reads, middle of read counts, etc to filter out false variants over and above chi-square statistic) 
// this was the main filter for the original CRISP prior to the genotype likelihood calculation approach 
//
int calculate_pool_statistics(READQUEUE* bq,struct  BAMFILE_data* bd,struct VARIANT* variant,int allele1,int allele2)
{
	double p=0; double pvstrand[variant->samples]; double pvallow[4] = {0.0,0.0,0.0,0.0};
	int spools =0; int variantpools =0; int a=0; int r=0,i=0; int varpool =0; int uniquess =0; int middle = 0;
	int C1,N1,C2,N2,C3,N3; // variables to store counts
	variant->strandpass = 0;
	//variant->hetpass = hetpass; variant->u2pass = u2pass; variant->u3pass = u3pass;

	for (a=0;a<variant->samples;a++)
	{
		pvstrand[a] = 0.0;
		if (variant->indcounts[a][allele2] + variant->indcounts[a][allele2+maxalleles]+variant->indcounts[a][allele2+2*maxalleles] < 2) continue;
		uniquess = 0; middle = 0; 
		calculate_uniquereads(bq,bd,variant,allele1,allele2,a,&uniquess,&middle);

		C1= variant->indcounts[a][allele2]; N1 = C1 + variant->indcounts[a][allele1];
                C2= variant->indcounts[a][allele2+maxalleles];  N2 = C2 + variant->indcounts[a][allele1+maxalleles];
                C3 = variant->indcounts[a][allele2+2*maxalleles];  N3 = C3 + variant->indcounts[a][allele1+2*maxalleles];

		pvallow[1] = minreads_pvalue(variant->indcounts[a][allele1]+variant->indcounts[a][allele2],variant->indcounts[a][allele2],1.0/variant->ploidy[a]);
		pvallow[2] = minreads_pvalue(variant->indcounts[a][allele1+maxalleles]+variant->indcounts[a][allele2+maxalleles],variant->indcounts[a][allele2+maxalleles],1.0/variant->ploidy[a]);
		pvallow[3] = minreads_pvalue(variant->indcounts[a][allele1+2*maxalleles]+variant->indcounts[a][allele2+2*maxalleles],variant->indcounts[a][allele2+2*maxalleles],1.0/variant->ploidy[a]);
		pvallow[0] =pvallow[1]+pvallow[2]+pvallow[3];

		if (spools ==0 && PFLAG >=1) 
		{
			fprintf(stdout,"\npoolstats %s:%d:%s:%s |",variant->chrom,variant->position,variant->itb[allele1],variant->itb[allele2]);
			spools++;
		}

		// impose filters on minimum number of reads and variant read position...
		if (variant->indcounts[a][allele2] + variant->indcounts[a][allele2+maxalleles] + variant->indcounts[a][allele2+2*maxalleles] >= MIN_READS && middle > 0 && uniquess >= MIN_READS-1)
varpool = 1; 
		else if (variant->indcounts[a][allele2] + variant->indcounts[a][allele2+maxalleles] + variant->indcounts[a][allele2+2*maxalleles] >= MIN_READS-1 && variant->indcounts[a][allele2]+variant->indcounts[a][allele2+2*maxalleles] > 0 && variant->indcounts[a][allele2+maxalleles]+variant->indcounts[a][allele2+2*maxalleles] > 0 && middle >0 && uniquess >= MIN_READS-1) varpool = 1; 
		else varpool = 0;
		variantpools += varpool;

		if ((pvallow[0] >= -3 || varpool ==1) && spools <10)
		{

			if (PFLAG >=1) fprintf(stdout," S:%d:%d %d:%d %d:%d %d:%d (%d,%d) %1.2f %1.2f |",a,varpool,N1,C1,N2,C2,N3,C3,uniquess,middle,pvallow[0],pvstrand[a]);
			spools++; // variant pools for strand pvalue threshold
		}
	}
	if (PFLAG >=1) fprintf(stdout,"varpools %d |",variantpools);
	/*
	for (a=0;a<variant->samples;a++)
	{
		if (pvstrand[a] < (log10(0.01) -log10(1.0+spools)) ) 
		{ 
			variant->strandpass++; 
			if (PFLAG >=1) fprintf(stdout,"strandbias:%d:%0.2f:%d ",a,pvstrand[a],spools); 
		}
	*/
	return variantpools;
}


