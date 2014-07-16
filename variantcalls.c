/////////////////////////////////////////////////////////////////////////////////////////////////////

int positions_processed =0;

// new function to update base-base error tables using non-variant sites... previous base is considered as covariate
int update_error_table(REFLIST* reflist,int current,int position,READQUEUE* bq,struct BAMFILE_data* bd,struct VARIANT* variant)
{
        int i=0,j=0,k=0;
        int index_fwd =-1,index_rev =-1; double s=0;
        if (position > 1 && reflist->sequences[current][position-1] != 'N' &&  reflist->sequences[current][position-1] != 'n' && reflist->sequences[current][position-2] != 'N')
        {
                index_fwd = BTI[reflist->sequences[current][position-2]] + 4*BTI[reflist->sequences[current][position-1]] + 16*BTI[reflist->sequences[current][position]];
        }
        if (position < reflist->lengths[current]-2 && reflist->sequences[current][position+1] != 'N' && reflist->sequences[current][position+1] != 'n' && reflist->sequences[current][position+2] != 'N')
        {
                index_rev = 3-BTI[reflist->sequences[current][position+2]] + 4*(3-BTI[reflist->sequences[current][position+1]]) + 16*(3-BTI[reflist->sequences[current][position]]);
        }

        for (i=0;i<variant->samples;i++)
        {
                if (index_fwd < 0 || index_rev < 0) continue;
                //if (index_fwd >=64 || index_rev >= 64) continue;
                for (j=0;j<4;j++) bd[i].counts[index_fwd][j] += variant->indcounts[i][j];
                for (j=0;j<4;j++) bd[i].counts[index_rev][3-j] += variant->indcounts[i][j+maxalleles];
                //for (j=0;j<4;j++) bd[i].counts[index_rev][3-j+4] += variant->indcounts[i][j+maxalleles]; 
        }
        if (positions_processed%20000 ==0 && positions_processed > 0)
        {
                for (i=0;i<variant->samples;i++)
                {
                        for (k=0;k<64;k++)
                        {
                                s = bd[i].counts[k][0]+bd[i].counts[k][1]+bd[i].counts[k][2]+bd[i].counts[k][3];
				fprintf(stdout,"%c%c ",variant->itb[(k%16)%4][0],variant->itb[(k%16)/4][0]);
                                fprintf(stdout,"%0.6f %0.6f %0.6f %0.6f %d | ",(float)bd[i].counts[k][0]/s,(float)bd[i].counts[k][1]/s,(float)bd[i].counts[k][2]/s,(float)bd[i].counts[k][3]/s,(int)s);
                                //fprintf(stdout,"%d,%d,%d,%d  %d,%d,%d,%d | ",bd[i].counts[k][0],bd[i].counts[k][1],bd[i].counts[k][2],bd[i].counts[k][3],bd[i].counts[k][4],bd[i].counts[k][5],bd[i].counts[k][6],bd[i].counts[k][7]);
                                if (k%16 ==15) fprintf(stdout,"ES:%d base %s\n ",i,variant->itb[k/16]);
                        }
                }
        }
}


int jointvariantcaller(REFLIST* reflist,int current,int position,READQUEUE* bq,struct BAMFILE_data* bamfiles_data,struct VARIANT* variant,struct OPTIONS* options,int lowcoverage)
{
	// no need of this filter since we can call variants at such positions
	if (reflist->sequences[current][position] == 'N' || reflist->sequences[current][position] == 'n') return 0;
	//printpileup(reflist,current,position,bq);
	initialize_variant(reflist,current,position,variant); // init variant counts
	if (ALLOW_AMBIGUOUS_BASES ==0 && variant->IUPAC ==1) 
	{
		fprintf(stderr,"reference base %c at position %s:%d is ambiguous, site will be ignored, provide reference fasta file with bases in [ACTGN] to call variants at these sites\n",variant->refb,variant->chrom,variant->position);
		return 0;
	}
	calculate_allelecounts(reflist,current,bq,bamfiles_data,variant);
	// function call to identify alleles at variant->position (pileup like approach) and their counts...

	ALLELE* maxvec = (ALLELE*)malloc(sizeof(ALLELE)*maxalleles);
	int i=0,j=0,v=0;
	for (i=0;i<maxalleles;i++){ maxvec[i].ct = 0; maxvec[i].al = i; }

	int allele1,allele2,allele3,allele4;
	int potentialvariant =0;

	// sorting algorithm can be different for different variant calling data (lowcov | pooled | diploid)
	if (PICALL ==1) sort_allelecounts_diploid(variant,maxvec,variant->refbase);  // refbase is always first allele in sorted order 
	else sort_allelecounts(variant,maxvec,variant->refbase);  // refbase is always either 0th or 1st allele in sorted order 

	allele1 = maxvec[0].al; allele2 = maxvec[1].al; allele3 = maxvec[2].al; allele4 = maxvec[3].al;
	variant->varalleles=0;
	if (PFLAG >=1) fprintf(stdout,"%s %d %c %s %s ",variant->chrom,variant->position+1,reflist->sequences[current][position],variant->itb[allele1],variant->itb[allele2]);
	if (PFLAG >=1) fprintf(stdout,"%5d %5d %5d %5d ",variant->counts[allele1],variant->counts[allele2],variant->counts[allele1+maxalleles],variant->counts[allele2+maxalleles]);
	if (variant->counts[allele3]+variant->counts[allele3+maxalleles]+variant->counts[allele3+2*maxalleles] >=4 && PFLAG >=1)  printf("%s:%d,%d ",variant->itb[allele3],variant->counts[allele3],variant->counts[allele3+maxalleles]);
	if (variant->counts[allele4]+variant->counts[allele4+maxalleles]+variant->counts[allele3+2*maxalleles] >=10 && PFLAG >=1 && allele4 >=4)  printf("%s:%d,%d,%d ",variant->itb[allele4],variant->counts[allele4],variant->counts[allele4+maxalleles],variant->counts[allele3+2*maxalleles]);
	//check if the variant can be a potential variant, added may 8 2012
	potentialvariant =0;
	for (i=0;i<maxalleles;i++)
	{
		if (i == variant->refbase) continue;
		// aug 19 2012 temporary filter change to use CRISP for lowcov
		if (variant->counts[i]+variant->counts[i+maxalleles]+variant->counts[i+2*maxalleles] >= 3 && PICALL == 2) potentialvariant =2;  // filter for low coverage
		for (j=0;j<variant->samples;j++)
		{
			// for CRISP this is the filter for potential variant
			if (variant->indcounts[j][i]+variant->indcounts[j][i+maxalleles]+variant->indcounts[j][i+2*maxalleles] >= 3 && (PICALL ==0 || PICALL ==3)) potentialvariant+=2;
			else if (variant->indcounts[j][i]+variant->indcounts[j][i+maxalleles] + variant->indcounts[j][i+2*maxalleles] >= 1 && PICALL ==3 && variant->ploidy[j] ==2) potentialvariant+=variant->indcounts[j][i]+variant->indcounts[j][i+maxalleles]+variant->indcounts[j][i+2*maxalleles];
			else if (variant->indcounts[j][i]+variant->indcounts[j][i+maxalleles]+variant->indcounts[j][i+2*maxalleles] >= 2 && PICALL ==1) potentialvariant++;
			else if (variant->indcounts[j][i]+variant->indcounts[j][i+maxalleles]+variant->indcounts[j][i+2*maxalleles] >= 3 && PICALL ==1) potentialvariant+=2;
		}
	}

	if (potentialvariant < 2 || variant->refb == 'N' || lowcoverage ==1) 
	{ 
		if (PFLAG>=1 && lowcoverage==1) fprintf(stdout," LC");
		if (PFLAG >=1) fprintf(stdout,"\n"); 
		if (lowcoverage != 5 && variant->refb != 'N') 
		{
			update_error_table(reflist,current,position,bq,bamfiles_data,variant);
			positions_processed++; 
		}
		free(maxvec); 
		return 0;
	}
	else 
	{
		if (PFLAG >=1 && OVERLAPPING_PE_READS==1 && variant->counts[allele1+2*maxalleles]+variant->counts[allele2+2*maxalleles] > 0) fprintf(stdout," | OPE=%d:%d:%d ",variant->counts[allele1+2*maxalleles],variant->counts[allele2+2*maxalleles],variant->filteredreads[4]);
		else if (PFLAG >=1) fprintf(stdout," | ");
		if (PFLAG>=1) fprintf(stdout,"\tMQ=%d,%d,%d,%d;FL=%d,%d;",variant->MQcounts[0],variant->MQcounts[1],variant->MQcounts[2],variant->MQcounts[3],variant->filteredreads[0],variant->filteredreads[1]);
		//if (PFLAG >=1 && variant->filteredreads[2] > 0) fprintf(stdout,"OPE=%d,%d,%d,%d",variant->filteredreads[2],variant->filteredreads[3],variant->filteredreads[4],variant->filteredreads[5]);
		if (PFLAG >=1) fprintf(stdout," | ");
		
		#if PICALL ==0 // old CRISP method that outputs allele frequencies and uses contingency table and qvalue test.. per pool
			VARIANTS_CALLED += CRISPcaller(reflist,current,position,bq,bamfiles_data,variant,maxvec,options->vfile);
		#elif PICALL ==3 // this is the new CRISP method that can estimate allele counts/genotypes
			v = newCRISPcaller(reflist,current,position,bq,bamfiles_data,variant,maxvec,options->vfile);
			VARIANTS_CALLED += v; 
			if (v ==0) 
			{
				update_error_table(reflist,current,position,bq,bamfiles_data,variant);
	                        positions_processed++;
			}
		#endif
	        //print_indel_lists(bamfiles_data,variant->samples,variant->position);
	}

	//if (PFLAG >=1)fprintf(stdout,"\n");
	free(maxvec);
	return 1;
}

//Zero-based index: Start and end positions are identified using a zero-based index. The end position is excluded

//if (reflist->sequences[current][k] =='N') continue; 
int callvariants(REFLIST* reflist,int current,int first,int last,READQUEUE* bq,struct BAMFILE_data* bamfiles_data,struct OPTIONS* options,struct VARIANT* variant)
{
	int i=0, k=0,reads=0; 

	if (bq->first == NULL) return 0;  else if ((bq->first)->position > first) first = (bq->first)->position; 

	// find first interval whose ending is after the current position
	if (targeted ==1 && reflist->cinterval < 0) reflist->cinterval = reflist->first_interval_chrom[current]; 

	// filter to ignore positions before start of target //fprintf(stdout,"target info %d %d \n",targettid,targetstart);
	if (options->targettid != -1 && first < options->targetstart) first= options->targetstart;
	//fprintf(stdout,"return first last %d %d %d %d\n",first,last,bq->first->position,bq->first->lastpos); 

	for (k=first;k<last;k++)
	{
		//for(i=0;i<variant->samples;i++) bamfiles_data[i].previous=NULL;
		// this loop can be shortened // june 19 2012
		bq->iter = bq->first; bq->previous = NULL; reads=0;
		while (bq->iter != NULL && bq->iter->position <= k) 
		{
			if (bq->iter->lastpos < k)  free_read_queue(bq,bamfiles_data); 
			else
			{
				//fprintf(stdout,"%d %d-%d | ",bq->iter->sampleid,bq->iter->position,bq->iter->lastpos);
				bamfiles_data[bq->iter->sampleid].previous = bq->iter;
				bq->previous = bq->iter; bq->iter = bq->iter->next;
				reads++;
			}
		}
		// continue;

		if (options->targettid != -1 && k >= options->targetend) 
		{
			fprintf(stderr,"position %d is outside target end %d\n",k,options->targetend);
			break;
		}
		if (CALL_VARIANTS ==0) continue;  // do not call variants, for evaluation purposes only

		if (targeted ==0) jointvariantcaller(reflist,current,k,bq,bamfiles_data,variant,options,0); // whole-genome variant calling 

		//if ((targeted ==0) || (targeted ==1 && (islower(reflist->sequences[current][k]) || reads >= 1000000)))  
		else if (reflist->cinterval >=0) //extra condition for chromosomes not in bedfile bug fixed sept 9 2012
		{
			while (reflist->intervallist[reflist->cinterval].end < k) 
			{
				if (reflist->intervallist[reflist->cinterval].chrom != current || reflist->cinterval >= reflist->intervals) break;
				reflist->cinterval++; 
			}
			if (k >= reflist->intervallist[reflist->cinterval].start && k <= reflist->intervallist[reflist->cinterval].end )
			{
				//fprintf(stdout,"BED:%d-%d %d| ",reflist->intervallist[reflist->cinterval].start,reflist->intervallist[reflist->cinterval].end,reads);
				if (reads >= HAPLOTYPES*MIN_COVERAGE_POOL || variant->ploidy[0] == 2) jointvariantcaller(reflist,current,k,bq,bamfiles_data,variant,options,0);
				//else jointvariantcaller(reflist,current,k,bq,bamfiles_data,variant,options,1);
			}
			else if (reflist->intervallist[reflist->cinterval].start-k > 0 && reflist->intervallist[reflist->cinterval].start-k < FLANKING_BASES &&reads >= MIN_COVERAGE_FLANKING )
			{
				if (PFLAG >=2) fprintf(stdout,"RF:%d-%d %d| ",reflist->intervallist[reflist->cinterval].start,reflist->intervallist[reflist->cinterval].end,reads);
				jointvariantcaller(reflist,current,k,bq,bamfiles_data,variant,options,0);
			}
			else if (reflist->cinterval-1 >=0 && k-reflist->intervallist[reflist->cinterval-1].end > 0 &&  k-reflist->intervallist[reflist->cinterval-1].end <FLANKING_BASES && reads >= MIN_COVERAGE_FLANKING)
			{
				if (PFLAG >=2) fprintf(stdout,"LF:%d-%d %d| ",reflist->intervallist[reflist->cinterval-1].start,reflist->intervallist[reflist->cinterval-1].end,reads);
				jointvariantcaller(reflist,current,k,bq,bamfiles_data,variant,options,0);
			}	
		}
	}
	return 1;
}

