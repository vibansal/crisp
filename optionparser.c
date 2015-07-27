
/* code for parsing arguments from command line for CRISP and for parsing the input file with information about BAM file and pool sizes 
 */

//void init_parser(struct OPTIONS* options);
int optparser(int argc, char* argv[],struct OPTIONS* options);

int parse_arguments(int argc,char* argv[],struct OPTIONS* options,char* bamfilepaths,char* phenotypefile)
{
	int i=0,bamfiles=0;
	// output command used to run  CRISP in stdout file 
	fprintf(stdout,"##CRISP_command "); for (i=0;i<argc;i+=1)fprintf(stdout," %s",argv[i]); fprintf(stdout,"\n"); 

	OUTPUT_ALLELE_COUNTS = NULL;
	for (i=1;i<argc-1;i+=2)
	{
		if (strcmp(argv[i],"--ref") ==0 || strcmp(argv[i],"--reference") ==0)        strcpy(options->fastafile,argv[i+1]);
		else if (strcmp(argv[i],"--phenotypes") ==0)        strcpy(phenotypefile,argv[i+1]); 
		else if (strcmp(argv[i],"--paths") ==0 || strcmp(argv[i],"--bams") ==0)        
		{
			strcpy(bamfilepaths,argv[i+1]);
		}
		else if (strcmp(argv[i],"-VCF") ==0 || strcmp(argv[i],"--VCF") ==0)        strcpy(options->vcffile,argv[i+1]);
		else if (strcmp(argv[i],"-indelfile") ==0 || strcmp(argv[i],"--indels") ==0)        strcpy(options->indelfile,argv[i+1]);
		else if (strcmp(argv[i],"--bam") ==0 || strcmp(argv[i],"--bamfile") ==0)    bamfiles++;
		else if (strcmp(argv[i],"--flankingbases") ==0 || strcmp(argv[i],"--fb") ==0 || strcmp(argv[i],"--flanking") ==0)  FLANKING_BASES =atoi(argv[i+1]);
		else if (strcmp(argv[i],"--pivotsample") ==0) { PIVOTSAMPLE = atoi(argv[i+1]); fprintf(stderr,"samples sequenced on two runs: 1...%d | %d+1... \n",PIVOTSAMPLE,PIVOTSAMPLE+1); } // SOLID reads bioscope 
		else if (strcmp(argv[i],"--solid") ==0) { SOLID = atoi(argv[i+1]); fprintf(stderr,"solid reads \n"); } // SOLID reads bioscope 
		else if (strcmp(argv[i],"--regions") ==0) 
		{
			options->regions = (char*)malloc(strlen(argv[i+1])+1); 
			strcpy(options->regions,argv[i+1]); // only look in certain regions/intervals for variant calling
		}
		else if (strcmp(argv[i],"--bedfile") ==0 || strcmp(argv[i],"--bed") ==0)        strcpy(options->bedfile,argv[i+1]);
		else if (strcmp(argv[i],"--QVoffset") ==0 || (strcmp(argv[i],"--qvoffset") ==0))      { QVoffset = atoi(argv[i+1]); QVset = 1; } 
		else if (strcmp(argv[i],"--weighted") ==0)      { USE_QV= atoi(argv[i+1]); fprintf(stderr,"USE_QV set to 1,weighted allele counts\n"); } 
		else if (strcmp(argv[i],"--mbq") ==0) MINQ = atoi(argv[i+1]);
		else if (strcmp(argv[i],"--minm") ==0 || strcmp(argv[i],"--mmq") ==0) MIN_M = atoi(argv[i+1]);// bug fixed jan 31 2012
		else if (strcmp(argv[i],"--maxm") ==0) MAX_MM = atoi(argv[i+1]);
		else if (strcmp(argv[i],"--filterreads") ==0) FILTER_READS_MISMATCHES = atoi(argv[i+1]);
		else if (strcmp(argv[i],"--mincov") ==0 || strcmp(argv[i],"--MC")==0) MIN_COVERAGE_POOL = atoi(argv[i+1]);
		else if (strcmp(argv[i],"--verbose") ==0 || strcmp(argv[i],"--pflag") ==0) PFLAG = atoi(argv[i+1]);
		else if (strcmp(argv[i],"--splitvariants") ==0 || strcmp(argv[i],"--sv") ==0) SPLIT_TRIALLELIC_VARS = atoi(argv[i+1]);
		else if (strcmp(argv[i],"--indelrealignment") ==0 || strcmp(argv[i],"--realignment") ==0) INDEL_REALIGNMENT = atoi(argv[i+1]);
		else if (strcmp(argv[i],"--callvariants") ==0) 
		{
			CALL_VARIANTS = atoi(argv[i+1]);
			if (CALL_VARIANTS ==0) fprintf(stderr,"\nno variants will be called, only for testing purposes \n"); 
		}
		else if (strcmp(argv[i],"--EM") ==0 || strcmp(argv[i],"--callgenotypes") ==0) 
		{
			EMflag = atoi(argv[i+1]);
			if (EMflag ==1) fprintf(stderr,"\nEM algorithm will be used for estimating genotypes and calling variants\n \n"); 
		}
		else if (strcmp(argv[i],"--callsnps") ==0) INDELSONLY = 1-atoi(argv[i+1]);
		else if (strcmp(argv[i],"--LLRthresh") ==0) LLRthresh = atof(argv[i+1]);
		else if (strcmp(argv[i],"--leftalignindels") ==0 || strcmp(argv[i],"--leftalign")==0 || strcmp(argv[i],"--LeftAlignIndels") ==0) 
		{
			LEFT_ALIGN_INDELS = atoi(argv[i+1]);
			if (LEFT_ALIGN_INDELS ==1) fprintf(stderr,"\ninsertions and deletions in reads will be left justified w.r.t. reference sequence (limited evaluation) \n"); 
		}
		else if (strcmp(argv[i],"--allowAB") ==0) 
		{
			ALLOW_AMBIGUOUS_BASES = atoi(argv[i+1]);
			if (ALLOW_AMBIGUOUS_BASES ==1) fprintf(stderr,"\nWARNING: CRISP will call variants at sites at which the reference base is ambiguous by choosing the most common allele as the reference, this may cause problems for downstream analysis. It is better to use a reference fasta file with no ambiguous bases (all bases should be A,C,T,G or N)\n"); 
		}
		else if (strcmp(argv[i],"--OPE") ==0) 
		{
			OVERLAPPING_PE_READS= atoi(argv[i+1]);
			if (OVERLAPPING_PE_READS ==1) fprintf(stderr,"\n overlapping paired-end reads will be detected and bases read from both ends will be considered once for variant calling, This can slow down the program in some cases, use \"--OPE 0\" to disable this option.\n"); 
			else fprintf(stderr,"overlapping paired-end reads will NOT be detected, use \"--OPE 1\" option to detect overlapping paired-end reads\n");
		}
		else if (strcmp(argv[i],"--refbias") ==0) AGILENT_REF_BIAS = atof(argv[i+1]); 
		else if (strcmp(argv[i],"--useduplicates") ==0) USE_DUPLICATES = atoi(argv[i+1]); 
		else if (strcmp(argv[i],"--fastfilter") ==0) FAST_FILTER = atoi(argv[i+1]); 
		else if (strcmp(argv[i],"--sse") ==0) CALCULATE_ERROR_RATES = atoi(argv[i+1]);  // allele counts for each sample to calculate context specific error rates 
		else if (strcmp(argv[i],"--outputAC") ==0) OUTPUT_ALLELE_COUNTS = argv[i+1];  // allele counts for each sample at indel sites will be output to the file specified

		else
		{
			#if (PICALL ==0 || PICALL ==3)
				if (strcmp(argv[i],"-e") ==0) SER = atof(argv[i+1]);
				else if (strcmp(argv[i],"--OUTPUTGENOTYPES") ==0) OUTPUTGENOTYPES = atoi(argv[i+1]);
				else if (strcmp(argv[i],"--genotypes") ==0) OUTPUTGENOTYPES = atoi(argv[i+1]);
				else if (strcmp(argv[i],"--minc") ==0) MIN_READS = atoi(argv[i+1]);
				else if (strcmp(argv[i],"--perms") ==0) MAXITER = atoi(argv[i+1]);
				else if (strcmp(argv[i],"-p") ==0 || strcmp(argv[i],"--poolsize") ==0) 
				{  
					options->POOLSIZE = atoi(argv[i+1]); options->varpoolsize =0;
				}
				else if (strcmp(argv[i],"--ctpval") ==0) thresh1 = atof(argv[i+1]);
				else if (strcmp(argv[i],"--qvpval") ==0) MINQpv = atof(argv[i+1]);
				else if (strcmp(argv[i],"--usebq") ==0 || strcmp(argv[i],"--useqv") ==0 ) USE_BASE_QVS = atoi(argv[i+1]);
				else if (strcmp(argv[i],"--chisquare") == 0) 
				{
					CHISQ_PERMUTATION = 1- atoi(argv[i+1]); 
					if (CHISQ_PERMUTATION ==0) fprintf(stderr,"asymptotic chisquare test will be used to calculate p-value for contingency tables \n");
					else fprintf(stderr,"permutation test will be used for calculating p-value of contingency tables \n");
				}
				else fprintf(stderr,"unrecognized option.... %s %s \n",argv[i],argv[i+1]);
			#endif

			#if PICALL ==1  // not used anymore 
				if (strcmp(argv[i],"--callsnps") ==0) INDELSONLY = 1-atoi(argv[i+1]);
				else if (strcmp(argv[i],"--LLRthresh") ==0) LLRthresh = atof(argv[i+1]);
				else if (strcmp(argv[i],"--HWE") ==0) 
				{
					HWE= atoi(argv[i+1]);
					if (HWE ==0) fprintf(stderr,"HWE set to 0, genotype prior will not favour heterozygotes \n"); 
				}
				else if (strcmp(argv[i],"--PFLAG") == 0 || strcmp(argv[i],"--pflag") == 0) 
				{
					PFLAG = atoi(argv[i+1]);
					fprintf(stderr,"verbose mode: program will output information from picall likelihood calculations \n");
				}
				else if (strcmp(argv[i],"--mincov") == 0) 
				{
					MINCOV = atoi(argv[i+1]); // minimum coverage used for picall 
					fprintf(stderr,"setting MINCOV for picall model to %d \n",MINCOV);
				}
				else if (strcmp(argv[i],"--alpha") == 0) 
				{
					alpha = atof(argv[i+1]); alpha0 = alpha; alpha1 = alpha; 
					fprintf(stderr,"setting alpha prior for picall model to %f \n",alpha);
				}
				else if (strcmp(argv[i],"--beta") == 0) 
				{
					beta = atof(argv[i+1]); beta0 = beta; beta1 = beta; 
					fprintf(stderr,"setting beta prior for picall model to %f \n",beta);
				}
				else fprintf(stderr,"unrecognized option.... %s %s \n",argv[i],argv[i+1]);
			#endif
		}
	}
	return bamfiles;
}

int parse_bamfilenames(char* bamfilepaths,struct OPTIONS* options)
{
	int i=0,j=0,k=0,s=0,e=0;
        char buffer[4096];
	FILE* fp = fopen(bamfilepaths,"r"); 
	options->bamfiles =0; options->samples=0;
	if(!fp) {fprintf(stderr,"could not open file %s, please provide a valid file with paths for bam files\n\n",bamfilepaths);return 0;}
	while (fgets(buffer,4096,fp) != NULL) { options->bamfiles++;} fclose(fp);
	options->bamfilelist = (char**)malloc(sizeof(char*)*options->bamfiles); 
	options->sampleids = (char**)malloc(sizeof(char*)*options->bamfiles); 
	options->ploidy = calloc(options->bamfiles,sizeof(int)); // allocate the ploidy array 
	fp = fopen(bamfilepaths,"r");
	for (i=0;i<options->bamfiles;i++)
	{
		options->ploidy[i]= 0;
		// SID=xx  in this file PS=40 (poolsize) PL=illumina etc... CaseControl=0/1/2/
		fgets(buffer,4096,fp); 
		j=0; while (buffer[j] == ' ' || buffer[j] == '\t') j++; s = j; // empty spaces 
		while (buffer[j]  != '\t' && buffer[j] != '\n' && buffer[j] !=' ') j++; 
		options->bamfilelist[i] = calloc(sizeof(char),(j-s+1));
		for (k=s;k<j;k++) options->bamfilelist[i][k-s] = buffer[k]; options->bamfilelist[i][k-s] ='\0';

		if (buffer[j] == '\n') continue; // nothing there

		while (buffer[j] == ' ' || buffer[j] == '\t') j++; s = j; // empty spaces 
		while (buffer[j]  != '\t' && buffer[j] != '\n' && buffer[j] !=' ') j++; 
		if (j-s >3 && buffer[s] == 'P' && buffer[s+1] == 'S' && buffer[s+2] == '=') 
		{
			s +=3; options->ploidy[i] = 0;
			while (s < j) 
			{ 
				options->ploidy[i] *=10; options->ploidy[i] += buffer[s]-48; s++; 
		
			} 
		}
		//fprintf(stdout,"%d ploidy %s\n",options->ploidy[i],buffer);
		if (buffer[j] == '\n') 
		{
			options->sampleids[i] = NULL;
			continue; 
		}

		// read the sample id for the bamfile
		while (buffer[j]  != '\t' && buffer[j] != '\n' && buffer[j] !=' ') j++; e = j; 
		options->sampleids[i] = calloc(sizeof(char),e-s+1);
		for (k=s;k<e;k++) options->sampleids[i][k-s] = buffer[k]; options->sampleids[i][k-s] ='\0';
		options->samples++;
		
	}
	options->varpoolsize = 1;
	for (i=0;i<options->bamfiles;i++) 
	{
		//fprintf(stderr,"buffer %d %d \n",i,options->ploidy[i]);
		if (options->ploidy[i] ==0) 
		{
			options->varpoolsize = 0; 
		}
	}
	fclose(fp);

	if (options->samples == options->bamfiles) // there are sampleids for each bamfile
	{
		options->BAM_TO_SAMPLE = calloc(sizeof(int),options->samples); k=0;
		for (i=0;i<options->bamfiles;i++)
		{
			s = 0;
			for (j=0;j<i;j++) 
			{
				if (strcmp(options->sampleids[i],options->sampleids[j]) ==0) 
				{  
					options->BAM_TO_SAMPLE[i] = options->BAM_TO_SAMPLE[j]; s = 1;
					break; 
				} 
			}
			if (s ==0) // new sample 
			{
				options->BAM_TO_SAMPLE[i] = k; k++; 
			}
		}
		options->samples = k; // number of unique samples 
		for (i=0;i<options->bamfiles;i++) fprintf(stderr,"bam |%s| |%s| %d\n",options->bamfilelist[i],options->sampleids[i],options->BAM_TO_SAMPLE[i]);
		fprintf(stderr,"no of unique sampleids %d \n",options->samples);
	}
	return 1;
}

// main function for parsing command line arguments to CRISP 
int optparser(int argc, char* argv[],struct OPTIONS* options)
{
	int i=0;
	char bamfilepaths[1024]; strcpy(bamfilepaths,"None"); char phenotypefile[1024]; strcpy(phenotypefile,"None");
	FILE* fp; 
	int bamfiles =0; int printoptions=0; options->association = 0;
	strcpy(options->fastafile,"None"); strcpy(options->bedfile,"None"); strcpy(options->vcffile,"None");
	strcpy(options->indelfile,"None");
	options->bamfiles =0; options->READLENGTH = 0; options->varpoolsize = 0; // variable pool size 0/1
	options->regions = NULL; options->samples =0; options->bamfiles =0;
	//else if (strcmp(argv[i],"--variantcallingmodel")==0 || strcmp(argv[i],"--model") == 0) PICALL = atoi(argv[i+1]); 
	// picall = 0 = pooled, picall = 1 = picall, picall = 2 = lowcoverage
	if (argc >= 2 && (strcmp(argv[1],"-h") ==0 || strcmp(argv[1],"--help") ==0))   
	{
		print_crispoptions(); return 0;
	}
	else if (argc >=2 && (strcmp(argv[1],"--helpfull") ==0 || strcmp(argv[1],"--full") ==0) )        
	{
		print_crispoptions(); print_crispoptions_additional(); return 0;
	}


	bamfiles = parse_arguments(argc,argv,options,bamfilepaths,phenotypefile);
	if (bamfiles > 0)
	{
		options->bamfilelist = (char**)malloc(sizeof(char*)*bamfiles); 
		for (i=0;i<bamfiles;i++) options->bamfilelist[i] = (char*)malloc(1024);
		bamfiles=0;
		for (i=1;i<argc;i+=2)
		{
			if (strcmp(argv[i],"--bam") ==0 || strcmp(argv[i],"--bamfile") ==0)
			{
				strcpy(options->bamfilelist[bamfiles],argv[i+1]); bamfiles++;
			}
		}
		options->bamfiles = bamfiles; 
	}
	else if (strcmp(bamfilepaths,"None") !=0) 
	{
		if (parse_bamfilenames(bamfilepaths,options) <=0) return 0;
	}

	
	if (strcmp(phenotypefile,"None") !=0) 
	{
		options->phenotypes = calloc(options->bamfiles,sizeof(int)); options->association = 1;
		fp = fopen(phenotypefile,"r");
		for (i=0;i<options->bamfiles;i++)
		{
			fscanf(fp,"%d\n",&options->phenotypes[i]); 
		}
		fprintf(stderr,"read phenotype file %s \n",phenotypefile);
		//for (i=0;i<options->bamfiles;i++) fprintf(stderr,"%d:%d \n",i,options->phenotypes[i]);
		fclose(fp);
	}

	// parse phenotypefile if it is provided by user

	#if (PICALL ==0 || PICALL == 3)
		if (options->bamfiles < 2 || options->POOLSIZE ==0  || printoptions ==1 || strcmp(options->fastafile,"None") ==0) 
		{
			print_crispoptions();
			return 0;
		}
		else fprintf(stderr,"\nCRISP options: QVoffset %d min_base_quality %d min_mapping_score %d max_permutations %d poolsize %d CT-pvalue-thresh %.1f QV-pvalue-thresh %.1f\n\n",QVoffset,MINQ,MIN_M,MAXITER,options->POOLSIZE,thresh1,MINQpv);
	#elif PICALL ==1 
		if (options->bamfiles < 1  || printoptions ==1 || strcmp(options->fastafile,"None") ==0) 
		{
			print_picalloptions(); 
			return 0;
		}
		else fprintf(stderr,"\noptions: readlength-for-indels %d QVoffset %d min_base_quality %d min_mapping_score %d \nmax_mismatches_read %d LLR-threshold %f \n\n",READLENGTH,QVoffset,MINQ,MIN_M,MAX_MM,LLRthresh);

	#endif
	return 1;

}
