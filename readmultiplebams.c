
/* MAIN CODE for running CRISP */

#include "parsebam/bamread.h"
#include "parsebam/bamsreader.h"
#include "parsebam/variant.h"
#include "parsebam/allelecounts.h"
#include "FET/contables.h"
#include "crisp/crispcaller.h"

#include "indels/generate_candidate_indels.c"
#include "indels/realignment.c"

// if this variable is set to 1, the bayesian-method is used for all variants SNPs/indels
//#define PICALL 1, picall variable is actually passed to the program at compile time using -D option (0 for crisp, 2 for low coverage, 3 for BFGS based crisp)

int USE_DUPLICATES = 0; // use PCR duplicate reads for variant calling or not 
uint32_t BAM_FILTER_MASK = 0;


int SOLID=0; // for reads on SOLiD platform
int targeted=0;

int QVset = 0; int USE_QV = 0; // added dec 21 2011 to allow for weighted table use in chi-square calculation
int QVoffset =33;
int MIN_M = 20, MAX_MM = 4; int MINQ = 13;  //changed to 13 as default
int maxalleles = 8; // maximum number of alleles at chrom, A,C,T,G -A -AA -AAA +A 
int READLENGTH = 100; // average read length used to calculate expected heterozygosity for indels 
int MAX_COV=60; // maximum coverage for lowcoverage variant calling per sample

/* parameters for log-likelihood ratio statistic */
double LLRthresh =4; double HWEprior = 1;
// alpha beta should be integers, can be doubles but then exact integration cannot be done
double alpha0 = 1; double beta0 =50; // prior error rate of 0.01 
double alpha = 1; double beta = 50; // prior error rate of 0.01 
double alpha1 = 1; double beta1 =50;
/* parameters for log-likelihood ratio statistic */

double theta = 0.0001;
int MINCOV = 2; int MINFLANK = 5; double MAXE=0.03;
int INDELSONLY = 1;
int OVERLAPPING_PE_READS = 1; 
int SPLIT_TRIALLELIC_VARS = 0;
int INDEL_REALIGNMENT = 0;
int PIVOTSAMPLE =0;

int PFLAG =2;
int VARIANTS_CALLED=0;

int MIN_COVERAGE_POOL=1, MIN_COVERAGE_FLANKING=0,HAPLOTYPES=0; // min. average coverage per haplotype for site to be considered for variant calling
int CALL_VARIANTS = 1;
int ALLOW_AMBIGUOUS_BASES = 0;

int CALCULATE_ERROR_RATES = 0;

#include "optionparser.c"

//  maintain for each pool: pointer to first (and last) read in shared queue, can be easily updated
// maintain for each read: pointer to next read of same pool in queue 

#include "variantcalls.c"

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// multi sample variant caller: CRISP, PICALL or low coverage method
int multisampleVC(struct OPTIONS* options,REFLIST* reflist,FILE* fp)
{
	if (USE_DUPLICATES ==1) BAM_FILTER_MASK = (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL); else BAM_FILTER_MASK = (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP);

	int bamfiles = options->bamfiles;

	int last=0; // last is the current position s.t. all reads have starting position > last
	int i=0; int h=0;
	unsigned long reads=0; int j=0; int prev_tid = -1;   int rf=0;
	int finishedfiles =0; 
	struct alignedread* pread = NULL;
	struct BAMFILE_data* bamfiles_data = calloc(bamfiles,sizeof(struct BAMFILE_data)); // added one extra to list to store indels for all samples combined

	READQUEUE* RQ = (READQUEUE*)malloc(sizeof(READQUEUE));  RQ->first = NULL; RQ->last = NULL; RQ->reads = 0; 
	int* fcigarlist = (int*)malloc(sizeof(int)*4096*128); // increased size to work for long reads, 07-25-2023

	// data structure for holding potential variants and read counts, etc 
	struct VARIANT variant;  variant.ploidy = calloc(options->bamfiles,sizeof(int)); 
	init_poolsizes(&variant,options,PICALL); 
	init_variant(&variant,options->bamfiles,options->bamfiles);
	variant.options = options;  // pointer to options

	BAMHEAP bheap; bheap.harray = (int*)malloc(sizeof(int)*bamfiles); bheap.length = bamfiles;
	for (i=0;i<bamfiles;i++) { bheap.harray[i] = i; bamfiles_data[i].finished= 0;}
	
	reflist->cinterval = -1; // first interval to the right of current base

        init_bamfiles(bamfiles_data,options->bamfilelist,bamfiles,options->regions,&options->targettid,&options->targetstart,&options->targetend);

	// error when reading indexed bam files probably due to lack of reads in some files resulting in heap error, fixed oct 17 2012
        j=0; for (i=0;i<bamfiles;i++) 
	{
		finishedfiles += bamfiles_data[i].finished; 
		if (bamfiles_data[i].finished ==0) bheap.harray[j++] = i; else bheap.length--; 	
	}
	buildminheap(&bheap,bamfiles_data); // initial minheap call
	//fprintf(stderr,"finishedfiles %d \n",finishedfiles);
	
	if (INDEL_REALIGNMENT >=1) allocate_mem_heap(bamfiles_data,bamfiles,100);

	
	HAPLOTYPES =0,MIN_COVERAGE_FLANKING =0;
	for (i=0;i<variant.samples;i++) 
	{
		MIN_COVERAGE_FLANKING += 2*variant.ploidy[i];  // enforced for regions outside the bedfile target
		HAPLOTYPES += variant.ploidy[i];
	}
	//int min_coverage_target = 1*variant->ploidy*variant->samples;  // enforced for regions outside the bedfile target
	int offset_readlength = 150;  // call variants in window (last,current_read_position-offset_readlength) to allow for indel analysis, set to 0 for original behavior of program
	// the value of offset should not affect the correctness or speed of the code
	int current_position =0;
	
	while (finishedfiles < bamfiles)
	{
		i = bheap.harray[0]; // take the top read off the heap
		if ( !(bamfiles_data[i].read->flag & BAM_FILTER_MASK))
		{
			if (bamfiles_data[i].read->tid != prev_tid) // read's chromosome is different from previousread 
			{
				if (prev_tid >=0)  // finish the processing of previous chromosome and cleanup
				{
					if (RQ->reads >0) 
					{
						fprintf(stderr,"processing %d reads left in queue for chrom %s...",RQ->reads,reflist->names[prev_tid]);
						callvariants(reflist,prev_tid,last,reflist->lengths[prev_tid],RQ,bamfiles_data,options,&variant);
						empty_queue(RQ,bamfiles_data); //clean thequeue
					}
					if (INDEL_REALIGNMENT >=1) clean_indel_lists(bamfiles_data,bamfiles,-1); current_position = 0; 
					for(j=0;j<bamfiles;j++) bamfiles_data[j].last=NULL; last =0; 
					free(reflist->sequences[prev_tid]); 
					fprintf(stderr,".....finished processing reads for chrom %s\n",reflist->names[prev_tid]);
					fprintf(stdout,".....finished processing reads for chrom %s\n",reflist->names[prev_tid]);
					reflist->cinterval = -1; // reset to -1 
				}
				read_chromosome(reflist,bamfiles_data[i].read->tid,fp); 
				prev_tid =bamfiles_data[i].read->tid;
			}

			if (bamfiles_data[i].read->position <last)
			{
				fprintf(stderr,"reads out of order i:%d h:%d pos: %d %d\n",i,h,bamfiles_data[i].read->position,last);
				fprintf(stderr,"the program will now exit, please sort the bamfiles\n");
				return 1;
			}

			if (INDEL_REALIGNMENT >=1 && bamfiles_data[i].read->position > current_position+offset_readlength) 
			{
				// need to clean up indel lists when we encounter a new chromosome... 
				print_indel_lists(bamfiles_data,bamfiles,current_position+offset_readlength); 
				clean_indel_lists(bamfiles_data,bamfiles,current_position);
				current_position = bamfiles_data[i].read->position;
			}
			// realign reads before calling variants, each read is realigned only once

			// small bug here, only call variants when last is less than current read position
			// bug fixed here, update last only when 'callvariants' is invoked, ???
			if (RQ->reads > 0 && bamfiles_data[i].read->position > last+offset_readlength) 
			{
				callvariants(reflist,bamfiles_data[i].read->tid,last,bamfiles_data[i].read->position-offset_readlength,RQ,bamfiles_data,options,&variant);  
			}
			last = bamfiles_data[i].read->position-offset_readlength; if (last < 0) last =0;

			bamfiles_data[i].read->cflag = 0; 
			// this function should only be called on reads inside/close_to targeted regions..
			parse_cigar(bamfiles_data[i].read,reflist,bamfiles_data[i].read->tid,fcigarlist); 

			if (INDEL_REALIGNMENT >=1 && bamfiles_data[i].read->gaps > 0 && bamfiles_data[i].read->mquality >= 20) extract_indel_reads(bamfiles_data[i].read,reflist,bamfiles_data[i].read->tid,i,bamfiles_data[i].ilist);
			
			//fprintf(stdout,"read s:%d IS:%d %s %d \n",i,bamfiles_data[i].read->IS,bamfiles_data[i].read->readid,bamfiles_data[i].read->position);
			if (RQ->last == NULL)
			{
				RQ->last = bamfiles_data[i].read; RQ->first = RQ->last; (RQ->last)->next = NULL;
				RQ->reads++;
			}
			else
			{
				(RQ->last)->next = bamfiles_data[i].read; RQ->last = bamfiles_data[i].read; 
				(RQ->last)->next = NULL;
				RQ->reads++;
			}
			if (bamfiles_data[i].last ==NULL) bamfiles_data[i].first = RQ->last;
			else bamfiles_data[i].last->nextread= RQ->last;
			bamfiles_data[i].last = RQ->last; (RQ->last)->nextread =NULL;
			// read that passes filters from 'i'th bam file is inserted in queue, should also add it to OPE queue 
			//if (bamfiles_data[i].read->position < bamfiles_data[i].read->mateposition && bamfiles_data[i].read->lastpos > bamfiles_data[i].read->mateposition) 
			//fprintf(stdout,"B %d %s %d %d %d \n",i,bamfiles_data[i].read->readid,bamfiles_data[i].read->position,bamfiles_data[i].read->mateposition,bamfiles_data[i].read->IS);
		}
		else free_read(bamfiles_data[i].read);
		//fprintf(stdout,"read from %d %d %s\n",i,bamfiles_data[i].read->position,bamfiles_data[i].read->readid);

		if (options->regions ==NULL) rf =samread(bamfiles_data[i].fp,bamfiles_data[i].b);
		else rf  = bam_iter_read(bamfiles_data[i].fp->x.bam,bamfiles_data[i].iter,bamfiles_data[i].b);
		if (rf >=0)
		{
			bamfiles_data[i].read = get_read_bamfile(bamfiles_data[i].b,bamfiles_data[i].fp,pread); 
			//if (options->samples ==0) bamfiles_data[i].read->sampleid = i;
			//else bamfiles_data[i].read->sampleid = options->BAM_TO_SAMPLE[i];  
			// bug here june 30 2013 commented out .... in 12 T2D pools 
			bamfiles_data[i].read->sampleid = i;
			if (!(bamfiles_data[i].read->flag & BAM_FILTER_MASK)) minHeapify(&bheap,0,bamfiles_data);
		}
		else // no more reads in file 'i' 
		{ 
			bamfiles_data[i].finished = 1; bamfiles_data[i].read= NULL; 
			bam_destroy1(bamfiles_data[i].b);
			h++; finishedfiles++; 
			//fprintf(stderr,"finished reading bam file %s \n",options->bamfilelist[i]); //return 1;
			bheap.harray[0] = bheap.harray[bheap.length-1]; bheap.length--;
			if (bheap.length > 0) minHeapify(&bheap,0,bamfiles_data);
			// call minheapify like function to push sample i off the heap, reduce heap size
		} 
		if ((++reads)%1000000 ==0 && RQ->reads >0) fprintf(stderr,".....processed %ld reads QSIZE:%d %s:%d:%d variants called %d\n",reads,RQ->reads,RQ->first->chrom,RQ->first->position,RQ->first->lastpos,VARIANTS_CALLED);
	}

	if (prev_tid >=0)  // finish the processing of last chromosome 
	{
		if (RQ->reads >0) 
		{
			fprintf(stderr,"processing %d reads left in queue for chrom %s.....",RQ->reads,reflist->names[prev_tid]);
			if (reflist->lengths[prev_tid] > last) callvariants(reflist,prev_tid,last,reflist->lengths[prev_tid],RQ,bamfiles_data,options,&variant);
			empty_queue(RQ,bamfiles_data); //clean thequeue
		}
		else fprintf(stderr,"queue for chrom %s is empty ",reflist->names[prev_tid]);
		free(reflist->sequences[prev_tid]); 
		fprintf(stderr,"finished processing reads for chrom %s \n\n",reflist->names[prev_tid]);
		if (INDEL_REALIGNMENT >=1) 
		{
			print_indel_lists(bamfiles_data,bamfiles,reflist->lengths[prev_tid]); 
			clean_indel_lists(bamfiles_data,bamfiles,reflist->lengths[prev_tid]);
		}
	}
	fprintf(stderr,"CRISP has finished processing bam files: total reads processed %ld total variants called %d \n\n",reads,VARIANTS_CALLED);

	//for (i=0;i<bamfiles;i++) bam_destroy1(bamfiles_data[i].b);
	free(bamfiles_data); free(bheap.harray); free(fcigarlist);
	//empty_queue(RQ); //clean thequeue
	//fprintf(stdout,"FILE %d %s %d %s %d %d %d mapped %d \n",i,read->readid,read->flag,read->chrom,read->position,read->mquality,read->IS,(read->flag &4));
	return 1;
}


int main(int argc, char* argv[])
{
	//fprintf(stderr,"size of indel element %d \n",sizeof(struct VCF_ALLELE));
	//fprintf(stdout,"%d %d %d\n",sizeof(uint8_t),BAM_PAIRED_READ1,BAM_PAIRED_READ2); return 1;
	int vflag=0;
	time_t now; time(&now);    unsigned int iseed = (unsigned int)time(NULL);  srand48(iseed);
	struct OPTIONS* options = (struct OPTIONS*)malloc(sizeof(struct OPTIONS)); options->POOLSIZE = 2;
	options->targettid =-1; options->targetstart=0; options->targetend=0;
	int flag = optparser(argc,argv,options); if (flag ==0) return 1;

	REFLIST reflist;
	if (read_fastaheader(options->fastafile,&reflist) == -1)  return -1;  strcpy(reflist.fastafile,options->fastafile);
	FILE* fp = fopen(options->fastafile,"r"); 

	if (read_bedfile(options->bedfile,&reflist) != -1) targeted = 1; else  targeted = 0; 

	if (strcmp(options->vcffile,"None") ==0) options->vfile = stdout; else options->vfile = fopen(options->vcffile,"w");

	// open indel file with candidate indels 
	if (strcmp(options->indelfile,"None") !=0) options->fp_indelfile = fopen(options->indelfile,"w"); else options->fp_indelfile = NULL;

	// print VCF to stdout as well if no VCFfile is specified 
	if (options->vfile != NULL && strcmp(options->vcffile,"None") !=0)
	{
		vflag =1; print_crispheader(options);
	}

	if (options->bamfiles >=2) fprintf(stderr,"processing %d bamfiles: %s ..... %s \n\n",options->bamfiles,options->bamfilelist[0],options->bamfilelist[options->bamfiles-1]);
	multisampleVC(options,&reflist,fp);
	if (vflag ==1) fclose(options->vfile); 
	fclose(fp); // close pointer to fasta reference file 
	return 1;

}

