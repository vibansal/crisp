// data structures and functions for reading multiple bam files to generate a list of base-calls that can be passed on to a variant caller..
#ifndef INC_bamsreader_H
#define INC_bamsreader_H
#include <stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include<string.h>
#include<ctype.h>
#include "sam.h" // samtools header file for bam definitions
#include "readfasta.h"
#include "bamread.h"

#define BAM_PAIRED_READ1 (BAM_FPAIRED | BAM_FPROPER_PAIR | BAM_FREAD1)
#define BAM_PAIRED_READ2 (BAM_FPAIRED | BAM_FPROPER_PAIR | BAM_FREAD2)

//extern int SOLID;

typedef struct 
{
	struct alignedread* first; struct alignedread* last;  int reads;
	struct alignedread* current; // reads before this are past the current position being considered for variant calling
	struct alignedread* iter; struct alignedread* previous;
} READQUEUE;


struct bamINDEL
{
        uint32_t position; int16_t length; int16_t counts[2]; char printed;
	int32_t* cigar; int32_t cigs;  // counts
        char* bases; // pointer to either read->sequence[offset] or reflist->sequences[chr][position] 
	char simple;

        // most indels will occur in small number of samples, bogus or rare true indels | common indels -> occur in many-many samples...
	// cigarlist 1M5I or 10I5D for complex indel uint32_t cigarlist[2];
	// also store best cigars (left and right) for each indel
} bamINDEL;

struct REALIGNMENT
{
        int varid; // index to variantlist, variant that is included in new cigar
        uint32_t cigarlist[256]; int32_t cigs; int newpos;
        int added; int mismatches; int delta;
};


struct INDELHEAP
{
        struct bamINDEL** harray;  // array of pointers, allow for easy sorting and removing elements from heap 
	unsigned int maxsize; unsigned int length;  // current length
        int recent; // most recent indel inserted/updated into INDELHEAP
	int fve; // first valid element, all elements before this need to be deleted
};

// store counts for alleles at non-variant sites to estimate error rates 
struct BASE_COUNTS
{
	int counts[64][8]; // read 1/2 separate  
	// previous base + current base = 16 possibilities, 4 bases x 2 = 8 
	// QV bins as well, 0-10,11-20,21-30,30-40, sequencing cycle... 
	// AA -> A C G T |  AC -> A C G T
	// our goal is to detect if for some samples, the average base specific error rate is different than other samples rather than to recalibrate individual base quality values.. 
};


//int parse_cigar(struct alignedread* read,REFLIST* reflist,int current,int* fcigarlist,int sample_id);
//struct alignedread* get_read_bamfile(const bam1_t *b, void *data,struct alignedread* read);

// data structures and bam pointers associated with each bam file 
struct BAMFILE_data
{
	bam_index_t *idx; // bam file index
	bam_iter_t iter; int ret; // for reading indexed bam files 

	struct alignedread* read;   samfile_t* fp;  bam1_t* b; 
	int finished;   int readflag; 
	struct alignedread* first; struct alignedread* last; //pointers to first and last read in queue
	struct alignedread* previous; //
	struct INDELHEAP* ilist; // list of candidate indels found in this bam file 10/24/13
	int** counts; 
	//struct alignedread* *OPEheap; int heaplength; 
	// declare an array of pointers to alignedread object as heap (nov 19 2012), int length of pqueue 
	// dynamic re-allocation of memory if the size goes above what is allocated initially...
} ;

int free_read_queue(READQUEUE* RQ,struct BAMFILE_data* bamfiles_data); //free read pointed to by RQ->iter, previous is RQ->previous

void empty_queue(READQUEUE* bq,struct BAMFILE_data* bamfiles_data);

int init_bamfiles(struct BAMFILE_data* bamfiles_data,char** bamfilelist,int bamfiles,char* regions,int* ref,int* start,int* end);

//void free_read(struct alignedread* read);

// binary heap that implements priority queue 
// heap basics: leftchild = 2i+1  rightchild = 2i + 2, parent = i-1/2
// 1 3 6 5 9 8 heap elements
// 0 1 2 3 4 5 array index 

typedef struct 
{
	int* harray; int length; //int maxlength; 
} BAMHEAP;

//void swap(BAMHEAP* heap,int i,int j);
void minHeapify(BAMHEAP* heap,int node,struct BAMFILE_data* bamfiles_data);
void buildminheap(BAMHEAP* heap,struct BAMFILE_data* bamfiles_data);

void mergesortedbams(char** bamfilelist,int bamfiles);

void printpileup(REFLIST* reflist,int current,int k,READQUEUE* bq);

int allocate_mem_heap(struct BAMFILE_data* bamfiles_data,int bamfiles,int heapsize);

#endif

