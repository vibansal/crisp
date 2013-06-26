// data structures and functions for reading multiple bam files to generate a list of base-calls that can be passed on to a variant caller..
#ifndef INC_bamsreader_H
#define INC_bamsreader_H
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include<string.h>
#include<ctype.h>
#include "sam.h" // samtools header file for bam definitions
#include "../readfasta.h"

#define BAM_PAIRED_READ1 (BAM_FPAIRED | BAM_FPROPER_PAIR | BAM_FREAD1)
#define BAM_PAIRED_READ2 (BAM_FPAIRED | BAM_FPROPER_PAIR | BAM_FREAD2)

extern int SOLID;
extern char INT_CIGAROP[];

// reduce memory for this datastructure by using uint32_t instead of 'int' 
struct pileupread
{
        // chrom, readid, matechrom are just pointers to bam structures
        char* readid; short flag; uint8_t strand; uint8_t bidir;
        char* chrom; int position; char* matechrom; int mateposition;
        int mquality;  int IS;
        int tid; //used for bamread, index of chromosome name in bamheader list
        char* sequence; char* quality; int readlength;
        int cigs; int* cigarlist; // 36 M 2 D 5 I 54 M
        int fcigs; int* fcigarlist;
        int mismatches; int gaps; int alignedbases; uint8_t CM;
        int cflag;
        int sampleid; // added may 25 2012 for multi-bam reading
        uint8_t* MDstring;
        struct pileupread* next; // next read in queue shared across all samples 
        struct pileupread* nextread; // next read from same BAM file as this read
        struct pileupread* matepair; // pointer to mate pair for handling overlapping paired-end reads...
        int lastpos; int lastpos_mate; // lastposition on reference sequence covered by read 

        short cigoffset; int l1,l2; // offset in cigar string and within the cigar.... partitioned by cigar blocks 30M 4D 50M
        short delta; char filter; short allele; short type;
        char OPE; // flag to mark reads which overlap with their mate, '0' default, '1' OPE, '2' IS is shorter than readlength

} pileupread;

typedef struct
{
        struct pileupread* first; struct pileupread* last;  int reads;
        struct pileupread* current; // reads before this are past the current position being considered for variant calling
        struct pileupread* iter; struct pileupread* previous;
} READQUEUE;


int parse_cigar(struct pileupread* read,REFLIST* reflist,int current,int* fcigarlist);

struct pileupread* get_read_bamfile(const bam1_t *b, void *data,struct pileupread* read);

struct BAMFILE_data
{
        bam_index_t *idx; // bam file index
        bam_iter_t iter; int ret; // for reading indexed bam files 

        struct pileupread* read;   samfile_t* fp;  bam1_t* b;
        int finished;   int readflag;
        struct pileupread* first; struct pileupread* last; //pointers to first and last read in queue
        struct pileupread* previous; //
        //struct pileupread* *OPEheap; int heaplength; 
        // declare an array of pointers to pileupread object as heap (nov 19 2012), int length of pqueue 
        // dynamic re-allocation of memory if the size goes above what is allocated initially...
} ;

int free_read_queue(READQUEUE* RQ,struct BAMFILE_data* bamfiles_data); //free read pointed to by RQ->iter, previous is RQ->previous

void empty_queue(READQUEUE* bq,struct BAMFILE_data* bamfiles_data);

int init_bamfiles(struct BAMFILE_data* bamfiles_data,char** bamfilelist,int bamfiles,char* regions,int* ref,int* start,int* end);

void free_read(struct pileupread* read);

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

#endif

~                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
~                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
~                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
~                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
~                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
~                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
~                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
~                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
~                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
~                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
~                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
~                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
~                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          1,1           All
