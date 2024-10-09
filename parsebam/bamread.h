// data structures and functions for reading multiple bam files to generate a list of base-calls that can be passed on to a variant caller..
#ifndef INC_bamread_H
#define INC_bamread_H
#include <stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include<string.h>
#include<ctype.h>
#include "sam.h" // samtools header file for bam definitions
#include "readfasta.h"

extern char INT_CIGAROP[];
extern int SOLID;

// reduce memory for this datastructure by using uint32_t instead of 'int' 
struct alignedread
{
	// chrom, readid, matechrom are just pointers to bam structures
	char* readid; short flag; uint8_t strand; uint8_t bidir; 
	char* chrom; int position; char* matechrom; int mateposition; 
	int mquality;  int IS;
	int tid; //used for bamread, index of chromosome name in bamheader list
	char* sequence; char* quality; int readlength;
	int cigs; int* cigarlist; /* 36 M 2 D 5 I 54 M*/	int fcigs; int* fcigarlist;

	// add extra fields for indel code compatibility, we can use same data structure for both codes
	short XC, cflag, clipped; char matestrand; int span; 

	int mismatches; int gaps; int alignedbases; uint8_t CM;
	int sampleid; // added may 25 2012 for multi-bam reading
	uint8_t* MDstring;
	struct alignedread* next; // next read in queue shared across all samples 
	struct alignedread* nextread; // next read from same BAM file as this read
	struct alignedread* matepair; // pointer to mate pair for handling overlapping paired-end reads...
	int lastpos; int lastpos_mate; // lastposition on reference sequence covered by read 

	short cigoffset; int l1,l2; int last1,last2; // offset in cigar string and within the cigar.... partitioned by cigar blocks 30M 4D 50M
	short delta; char filter; short allele; short type;
	char OPE; // flag to mark reads which overlap with their mate, '0' default, '1' OPE, '2' IS is shorter than readlength
	uint8_t realigned;

};

struct alignedread* get_read_bamfile(const bam1_t *b, void *data,struct alignedread* read);

int parse_cigar(struct alignedread* read,REFLIST* reflist,int current,int* fcigarlist);//int sample_id);

// basic function to left align indels in fcigarlist generated after parsing read from bam file 
// read->cigarlist will not match fcigarlist after this function is called 
int left_align_indels(struct alignedread* read,REFLIST* reflist,int current);

void free_read(struct alignedread* read);

#endif

