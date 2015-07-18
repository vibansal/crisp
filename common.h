#ifndef INC_common_H
#define INC_common_H

#include <stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>

extern int EMflag; // if this is 1, then we use new crisp EM method, otherwise old method

extern int PIVOTSAMPLE; // sample that divides list of pools into two groups
extern int QVoffset;
extern int MAX_MM;
extern int MIN_M;
extern int MINQ;
extern int maxalleles;
extern int readlength, READLENGTH;
//extern int PICALL; // model for variant calling, CRISP / PICALL / LOWCOVERAGE 
extern int PFLAG;
extern int USE_QV, QVset;
extern int FILTER_READS_MISMATCHES;
extern int FLANKING_BASES;

extern int MINFLANK; // for indels
extern int CLIP_START,CLIP_END; // discard bases at beginning of reads
extern int OVERLAPPING_PE_READS; // if set to 1 we determine these reads and count them properly, if not we just ignore this
extern int ALLOW_AMBIGUOUS_BASES; // 0 = positions with ambiguous base will not be used for variant calling, 1 will be used 
extern int LEFT_ALIGN_INDELS; 
extern int SPLIT_TRIALLELIC_VARS;
extern int INDEL_REALIGNMENT;
extern double AGILENT_REF_BIAS; 
extern int USE_DUPLICATES;
extern int VARIANTS_CALLED;

extern int FAST_FILTER; // filter false variants using contingency table or chi-square method

extern unsigned int BTI[]; // array that stores mapping of bases to integers 'A' -> 0 , 'C' -> 1, 'G' -> 2

char* OUTPUT_ALLELE_COUNTS; 

struct OPTIONS
{
        char** bamfilelist; // list of filenames
	char** sampleids; // list of sampleids corresponding to the bamfilelist 
	int* BAM_TO_SAMPLE; // which bamfile corresponds to which SAMPLE
        int bamfiles; int samples; 
	int multiplebams; // 0 => one bam for each sample, 1 -> some samples have multiple bam files 

        int READLENGTH; int qvoffset;
        int MIN_Q; int MIN_M; int MAX_MM; int MIN_POS; int MAX_POS; int MIN_READS; int POOLSIZE;
        double ctpval; double qvpval; double binpval; double SER;
        int MAXITER; char fastafile[1024]; char vcffile[1024]; char bedfile[1024];
        int VCmodel; // variant calling model crisp=0, picall=1, lowcoverage=2
        FILE* vfile; // vcffile
	char* regions; int targettid,targetstart,targetend; // targeted interval using bamindex
	char indelfile[1024]; FILE* fp_indelfile; // VCF file with indel candidates 
	int *ploidy; int varpoolsize;
	int* phenotypes; int association; // 0/1/2 for each pool to do case control association testing

};

typedef struct
{
	int al; int ct; double vars;
} ALLELE;

#ifndef max
        #define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif
#ifndef min
        #define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif


#endif
