#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include "variant.h"
//#include "allelecounts.h"

int FILTER_READS_MISMATCHES =1; // if set to 0, reads will not be filtered based on # of mismatches| use --filterreads 0 to set to 0
int FLANKING_BASES =0;
int CLIP_START=0,CLIP_END=0;

int OPE_OPS =0;

// The VCF spec states that only A|C|G|T|N bases are allowed in the ref position
// changed to reflect IUPAC ambiguity codes 'S' = C/G (83) before 'T' | ... and others that are not 'A'
unsigned int BTI[] = {
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 1, 1,  0, 0, 0, 2,  0, 0, 0, 2,  0, 0, 0, 0, 0, 0, 0, 1,  3, 0, 0, 0,  0, 1, 0, 0,  0, 0, 0, 0,  // A | C | G | T
	0, 0, 1, 1,  0, 0, 0, 2,  0, 0, 0, 2,  0, 0, 0, 0, 0, 0, 0, 1,  3, 0, 0, 0,  0, 1, 0, 0,  0, 0, 0, 0,  // a | c | g | t 
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0
};


#include "process_indel_variant.c"
#include "init_variant.c"

//for (i=0;i<bcall->fcigs;++i) fprintf(stdout,"%d:%d ",bcall->fcigarlist[i]>>4,bcall->fcigarlist[i]&0xf); fprintf(stdout," fc:%d %s %d\t",bcall->cigoffset,bcall->readid,bcall->position); fprintf(stdout,"mismatch %d %d %d \n",variant->position,bcall->l2,bcall->l1);

void print_read(struct alignedread* bcall, struct alignedread* bcall_mate,int j)
{
	fprintf(stdout,"OPE_INDEL | ");
	int i=0;
	for (i=0;i<bcall->fcigs;i++) fprintf(stdout,"%d%c",bcall->fcigarlist[i]>>4,INT_CIGAROP[bcall->fcigarlist[i]&0xf]); 
	fprintf(stdout," ");
   	for (i=0;i<bcall_mate->fcigs;i++) fprintf(stdout,"%d%c",bcall_mate->fcigarlist[i]>>4,INT_CIGAROP[bcall_mate->fcigarlist[i]&0xf]); 
	fprintf(stdout," | %s %d \nread %s %d OPE_INDEL\n",bcall->readid,bcall->IS,bcall->sequence,bcall->position);
	for (i=0;i<bcall_mate->position-bcall->position;i++) fprintf(stdout,"."); 
	fprintf(stdout,"mate %s %d OPE_INDEL\n",bcall_mate->sequence,bcall_mate->position);
}

// advance read function should be changed to call in an interval (start-end) and a set of haplotypes | return haplotype that the read matches....
// nov 18 2013 
// haplotypes could be single base or multi-base....  for SNP: A|G|C|T   for deletion:  CAAAAT | CAAAT  
// read at end could match multiple haplotypes or partial match
// read at start ??? realignment needs to be done before
// for complex variant TATT | GGTAATTTCC 

// if advance_read has been called once, if we call it again with same 'variant->position' it will not make a difference
// bcall->delta is a variable used for storing offset in a string of base matches (no gaps) 
// find the position in read that matches current position for variant analysis, advance bcall->cigoffset, l1, l2 values
int advance_read(struct alignedread* bcall,struct VARIANT* variant)
{
	//if (bcall->l2 >variant->position && bcall->cigoffset > 0 && (bcall->fcigarlist[bcall->cigoffset-1]&0xf) == BAM_CDEL) fprintf(stdout,"error %d %d %d\n",bcall->l2,variant->position,bcall->fcigarlist[bcall->cigoffset-1]>>4);
	
	//if (bcall->l2 > variant->position && bcall->cigoffset > 0 && (bcall->fcigarlist[bcall->cigoffset-1]&0xf) == BAM_CDEL) return 2;
	int op=0,l=0,i=0;
	while (bcall->l2 <= variant->position && bcall->cigoffset < bcall->fcigs)
	{
		op = bcall->fcigarlist[bcall->cigoffset]&0xf; l = bcall->fcigarlist[bcall->cigoffset]>>4;
		if (op == BAM_CDIFF) // single base mismatch 
		{
			if (bcall->l2 == variant->position)
			{
				bcall->delta  = 0; return 1;
			}
			else if (bcall->l2 < variant->position)	
			{
				bcall->l2++; bcall->cigoffset++; bcall->l1++;
			}
			else
			{
				if (PFLAG >=2) fprintf(stdout,"error, code should not be here\n");
				return 0;
			}
		}
		else if (op == BAM_CMATCH)
		{
			if (variant->position >= bcall->l2 && variant->position < bcall->l2 + l) // position is between the start and end
			{
				// ignore base if first base in 'M' and previous base was insertion allele
				if (bcall->l2 == variant->position  && bcall->cigoffset > 0 && (bcall->fcigarlist[bcall->cigoffset-1]&0xf) == BAM_CINS) return 0;
				bcall->delta = variant->position- bcall->l2;  return 1;
			}
			else //if (variant->position < bcall->l2)  
			{
				//fprintf(stderr,"here %d %d %d\n",variant->position,bcall->l2,l);
				bcall->l2 += l; bcall->l1 += l; bcall->cigoffset++;
			}
		}
		else if (op == BAM_CDEL || op == BAM_CREF_SKIP ) // allowing 'N' in cigar for RNA-seq
		{
			if (bcall->l2 == variant->position)
			{
				bcall->type =-1*l; 
				if (op == BAM_CDEL) return 1; else return 0;
			}
			// 11/13/13 this should be changed to add deletion allele if we are still in the bases covered by the deletion 
			else if (bcall->l2 < variant->position && bcall->l2 + l > variant->position)
			{
				return 2;
			}
			else 
			{
				bcall->l2 += l; bcall->cigoffset++; 
			}
		}
		else if (op == BAM_CINS)
		{
			if (bcall->l2 == variant->position)
			{
				bcall->type =l; return 1;
			}
			else if (bcall->l2 < variant->position)
			{
				bcall->l1 += l; bcall->cigoffset++;
			}
			else return 0;
		}
		else if (op == BAM_CSOFT_CLIP)
		{
			bcall->l1 += l; bcall->cigoffset++; 
		}
		else if (op == BAM_CHARD_CLIP) bcall->cigoffset++;
	}		
	return 0;
}
 

// this creates strand bias (since Forward reads are filtered) for overlapping paired-end reads, june 18 2012
// this is slowing down things for high depth of coverage....need to link read and mate once for each read pair 

int evaluate_OPE_readpair(struct alignedread* bcall, struct alignedread* bcall_mate,struct VARIANT* variant,uint8_t* strand,uint8_t* randomstrand,int* quality)
{
	int matefound = 0; int quality_m=0;	OPE_OPS=0;

	if (bcall->matepair != NULL && (bcall->matepair)->position == bcall->mateposition) bcall_mate = bcall->matepair; 
	else  // try to find the mate by moving in linked list
	{
		bcall_mate = bcall->nextread;
		while (bcall_mate != NULL && bcall->mateposition >= bcall_mate->position) 
		{
			OPE_OPS++;
			if (bcall->mateposition == bcall_mate->position && bcall->IS == -1*bcall_mate->IS && strcmp(bcall_mate->readid,bcall->readid) ==0) 
			{ 
				bcall->matepair = bcall_mate; bcall_mate->matepair = bcall; 
				bcall->lastpos_mate = bcall_mate->lastpos; bcall_mate->lastpos_mate = bcall->lastpos; 
				matefound = 1; break; 
			} 
			bcall_mate = bcall_mate->nextread; 

		}
		if (matefound ==0) return 0;
	}
	//if (matefound ==0 && PFLAG >=20) fprintf(stdout,"Error...matenotfound %d j %d OPE_OPS %d found:%d read %s %d %d %d %d\n",variant->position,j+1,OPE_OPS,matefound,bcall->readid,bcall->position,bcall->mateposition,bcall->IS,bcall->lastpos);

	bcall_mate->type = 0; 	
	if  (advance_read(bcall_mate,variant) ==0) bcall_mate->filter = '1'; 
	else if (bcall->type != bcall_mate->type)   // one of them is indel and not same as other one
	{
		if (OVERLAPPING_PE_READS ==2 && bcall->IS > 100) print_read(bcall,bcall_mate,0); // last arg is sample-id  05/07/2015 remove comment to print such reads 
		// this may happen if one end is not aligned with indel while other end is ?? 
		// we need to decide if we want to filter all such pairs or not ??
		// indel position ambiguity needs to be accounted for 
		if (bcall->type ==0) 
		{ 
			bcall->filter = '1'; variant->filteredreads[(bcall->strand==0) ? 5:6]++; 
		}
		else if (bcall_mate->type ==0)
		{
			bcall_mate->filter = '1'; variant->filteredreads[(bcall_mate->strand==0) ? 5:6]++;
		}
		else
		{
			bcall->filter = '1'; bcall_mate->filter = '1';
			variant->filteredreads[7]++;
			//fprintf(stdout,"filt %d:%d %d:%d \n",bcall->strand,bcall_mate->strand,bcall->type,bcall_mate->type);
		}
		// these reads are not counted ... important for indels in homopolymer runs...
	}
	else if (bcall_mate->type ==0)  // match for SNPs
	{
		quality_m = (int)bcall_mate->quality[bcall_mate->l1+bcall_mate->delta]-QVoffset;
		 
		if (bcall->sequence[bcall->l1+bcall->delta] == bcall_mate->sequence[bcall_mate->l1+bcall_mate->delta] || (*quality-quality_m) >=MINQ) 
		{
			bcall_mate->filter = '1'; bcall->bidir = 1; 
			if (*randomstrand == 0) *randomstrand = 1; else *randomstrand = 0;
			*strand = *randomstrand + 2; // bidirectional marked 
			variant->filteredreads[(*randomstrand ==0) ? 2:3]++;
			if (bcall_mate->mquality < quality_m) quality_m = bcall_mate->mquality;
			if (*quality < quality_m)
			{	
				*quality =quality_m;
				bcall->quality[bcall->l1+bcall->delta]= quality_m+QVoffset;
			}
			//if (bcall->sequence[bcall->l1+bcall->delta] != variant->refb) fprintf(stdout,"%d read %s %d %c %c %d:%d \n",j+1,bcall->readid,bcall->position,base,bcall_mate->sequence[bcall_mate->l1+bcall_mate->delta],bcall->strand,bcall_mate->strand);
		}
		else
		{ 
			// single base  mismatches
			bcall->filter = '1'; bcall_mate->filter = '1'; variant->filteredreads[4]++;
			if (bcall->sequence[bcall->l1+bcall->delta] == variant->refb) variant->filteredreads[(bcall->strand==0) ? 5:6]++;
			else if (bcall_mate->sequence[bcall_mate->l1+bcall_mate->delta] == variant->refb) variant->filteredreads[(bcall_mate->strand==0) ? 5:6]++;
		}	
	}
	else if (bcall->mismatches < bcall_mate->mismatches)  // case: indels matching from both ends
	{
		bcall_mate->filter ='1'; bcall->bidir = 1;  // counted in bidirection read count 
	}
	else 
	{
		bcall->filter = '1'; bcall_mate->bidir =1;
	}
	return 1;
}


// TOFIX for deletions, when we move to next base deleted allele should be counted... june 24 2012 
// iterate over all reads covering 'variant->position' and calculate #ofreads supporting A,C,T,G,+A,-A.....
void calculate_allelecounts(REFLIST* reflist,int current,READQUEUE* bq,struct BAMFILE_data* bamfiles_data,struct VARIANT* variant)
{
	int insertions =0,deletions =0,indelalleles=0;
	char base,base_m; uint8_t strand,randomstrand=0;
	int quality;
	int ibase;
	int offset=0; double ep; int allele; int allelematch=0; int bin=0;
	struct alignedread* bcall; struct alignedread* bcall_mate; // bcall_mate is matepair of the read pointed to by bcall
	int i=0; int reads=0; int j=0; int ar=0;

	// before we calculate reads supporting each allele -> merge candidate indels from all samples and make a single sorted list of candidate indels...
	//if (INDEL_REALIGNMENT ==1) 
	{
		//for (j=0;j<variant->samples;j++) print_indels_interval(bamfiles_data,j,variant->position,variant->position);
	}

	for (j=0;j<variant->samples;j++)
	{
		//if (variant->samples == variant->bamfiles) j = j; else j = variant->options->BAM_TO_SAMPLE[j]; 
		for (bcall=bamfiles_data[j].first; bcall != NULL && bcall->position <= variant->position;bcall=bcall->nextread) 
		{
			// reads could be realigned here: if there is an indel close by or current position has an indel allele 
			bcall->filter = '0'; bcall->bidir = 0; // initialize the two variables for each read/base-call
			if (INDEL_REALIGNMENT ==2 && bcall->realigned ==0 && bcall->mismatches > 0) realign_read(bamfiles_data,j,bcall,reflist);
		}
		// even if there are additional reads in the bcall queue at the end, they will not affect performance of this loop..
		for (bcall=bamfiles_data[j].first; bcall != NULL && bcall->position <= variant->position;bcall=bcall->nextread)
		{
			// reads in queue sorted by position so lastpos not correct to break out of loop if there are reads of different lengths....
			if (bcall->lastpos < variant->position) continue; // changing this to <= caues BUG | june 20 2012
			bcall->allele = -1; bcall->type = 0;  
			ar = advance_read(bcall,variant);
			if (ar ==0) 
			{
				bcall->filter = '1'; continue; // extrafilter
			}
			
			if (bcall->type ==0 && ar ==1)  // advance_read can have 3 values: 0/1/2
			{
				quality = (int)bcall->quality[bcall->l1+bcall->delta]-QVoffset; 
				base = bcall->sequence[bcall->l1+bcall->delta];
				// extra check min of base and mapping quality | if (bcall->mquality < quality) quality = bcall->mquality;
			}
			strand = bcall->strand;

			// special filter for two cases: 1st read starts after 2nd read, 2nd read ends before 1st read in pair 
			if ((bcall->IS > 0 && bcall->lastpos_mate > 0 && variant->position >= bcall->lastpos_mate) || (bcall->IS < 0 &&  variant->position < bcall->mateposition)) bcall->filter = '1';

			// calculation of OPE for indel reads does not account for ambiguous reads marked as reference !!
			if (OVERLAPPING_PE_READS >=1 && bcall->position < bcall->mateposition && variant->position > bcall->mateposition && bcall->filter == '0' && ar ==1)
			{
				evaluate_OPE_readpair(bcall,bcall_mate,variant,&strand,&randomstrand,&quality);
			}
			if (bcall->filter == '1') continue;

			variant->readdepths[j]++; 
			// 0-9, 10-19,20-39 40+ 37,60 are BWA mapping qualities for single-reads and paired-end reads
			if (bcall->mquality <10) variant->MQcounts[0]++;
			else variant->MQcounts[(bcall->mquality >= 40) ? 3: bcall->mquality/20+1]++; 
			// filter for MAX_MM mismatches allowed per read added june 12 2012
			// this should be local filter so that mismatches at one end of read don't affect full read
			MAX_MM = bcall->alignedbases/40 + 2; 
			//if (bcall->alignedbases < 40) MAX_MM = 2; else if (bcall->alignedbases < 75) MAX_MM = 3; 	else if (bcall->alignedbases < 100) MAX_MM = 4; 
			if (bcall->mismatches+bcall->gaps > MAX_MM && bcall->mquality >= MIN_M) variant->filteredreads[0]++; // # of mismatches greater than expected
			else if (quality < MINQ && bcall->mquality >= MIN_M) variant->filteredreads[1]++;

			bcall->filter = '1'; // default here is 1, if read passes the filters below it is set to 0

			if (ar ==2) // special tag for reads that span deletion and variant->position is in middle of deletion, artificially add reads supporting reference allele, hack 13/11/13  
			{
				base = variant->refb;	bcall->type =0; quality= 30; bcall->bidir = 0;
				// problem will arise in code for calculating likelihoods that iterates over bcall
			}

			if (quality >= MINQ && bcall->mquality >= MIN_M && (bcall->mismatches+bcall->gaps <= MAX_MM || FILTER_READS_MISMATCHES ==0) && bcall->type == 0)
			{
				if (strand == 1 || strand ==3) offset = maxalleles; else offset = 0;
				if ( (base == variant->refb && (bcall->mismatches+bcall->gaps <= MAX_MM-1 || FILTER_READS_MISMATCHES ==0) ) || base != variant->refb) 
				{
					ibase = BTI[base]; bcall->allele = ibase+offset;
					if (bcall->bidir ==1) 
					{
						variant->counts[ibase+2*maxalleles]++; variant->indcounts[j][ibase+2*maxalleles]++; 
					}
					else
					{
						variant->counts[ibase+offset]++; variant->indcounts[j][ibase+offset]++; 
					}
					bin = quality/10-1;
					if (bin >=0) variant->indcounts_binned[j][ibase+offset][bin > 2 ? 2: bin]++;
					if ( (bcall->flag & BAM_PAIRED_READ1) == BAM_PAIRED_READ1) variant->indcounts_binned[j][ibase+offset][3]++;
					else if ((bcall->flag & BAM_PAIRED_READ2) == BAM_PAIRED_READ2 ) variant->indcounts_binned[j][ibase+offset][4]++;
					else variant->indcounts_binned[j][ibase+offset][5]++;
					ep = pow(0.1,(double)quality/10); 
					variant->stats[j][ibase+offset] +=ep; variant->tstats[ibase+offset] += ep; 
					if (USE_QV ==1) variant->Qhighcounts[j][ibase+offset] += 1.0-ep; 
					else variant->Qhighcounts[j][ibase+offset] +=1;
					bcall->filter = '0'; //bcall->type =type;
					if (ar ==2) bcall->filter = '1'; // this ensures that these bases will not be used for likelihood calculation...
				}
			}
			else if (bcall->mquality >= MIN_M && (bcall->mismatches+bcall->gaps <= MAX_MM || FILTER_READS_MISMATCHES ==0)  && bcall->type != 0) 
			{
				if (bcall->type < 0) deletions++; else insertions++;
				bcall->filter = '0'; //bcall->type = type;
			}
			reads++;
		}
	}

	// this function can be merged
	indelalleles = 0; 
	if (insertions >= 3 || deletions >= 3) 
	{
		//fprintf(stdout," allele counts for pos %d reads %d ireads %d ",variant->position,reads,insertions+deletions);
		indelalleles = identify_indelalleles(reflist,current,bq,variant);
	}
}




/*
   for (i=0;i<bcall->cigs;i++) fprintf(stdout,"%d%c",bcall->cigarlist[i]>>4,INT_CIGAROP[bcall->cigarlist[i]&0xf]); 
   fprintf(stdout," j:%d ",j);
   for (i=0;i<bcall_mate->cigs;i++) fprintf(stdout,"%d%c",bcall_mate->cigarlist[i]>>4,INT_CIGAROP[bcall_mate->cigarlist[i]&0xf]); 
   fprintf(stdout," filter:%c:%c ",bcall->filter,bcall_mate->filter);
   fprintf(stdout,"strand %d read %s %d IS %d matepos %d variant %d %d %d\n",bcall->flag,bcall->readid,bcall->position,bcall->IS,bcall->mateposition,variant->position,bcall->sampleid,bcall->lastpos);
   fprintf(stdout,"matefound %d type %d %d ",matefound,bcall->type,bcall_mate->type);
   fprintf(stdout,"%c:%c %c:%c ",base,quality,base_m,quality_m);
   if (bcall->type != bcall_mate->type) fprintf(stdout,"diffalleles ");
   bcall_mate->filter = '1'; variant->filteredreads[2]++;  
 */
//if (j==6 && strcmp(bcall->readid,"HISEQ:73:B00NUACXX:1:1304:13523:102166")==0) fprintf(stderr,"testread %d j %d OPE_OPS %d found:%d read %s %d %d %d\n",variant->position,j+1,OPE_OPS,matefound,bcall->readid,bcall->position,bcall->mateposition,bcall->IS);


//fprintf(stdout,"filter %d j %d read %s %d %d %d %d\n",variant->position,j+1,bcall->readid,bcall->position,bcall->mateposition,bcall->IS,bcall->lastpos_mate);
