//#ifndef INC_generatecandidates_H
//#define INC_generatecandidates_H
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include<string.h>
#include<ctype.h>

// store array of indels and their counts for each sample separately: for realignment, merge and sort arrays...
// only stores indels from a single chromosome at a time, so no need for chromosome identifier

// current approach is to realign using indels for each sample individually... delete indels that start after current read position
// this could also be useful for finding candidate indels at each locus for variant calling, this is currently implemented using a different approach
// easier way of finding complex indels on single read, compared to other approach

//struct bamINDEL: uint32_t position; int16_t length; int16_t counts[2];
//struct INDELHEAP: struct bamINDEL* harray; int maxsize;  int length; int recent; 

struct bamINDEL** tp;
struct bamINDEL* swap_pointer; 
struct INDELHEAP* heap;

int compare_elements(struct INDELHEAP* heap,int lc,int rc)
{
	if (heap->harray[lc]->position == heap->harray[rc]->position) return heap->harray[lc]->length -heap->harray[rc]->length;
	else return heap->harray[lc]->position - heap->harray[rc]->position; 
}

// we can keep the list sorted by position, if adding new element that is less than current last, swap pointers
int addElement(struct INDELHEAP* heap,int position,int length,int strand)
{
	int i = 0,addnewelement = -1;
	// first check if the indel is present in heap or not | if not create a new element 
	if (heap->length ==0) addnewelement = 1;
	else if (heap->harray[heap->recent]->position == position && heap->harray[heap->recent]->length == length) addnewelement = 0;
	else // check list in reverse order to see if indel is there 
	{
		i=heap->length-1; addnewelement = 1; 
		while (i >= 0 && position <= heap->harray[i]->position)
		{
			if (heap->harray[i]->position == position && heap->harray[i]->length == length) 
			{ 
				heap->recent = i; addnewelement =0; break;
			} 
			i--;
		}
	}
	if (addnewelement ==1)
	{
		if (heap->length >= heap->maxsize) 
		{
			tp = (struct bamINDEL**)realloc(heap->harray,2*heap->maxsize*sizeof(struct bamINDEL*)); 
			if (tp == NULL) { fprintf(stdout,"error could not re-allocate heap->harray problem.... \n"); exit(0); }
			else heap->harray = tp; //fprintf(stdout,"increased size of heap to %d \n",heap->maxsize);
			for (i=heap->maxsize;i<2*heap->maxsize;i++) heap->harray[i] = NULL; 
			heap->maxsize *=2;
		}
		if (heap->harray[heap->length] == NULL) heap->harray[heap->length] = calloc(1,sizeof(struct bamINDEL)); 
		heap->harray[heap->length]->position = position; heap->harray[heap->length]->length = length; 
		heap->harray[heap->length]->counts[0] = 0; heap->harray[heap->length]->counts[1] =0; heap->harray[heap->length]->printed = '0';
		if (strand ==0) heap->harray[heap->length]->counts[0] = 1; else heap->harray[heap->length]->counts[1] =1; 

		// if new element's position is less than previous indel in list, bubble this new indel upwards in list to maintain sorted order..
		heap->recent = heap->length; 

		i =heap->length; while (i > 0 && compare_elements(heap,i-1,i) > 0) 
		{
			//fprintf(stderr,"swapping pointers %d %d \n",heap->harray[i-1]->position,heap->harray[i]->position);
			swap_pointer = heap->harray[i-1]; heap->harray[i-1] = heap->harray[i]; heap->harray[i] = swap_pointer;					
			heap->recent = i-1;
			i--;
		}
		heap->length++; 
	}
	else if (addnewelement ==0)
	{
		if (strand ==0) heap->harray[heap->recent]->counts[0] += 1; else heap->harray[heap->recent]->counts[1] +=1;
	}
	return 1;
}

// remove all indels in list that start before position | position = -1, -> remove all indels in heaps for new chromosome
void clean_indel_lists(struct BAMFILE_data* bamfiles_data,int bamfiles,int position)
{
	int i=0,b=0,k=0,t=0;
	for (b=0;b<bamfiles;b++)
	{
		heap = bamfiles_data[b].ilist;  
		if (position == -1) { heap->length = 0; heap->recent = 0; continue; }
		if (heap->length ==0) continue;
 
		i=0; while (i < heap->length && heap->harray[i]->position < position) i++;  
		if (i==0) continue;  // do nothing
		for (t=0,k=i;k<heap->length;k++) // swap 'k' element with t-th element
		{
			swap_pointer = heap->harray[t]; heap->harray[t] = heap->harray[k]; heap->harray[k] = swap_pointer; t++;
		}
		heap->length = t; heap->recent = 0; 
	}
	
}

void print_indels_interval(struct BAMFILE_data* bamfiles_data,int b,int start,int end)
{
	int i=0,k=0;
	heap = bamfiles_data[b].ilist; k=0;
	for (i=0;i<heap->length;i++)
	{
		if (heap->harray[i]->position < start || heap->harray[i]->position > end) continue;
		if (heap->harray[i]->counts[0] + heap->harray[i]->counts[1] < 2) continue;
		heap->harray[i]->printed = '1';
		fprintf(stdout,"%d:%d:%d,%d | ",heap->harray[i]->position,heap->harray[i]->length,heap->harray[i]->counts[0],heap->harray[i]->counts[1]);
		k++;
	}
	if (k > 0) fprintf(stdout,"sample:%d indels %d good_indels %d \n",b,heap->length,k);
} 

void print_indel_lists(struct BAMFILE_data* bamfiles_data,int bamfiles,int current_position)
{
	int i=0,b=0,k=0;
	for (b=0;b<bamfiles;b++)
	{
		heap = bamfiles_data[b].ilist; k=0;
		for (i=0;i<heap->length;i++)
		{
			if (heap->harray[i]->position > current_position) break; 
			if (heap->harray[i]->counts[0] + heap->harray[i]->counts[1] < 2) continue;
			heap->harray[i]->printed = '1';
			fprintf(stdout,"%d:%d:%d,%d | ",heap->harray[i]->position,heap->harray[i]->length,heap->harray[i]->counts[0],heap->harray[i]->counts[1]);
			k++;
		}
		if (k > 0) fprintf(stdout,"sample:%d indels %d pos %d good_indels %d \n",b,heap->length,current_position,k);
	}
	//fprintf(stdout,"\n");
} 

// if we have list of indels/complex variants, we can realign on the fly + output haplotypes instead of single variants 
// parse using fcigarlist to also get SNPs close to indels..
// this function is only called if the read has an indel in it
int extract_indel_reads(struct alignedread* read,REFLIST* reflist,int current,int sample,struct INDELHEAP* heap)
{
        int i=0,j=0,op=0,l=0,flag=0, l1=0,l2=0;
	int indel_added =0;
        for (i=0;i<read->fcigs;i++)
        {
                op = read->fcigarlist[i]&0xf; l = read->fcigarlist[i]>>4;
                if (op == BAM_CMATCH) { l1 += l; l2 += l; }
                else if (op == 8) 
		{
			//addElement(heap,read->position+l2,0,read->flag & 16);
			l1 += l; l2 += l; 
		} // substitution
                else if (op == BAM_CSOFT_CLIP) l1 += l;
                else if (op == BAM_CDEL)
                {
			addElement(heap,read->position+l2,-1*l,read->flag & 16); indel_added++;
                        l2 += l;
                }
                else if (op == BAM_CINS)
                {
                        flag =0;
                        if ((i ==1 && (read->fcigarlist[0]>>4) < 5) || (i ==read->fcigs-2 && (read->fcigarlist[read->fcigs-1]>>4) < 5) ) flag = 1;
                        if (flag ==0) 
			{	addElement(heap,read->position+l2,l,read->flag & 16); indel_added++;
	                        //for (j=0;j<read->fcigs;++j) fprintf(stdout,"%d%c:",read->fcigarlist[j]>>4,INT_CIGAROP[read->fcigarlist[j]&0xf]);
        	                //fprintf(stdout,"CAN_INS %d %dI S:%d flag %d i %d\n",read->position+l2,l,sample,flag,i);
			}
                        l1 += l;
                }
                else if (op == BAM_CHARD_CLIP) {}
                else if (op == BAM_CREF_SKIP) l2 += l;
        }
	return 1;
}


