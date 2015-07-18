
#include "bamsreader.h"
int seqplatform = 0; // 0 is for Illumina, 1 is for 454 reads

void removefromHeap(BAMHEAP* heap,struct BAMFILE_data* bamfile_data)
{
	// element is at position '0' in heap 
	heap->harray[0] = heap->harray[heap->length-1]; heap->length--; 
	minHeapify(heap,0,bamfile_data);
}

// trickledown can also be used to update the BAMHEAP if the score of a node is decreased via an update 
// (chromosome,position) of a node is increased, we need to trickle 
void minHeapify(BAMHEAP* heap,int node,struct BAMFILE_data* bamfiles_data)
{
	int lc = 2*node+1; int rc = 2*node + 2; int minindex = node; int temp; 
	if (rc >= heap->length)
	{
		if (lc < heap->length) minindex = lc; 
		// else case is equal to minindex = node
	}
	else
	{
		if (bamfiles_data[heap->harray[lc]].read->tid == bamfiles_data[heap->harray[rc]].read->tid)
		{
			if (bamfiles_data[heap->harray[lc]].read->position <= bamfiles_data[heap->harray[rc]].read->position) minindex = lc; 
			else minindex = rc; 
		}
		else if (bamfiles_data[heap->harray[lc]].read->tid < bamfiles_data[heap->harray[rc]].read->tid) minindex = lc;
		else minindex = rc;
	}
	if (bamfiles_data[heap->harray[node]].read->tid > bamfiles_data[heap->harray[minindex]].read->tid || (bamfiles_data[heap->harray[node]].read->tid == bamfiles_data[heap->harray[minindex]].read->tid && bamfiles_data[heap->harray[node]].read->position > bamfiles_data[heap->harray[minindex]].read->position ))
	{
		temp = heap->harray[node]; heap->harray[node] = heap->harray[minindex]; heap->harray[minindex] = temp;
		minHeapify(heap,minindex,bamfiles_data);
	}
}

void buildminheap(BAMHEAP* heap,struct BAMFILE_data* bamfiles_data)
{
	int i=0; 
	for (i=heap->length/2-1;i>=0;i--)  minHeapify(heap,i,bamfiles_data);
	//fprintf(stdout,"heapify %d hl %d\n",i,heap->length);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////

void printpileup(REFLIST* reflist,int current,int k,READQUEUE* bq)
{
	printf("\n");
}

int free_read_queue(READQUEUE* bq,struct BAMFILE_data* bamfiles_data) //free read pointed to by RQ->iter, previous is RQ->previous
{
	// update queue pointers for the pool 'i'
	//fprintf(stdout,"free %d\n",bq->reads);
	int i=bq->iter->sampleid;
	//adding the line below
	if(bq->iter==bamfiles_data[i].previous)bamfiles_data[i].previous=NULL;
	if(bq->iter == bamfiles_data[i].last) bamfiles_data[i].last = bamfiles_data[i].previous;
	if(bq->iter == bamfiles_data[i].first) // first read in queue for pool 'i'
	{
		bamfiles_data[i].first =bq->iter->nextread;
		//andthislineworks
		bamfiles_data[i].previous=NULL;
	}
	else
	{
		//if(bamfiles_data[i].previous==NULL) fprintf(stdout,"error bq-iter\n");
		(bamfiles_data[i].previous)->nextread =(bq->iter)->nextread;
	}

	if (bq->iter == bq->last) bq->last= bq->previous; // special case 
	if (bq->iter == bq->first)
	{
		bq->first = bq->iter->next; // set the first read in queue 
		// add NULL check here but it would mask BUGS in code
		//for (i=0;i<bq->iter->fcigs;i++) fprintf(stdout,"%d%c ",bq->iter->fcigarlist[i]>>4,INT_CIGAROP[bq->iter->fcigarlist[i]&0xf]);
		//fprintf(stdout,"%s ",bq->iter->sequence);
                //fprintf(stdout,"readqueuef s:%d IS:%d %s %d r:%d \n",bq->iter->sampleid,bq->iter->IS,bq->iter->readid,bq->iter->position,bq->iter->mismatches);
		if ((bq->iter)->matepair != NULL) ((bq->iter)->matepair)->matepair = NULL; 
		free(bq->iter->sequence); free(bq->iter->quality); free(bq->iter->cigarlist);  free(bq->iter->readid);
		if (bq->iter->fcigs >0) free(bq->iter->fcigarlist);
		free(bq->iter);  // free current read
		bq->iter = bq->first; 
	}
	else // handles all other cases
	{
		bq->previous->next = bq->iter->next;	
		if ((bq->iter)->matepair != NULL) ((bq->iter)->matepair)->matepair = NULL; 
		free(bq->iter->sequence); free(bq->iter->quality); free(bq->iter->cigarlist); free(bq->iter->readid);
		if (bq->iter->fcigs >0) free(bq->iter->fcigarlist);
		free(bq->iter);
		bq->iter = bq->previous->next;  
	}
	bq->reads--;
}

void empty_queue(READQUEUE* bq,struct BAMFILE_data*  bamfiles_data)
{
        fprintf(stderr,"cleaning read queue from prev chrom\n");
        bq->iter = bq->first; bq->previous = NULL;
        while (bq->iter != NULL) 
        {
                free_read_queue(bq,bamfiles_data);
        }
	bq->last = NULL; bq->first = NULL;
}

int allocate_mem_heap(struct BAMFILE_data* bamfiles_data,int bamfiles,int heapsize)
{
        int i=0,j=0;
        for (i=0;i<bamfiles;i++)
        {
		bamfiles_data[i].ilist = calloc(1,sizeof(struct INDELHEAP));
                (bamfiles_data[i].ilist)->harray = calloc(heapsize,sizeof(struct bamINDEL*));
                (bamfiles_data[i].ilist)->maxsize = heapsize; (bamfiles_data[i].ilist)->length = 0;
		(bamfiles_data[i].ilist)->fve = 0;
		for (j=0;j<heapsize;j++) (bamfiles_data[i].ilist)->harray[j] = NULL; //calloc(1,sizeof(struct bamINDEL));
        }
}


int init_bamfiles(struct BAMFILE_data* bamfiles_data,char** bamfilelist,int bamfiles,char* regions,int* reference,int* startpos,int* endpos)
{
	int ref=-1,beg=0,end=0,i=0,finishedfiles=0;
        for (i=0;i<bamfiles;i++)
        {
                bamfiles_data[i].first=NULL; bamfiles_data[i].last=NULL;bamfiles_data[i].previous=NULL;
                if ((bamfiles_data[i].fp = samopen(bamfilelist[i], "rb", 0)) == 0) 
                {
                        fprintf(stderr, "unable to open the BAM file %s\n", bamfilelist[i]);  exit(0);
                }
                if (regions == NULL) 
                {
                        bamfiles_data[i].b = bam_init1();
                        if (samread(bamfiles_data[i].fp,bamfiles_data[i].b) >=0)
                        {
                                bamfiles_data[i].read = get_read_bamfile(bamfiles_data[i].b,bamfiles_data[i].fp,bamfiles_data[i].read);
                                bamfiles_data[i].read->sampleid = i;
                        }
                        else
                        {
                                bamfiles_data[i].finished = 1; finishedfiles++; bamfiles_data[i].read= NULL;
                        }
                }
		else // load BAM index and find region 
                {
                        if ( (bamfiles_data[i].idx  = bam_index_load(bamfilelist[i])) ==0) 
                        {
                                fprintf(stderr,"unable to load bam index for file %s\n",bamfilelist[i]); exit(0);
                        }
                        bam_parse_region(bamfiles_data[i].fp->header,regions,&ref,&beg,&end); 
                        if (ref < 0)
                        {
                                fprintf(stderr,"invalid region for bam file %s \n",regions); exit(0);
                        }
                        bamfiles_data[i].b = bam_init1();
                        bamfiles_data[i].iter = bam_iter_query(bamfiles_data[i].idx,ref,beg,end);
			//if (i==0) fprintf(stderr,"region argument %s %d:%d:%d\n",regions,ref,beg,end);
                        if ( (bamfiles_data[i].ret = bam_iter_read(bamfiles_data[i].fp->x.bam,bamfiles_data[i].iter,bamfiles_data[i].b)) >=0)
                        {
                                if (i==0 || i+1 == bamfiles) fprintf(stderr,"reading bam file index: %s.bai region %d:%d-%d\n",bamfilelist[i],ref,beg,end);
                                bamfiles_data[i].read = get_read_bamfile(bamfiles_data[i].b,bamfiles_data[i].fp,bamfiles_data[i].read);
                                bamfiles_data[i].read->sampleid = i;

                        }
                        else
                        {
                                bamfiles_data[i].finished = 1; finishedfiles++; bamfiles_data[i].read= NULL;
                                bam_iter_destroy(bamfiles_data[i].iter); bam_destroy1(bamfiles_data[i].b); bam_index_destroy(bamfiles_data[i].idx);
                        }
                }
        }
	*reference = ref; *startpos = beg; *endpos = end;
	fprintf(stderr,"finished reading bam file indexes \n");
	int j=0,k=0;
	for (i=0;i<bamfiles;i++)	
	{
		bamfiles_data[i].counts = calloc(64,sizeof(int*)); for (j=0;j<64;j++) bamfiles_data[i].counts[j] = calloc(8,sizeof(int));
		for (j=0;j<64;j++) { for (k=0;k<8;k++) bamfiles_data[i].counts[j][k] = 0; } 
	}

	return 1;
}


	/*
	int k=0; int flag =0;
        while (read->MDstring[k] != '\0')
        {
                if (read->MDstring[k] == '^') flag =1;
		else if(flag ==1 && (int)read->MDstring[k] >=65 && (int)read->MDstring[k] <=84) flag=1;
		else if (flag ==0 && (int)read->MDstring[k] >=65 && (int)read->MDstring[k] <=84) mismatches++;
                else flag =0;
                k +=1;
        }
//	fprintf(stdout,"MDstring %d %s %d\n",read->position,read->MDstring,k);
	*/
