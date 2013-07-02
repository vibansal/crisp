#include "bamsreader.h"
int seqplatform = 0; // 0 is for Illumina, 1 is for 454 reads

// june 1 2012, we can use MDstring and cigar to speed up this function, no need of reference base checking 
// this parse_cigar is used for bamtopileup function for calcuating the number of mismatches in an aligned read, maybe we should rename it as such 

// character array: 0->BAM_CMATCH, 1-> BAM_CINS ....
char INT_CIGAROP[] = {'M','I','D','N','S','H','P','E','X'};

int parse_cigar(struct pileupread* read,REFLIST* reflist,int current,int* fcigarlist)
{
	int i=0,t=0, l1=0,l2=0; int l=0;
	int f=0,m=0; int op;
	//int* fcigarlist = (int*)malloc(sizeof(int)*read->cigs+sizeof(int)*(int)(read->readlength)); 
	read->fcigs =0;
	read->alignedbases = 0; read->lastpos = read->position; read->cflag =0; read->mismatches=0; read->gaps = 0;
	read->cigoffset = 0; read->l1 =0; read->l2 = read->position; // l1 and l2 are indexes on read and reference sequence
	read->lastpos_mate = -1; 
	for (i=0;i<read->cigs;i++)
	{
		op = read->cigarlist[i]&0xf; l = read->cigarlist[i]>>4; 
		// ignore hard clips as well (SOLID data can cause BUG)
		if (op != BAM_CMATCH && op != BAM_CHARD_CLIP) fcigarlist[f++] = read->cigarlist[i]; 
                if (op == BAM_CMATCH)
		{
			m=0;
			for (t=0;t<l;t++)
                        {
                                if (read->sequence[l1+t]  != reflist->sequences[current][read->position+l2+t] && read->sequence[l1+t] != reflist->sequences[current][read->position+l2+t]-32 && read->sequence[l1+t]  !='N') 
				{
					read->mismatches++;
					if (m > 0) fcigarlist[f++] = m<<4; fcigarlist[f++] = 24; m=0;
					//if (m > 0) fcigarlist[f++] = m<<4; fcigarlist[f++] = (read->sequence[l1+t]<<4)+8; m=0;
				}
				else m++; 
                        }
			if (m > 0) fcigarlist[f++] = m<<4; 

			l1 += l; l2 +=l; read->alignedbases += l; read->lastpos += l;
		}
		else if (op == BAM_CDEL || op == BAM_CREF_SKIP) // allow 'N' in cigar  
		{
			read->gaps++; l2 += l; read->lastpos += l; 
		}
		else if (op == BAM_CINS) 
		{ 
			read->gaps++; l1 +=l; read->alignedbases += l; 
			//if (i==0) { read->l1 += l; read->cigoffset +=1; }  // if insertion or softclip is first op, ignore it
		}
		else if (op == BAM_CSOFT_CLIP) 
		{
			l1 += l;  
			//if (i==0) { read->l1 += l; read->cigoffset +=1; } 
		}
		else if (op == BAM_CHARD_CLIP) {}
		else { read->cflag =1; break; } 
	}
	char indelhaplotype_L[120]; strcpy(indelhaplotype_L,"TATGCAAATACAAAATACATATGACAAAAATACA");
	char* strfound;
	
	if (f > 0) 
	{ 
		read->fcigs = f;
		read->fcigarlist = (int*)malloc(sizeof(int)*read->fcigs); 
		for (i=0;i<read->fcigs;i++) read->fcigarlist[i] = fcigarlist[i];  
		//for (i=0;i<read->cigs;++i) fprintf(stdout,"%d:%d ",read->cigarlist[i]>>4,read->cigarlist[i]&0xf); 
		if (f >= 2000)  
		{
			strfound = strstr(read->sequence,indelhaplotype_L);
			if (strfound != NULL) fprintf(stdout,"%s %d %s FOUND %ld ",read->sequence,read->position,read->readid,strfound-read->sequence);
			else fprintf(stdout,"%s %d %s NULL ",read->sequence,read->position,read->readid);
			for (i=0;i<f;++i) fprintf(stdout,"%d%c ",read->fcigarlist[i]>>4,INT_CIGAROP[read->fcigarlist[i]&0xf]); 
			fprintf(stdout,"\n");
		}
		/*
		*/
	}
	//free(fcigarlist); 
	return 0;
	//fprintf(stdout,"mismatches %d %d %s %d %d %s\n",read->mismatches,read->NM,read->sequence,current,read->position,read->cigar);
}

// b is full bam record, c is core alignment record
// data is pointer to bamfile
struct pileupread* get_read_bamfile(const bam1_t *b, void *data,struct pileupread* read)
{
	samfile_t *fp = (samfile_t*)data; uint32_t *cigar = bam1_cigar(b);  const bam1_core_t *c = &b->core;
	//	if (c->flag & (BAM_FUNMAP|BAM_FSECONDARY|BAM_FQCFAIL|BAM_FDUP)) return 0; /* skip unmapped reads  remove this filter */
	// CIGAR array: lower 4 bits encode CIGAR operation and higher 28 bits keep length of cigar 
	read = (struct pileupread*)malloc(sizeof(struct pileupread)); 
	//if (read == NULL) fprintf(stderr,"malloc error \n");
	read->readlength = b->core.l_qseq; 
	read->sequence = (char*)malloc(b->core.l_qseq+1); read->quality = (char*)malloc(b->core.l_qseq+1);
	read->matepair = NULL;

	// pos is 0-based leftmost coordinate 
	read->flag = c->flag; read->mquality= c->qual; read->position = c->pos; read->mateposition = c->mpos; read->IS = c->isize;
	read->strand = 0; if ((read->flag & 16) == 16) read->strand = 1; 

	read->cigarlist = (int*)malloc(sizeof(int)*c->n_cigar);
	read->cigs =c->n_cigar; read->fcigs =0;
	int i=0;
	for (i =0; i < c->n_cigar; ++i) read->cigarlist[i] = cigar[i]; 

	read->readid=(char*)malloc(c->l_qname+1); 
	for (i=0;i<c->l_qname;i++) read->readid[i] = b->data[i]; read->readid[i]= '\0';

	// added this special case for  HLA bams due to readid having extra /3 but not good  ingeneral
	if (read->readid[c->l_qname-3] =='/') read->readid[i-3]='\0';
	//fprintf(stdout,"ql %d %s \n",c->l_qname,read->readid);

	if (c->tid >= 0) read->chrom = fp->header->target_name[c->tid]; else read->chrom = NULL; 
	read->tid = c->tid;
	if (c->mtid >= 0) read->matechrom = fp->header->target_name[c->mtid]; else read->matechrom = NULL;

	// sequence, quality XM, NM 	
	uint8_t* sequence = bam1_seq(b); uint8_t* quality = bam1_qual(b);
	for (i=0;i<b->core.l_qseq;i++) read->sequence[i] = bam_nt16_rev_table[bam1_seqi(sequence,i)]; read->sequence[i] = '\0';
	for (i=0;i<b->core.l_qseq;i++) read->quality[i] = (char)(quality[i]+33); read->quality[i] = '\0';

	// c->tid is integer pointing to reference fasta to which read is mapped 
	//fp->header is header of bam file 
	// fp->header->target_name gives actual name 'chrom6' of read mapping 
	//uint8_t* auxp = bam1_aux(b); 
	if ( !(read->flag & 4) && SOLID ==1)
	{
		//uint8_t* auxp = bam_aux_get(b,"MD"); 
		//if (auxp != NULL) read->MDstring = auxp+1; else read->MDstring = NULL; 
		uint8_t* cmstring = bam_aux_get(b,"CM"); 
		if (cmstring != NULL)
		{
			read->CM = *(cmstring+1); 
			//for (i=0;i<read->cigs;i++) fprintf(stdout,"%d%c",read->cigarlist[i]>>4,INT_CIGAROP[read->cigarlist[i]&0xf]);
			//fprintf(stdout," read %s %d %d\n",read->readid,read->position,read->CM);
		}
		else read->CM =0;
		//fprintf(stdout," read %s %d %s\n",read->readid,read->position,auxp+1);
	}
	/*
	*/
	//fprintf(stdout,"read %s %d \n",read->readid,read->position);
	//	printf(" qname %s ref %s position %d mate %s matepos %d IS %d mappingq %d flag %d\n",bam1_qname(b),fp->header->target_name[c->tid],c->pos+1,fp->header->target_name[c->mtid],c->mpos+1,c->isize,c->qual,c->flag);
	//	free(read->sequence);free(read->quality); free(read->cigarlist);
	return read;
}

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

void free_read(struct pileupread* read)
{
	if (read->matepair != NULL) (read->matepair)->matepair = NULL; 
	free(read->sequence);free(read->quality);free(read->readid); free(read->cigarlist); free(read);
	if (read->fcigs > 0) free(read->fcigarlist);
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
