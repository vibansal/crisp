#include "bamread.h"
// character array: 0->BAM_CMATCH, 1-> BAM_CINS .... 2 -> BAM_CDEL  3 -> BAM_CREF_SKIP 4 -> BAM_CSOFT_CLIP
char INT_CIGAROP[] = {'M','I','D','N','S','H','P','E','X'};

int LEFT_ALIGN_INDELS =0;
#include "left_align_indels.c"
//#include "../indels_shared/left_align_indels.c"

// we can use MDstring or XM tag to speed up this function, no need of reference base checking but not big difference in speed compared to other parts of code
int parse_cigar(struct alignedread* read,REFLIST* reflist,int current,int* fcigarlist)
{
	int i=0,t=0, l1=0,l2=0; int l=0,j=0;
	int f=0,m=0; int op;
	read->fcigs =0;
	read->alignedbases = 0; read->lastpos = read->position; read->cflag =0; read->mismatches=0; read->gaps = 0;
	read->cigoffset = 0; read->l1 =0; read->l2 = read->position; // l1 and l2 are indexes on read and reference sequence
	read->lastpos_mate = -1; 
	for (i=0;i<read->cigs;i++)
	{
		op = read->cigarlist[i]&0xf; l = read->cigarlist[i]>>4; 
		//if ((i ==0 || i == read->cigs-1) && op == BAM_CINS) { read->cigarlist[i] += 3; op = BAM_CSOFT_CLIP; } 
		// ignore hard clips as well (SOLID data can cause BUG)
		if (op != BAM_CMATCH && op != BAM_CHARD_CLIP) fcigarlist[f++] = read->cigarlist[i];  // copy as it is 
                if (op == BAM_CMATCH)
		{
			m=0;
			for (t=0;t<l;t++)
                        {
				// what if reference sequence is non-[ACTGN], ambiguous bases -> how to handle..
                                if (read->sequence[l1+t]  != reflist->sequences[current][read->position+l2+t] && read->sequence[l1+t] != reflist->sequences[current][read->position+l2+t]-32 && read->sequence[l1+t]  !='N') 
				{
					read->mismatches++;
					if (m > 0) fcigarlist[f++] = m<<4; fcigarlist[f++] = 24; m=0; // 24 = 1X in cigar code 
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
			//for (j=0;j<read->cigs;++j) fprintf(stdout,"%d%c:",read->cigarlist[j]>>4,INT_CIGAROP[read->cigarlist[j]&0xf]); 
			//fprintf(stdout,"INS %d %dI S:%d\n",read->position+l2,l,sample); 
			read->gaps++; l1 +=l; read->alignedbases += l; 
			//if (i==0) { read->l1 += l; read->cigoffset +=1; }  // if insertion or softclip is first op, ignore it
		}
		else if (op == BAM_CSOFT_CLIP)  l1 += l;  
		else if (op == BAM_CHARD_CLIP) {}
		else 
		{ 
			read->cflag =1; break;  // do not handle B and P cigars
		} 
	}
        read->last1 = l1; read->last2 = l2;
	
	if (f > 0) 
	{ 
		read->fcigs = f; read->fcigarlist = (int*)malloc(sizeof(int)*read->fcigs); 
		for (i=0;i<read->fcigs;i++) read->fcigarlist[i] = fcigarlist[i];  
		//for (i=0;i<read->cigs;++i) fprintf(stdout,"%d:%d ",read->cigarlist[i]>>4,read->cigarlist[i]&0xf); 
	}
	if (read->gaps > 0 && LEFT_ALIGN_INDELS ==1) left_align_indels(read,reflist, current); 

	//if (read->gaps > 0 && read->mquality >= 20) extract_indel_reads(read, reflist,current,sample);
	//free(fcigarlist); 
	return 0;
	//fprintf(stdout,"mismatches %d %d %s %d %d %s\n",read->mismatches,read->NM,read->sequence,current,read->position,read->cigar);
}

// b is full bam record, c is core alignment record
// data is pointer to bamfile
struct alignedread* get_read_bamfile(const bam1_t *b, void *data,struct alignedread* read)
{
	samfile_t *fp = (samfile_t*)data; uint32_t *cigar = bam1_cigar(b);  const bam1_core_t *c = &b->core;
	//	if (c->flag & (BAM_FUNMAP|BAM_FSECONDARY|BAM_FQCFAIL|BAM_FDUP)) return 0; /* skip unmapped reads  remove this filter */
	// CIGAR array: lower 4 bits encode CIGAR operation and higher 28 bits keep length of cigar 
	read = (struct alignedread*)malloc(sizeof(struct alignedread)); 
	//if (read == NULL) fprintf(stderr,"malloc error \n");
	read->readlength = b->core.l_qseq; 
	read->sequence = (char*)malloc(b->core.l_qseq+1); read->quality = (char*)malloc(b->core.l_qseq+1);
	read->matepair = NULL;

	// pos is 0-based leftmost coordinate 
	read->flag = c->flag; read->mquality= c->qual; read->position = c->pos; read->mateposition = c->mpos; read->IS = c->isize;
	read->strand = 0; if ((read->flag & 16) == 16) read->strand = 1; 
	read->fcigs =0;

	// ignore hard clip in read->cigarlist 
	read->cigs = c->n_cigar; int s=0,e=c->n_cigar,i=0;
	//if ((cigar[0]&0xf) == BAM_CHARD_CLIP) { read->cigs--; s++; }
	//if ((cigar[c->n_cigar-1]&0xf) == BAM_CHARD_CLIP) { read->cigs--; e--; } 
	read->cigarlist = calloc(read->cigs,sizeof(int));
	for (i =s; i <e;i++) read->cigarlist[i-s] = cigar[i]; 

	// change the first 'I' or last 'I' to soft clip 07/03/13 | this fixed a segfault using novoalign (12 pools)
	// move this to parse_cigar function 
	if ((read->cigarlist[0]&0xf) == BAM_CINS) read->cigarlist[0] += 3; 
	else if ((read->cigarlist[0]&0xf) == BAM_CHARD_CLIP && (read->cigarlist[1]&0xf) == BAM_CINS) read->cigarlist[1] += 3; 
	if ((read->cigarlist[read->cigs-1]&0xf) == BAM_CINS) read->cigarlist[read->cigs-1] += 3; 
	else if ((read->cigarlist[read->cigs-1]&0xf) == BAM_CHARD_CLIP && (read->cigarlist[read->cigs-2]&0xf) == BAM_CINS) read->cigarlist[read->cigs-2] += 3; 

	read->readid=(char*)malloc(c->l_qname+1); 
	for (i=0;i<c->l_qname;i++) read->readid[i] = b->data[i]; read->readid[i]= '\0';

	// added this special case for  HLA bams due to readid having extra /3 but not good  ingeneral
	//if (read->readid[c->l_qname-3] =='/') read->readid[i-3]='\0';
	//fprintf(stdout,"ql %d %s \n",c->l_qname,read->readid);

	if (c->tid >= 0) read->chrom = fp->header->target_name[c->tid]; else read->chrom = NULL; 
	read->tid = c->tid;
	if (c->mtid >= 0) read->matechrom = fp->header->target_name[c->mtid]; else read->matechrom = NULL;

	// sequence, quality XM, NM 	
	uint8_t* sequence = bam1_seq(b); uint8_t* quality = bam1_qual(b);
	for (i=0;i<b->core.l_qseq;i++) read->sequence[i] = bam_nt16_rev_table[bam1_seqi(sequence,i)]; read->sequence[i] = '\0';
	for (i=0;i<b->core.l_qseq;i++) read->quality[i] = (char)(quality[i]+33); read->quality[i] = '\0';

	read->XC = 0;
	// only do this for BWA aligned reads and those that have soft clip in cigar  
        //uint8_t* cmstring = bam_aux_get(b,"XC"); if (cmstring != NULL) read->XC = read->readlength - *(cmstring+1);


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
	//fprintf(stdout,"read %s %d \n",read->readid,read->position);
	//	printf(" qname %s ref %s position %d mate %s matepos %d IS %d mappingq %d flag %d\n",bam1_qname(b),fp->header->target_name[c->tid],c->pos+1,fp->header->target_name[c->mtid],c->mpos+1,c->isize,c->qual,c->flag);
	return read; // return pointer containing address of new memory allocated 
}

void free_read(struct alignedread* read)
{
	if (read->matepair != NULL) (read->matepair)->matepair = NULL; 
	free(read->sequence);free(read->quality);free(read->readid); free(read->cigarlist); free(read);
	if (read->fcigs > 0) free(read->fcigarlist);
}
