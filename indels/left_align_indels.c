
// basic function to left align indels in fcigarlist generated after parsing read from bam file 
//
// works for single del or ins... implemented 07/15/13
// read->cigarlist will not match fcigarlist after this function is called 
int left_align_indels(struct alignedread* read,REFLIST* reflist,int current)
{
	int i=0,j=0,op=0,l=0,l1=0,l2=0,j1=0,j2=0,op_prev = 0,l_prev =0;
	//for (j=0;j<read->fcigs;++j) fprintf(stdout,"%d%c:",read->fcigarlist[j]>>4,INT_CIGAROP[read->fcigarlist[j]&0xf]); 
	
	// iterate over fcigarlist 
	for (i=0;i<read->fcigs;i++)
	{
		op = read->fcigarlist[i]&0xf; l = read->fcigarlist[i]>>4; 

		if (op == BAM_CMATCH || op == 8) { l1 += l; l2 += l;  }
		else if (op == BAM_CSOFT_CLIP) l1 += l;
		else if (op == BAM_CREF_SKIP) l2 += l;
		else if (op == BAM_CHARD_CLIP) { }

		else if (op == BAM_CDEL)
		{
			j1=0; j2 = l-1;  // last base of deletion 
			while (read->sequence[l1-1-j1] == reflist->sequences[current][read->position+l2+j2] && j1 < l_prev-1) 
			{ 
				//fprintf(stdout,"j1 %d %d %d %c:%c %d\t",j1,l1-1-j1,read->position+j2,read->sequence[l1-1-j1],reflist->sequences[current][read->position+j2],l_prev-1);
				j1++; j2--; 
			}
			// previous cigar has to BAM_CMATCH, next cigar need not be BAM_CMATCH ! if it is something else -> change fcigarlist
			if (j1 > 0 && op_prev == BAM_CMATCH && i < read->fcigs-1 && (read->fcigarlist[i+1]&0xf) == BAM_CMATCH) 
			{
				read->fcigarlist[i-1] -= j1<<4; read->fcigarlist[i+1] += j1<<4; 
				l2 += l-j1;		
			}
			else l2 += l;
		}
		else if (op == BAM_CINS)
		{
			j1=0; j2 = l-1;  // last base of deletion 
			while (read->sequence[l1-1-j1] == read->sequence[l1+j2] && j1 < l_prev-1) 
			{ 
				j1++; j2--; 
			}
			if (j1 > 0 && op_prev == BAM_CMATCH && i < read->fcigs-1 && (read->fcigarlist[i+1]&0xf) == BAM_CMATCH) 
			{
				read->fcigarlist[i-1] -= j1<<4; read->fcigarlist[i+1] += j1<<4; 
				l2 += l-j1;		
			}
			else l2 += l;
		}
		op_prev = op; l_prev = l; 
	}
	//for (j=0;j<read->fcigs;++j) fprintf(stdout,"%d%c:",read->fcigarlist[j]>>4,INT_CIGAROP[read->fcigarlist[j]&0xf]); fprintf(stdout," %s %d \n",read->readid,read->position);
	
}

