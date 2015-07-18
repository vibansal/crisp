/*
find best alignment of read upto indel position and then align remainder to read to indel haplotype...
best alignment of read upto indel may not be fully consistent with cigar of read...

if we extend alignment into 'S' region, we should also original alignment and compare... | calculate XC tag for BWA 

1. to fix, no 'H' cigar in read->cigarlist... since we assume last op to be 'S' only...
 variables read->last1 and read->last2 need to be set before this function is called...
read->XC value shhould be set before realignment...

code for generating new cigar using indel: algorithm description 
1. find all potential indels that overlap the aligned/partially aligned read (variant index)
2. for each indel, evaluate if adding the indel into the existing cigar of the read reduces the # of mismatches or extends into the clipped part of the read 
3. the above step is done for both directions where indel is inserted L->R and R->L 
*/

//char INT_CIGAROP[] = {'M','I','D','N','S','H','P','E','X'};


int realignread_LR(struct alignedread* read,REFLIST* reflist,struct bamINDEL* bam_indel,int var,struct REALIGNMENT* rl)
{
	rl->cigs =0; //initilze
	// extend partially mapped reads to calculate total score without clips, to avoid making spurios indel realignments...
	int i=0,j=0,l1=0,l2=0,t=0,matches =0, carryover = 0,dl=0,il=0;
	int l11=0,l21=0; int added = -1; int MM[3] = {0,0,0};  // pair of new mismatches,old mismatches before indel, new mismatches after indel event
	int op=0,ol=0,p = 0; int indeladded =0; int b =0;
	//int delta =0; if (bam_indel->length < 0) delta = -1*bam_indel->length; 
	int lastposition = read->readlength; // lastposition in read, if read is clipped due to poor quality bases, lastposition < readlength
	op = read->cigarlist[read->cigs-1]&0xf;  ol = read->cigarlist[read->cigs-1]>>4;
	if ( (read->flag &16) ==0 && op == BAM_CSOFT_CLIP) lastposition = read->readlength-read->XC; 
	bam_indel->simple = '1';
	//if (read->readlength > 90) fprintf(stdout,"%c char \n",read->sequence[64]);
	
	// check read going from left to right 
	for (i=0;i<read->cigs;i++)
	{
	        op = read->cigarlist[i]&0xf; ol = read->cigarlist[i]>>4; p = ol;
		// only move into soft clip region if previous cigar is not 'M' extra condition added jan 2 2013 
		if (op == BAM_CSOFT_CLIP && i  > 0 && rl->cigs > 0 && (rl->cigarlist[rl->cigs-1]&0xf) ==0 && (rl->cigarlist[rl->cigs-1]>>4) >=10 && indeladded > 0)  
		{
			//fprintf(stdout,"ignoring soft clip region \t");
			rl->cigarlist[rl->cigs++] = read->cigarlist[i]; l1 += ol; l11 += ol; continue;
		}
		else if (op == BAM_CMATCH || (op == BAM_CSOFT_CLIP && i > 0) || (op == BAM_CINS && indeladded > 0) || (op == BAM_CINS  && read->position + l2 == bam_indel->position && bam_indel->length > 0 && indeladded ==0) )
		{
			//if (op == BAM_CINS && indeladded > 0) MM[1]++;  
			for (t=0;t<ol;t++)
			{
				if (read->position + l2 == bam_indel->position && indeladded ==0)  // indel insertion match in 'M' stretch
				{
					// here 't > 0' is necessary otherwise p=0, indel cannot happen | maybe not check this 06/14/13
					if (bam_indel->length < 0 && bam_indel->simple == '1' && (rl->cigs > 0 || t > 0))
					// position in read matches the known position of deletion event 
					{
						p = ol-t; 
						if (t > 0 && rl->cigs ==0) rl->cigarlist[rl->cigs++] = t<<4; 
						else if (t > 0)
						{
							b = rl->cigarlist[rl->cigs-1]&0xf;
							if (b ==0) rl->cigarlist[rl->cigs-1] += t<<4; else rl->cigarlist[rl->cigs++] = t<<4;
						}
						rl->cigarlist[rl->cigs++] = ((-1*bam_indel->length)<<4) + 2; 
						l2 += -1*bam_indel->length; added = 0; carryover = 0; indeladded = 1;
						printf("LR_del->%d %d ",l1,p);
					}
					else if (bam_indel->length > 0 && bam_indel->simple == '1' && (rl->cigs > 0 || t > 0))
					{ 
						p = ol-t - bam_indel->length; 
						if (p < 0) carryover = -1*p; else carryover = 0;
						if (t > 0 && rl->cigs== 0) rl->cigarlist[rl->cigs++] = t<<4;
						else if (t > 0)
						{
							b = rl->cigarlist[rl->cigs-1]&0xf;
							if (b ==0) rl->cigarlist[rl->cigs-1] += t<<4; else rl->cigarlist[rl->cigs++] = t<<4;
						}
						rl->cigarlist[rl->cigs++] = (bam_indel->length<<4) + 1;
						for (j=0;j<bam_indel->length && j < ol-t;j++)  // count mismatches of inserted bases to known insertion allele
						{
							//if (read->sequence[l1+j]  != bam_indel->bases[j+1]) MM[2]++;  
						}
						printf("ins->%d %d %d %c carryover %d ",l1,p,t,INT_CIGAROP[op],carryover);
						l1 += bam_indel->length; t += bam_indel->length;
						// when t is incremented by more than is covered by this segment, then something needs to subtracted from next segment.....
						added = 0; indeladded =1;
					}
					else if (bam_indel->simple == '0' && (rl->cigs > 0 || t > 0)) // complex indel still not fully fixed 
					{ 
						//dl = strlen(bam_indel[var].RA)-1; il = strlen(bam_indel[var].AA)-1;
						p = ol-t; p -= il; 
						if (p < 0) carryover = -1*p; else carryover = 0;
						if (t > 0 && rl->cigs== 0) rl->cigarlist[rl->cigs++] = t<<4;
						else if (t > 0)
						{
							b = rl->cigarlist[rl->cigs-1]&0xf;
							if (b ==0) rl->cigarlist[rl->cigs-1] += t<<4; else rl->cigarlist[rl->cigs++] = t<<4;
						}
						rl->cigarlist[rl->cigs++] = (dl<<4) + 2; rl->cigarlist[rl->cigs++] = (il<<4) + 1;
						for (j=0;j<il && j < ol-t;j++)  // count mismatches of inserted bases to known insertion allele
						{
							if (read->sequence[l1+j]  != bam_indel->bases[j+1]) MM[2]++;  
						}
						printf("complexindel->%d %d %d %c carryover %d ",l1,p,t,INT_CIGAROP[op],carryover);
						l1 += il; t += il; l2 += dl;
						// when t is incremented by more than is covered by this segment, then something needs to subtracted from next segment.....
						added = 0; indeladded =1;
					}
				}

				//fprintf(stdout,"%c:%c:%d:%d %s \n",read->sequence[l11+t],reflist->sequences[read->tid][read->position+l21+t],t,l11,read->readid);
				if (t >=ol) break;
				if ( (op == BAM_CMATCH || op == BAM_CSOFT_CLIP || op == BAM_CINS) && carryover > 0) { carryover--; p--; continue; } 
				if (l1 >= lastposition) { break;} // important check
				if (added >= 0 && op == BAM_CSOFT_CLIP) added++; // variant has been added and now into clipped region 

				// if t >= ol because of insertion -> these won't be counted jan 1 2013 
				//if (read->sequence[l1]  != reflist->sequences[read->tid][read->position+l2] && indeladded ==0) MM[1]++;
				if (read->sequence[l11+t]  != reflist->sequences[read->tid][read->position+l21+t] && indeladded ==0 && op != BAM_CINS) MM[0]++;

				if (read->sequence[l1]  != reflist->sequences[read->tid][read->position+l2]) 
				{
					if (op == BAM_CSOFT_CLIP && added >= 5 && carryover <=0 && bam_indel->length < 0) 
					{ 
						t = ol; added--;   // we should break here rather than continue since cigar 
					} // get out
					else if (indeladded > 0 || op != BAM_CMATCH) MM[2]++; 
				}
				else matches++; 
				l1++; l2++; 
			}
			if (added < 0)  // no indel was inserted into newcigarlist
			{ 
				if (rl->cigs > 0 && (rl->cigarlist[rl->cigs-1]&0xf) ==0) rl->cigarlist[rl->cigs-1] += ol<<4;
				else rl->cigarlist[rl->cigs++] = read->cigarlist[i]; 
			}
			else if (p > 0 && (op == BAM_CMATCH || op == BAM_CINS)) // left over 'M' after indel insertion  
			{
				if (rl->cigs > 0 && (rl->cigarlist[rl->cigs-1]&0xf) ==0) rl->cigarlist[rl->cigs-1] += p<<4;
                                else rl->cigarlist[rl->cigs++] = p<<4; 
			}
			else if (carryover ==0 && op == BAM_CSOFT_CLIP) 
			{ 
				b = rl->cigarlist[rl->cigs-1]&0xf;
				if (added > 0 && b != 0) rl->cigarlist[rl->cigs++] = added<<4; 
				else if (added > 0) rl->cigarlist[rl->cigs-1] += added<<4; 
				if (p-added > 0) rl->cigarlist[rl->cigs++] = ((p-added)<<4) + 4; // remaining portion of soft clip
			}
			else if (carryover > 0 && op == BAM_CSOFT_CLIP) 
			{ 
				rl->cigarlist[rl->cigs++] = (ol-carryover)<<4; carryover =0;
			}
			if (op == BAM_CINS) l11 += ol; else { l11 += ol; l21 += ol;} 
		}
		else if (op == BAM_CDEL) 
		{ 
			if (indeladded > 0) { MM[1]++; } // already added the indel event we are evaluating, so ignore extra deletion event
			else if (read->position + l2 + ol < bam_indel->position ) // added 07/17/13, deletion is ignored if the new variant being evaluated falls within the interval deleted on reference 
			{
				if (carryover > 0) {  rl->cigarlist[rl->cigs-1] -= carryover<<4; carryover = 0; } // probably never invoked
				rl->cigarlist[rl->cigs++] = read->cigarlist[i]; 
				l2 += ol;  
			}
			l21 += ol; 
		}
		else if (op == BAM_CINS) // add filter on insertion being close to variant being evaluated...  07/17/13
		{ 
			rl->cigarlist[rl->cigs++] = read->cigarlist[i]; l1 += ol; l11 += ol; 
		}
		else if (op == BAM_CSOFT_CLIP) 
		{ 
			rl->cigarlist[rl->cigs++] = read->cigarlist[i]; l1 += ol; l11 += ol; 
		}
	}
	// if the previous cigarlength is less than carryover -> negative BUG 
	if (carryover > 0)  
	{
	       	op = rl->cigarlist[rl->cigs-1]&0xf; ol = rl->cigarlist[rl->cigs-1]>>4; 
 		if (op != BAM_CDEL && ol > carryover)  rl->cigarlist[rl->cigs-1] -= carryover<<4; 
		else return 0;
	}

	// MM[0] is number of mismatches prior to indel event, MM[1] -> gaps
	// MM[2] is number of mismatches introduced after indel event is added to new cigar
	//fprintf(stdout,"newciglength %d \n",rl->cigs);
	if ( (added ==0 && MM[2] < read->mismatches-MM[0]+MM[1]-1 && MM[2] <2) || (added == 0 && MM[2] < read->mismatches-MM[0] && MM[2] < 3) || (added >= 2 && MM[2] < 1) || (added >= 5 && MM[2] <2)) 
	{
		rl->newpos = read->position; rl->added = added; rl->mismatches = MM[2]; rl->delta = read->mismatches-MM[0]-MM[2]; 
		if (MM[1] > 0) fprintf(stdout,"newflag ");
		printf("newms:%d:%d:%d MM %d:%d added %d | NEWCIGAR_LR: ",MM[2],read->mismatches-MM[0],MM[1],MM[2],MM[0],added);
		for (i=0;i<rl->cigs;i++) fprintf(stdout,"%d%c ",rl->cigarlist[i]>>4,INT_CIGAROP[rl->cigarlist[i]&0xf]); //fprintf(stdout,"NEW "); 
		for (i=0;i<read->fcigs;i++) fprintf(stdout,"%d%c",read->fcigarlist[i]>>4,INT_CIGAROP[read->fcigarlist[i]&0xf]);
		fprintf(stdout," %d %s %s\n",read->position,read->sequence,read->readid);
		return 1;
	}
	return 0;
}

// what about 'D' events that are close to the variant we are going to evaluate, try with ignoring all 'deletions' separaterly..
// once indel event is added to new cigar, insertions are treated as 'M' and 'D' are ignored 
// indel realignment when indel is near beginning of read 
int realignread_RL(struct alignedread* read,REFLIST* reflist,struct bamINDEL* bam_indel,int var,struct REALIGNMENT* rl)
{
	rl->cigs =0; //initilze
	int i=0,j=0,l1=0,l2=0,t=0,matches =0,carryover = 0,il=0,dl=0;
	int l11=read->last1-1,l21=read->last2-1; int added = -1; int MM[3] = { 0,0,0};
	int op=0,ol=0,p = 0; int indeladded =0; int b =0;  
	l1 = read->last1-1; l2 = read->last2-1; 
	int delta =0; if (bam_indel->length < 0) delta = -1*bam_indel->length;  //if (bam_indel[var].simple == '0') delta = strlen(bam_indel[var].RA)-1;

	//for (i=0;i<read->fcigs;i++) fprintf(stdout,"%d%c",read->fcigarlist[i]>>4,INT_CIGAROP[read->fcigarlist[i]&0xf]); printf( " %s ",read->sequence);
	//printf(" %s readpos %d indel %d length %d ct %d end %d \n",read->readid,read->position,bam_indel->position,bam_indel->length,bam_indel->counts[0],read->span+read->position);
	bam_indel->simple = '1';
	
	int firstposition = 0; op = read->cigarlist[0]&0xf;  ol = 0>>4;
	if ( (read->flag &16) ==16 && op == BAM_CSOFT_CLIP) firstposition = read->XC+1;  
	// if read->XC is non zero then we cannot extend the alignment into the clipped region (XC is # low quality bases clipped by BWA)

	//check from other direction right to left so we move over cigar string in that order
	for (i=read->cigs-1;i>=0;i--)
	{
	       	op = read->cigarlist[i]&0xf; ol = read->cigarlist[i]>>4; p = ol; 
		// only move into soft clip region if previous cigar is not 'M' extra condition added jan 2 2013 
		if (op == BAM_CSOFT_CLIP && i ==0 && rl->cigs > 0 && (rl->cigarlist[rl->cigs-1]&0xf) ==0 && (rl->cigarlist[rl->cigs-1]>>4) >=10 && indeladded > 0)  
		{
			//fprintf(stdout,"ignoring soft clip region \t");
			rl->cigarlist[rl->cigs] = read->cigarlist[i]; rl->cigs +=1; l1 -= ol; l11 -= ol; continue;
		}
		else if (op == BAM_CMATCH || (op == BAM_CSOFT_CLIP && i ==0) || (op == BAM_CINS && indeladded > 0) || (op == BAM_CINS  && read->position + l2 == bam_indel->position + delta-1 && bam_indel->length > 0 && indeladded ==0) )
		{
			//if (op == BAM_CINS && indeladded > 0) MM[1]++;  
			for (t=0;t<ol;t++) // going from left to right in this case
			{
				//fprintf(stdout,"here t %d %d %d \n",read->position+l2,bam_indel->position,t);
				if (read->position + l2 == bam_indel->position + delta-1 && indeladded ==0)  // indel insertion match in 'M' stretch
				{
					// changed 't < ol-1' -> 't < ol' 06/14/13
					if (bam_indel->length < 0 && bam_indel->simple == '1' && (rl->cigs > 0 || t < ol )) // do not add prior to start of 'M' 
					{
						fprintf(stdout,"cigs %d %d\n",rl->cigs,op);
						p = ol-t; carryover=0;
						if (t > 0 && rl->cigs ==0) 
						{ 
							rl->cigarlist[rl->cigs] = t<<4; rl->cigs +=1; 
						} 
						else if (t > 0)
						{
							b = rl->cigarlist[rl->cigs-1]&0xf;
							if (b ==0) rl->cigarlist[rl->cigs-1] += t<<4;
							else { rl->cigarlist[rl->cigs] = t<<4; rl->cigs +=1; }
						}
						rl->cigarlist[rl->cigs] = ((-1*bam_indel->length)<<4) + 2; rl->cigs +=1; 
						printf("RL_deletion:%d:%d:%d %d:%d ",l1,t,ol,read->position,l2);
						l2 -= -1*bam_indel->length; added = 0; indeladded =1;
					}
					else if (bam_indel->length > 0 && bam_indel->simple == '1' &&  (rl->cigs > 0 || t < ol))
					{
						p = ol-t - bam_indel->length; carryover =0; if (p < 0) carryover = -1*p; 
						// carryover is used for storing # inserted bases over to next cigar
						if (t > 0 && rl->cigs == 0) rl->cigarlist[rl->cigs++] = t<<4; // 'M' to cigar 
						else if (t > 0)
						{
							b = rl->cigarlist[rl->cigs-1]&0xf;
							if (b ==0) rl->cigarlist[rl->cigs-1] += t<<4; else rl->cigarlist[rl->cigs++] = t<<4;
						}
						rl->cigarlist[rl->cigs++] = (bam_indel->length<<4) + 1;
						for (j=bam_indel->length-1;j>=0 && l1-(bam_indel->length-j)+1 >= firstposition;j--) 
						{
							//if (read->sequence[l1-(bam_indel->length-j)+1]  != bam_indel->bases[j+1]) MM[2]++;  
						}
						printf("ins<- %d carry %d ",l1,carryover);
						l1 -= bam_indel->length; t += bam_indel->length;
						added =0; indeladded = 1;
					}
					else if (bam_indel->simple == '0' && (rl->cigs > 0 || t < ol-1)) // complex indel
					{ 
						//dl = strlen(bam_indel[var].RA)-1; il = strlen(bam_indel[var].AA)-1;
						p = ol-t; p -= il; carryover =0;if (p < 0) carryover = -1*p; 

						if (t > 0 && rl->cigs== 0) rl->cigarlist[rl->cigs++] = t<<4;
						else if (t > 0)
						{
							b = rl->cigarlist[rl->cigs-1]&0xf;
							if (b ==0) rl->cigarlist[rl->cigs-1] += t<<4; else rl->cigarlist[rl->cigs++] = t<<4;
						}
						rl->cigarlist[rl->cigs++] = (il<<4) + 1; rl->cigarlist[rl->cigs++] = (dl<<4) + 2;
						for (j=il-1;j>=0 && l1-(il-j)+1 >= firstposition;j--) 
						{
							if (read->sequence[l1-(il-j)+1]  != bam_indel->bases[j+1]) MM[2]++;  
						}
						printf("RLcomplexindel->%d p %d  %c carryover %d %d %d dl:%d il:%d ",l1,p,INT_CIGAROP[op],carryover,MM[2],t,dl,il);
						l1 -= il; t += il; l2 -= dl;
						added = 0; indeladded =1;
					}
				}
				if (t >=ol) break;
				if ( (op == BAM_CMATCH || op == BAM_CSOFT_CLIP || op == BAM_CINS) && carryover > 0) 
				{ 
					carryover--; p--; continue; 
				}
				if (l1 < firstposition) 
				{
					//fprintf(stdout,"| matches %d l1 %d first %d %d | ",matches,l1,firstposition,t);
					break;  // as soon as the read position l1 is smaller than 0/firstnonclipped base, break
				}
				if (added >= 0 && op == BAM_CSOFT_CLIP) added++;  // BUG FIXED 06/14/13 moved this after previous if condition

				// if t >= ol because of insertion -> these won't be counted jan 1 2013 
				//if (read->sequence[l1]  != reflist->sequences[read->tid][read->position+l2] && indeladded ==0) MM[1]++;
				if (read->sequence[l11-t]  != reflist->sequences[read->tid][read->position+l21-t] && indeladded ==0 && op != BAM_CINS)	MM[0]++;

				if (read->sequence[l1]  != reflist->sequences[read->tid][read->position+ l2]) 
				{
					if (op == BAM_CSOFT_CLIP && added >= 5 && carryover <=0 && bam_indel->length < 0) 
					{
						t = ol; added--; // subtract one since we don't count the last base
					}
					else if (indeladded > 0 || op != BAM_CMATCH) MM[2]++; 
				}
				else matches++; 
				//if (added >= 0) fprintf(stdout,"bbb added %d %c:%c",added,read->sequence[l1], reflist->sequences[read->tid][read->position+ l2]);
				if (t < ol) { l1--; l2--; } // BUG fixed, the start position of new cigar is affected ERR024183.19750915 
			}

			if (added < 0) 
			{ 
				if (rl->cigs > 0 && (rl->cigarlist[rl->cigs-1]&0xf) ==0) rl->cigarlist[rl->cigs-1] += ol<<4;
				else { 	rl->cigarlist[rl->cigs] = read->cigarlist[i];  rl->cigs +=1; } 
			}
			else if (p > 0 && (op == BAM_CMATCH || op == BAM_CINS)) 
			{ 
				if (rl->cigs > 0 && (rl->cigarlist[rl->cigs-1]&0xf) ==0) rl->cigarlist[rl->cigs-1] += p<<4;
                                else { rl->cigarlist[rl->cigs] = p<<4; rl->cigs +=1; } 
			}
			else if (carryover ==0 && op == BAM_CSOFT_CLIP) 
			{ 
				printf("no carry S....%d %d %d ",p,carryover,added);
				b = rl->cigarlist[rl->cigs-1]&0xf;
                                if (added > 0 && b != 0)
                                {
                                        rl->cigarlist[rl->cigs] = added<<4; rl->cigs +=1;
                                }
                                else if (added > 0) rl->cigarlist[rl->cigs-1] += added<<4;

                                if (p-added > 0)
                                {
                                        rl->cigarlist[rl->cigs] = ((p-added)<<4) + 4; rl->cigs +=1; // remaining portion of soft clip
                                }
			}
			else if (carryover > 0 && op == BAM_CSOFT_CLIP && indeladded ==0) 
			{ 
				rl->cigarlist[rl->cigs] = (ol-carryover)<<4; rl->cigs +=1; carryover =0;
				// need to take care when carrayover is not from last iteration but current iteration...
			}
			if (op == BAM_CINS) l11 -= ol; else { l11 -= ol; l21 -= ol;} 
		}
		else if (op == BAM_CDEL) 
		{ 
			if (indeladded > 0) MM[1]++; // already added the indel event we are evaluating, so ignore extra deletion event
			else if (bam_indel->position < read->position + l2 - ol ) // added 07/17/13, deletion is ignored if the new variant being evaluated falls within the interval deleted on reference 
			{
				//fprintf(stdout,"del var %d %d \t",bam_indel->position,read->position+l2,ol);
				if (carryover > 0) {  rl->cigarlist[rl->cigs-1] -= carryover<<4; carryover = 0; } // probably never invoked
				rl->cigarlist[rl->cigs] = read->cigarlist[i]; rl->cigs +=1;
				l2 -= ol; 
			}
			l21 -= ol;
		}
		else if (op == BAM_CINS) // only invoked when insertion in old cigar is prior to indeladded event 
		{
			rl->cigarlist[rl->cigs] = read->cigarlist[i]; rl->cigs +=1; l1 -= ol; l11 -= ol;
		}
		else if (op == BAM_CSOFT_CLIP && i > 0)  // soft clip at end of read, copy as it is
		{ 
			rl->cigarlist[rl->cigs] = read->cigarlist[i]; rl->cigs +=1; l1 -= ol; l11 -= ol;
		}

	}
	if (carryover > 0)
	{
	       	op = rl->cigarlist[rl->cigs-1]&0xf; ol = rl->cigarlist[rl->cigs-1]>>4; 
 		if (op != BAM_CDEL && ol > carryover)  rl->cigarlist[rl->cigs-1] -= carryover<<4; 
		else return 0;
		//fprintf(stdout,"carry %d %d \t",carryover,rl->cigarlist[rl->cigs-1]);
	}
	
	// condition for passing read as better aligned as original read 
	int flag = 0;
	if (added ==0 && MM[2] < read->mismatches-MM[0]+MM[1]-1 && MM[2] <2) flag =1;
	if (added == 0 && MM[2] < read->mismatches-MM[0] && MM[2] < 3) flag = 2; 
	if ((added >= 2 && MM[2] < 1) || (added >= 5 && MM[2] <2)) flag = 3; 
	//fprintf(stdout,"flag %d MM %d:%d:%d %d %d\n",flag,MM[0],MM[1],MM[2],read->mismatches,added);
	if (flag > 0)
	{
		rl->newpos = read->position+l2+1; rl->added = added; rl->mismatches = MM[2]; rl->delta = read->mismatches-MM[0]-MM[2]; 
		if (MM[1] > 0) fprintf(stdout,"newflag ");
		printf("newms:%d:%d:%d MM %d:%d added %d newpos %d l2 %d NEWCIGAR_RL: ",MM[2],read->mismatches-MM[0],MM[1],MM[2],MM[0],added,read->position+l2+1,l2);
		for (i=0;i<rl->cigs/2;i++)  // reverse the order of cigar list
		{
			j = rl->cigarlist[i]; rl->cigarlist[i] = rl->cigarlist[rl->cigs-1-i]; rl->cigarlist[rl->cigs-1-i] = j; 
		}
		for (i=0;i<rl->cigs;i++) fprintf(stdout,"%d%c ",rl->cigarlist[i]>>4,INT_CIGAROP[rl->cigarlist[i]&0xf]); 
		for (i=0;i<read->fcigs;i++) fprintf(stdout,"%d%c",read->fcigarlist[i]>>4,INT_CIGAROP[read->fcigarlist[i]&0xf]);
		fprintf(stdout," %d %s\n",read->position,read->sequence);
		//for (i=rl->cigs-1;i>=0;i--) fprintf(stdout,"%d%c ",rl->cigarlist[i]>>4,INT_CIGAROP[rl->cigarlist[i]&0xf]); 
		return 1;
	}
	return 0;
}

//fprintf(stdout,"l1 %d l2 %d %c %c \n",l1,read->position+l2,read->sequence[l1],reflist->sequences[read->tid][read->position+l2]);
//if ((bamfiles_data[j].ilist)->length > 0 && (bamfiles_data[j].ilist)->harray[0]->position < bcall->lastpos && strcmp(bcall->readid,"HISEQ:69:C00NKACXX:5:2305:19838:188000") ==0) 
//fprintf(stderr,"indel list size %d %d %d\n",(bamfiles_data[j].ilist)->length,(bamfiles_data[j].ilist)->harray[0]->position,bcall->lastpos);
int realign_read(struct BAMFILE_data* bamfiles_data,int sample,struct alignedread* bcall,REFLIST* reflist)
{
	struct REALIGNMENT realignments[100];
	int alignments =0; int rlflag =0; int best = -1; int bestscore[2] = {-100,100};  int ss=0;
	int i=0,j=0; 
	
	/*
	fprintf(stdout,"%s %s %s %s:%d ",read->sequence,read->quality,read->readid,read->chrom,read->position);
	for (i=0;i<read->fcigs;++i) fprintf(stdout,"%d%c:",read->fcigarlist[i]>>4,INT_CIGAROP[read->fcigarlist[i]&0xf]);
	fprintf(stdout," XM:%d XC:%d\n",read->mismatches,read->XC);
	*/

	for (i=0;i<(bamfiles_data[sample].ilist)->length;i++) 
	{
		if ((bamfiles_data[sample].ilist)->harray[i]->position > bcall->lastpos) continue;
		if ((bamfiles_data[sample].ilist)->harray[i]->counts[0] + (bamfiles_data[sample].ilist)->harray[i]->counts[1] < 2) continue;

		realignments[0].cigs = 0; rlflag = 0;
		rlflag = realignread_LR(bcall,reflist,(bamfiles_data[sample].ilist)->harray[i],0,&realignments[0]);
		if (rlflag ==1)
		{
			realignments[alignments+1].cigs = realignments[0].cigs; 
			for (j=0;j<realignments[0].cigs;j++) realignments[alignments+1].cigarlist[j] = realignments[0].cigarlist[j]; 
			realignments[alignments+1].added = realignments[0].added;  
			realignments[alignments+1].newpos = realignments[0].newpos;  
			realignments[alignments+1].delta = realignments[0].delta; realignments[alignments+1].mismatches = realignments[0].mismatches; 
			alignments++;
		}

		//fprintf(stdout," === ");
		realignments[0].cigs = 0; rlflag = 0;
		rlflag = realignread_RL(bcall,reflist,(bamfiles_data[sample].ilist)->harray[i],0,&realignments[0]);
		if (rlflag ==1)
		{
			realignments[alignments+1].cigs = realignments[0].cigs; 
			for (j=0;j<realignments[0].cigs;j++) realignments[alignments+1].cigarlist[j] = realignments[0].cigarlist[j]; 
			realignments[alignments+1].added = realignments[0].added; 
			realignments[alignments+1].newpos = realignments[0].newpos;  
			realignments[alignments+1].delta = realignments[0].delta; realignments[alignments+1].mismatches = realignments[0].mismatches; 
			// check if cigars are identical or not...to avoid duplicate LR and RL alignments
			if (realignments[0].cigs != realignments[alignments].cigs || alignments ==0) alignments++; 
			else
			{
				for (j=0;j<realignments[0].cigs;j++) 
				{
					if (realignments[alignments].cigarlist[j] != realignments[0].cigarlist[j]) { alignments++; break; } 
				}
			}
		}
		//fprintf(stdout,"\n");
	}
	// find best realignment if there are multiple ones, if equally good return -1, else return index of best realignment...
	// having two identical alignments (RL and LR) is possible -> treat as one...
	for (i=0;i<alignments;i++) 
	{
		if (realignments[i+1].delta > bestscore[0] && realignments[i+1].mismatches < bestscore[1]) 
		{
			bestscore[0] = realignments[i+1].delta; bestscore[1] = realignments[i+1].mismatches; best = i+1;
		}
		else if (realignments[i+1].delta == bestscore[0]) 
		{ 
			//fprintf(stdout," %d:%d:%d:%d ",realignments[i+1].delta,realignments[i+1].mismatches,bestscore[0],i); 
			best = -1; break; 
		} 
	}
	//if (alignments > 1 && best >= 0 && realignments[best].delta < 2) best = -1; // for multiple realignment case, require delta >=2

	if (best >= 0) 
	{
		fprintf(stdout,"best %d %d %d \t",best,realignments[best].delta,realignments[best].mismatches);
		fprintf(stdout,"realignments %d\n\n",alignments);
	}
	//else fprintf(stdout,"best -1 ");

	bcall->realigned =1;
	return best; 
}
