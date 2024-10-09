/* FUNCTIONS to extract indel alleles from reads, calculate counts and reduce counts for ambigous reads that do not span the entire indel homopolymer or repeat tract , created 07/04/11 */

#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>

int indel_cmp(const void *a,const void *b)
{
	const struct INDEL *ia = (const struct INDEL*)a; const struct INDEL *ib = (const struct INDEL*)b;
	return ib->reads-ia->reads;// sort in descending order
	//return strcmp(ia->bases,ib->bases);
}

// function to reduce count in favor of reference allele due to reads that do not span the indel homopolymer stretch at the end of a read | reads with the homopolymer stretch at beginning of read will automatically not be counted towards reference allele count 
//HPlength is ambiguity in positioning of indel 
// this function can reduce the allele counts at a given site substantially .... 
// if delta = +1 counts are reduced | if delta = -1 counts are set back to original values
int remove_ambiguous_bases(READQUEUE* bq,struct VARIANT* variant,int HPlength,int delta) 
{
	int base; int offset=0; double ep; int op = 0; int ol=0; int bin=0; int bases_removed0 =0,bases_removed1=0;
	struct alignedread* bcall = bq->first; 
	//if (delta < 0) fprintf(stdout,"resetting counts\n");  
	for (bcall = bq->first; bcall != NULL && bcall->position <= variant->position; bcall = bcall->next)
	{
		if (bcall->lastpos <= variant->position) continue; 
		//if (bcall->position > variant->position || bcall->lastpos < variant->position) break; 
		if (bcall->strand == 1) offset = maxalleles; else offset = 0;

		if (bcall->filter =='0' && bcall->allele-offset == variant->refbase)
		{
			op = bcall->fcigarlist[bcall->cigoffset]&0xf; ol = bcall->fcigarlist[bcall->cigoffset]>>4;
			//if (delta < 0) fprintf(stdout,"bcall %s position %d %d op %d ol %d %d %d\n",bcall->readid,bcall->position,variant->position,op,ol,bcall->delta,HPlength);
			if (op != BAM_CMATCH || (ol-bcall->delta <= HPlength)) // condition for removing
			{
				base = variant->refbase; // 0/1/2/3
				//fprintf(stdout,"pos %d base %d sample %d %d %d counts %d %d\n",variant->position,base,bcall->sampleid,HPlength,bcall->dfe,variant->counts[base+offset],bcall->pos); 
				if (bcall->bidir == 1) 
				{
					// this is unlikely to happen, since indel will be start of read from one end and not for other read
					variant->indcounts[bcall->sampleid][base+2*maxalleles] -= delta; variant->counts[base+2*maxalleles] -=delta; 
				}
				else 
				{
					variant->indcounts[bcall->sampleid][base+offset] -=delta; variant->counts[base+offset] -=delta; 
					if (offset ==0) bases_removed0++; else bases_removed1++;
				}
				// what happens to bidirectional reads in code below | are they binned june 13 2013

				bin = (int)(bcall->quality[bcall->l1+bcall->delta]-QVoffset)/10-1;
				if (bin >=0) variant->indcounts_binned[bcall->sampleid][base+offset][bin > 2 ? 2: bin] -=delta;
				if ( (bcall->flag & BAM_PAIRED_READ1) == BAM_PAIRED_READ1) variant->indcounts_binned[bcall->sampleid][base+offset][3] -=delta;
				else if ((bcall->flag & BAM_PAIRED_READ2) ==BAM_PAIRED_READ2 ) variant->indcounts_binned[bcall->sampleid][base+offset][4] -=delta;
				else variant->indcounts_binned[bcall->sampleid][base+offset][5] -=delta;

				ep = pow(0.1,((double)bcall->quality[bcall->l1+bcall->delta]-QVoffset)/10); 
				variant->stats[bcall->sampleid][base+offset] -= delta*ep; variant->tstats[base+offset] -= (double)delta*ep; 
				if (USE_QV ==1) variant->Qhighcounts[bcall->sampleid][base+offset] -= (double)delta*(1.0-ep); 
				else variant->Qhighcounts[bcall->sampleid][base+offset] -= delta;
			}
		}
	}
        if (PFLAG >=1) fprintf(stdout,"%d:%d_removed HP=%d,%d\t",bases_removed0,bases_removed1,HPlength,delta);
	return bases_removed0+bases_removed1; 
}

// calculate the right ambiguity in positioning of indel, we will assume that indels are left justified in the alignments.
int calculate_indel_ambiguity(struct VARIANT* variant,REFLIST* reflist,int current,int allele)
{
	int j=0,k=0,rightamb=0;
	int il = strlen(variant->itb[allele])-1;
	if (il ==0) 
	{
		fprintf(stderr,"indel of length 0 position %d allele %d %s \n",variant->position,allele,variant->itb[allele]);
		variant->HPlength[allele] = 0; return 1;
	}
	if (variant->itb[allele][0]=='-')  // deletion allele
	{
		j = variant->position; k = variant->position+il;  // k starts after the deleted bases assuming left justified
		while (j < reflist->lengths[current] && toupper(reflist->sequences[current][j]) == toupper(reflist->sequences[current][k])) 
		{ 
			j++; k++; rightamb++; 
		}  
	}
	else 
	{
		j = variant->position;
		while (j < reflist->lengths[current] && toupper(reflist->sequences[current][j]) == variant->itb[allele][rightamb%il+1]) { j++; rightamb++;}
	}
	variant->HPlength[allele] = rightamb;
}

// we need to go through the baselist and identify the three most common indel alleles and instantiate them in variant->itb
int identify_indelalleles(REFLIST* reflist,int current,READQUEUE* bq, struct VARIANT* variant)
{
	struct INDEL* indelalleles = (struct INDEL*)malloc(sizeof(struct INDEL)*10); 
	int alleles =0; // number of indel alleles
	int offset=0; int i=0,j=0,k=0;
	int flag =0; 
	int varalleles = 0; // final number of indel alleles identified
	int discard=0; int match=0;
	struct alignedread* bcall = bq->first; 
	for (bcall = bq->first; bcall != NULL && bcall->position <= variant->position; bcall = bcall->next)
	{
		if (bcall->lastpos <= variant->position) continue; 
		if (bcall->filter == '0' && (bcall->type != 0)) 
		{
			flag = 0;
			for (i=0;i<alleles && flag ==0;i++)
			{
				if (indelalleles[i].length == bcall->type && bcall->type < 0) {indelalleles[i].reads++; flag = 1;}
				else if (indelalleles[i].length == bcall->type)
				{
					match =1;
					for (j=0;j<bcall->type;j++)
					{
						if (indelalleles[i].bases[j+1] != bcall->sequence[bcall->l1+j]) match = 0; 
					}
					if (match ==1) {indelalleles[i].reads++; flag = 1;}
				}
				//if (indelalleles[i].length == bcall->type) {indelalleles[i].reads++; flag = 1; break;}
			}
			if (flag ==0 && alleles < 10)
			{
				indelalleles[alleles].length = bcall->type; 
				if (bcall->type < 0) 
				{
					indelalleles[alleles].bases = (char*)malloc(-1*bcall->type+2); 
					indelalleles[alleles].bases[0] = '-';
					for (i=0;i<-1*bcall->type;i++) indelalleles[alleles].bases[i+1] = toupper(reflist->sequences[current][bcall->l2+i]);
					indelalleles[alleles].bases[-1*bcall->type+1] = '\0';
					indelalleles[alleles].reads = 1; 
					alleles++; 

				}
				// discard insertion alleles with 'N' in them
				else if (bcall->type >0)
				{
					discard =0;
					for (i=0;i<bcall->type;i++) { if (bcall->sequence[bcall->l1+i] == 'N') discard=1; } 
					if (discard ==0) 
					{
						indelalleles[alleles].bases = (char*)malloc(bcall->type+2); 
						indelalleles[alleles].bases[0] = '+'; 
						for (i=0;i<bcall->type;i++) indelalleles[alleles].bases[i+1] = bcall->sequence[bcall->l1+i];
						indelalleles[alleles].bases[bcall->type+1] = '\0';
						indelalleles[alleles].reads = 1; 
						alleles++; 
					}
				}
			}
		}
	}
	if (alleles ==0)  return 0; // no indel alleles detected, nothing to deallocate, fix 10-09-2024
	if (alleles > 1) qsort(indelalleles,alleles,sizeof(struct INDEL),indel_cmp); 
	if (indelalleles[0].reads < 3) // first indel allele has less than three reads, ignore
	{
		for (i=0;i<alleles;i++) free(indelalleles[i].bases); free(indelalleles);
		return 0;
	}

	varalleles =1; 
	strcpy(variant->itb[4],indelalleles[0].bases); 
	if (alleles >=2  && indelalleles[1].reads >=3) 
	{ 
		strcpy(variant->itb[5],indelalleles[1].bases);	
		varalleles++; 
	}
	if (alleles >=3  && indelalleles[2].reads >=3) 
	{ 
		strcpy(variant->itb[6],indelalleles[2].bases); 
		varalleles++; 
	}

	// add counts for indel alleles to variant->counts and variant->indcounts
	for (bcall = bq->first; bcall != NULL && bcall->position <= variant->position; bcall = bcall->next)
	{
		if (bcall->lastpos < variant->position || bcall->filter =='1' || bcall->type ==0) continue; 
		if (bcall->strand == 1) offset = maxalleles; else offset =0;
		for (k=0;k<varalleles;k++) 
		{
			if (indelalleles[k].length == bcall->type) 
			{
				bcall->allele = 4+offset+k;
				if (bcall->bidir ==1)
				{
					variant->counts[4+2*maxalleles+k]++; variant->indcounts[bcall->sampleid][4+2*maxalleles+k]++; 
				}
				else
				{
					variant->counts[4+offset+k]++; variant->indcounts[bcall->sampleid][4+offset+k]++; 
				}
				if ((bcall->flag & BAM_PAIRED_READ1) == BAM_PAIRED_READ1) variant->indcounts_binned[bcall->sampleid][4+offset+k][3]++;
				else if ((bcall->flag & BAM_PAIRED_READ2)== BAM_PAIRED_READ2 ) variant->indcounts_binned[bcall->sampleid][4+offset+k][4]++;
				else variant->indcounts_binned[bcall->sampleid][4+offset+k][5]++;
				variant->Qhighcounts[bcall->sampleid][4+offset+k]++;
				k = varalleles; 
			}
		}
	}
	for (i=0;i<alleles;i++) free(indelalleles[i].bases); free(indelalleles);
	return varalleles;
	// +AA, +A, -AC determine counts for each allele 

}

// this can also work for false SNPs close to indels...
// print indel haplotypes in order to identify neighboring variants on same haplotypes | goal is to merge complex indels or SNPs in phase -> single variant
// cigar of haplotype 0- +20 around variant and find frequent cigars.... we are looking for case where single cigar is consensus | if read is shorter than 20 bp after variant -> partial cigar.. that is not considered
// assumption is that bcall->type and bcall->cigoffset are correct 
// current code is wrong, we need to first find the position in read that matches the variant position and then check it ??
int print_indel_haplotypes(READQUEUE* bq,struct VARIANT* variant,int allele2) 
{
	int indel_length = 0; 
	if (variant->itb[allele2][0] == '-') indel_length = 1-strlen(variant->itb[allele2]); 
	else if (variant->itb[allele2][0] == '+') indel_length = strlen(variant->itb[allele2])-1; 

	int base; int offset=0; double ep; int op = 0; int ol=0; int i=0;

	struct alignedread* bcall = bq->first; 
	for (bcall = bq->first; bcall != NULL && bcall->position <= variant->position; bcall = bcall->next)
	{
		if (bcall->type ==0 && bcall->sequence[bcall->l1+bcall->delta] != variant->itb[allele2][0]) continue; 
		if (bcall->filter =='0' && bcall->type == indel_length && bcall->lastpos > variant->position )
		{
			//for (i=bcall->cigoffset;i<bcall->fcigs;i++) fprintf(stdout,"%d%c:",bcall->fcigarlist[i]>>4,INT_CIGAROP[bcall->fcigarlist[i]&0xf]);
			for (i=0;i<bcall->fcigs;i++) fprintf(stdout,"%d%c:",bcall->fcigarlist[i]>>4,INT_CIGAROP[bcall->fcigarlist[i]&0xf]);
			fprintf(stdout,"%s %d\n",bcall->sequence,bcall->position);
			fprintf(stdout," %d %d %d:%d\n",bcall->position,bcall->cigoffset,bcall->type,bcall->l1);
		}
	}
	//int frequent_cigs[10][16]; int freq[10]; int cig[16]; 
}
