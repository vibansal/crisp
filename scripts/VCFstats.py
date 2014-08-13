## analyze VCF file to calculate statistics on Ti/Tv and indel lengths and indels in homopolymers

#! /usr/bin/env python
import sys, os, glob, string

##https://github.com/samtools/htslib/blob/master/vcfcheck.c

###################### TODO #############
"""
1. length distribution of indels 
2. use bed file of exons to see how many indels overlap 
3. analyze genotypes: number of individuals, singleton SNPs etc...
"""


## different variant types 
def variantstats(vfile1):
        print >>sys.stderr, "\ncalculating variant statistics for VCF file",vfile1;
        SNPs=0; indels=0; indels_3n=0; deletions=0; insertions=0; MNPs=0; varsites=0;
        Ts=0; Tv=0; AG =0; CT =0; AC =0; AT = 0; CG =0; GT = 0;
        if vfile1 != '-' and vfile1 != 'sys.stdin': File1 = open(vfile1,'r');
        else: File1 = sys.stdin;


        for line in File1:
                if line[0] == '#': continue;
                varsites +=1;
                if varsites > 0 and varsites%500000==0: print >>sys.stderr, 'processed',varsites,'records';
                var = line.strip().split();  # this should actually be split by '\t' but some CRISP vcf files are space separated
                refallele = var[3]; altalleles =var[4].split(',');
                for i in xrange(len(altalleles)):
                        rl = len(refallele); al = len(altalleles[i]);
                        if rl == al:
                                if len(refallele) ==1: SNPs +=1;
                                else: MNPs +=1;
                        elif rl < al:
                                insertions +=1; indels +=1;
                                if (al-rl)%3 ==0: indels_3n +=1;
                        elif rl > al:
                                deletions +=1; indels +=1;
                                if (rl-al)%3 ==0: indels_3n +=1;

                        if rl != al: continue;
                        if (refallele == 'A' and altalleles[i] == 'G') or (refallele == 'G' and altalleles[i] == 'A'): Ts +=1; AG +=1;
                        elif (refallele == 'C' and altalleles[i] == 'T') or (refallele == 'T' and altalleles[i] == 'C'): Ts +=1; CT +=1;

                        elif (refallele == 'A' and altalleles[i] == 'C') or (refallele == 'C' and altalleles[i] == 'A'): Tv +=1; AC +=1;
                        elif (refallele == 'A' and altalleles[i] == 'T') or (refallele == 'T' and altalleles[i] == 'A'): Tv +=1; AT +=1;
                        elif (refallele == 'C' and altalleles[i] == 'G') or (refallele == 'G' and altalleles[i] == 'C'): Tv +=1; CG +=1;
                        elif (refallele == 'G' and altalleles[i] == 'T') or (refallele == 'T' and altalleles[i] == 'G'): Tv +=1; GT +=1;
                        else: Tv +=1;
        if vfile1 != '-' and vfile1 != 'sys.stdin': File1.close();

        print 'VARSITES:%d SNPs:%d INDELS:%d insertions:%d deletions:%d 3n-indels:%d:%0.3f MNPs:%d' %(varsites,SNPs,indels,insertions,deletions,indels_3n,float(indels_3n)/(indels+0.001),MNPs);
        print 'Ts/Tv: %0.3f %d/%d Transitions: A-G:%d C-T:%d Transversions: A-C:%d A-T:%d C-G:%d G-T:%d\n' %(float(Ts)/Tv,Ts,Tv,AG,CT,AC,AT,CG,GT);


variantstats(sys.argv[1]);
