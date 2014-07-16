
void print_sampleids(struct OPTIONS* options)
{
	// need to extract the bam file path (strip the path and .bam from file name)...
	int i=0,s=0, e=0,j=0;
	for (i=0;i<options->bamfiles;i++) 
	{
		s = strlen(options->bamfilelist[i])-1; j= s+1; 
		while (options->bamfilelist[i][s] != '/' && s >0) s--; s++;  // go backwards to find first '/' character 
		e = s; while (options->bamfilelist[i][e] != '.' && e <j ) e++; // go forward to find first '.' character 
		if (e <= s) e = j;
		if (e > s) 
		{
			for (j=s;j<e;j++) fprintf(options->vfile,"%c",options->bamfilelist[i][j]);  // print string between them
		}
		else fprintf(options->vfile,"%s",options->bamfilelist[i]);
		if (i < options->bamfiles-1) fprintf(options->vfile,"\t"); else fprintf(options->vfile,"\n");
	}
}

void print_crispheader(struct OPTIONS* options)
{
	FILE* vfile = options->vfile;
	// also print program options, poolsize, qvoffset, etc 
	time_t now = time(NULL);
	fprintf(options->vfile,"##fileformat=VCFv4.0\n##fileDate=%s",ctime(&now));
	fprintf(vfile,"##source=CRISP_V0.1\n");
	fprintf(vfile,"##reference=%s\n",options->fastafile);
	fprintf(vfile,"##options: poolsize=%d,bamfiles=%d,qvoffset=%d,bedfile=%s,min_base_quality=%d\n",options->POOLSIZE,options->bamfiles,QVoffset,options->bedfile,MINQ);
	fprintf(vfile,"##INFO=<ID=NP,Number=1,Type=Integer,Description=\"Number of Pools With Data\">\n");
	fprintf(vfile,"##INFO=<ID=DP,Number=2,Type=Integer,Description=\"Total number of reads (+strand,-strand) across all pools (filtered reads only)\">\n");
	fprintf(vfile,"##INFO=<ID=CT,Number=.,Type=Float,Description=\"contingency table p-value for each variant allele in same order as listed in column 5\">\n");
	//fprintf(vfile,"##INFO=<ID=QVpf,Number=.,Type=Float,Description=\"quality values based p-value for each variant allele using forward strand reads\">\n");
	//fprintf(vfile,"##INFO=<ID=QVpr,Number=.,Type=Float,Description=\"quality values based p-value for each variant allele using reverse strand reads\">\n");
	fprintf(vfile,"##INFO=<ID=VP,Number=.,Type=Integer,Description=\"Number of Pools with variant allele(s)\">\n");
	fprintf(vfile,"##INFO=<ID=HP,Number=.,Type=Integer,Description=\"Ambiguity in positioning of indel (homopolymer run length or microsatellite length)\">\n");
	fprintf(vfile,"##INFO=<ID=MQ,Number=.,Type=Integer,Description=\"# of reads with mapping qualities 0-9,10-19,20-39,40-255\">\n");
	fprintf(vfile,"##FILTER=<ID=StrandBias,Description=\"strand bias in distribution of reads between reference and variant alleles on two strands\">\n");
	fprintf(vfile,"##FILTER=<ID=LowDepth,Description=\"low average coverage: less than 1 (filtered) read per haplotype across all samples \">\n");
	fprintf(vfile,"##FILTER=<ID=LowMQ20,Description=\" >20 percent of reads have mapping quality score less than 20\">\n");
	fprintf(vfile,"##FILTER=<ID=LowMQ10,Description=\" >10 percent of reads have mapping quality score less than 20\">\n");
	fprintf(vfile,"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
	//fprintf(vfile,"##FORMAT=<ID=PP,Number=.,Type=Float,Description=\"Probability that the pool contains the variant allele(s), one value for each variant allele listed in column 5\">\n");
	fprintf(vfile,"##FORMAT=<ID=MLAC,Number=.,Type=Integer,Description=\"Maximum likelihood estimate for the allele counts for the ALT allele(s), in the same order as listed in column 5\">\n");
	fprintf(vfile,"##FORMAT=<ID=AF,Number=.,Type=Float,Description=\"variant allele frequency in pool, one value for each variant allele listed in column 5\">\n");
	fprintf(vfile,"##FORMAT=<ID=ADf,Number=.,Type=Integer,Description=\"Number of reads aligned to the forward strand of the genome supporting reference allele and the alternate alleles in the order listed\">\n");
	fprintf(vfile,"##FORMAT=<ID=ADr,Number=.,Type=Integer,Description=\"Number of reads aligned to the reverse strand of the genome supporting reference allele and the alternate alleles in the order listed\">\n");
	fprintf(vfile,"##FORMAT=<ID=ADb,Number=.,Type=Integer,Description=\"Number of overlapping paired-end reads (read from both strands) supporting reference allele and the alternate alleles in the order listed\">\n");
        fprintf(options->vfile,"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t");
	print_sampleids(options);
}


// add option for markduplicates (default is set to 1), --verbosity of output, logfilename, default is stdout 
// change --minc to --minReads 
// add --allele bias for diploid samples... filter reads with small insert size ... (IS < read-length)

void print_crispoptions()
{
	fprintf(stderr,"\nCRISP: statistical method to identify SNVs and indels from pooled DNA sequencing data (requires multiple samples (pools), ideally >=5 samples)\n\n");
	fprintf(stderr,"./CRISP [options] --bams file_bam_paths --ref reference.fasta --VCF variantcalls.VCF -p poolsize > variantcalls.log\n\n");
	fprintf(stderr,"Options:\n");
	fprintf(stderr,"         --bams         	textfile with list of bam file paths (one for each pool)\n");
	fprintf(stderr,"         --bam          	bam file for one pool, specify file for each pool using --bam pool1.bam --bam pool2.bam .... --bam pooln.bam\n");
	fprintf(stderr,"         --ref       		Indexed Reference Sequence file (fasta)\n");
	fprintf(stderr,"         --bed       		bed file for list of regions in which variants should be called (format is chrom start end on each line)\n");
	fprintf(stderr,"         -p/--poolsize <int>    poolsize (number of haploid genomes in each pool), for diploid genomes: 2 x # individuals\n");
	fprintf(stderr,"         --VCF       		VCF file to which the variant calls will be output \n");
	fprintf(stderr,"         --qvoffset <int> 	quality value offset, 33 for Sanger format (default)\n");
	// change this to --phred33 or phred 64 
	fprintf(stderr,"         --mbq     <int>  	minimum base quality to consider a base for variant calling, default 10\n");
	fprintf(stderr,"         --mmq     <int>  	minimum read mapping quality to consider a read for variant calling, default 20\n");
	fprintf(stderr,"         --regions      	region(s) in which variants will be called, e.g chr1:654432-763332. BAM files should be indexed for using this option.\n");
	fprintf(stderr,"         --minc    <int>  	minimum number of reads with alternate allele required for calling a variant, default 4\n");
	fprintf(stderr,"         --ctpval  <float> 	threshold on the contingency table p-value for calling position as variant (specified as log10), default is -3.5\n");
	fprintf(stderr,"         --qvpval  <float> 	threshold on the quality values based p-value for calling position as variant (specified as log10), default is -5\n");
	fprintf(stderr,"         --perms   <int> 	maximum number of permutations for calculating contingency table p-value, default 20000\n");
	fprintf(stderr,"         --filterreads <0/1>	filter reads with excessive number of mismatches (and gaps) compared to the reference sequence, default is 1. Set to 0 to disable filtering\n");
	fprintf(stderr,"         --verbose  <0/1/2> 	amount of information to output to log file, 0: no output, 1: medium (default), 2: detailed\n");
	fprintf(stderr,"         --OPE  <0/1> 		identify overlapping paired-end reads and treat as single read in the overlapping region, default is 1, set to 0 to disable this (can be slow for high-coverage datasets)\n");
	fprintf(stderr,"         --refbias  <float> 	reference allele bias for targeted sequencing data, default is 0.5, use 0.52-0.54 for Agilent SureSelect targeted sequencing experiments\n");
	fprintf(stderr,"	 --EM		<0/1>	0 = old CRISP method with allele frequencies, 1 = EM algorithm will be used for estimating pooled genotypes and calling variants, default is 1\n");
	fprintf(stderr,"	 --flankingbases <int>	call variants in regions that flank target regions in bed file, use 50 or 100 for targeted sequencing, default value is 0 \n");

	//fprintf(stderr,"         --seqerate  	estimate of average sequencing error rate. Default value is -1 which uses individual base qualities to estimate p-values\n");
	//fprintf(stderr,"         --maxm      	maximum number of mismatches allowed in a read to be considered for variant calling (program automatically calculates this based on readlength) set to -1 to disable filtering\n");

	fprintf(stderr,"\nNotes:\n\n");
	fprintf(stderr," 1. CRISP requires poolsize and reference fasta file for making variant calls\n");
	fprintf(stderr," 2. CRISP requires at least two pools to make variant calls, but at least 5 pools are ideal\n");
	fprintf(stderr," 3. The reference sequence file should be indexed using 'samtools faidx' or a similar program and placed in same directory as fasta file with extension .fai\n");
	fprintf(stderr," 4. The ploidy of each pool is assumed to be the same, this will be changed in future releases\n");
	fprintf(stderr," 5. For human re-sequencing, if a bedfile is not specified, the program will evaluate each base for variant calling and the output log file can be huge\n");
	fprintf(stderr," 6. Please set the quality value offset correctly: 33 for quality values encoded in Sanger format and 64 for Illumina 1.3+ format. Setting this value incorrectly can result in significant overcalling/undercalling of variants. Default is 33\n");
	fprintf(stderr," 7. Please make sure that the reference sequence file is the same as the one used to align the reads in the BAM files and that the BAM files are coordinate sorted\n");
	fprintf(stderr," 8. For indel analysis, CRISP assumes that indels are left justified, --leftalign 1 option can be used to left justify gaps in aligned reads\n");
	fprintf(stderr," 9. the program uses the samtools API for reading bam files \n");
	fprintf(stderr," 10. BAM files should be indexed in order to use the --regions option with the indexed bam file as pooln.bam.bai\n\n");
        fprintf(stderr," 11. The aligned reads for each pool should be in a single bam file that is sorted by chromosomal coordinates\n");
        fprintf(stderr," 12. The ctpval and qvpval thresholds are only used with EM = 0 option (older version of CRISP)\n");

	//fprintf(stderr," 11. The bed file can be used to re-call variants at specified positions in the genome\n\n");
	//fprintf(stderr," 7. The value of the parameter --maxm should be choosen based on average length of reads (3 for < 75-bp reads, 4 for < 100 bp reads...)\n\n");

}


void print_crispoptions_additional()
{
	fprintf(stderr,"\n ################### Additional Options for CRISP ##################### \n\n");
	fprintf(stderr,"         --pivotsample  <int>	samples sequenced on two runs, analyze samples 1...PIVOTSAMPLE and PIVOTSAMPLE+1.... separately\n");
	fprintf(stderr,"	 --leftalign	<0/1>	1 = insertions and deletions in reads will be left justified w.r.t. reference sequence (limited evaluation) \n");
	fprintf(stderr,"	 --allowAB 	<0/1>	1 = CRISP will call variants at sites at which the reference base is ambiguous by choosing the most common allele as the reference, this may cause problems for downstream analysis. It is better to use a reference fasta file with no ambiguous bases (all bases should be A,C,T,G or N)\n\n\n");
	fprintf(stderr,"	 --callvariants	<0/1>	0/1 variants called or not\n");
	fprintf(stderr,"	 --mincov	<int>	0/1 minimum average coverage per haploid genome to call variants at a site\n");

}
