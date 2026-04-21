
void print_sampleids(struct OPTIONS* options)
{
	int i=0,s=0,e=0,j=0;
	for (i=0;i<options->bamfiles;i++)
	{
		if (options->sampleids != NULL && options->sampleids[i] != NULL)
		{
			fprintf(options->vfile,"%s",options->sampleids[i]);
		}
		else
		{
			// derive from BAM file path: strip directory and extension
			s = strlen(options->bamfilelist[i])-1; j= s+1;
			while (options->bamfilelist[i][s] != '/' && s >0) s--; s++;
			e = s; while (options->bamfilelist[i][e] != '.' && e <j ) e++;
			if (e <= s) e = j;
			if (e > s) for (j=s;j<e;j++) fprintf(options->vfile,"%c",options->bamfilelist[i][j]);
			else fprintf(options->vfile,"%s",options->bamfilelist[i]);
		}
		if (i < options->bamfiles-1) fprintf(options->vfile,"\t"); else fprintf(options->vfile,"\n");
	}
}

// NP=10;DP=6381,4908,1591;VT=SNV;CT=-78.9;VP=1;VF=EMpass;AC=1;AF=0.00258;EMstats=56.52:-31.26;HWEstats=-0.0;MQS=1,0,0,13106;FLANKSEQ=cccatgcagg:C:atctggcagc
void print_crispheader(struct OPTIONS* options)
{
	FILE* vfile = options->vfile;
	// also print program options, poolsize, qvoffset, etc 
	time_t now = time(NULL);
	fprintf(options->vfile,"##fileformat=VCFv4.2\n##fileDate=%s",ctime(&now));
	fprintf(vfile,"##source=CRISP_V0.1\n");
	fprintf(vfile,"##reference=%s\n",options->fastafile);
	fprintf(vfile,"##options: poolsize=%d,bamfiles=%d,qvoffset=%d,bedfile=%s,min_base_quality=%d\n",options->POOLSIZE,options->bamfiles,QVoffset,options->bedfile,MINQ);
	fprintf(vfile,"##INFO=<ID=NP,Number=1,Type=Integer,Description=\"Number of Pools With Data\">\n");
	fprintf(vfile,"##INFO=<ID=DP,Number=.,Type=Integer,Description=\"Total number of reads (+strand,-strand,bidirectional) across all pools (filtered reads only)\">\n");
	fprintf(vfile,"##INFO=<ID=CT,Number=.,Type=Float,Description=\"contingency table p-value for each variant allele in same order as listed in column 5\">\n");
	fprintf(vfile,"##INFO=<ID=AC,Number=.,Type=Integer,Description=\"estimated allele count (across all pools) for each variant allele in same order as listed in column 5\">\n");
	fprintf(vfile,"##INFO=<ID=AF,Number=.,Type=Float,Description=\"estimated allele frequency for (estimated by EM algorithm jointly across all pools) each variant allele in same order as listed in column 5\">\n");
	fprintf(vfile,"##INFO=<ID=FLANKSEQ,Number=1,Type=String,Description=\"sequence flanking the variant\">\n");
	fprintf(vfile,"##INFO=<ID=VT,Number=1,Type=String,Description=\"variant type, SNV, INSERTION or DELETION\">\n");
	fprintf(vfile,"##INFO=<ID=EMstats,Number=.,Type=String,Description=\"statistics from EM algorithm\">\n");
	fprintf(vfile,"##INFO=<ID=HWEstats,Number=.,Type=String,Description=\"p-value for HWE test across all pools\">\n");
	fprintf(vfile,"##INFO=<ID=VF,Number=.,Type=String,Description=\"EMpass or EMfail, EMfail variants should be discarded\">\n");
	//fprintf(vfile,"##INFO=<ID=QVpf,Number=.,Type=Float,Description=\"quality values based p-value for each variant allele using forward strand reads\">\n");
	//fprintf(vfile,"##INFO=<ID=QVpr,Number=.,Type=Float,Description=\"quality values based p-value for each variant allele using reverse strand reads\">\n");
	fprintf(vfile,"##INFO=<ID=VP,Number=.,Type=Integer,Description=\"Number of Pools with variant allele(s)\">\n");
	fprintf(vfile,"##INFO=<ID=HP,Number=.,Type=Integer,Description=\"Ambiguity in positioning of indel (homopolymer run length or microsatellite length)\">\n");
	fprintf(vfile,"##INFO=<ID=MQS,Number=.,Type=Integer,Description=\"# of reads with mapping qualities 0-9,10-19,20-39,40-255\">\n");
	fprintf(vfile,"##FILTER=<ID=StrandBias,Description=\"strand bias in distribution of reads between reference and variant alleles on two strands\">\n");
	fprintf(vfile,"##FILTER=<ID=LowDepth,Description=\"low average coverage: less than 1 (filtered) read per haplotype across all samples \">\n");
	fprintf(vfile,"##FILTER=<ID=LowMQ20,Description=\" >20 percent of reads have mapping quality score less than 20\">\n");
	fprintf(vfile,"##FILTER=<ID=LowMQ10,Description=\" >10 percent of reads have mapping quality score less than 20\">\n");
	fprintf(vfile,"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
	fprintf(vfile,"##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth for each pool/sample (includes all reads)\">\n");
	fprintf(vfile,"##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n");
	//fprintf(vfile,"##FORMAT=<ID=PP,Number=.,Type=Float,Description=\"Probability that the pool contains the variant allele(s), one value for each variant allele listed in column 5\">\n");
	fprintf(vfile,"##FORMAT=<ID=AC,Number=.,Type=Integer,Description=\"Maximum likelihood estimate for the allele counts for the ALT allele(s), in the same order as listed in column 5\">\n");
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
	fprintf(stderr,"\nCRISP: statistical method to identify SNVs and indels from pooled DNA sequencing data (requires multiple pools/samples, ideally >=5)\n\n Please provide a bed-file for calling variants on human samples to avoid large output files\n\n");
	fprintf(stderr,"./CRISP [options] --bams file_bam_paths --ref reference.fasta --VCF variantcalls.VCF -p poolsize --bed targets.bed > variantcalls.log\n\n");
	fprintf(stderr,"Options:\n");
	fprintf(stderr,"         --bams         \ttextfile with list of BAM/CRAM file paths (one per pool)\n");
	fprintf(stderr,"         --bam          \tBAM/CRAM file for one pool; repeat for each pool: --bam pool1.bam --bam pool2.bam ...\n");
	fprintf(stderr,"         --ref          \tindexed reference FASTA file (required; mandatory for CRAM input)\n");
	fprintf(stderr,"         --bed          \tBED file with intervals for variant calling (chrom start end per line)\n");
	fprintf(stderr,"         -p/--poolsize <int>\tpoolsize (number of haploid genomes per pool); for diploid: 2 x number of individuals\n");
	fprintf(stderr,"         --VCF          \toutput VCF file for variant calls\n");
	fprintf(stderr,"         --sm <0/1>     \tuse @RG SM tag from BAM/CRAM header as VCF sample name (default 0: derive from filename)\n");
	fprintf(stderr,"         --qvoffset <int>\tquality value offset: 33 (Sanger, default) or 64 (Illumina 1.3+)\n");
	fprintf(stderr,"         --mbq <int>    \tminimum base quality for variant calling, default 10\n");
	fprintf(stderr,"         --mmq <int>    \tminimum read mapping quality for variant calling, default 20\n");
	fprintf(stderr,"         --regions      \tregion for variant calling, e.g. chr1:654432-763332 (requires indexed BAM/CRAM)\n");
	fprintf(stderr,"         --minc <int>   \tminimum reads with alternate allele to call a variant, default 4\n");
	fprintf(stderr,"         --ctpval <float>\tcontingency table p-value threshold (log10), default -3.5; increase as pool count grows\n");
	fprintf(stderr,"         --perms <int>  \tmaximum permutations for contingency table p-value, default 20000\n");
	fprintf(stderr,"         --filterreads <0/1>\tfilter reads with excessive mismatches vs reference, default 1 (0 to disable)\n");
	fprintf(stderr,"         --verbose <0/1/2>\tlog verbosity: 0 = none, 1 = medium (default), 2 = detailed\n");
	fprintf(stderr,"         --OPE <0/1>    \tidentify overlapping paired-end reads, treat as single read in overlap, default 1\n");
	fprintf(stderr,"         --refbias <float>\treference allele bias for targeted sequencing, default 0.5 (use 0.52-0.54 for Agilent SureSelect)\n");
	fprintf(stderr,"         --EM <0/1>     \t0 = allele frequency method, 1 = EM algorithm for pooled genotype calling (default)\n");
	fprintf(stderr,"         --mincov <int> \tminimum average coverage per haploid genome to call variants at a site\n");
	fprintf(stderr,"         --flankingbases <int>\tcall variants in flanking regions around bed targets (use 50-100 for targeted sequencing), default 0\n");

	fprintf(stderr,"\nNotes:\n\n");
	fprintf(stderr," 1. CRISP requires poolsize (-p), reference FASTA (--ref), BAM/CRAM file list (--bams), and output VCF (--VCF)\n");
	fprintf(stderr," 2. Both BAM and CRAM input files are supported; CRAM files require --ref\n");
	fprintf(stderr," 3. CRISP uses htslib for reading BAM/CRAM files (install: apt-get install libhts-dev or brew install htslib)\n");
	fprintf(stderr," 4. The reference FASTA must be indexed (samtools faidx reference.fasta); the .fai index must be in the same directory\n");
	fprintf(stderr," 5. CRISP validates that contig names and lengths in the BAM/CRAM header match the reference FASTA\n");
	fprintf(stderr," 6. CRISP requires at least two pools; at least 5 pools are recommended\n");
	fprintf(stderr," 7. Pool sizes are assumed equal across pools; variable pool sizes can be specified in the --bams input file\n");
	fprintf(stderr," 8. For human re-sequencing, use --bed to restrict variant calling; without it, every position is evaluated\n");
	fprintf(stderr," 9. BAM/CRAM files must be indexed for the --regions option (pooln.bam.bai or pooln.cram.crai)\n");
	fprintf(stderr," 10. Each pool should be in a single coordinate-sorted BAM/CRAM file\n");
	fprintf(stderr," 11. For indel calling, --leftalign 1 left-justifies gaps in aligned reads\n\n");

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
