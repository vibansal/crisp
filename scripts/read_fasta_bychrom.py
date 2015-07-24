import sys,os

#http://stackoverflow.com/questions/620367/python-how-to-jump-to-a-particular-line-in-a-huge-text-file
## read fasta file, build simple index, next time reuse it by writing to file, allows reading of one chromosome at a time... good enough for processing data efficiently... 

## CODED dec 19 2014, simple code to read fasta file, create index allowing to read one chromosome at a time 

## make index of start positions for each chromosome...
def make_fasta_index(fastafile):
	if os.path.isfile(fastafile + '.python.fai'): 
		print >>sys.stderr, 'fasta index exists, reading from file';
		offset_index = {};
		File = open(fastafile + '.python.fai','r');
		for line in File: chrom = line.strip().split(); offset_index[chrom[0]] = int(chrom[1]); 
		File.close(); 
		return offset_index; 

	print >>sys.stderr, 'generating fasta index';
	File= open(fastafile,'r');
	line_offset = []
	offset = 0;
	for line in File:
		if line[0] == '>': line_offset.append([offset,line.strip().split()[0]]); print >>sys.stderr,line, 
		offset += len(line)
	print >>sys.stderr, 'finished reading fasta file'
	File.close();
	#File.seek(16909); print 'line',File.readline(); 
	offset_index = {}; 
	for a in line_offset: offset_index[a[1][1:]] = a[0]; 
	## write index to file 'fastafile' + '.pythonfai'
	File = open(fastafile + '.python.fai','w');
	for a in line_offset: File.write(a[1][1:] + '\t'+ `a[0]` + '\n'); 
	File.close(); 

	return offset_index

## read a single chromosome from fasta file, all lines from a particular seek
def read_chromosome(fastafile,offset_index,chrom_name):
	File= open(fastafile,'r');
	File.seek(offset_index[chrom_name]); 
	sequences = []; length = 0;
	line = File.readline(); print >>sys.stderr, chrom_name,offset_index[chrom_name],line[1:],
	while 1: 
		line = File.readline()
		if line[0] == '>': break; 
		else: sequences.append(line.strip('\n')); length += len(sequences[-1]);
	File.close(); 
	print >>sys.stderr, 'chromosome stats for',chrom_name,'length:',len(sequences),length;
	sequence = ''.join(sequences); 
	return [sequence,length]

#offset_index = make_fasta_index(sys.argv[1]); current_chrom = '-';

