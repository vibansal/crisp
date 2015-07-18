
CC=gcc -Wall -D_GNU_SOURCE
CFLAGS= -D_GNU_SOURCE -O2 -Wall -std=gnu99 -std=c99 -lm

## modify this variable to point to the location of the samtools directory in order to compile CRISP 
SAMTOOLS=samtools

## new version of CRISP that calls genotypes using likelihood and MLE method 
all: 
	$(MAKE) -C samtools all 
	$(MAKE) -C . new 
	#$(MAKE) -C . crisp
	$(MAKE) -C . indel
	cp bin/CRISP-genotypes CRISP

new: include/readfasta.o include/bamsreader.o include/variant.o include/allelecounts.o include/chisquare.o include/contables.o include/newcrispcaller.o optionparser.c FET/lowcovFET.c include/bamread.o
	$(CC) -DPICALL=3 -I$(SAMTOOLS) -g -O2 include/readfasta.o include/variant.o include/allelecounts.o include/bamsreader.o include/newcrispcaller.o include/chisquare.o include/contables.o include/bamread.o -lm readmultiplebams.c -o bin/CRISP-genotypes -L$(SAMTOOLS) -lbam -lz
	

## this is the older version of CRISP that does not call genotypes 
crisp: include/readfasta.o include/bamsreader.o include/variant.o include/allelecounts.o include/chisquare.o include/contables.o include/crispcaller.o include/bamread.o optionparser.c
	$(CC) -DPICALL=0 -I$(SAMTOOLS) -g -O2 include/readfasta.o include/variant.o include/allelecounts.o include/bamsreader.o include/crispcaller.o include/chisquare.o include/contables.o include/bamread.o -lm readmultiplebams.c -o bin/CRISP-AF -L$(SAMTOOLS) -lbam -lz

include/crispcaller.o:	crisp/crispcaller.c crisp/crispcaller.h FET/contables.h variant.h bamsreader.h FET/pooledFET.c crisp/crispprint.c
	$(CC) -I$(SAMTOOLS) -c crisp/crispcaller.c -o include/crispcaller.o

## picall replacement to call indels only 
indel: include/readfasta.o include/bamsreader.o include/variant.o include/allelecounts.o include/chisquare.o include/contables.o include/newcrispcaller.o optionparser.c FET/lowcovFET.c include/bamread.o
	$(CC) -DPICALL=3 -I$(SAMTOOLS) -g -O2 include/readfasta.o include/variant.o include/allelecounts.o include/bamsreader.o include/newcrispcaller.o include/chisquare.o include/contables.o include/bamread.o -lm readmultiplebams.c -o bin/CRISPindel -L$(SAMTOOLS) -lbam -lz

### shared objects ### 

include/newcrispcaller.o:	crisp/newcrispcaller.c crisp/crispcaller.h FET/contables.h variant.h bamsreader.h bamread.h crisp/crispEM.c FET/pooledFET.c crisp/newcrispprint.c
	$(CC) -I$(SAMTOOLS) -c crisp/newcrispcaller.c -o include/newcrispcaller.o

include/variant.o: variant.c variant.h bamsreader.h bamsreader.c readfasta.h readfasta.c process_indel_variant.c
	$(CC) -I$(SAMTOOLS) -c variant.c -o include/variant.o

include/bamsreader.o: readfasta.h readfasta.c bamsreader.h bamsreader.c bamread.h bamread.c
	$(CC) -I$(SAMTOOLS) -c bamsreader.c -o include/bamsreader.o

include/bamread.o: readfasta.h readfasta.c bamread.h bamread.c 
	$(CC) -I$(SAMTOOLS) -c bamread.c -o include/bamread.o

include/allelecounts.o: readfasta.h allelecounts.h allelecounts.c variant.h
	$(CC) -I$(SAMTOOLS) -c allelecounts.c -o include/allelecounts.o

include/readfasta.o: readfasta.c readfasta.h 
	$(CC) -c readfasta.c -o include/readfasta.o

include/chisquare.o: FET/chisquare.h FET/chisquare.c
	$(CC) -I$(SAMTOOLS) -c FET/chisquare.c -o include/chisquare.o  

include/contables.o: FET/contables.c FET/contables.h FET/chisquare.h FET/chisquare.c
	$(CC) -c FET/contables.c -o include/contables.o

clean:
	rm -f include/*.o CRISP bin/CRISP 
