
CC=gcc -Wall -D_GNU_SOURCE
CFLAGS= -D_GNU_SOURCE -O2 -Wall -std=gnu99 -std=c99 -lm

## modify this variable to point to the location of the samtools directory in order to compile CRISP 
SAMTOOLS=/home/vikas/Public/WorkingProjects/samtools-0.1.18

all:
	crispdip

## this is the older version of CRISP that does not call genotypes 
crisp: include/readfasta.o include/bamsreader.o include/variant.o include/allelecounts.o include/chisquare.o include/contables.o include/crispcaller.o include/bamread.o optionparser.c
	$(CC) -DPICALL=0 -I$(SAMTOOLS) -g -O2 include/readfasta.o include/variant.o include/allelecounts.o include/bamsreader.o include/crispcaller.o include/chisquare.o include/contables.o include/bamread.o -lm readmultiplebams.c -o bin/CRISP -L$(SAMTOOLS) -lbam -lz

include/crispcaller.o:	crisp/crispcaller.c crisp/crispcaller.h FET/contables.h variant.h bamsreader.h FET/pooledFET.c crisp/crispprint.c
	$(CC) -I$(SAMTOOLS) -c crisp/crispcaller.c -o include/crispcaller.o


## picall replacement to call indels only 
crispindel: include/readfasta.o include/bamsreader.o include/variant.o include/allelecounts.o include/chisquare.o include/contables.o include/newcrispcaller.o optionparser.c FET/lowcovFET.c include/bamread.o
	$(CC) -DPICALL=3 -I$(SAMTOOLS) -g -O2 include/readfasta.o include/variant.o include/allelecounts.o include/bamsreader.o include/newcrispcaller.o include/chisquare.o include/contables.o include/bamread.o -lm readmultiplebams.c -o bin/CRISPindel -L$(SAMTOOLS) -lbam -lz



## new version of CRISP that calls genotypes using likelihood and MLE method 
crispdip: include/readfasta.o include/bamsreader.o include/variant.o include/allelecounts.o include/chisquare.o include/contables.o include/newcrispcaller.o optionparser.c FET/lowcovFET.c include/bamread.o
	$(CC) -DPICALL=3 -I$(SAMTOOLS) -g -O2 include/readfasta.o include/variant.o include/allelecounts.o include/bamsreader.o include/newcrispcaller.o include/chisquare.o include/contables.o include/bamread.o -lm readmultiplebams.c -o bin/CRISP-diploid -L$(SAMTOOLS) -lbam -lz

include/newcrispcaller.o:	crisp/newcrispcaller.c crisp/crispcaller.h FET/contables.h variant.h bamsreader.h bamread.h crisp/crispEM.c FET/pooledFET.c crisp/newcrispprint.c
	$(CC) -I$(SAMTOOLS) -c crisp/newcrispcaller.c -o include/newcrispcaller.o



### shared objects ### 

include/variant.o: variant.c variant.h bamsreader.h bamsreader.c ../readfasta.h ../readfasta.c process_indel_variant.c
	$(CC) -I$(SAMTOOLS) -c variant.c -o include/variant.o

include/bamsreader.o: ../readfasta.h ../readfasta.c bamsreader.h bamsreader.c bamread.h bamread.c
	$(CC) -I$(SAMTOOLS) -c bamsreader.c -o include/bamsreader.o

include/bamread.o: ../readfasta.h ../readfasta.c bamread.h bamread.c 
	$(CC) -I$(SAMTOOLS) -c bamread.c -o include/bamread.o

include/allelecounts.o: ../readfasta.h allelecounts.h allelecounts.c variant.h
	$(CC) -I$(SAMTOOLS) -c allelecounts.c -o include/allelecounts.o

include/readfasta.o: ../readfasta.c ../readfasta.h 
	$(CC) -c ../readfasta.c -o include/readfasta.o

include/chisquare.o: FET/chisquare.h FET/chisquare.c
	$(CC) -I$(SAMTOOLS) -c FET/chisquare.c -o include/chisquare.o  

include/contables.o: FET/contables.c FET/contables.h FET/chisquare.h FET/chisquare.c
	$(CC) -c FET/contables.c -o include/contables.o

clean:
	rm -f include/*.o 
	#readfasta.o readbams bamsreader.o allelecounts.o variant.o chisquare.o contables.o crispcaller.o newcrispcaller.o




##########################################################################################################

bfgs.o:	BFGSmethod/bfgs.c BFGSmethod/bfgs.h
	$(CC) -c BFGSmethod/bfgs.c 

indels: bfgs.o readfasta.o bamsreader.o variant.o allelecounts.o chisquare.o contables.o crispcaller.o optionparser.c crisp/crispEM.c
	$(CC) -DPICALL=0 -I$(SAMTOOLS) -g -O2 bfgs.o readfasta.o variant.o allelecounts.o bamsreader.o crispcaller.o chisquare.o contables.o -lm extractindelreads.c -o INDELS -lz -L$(SAMTOOLS) -lbam


picall: readfasta.o bamsreader.o variant.o allelecounts.o chisquare.o contables.o optionparser.c picall.o 
	$(CC) -DPICALL=1 -I$(SAMTOOLS) -g -O2 readfasta.o variant.o allelecounts.o bamsreader.o chisquare.o contables.o picall.o -lm readmultiplebams.c -o PICALL -lz -L$(SAMTOOLS) -lbam

## new version of CRISP as well but this one uses the BFGScode adapted (open source) rather than the BMC code 
crispbfgs: readfasta.o bamsreader.o variant.o allelecounts.o chisquare.o contables.o newcrispcaller.o optionparser.c ../Lbfgsb.3.0/lbfgsb.o
	$(CC) -DPICALL=3 -I$(SAMTOOLS) -g -O2 ../Lbfgsb.3.0/lbfgsb.o readfasta.o variant.o allelecounts.o bamsreader.o newcrispcaller.o chisquare.o contables.o -lm readmultiplebams.c -o bin/CRISP-diploid -lz -L$(SAMTOOLS) -lbam
lowcov: readfasta.o bamsreader.o variant.o allelecounts.o chisquare.o contables.o optionparser.c picall.o lowcoveragecaller.o
	$(CC) -DPICALL=2 -I$(SAMTOOLS) -g -O2 readfasta.o variant.o allelecounts.o bamsreader.o chisquare.o contables.o picall.o lowcoveragecaller.o -lm readmultiplebams.c -o bin/LCvarcaller -lz -L$(SAMTOOLS) -lbam
lowcoveragecaller.o:	lowcov/lowcoveragecaller.c lowcov/lowcoveragecaller.h contables.h variant.h bamsreader.h allelecounts.h
	$(CC) -I$(SAMTOOLS) -c lowcov/lowcoveragecaller.c

picall.o:	picall/picall.c picall/picall.h variant.h crisp/crispcaller.h allelecounts.h contables.h bamsreader.h picall/integrals.c picall/popgenotypes.c
	$(CC) -I$(SAMTOOLS) -c picall/picall.c
