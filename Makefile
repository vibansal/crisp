
CC=gcc -Wno-all -D_GNU_SOURCE
CFLAGS= -D_GNU_SOURCE -O2 -Wno-all -std=gnu99 -std=c99 -lm

## Use system-installed htslib (provides BAM/CRAM support)
## Install htslib: apt-get install libhts-dev (Debian/Ubuntu) or brew install htslib (macOS)
HTSLIB_CFLAGS=$(shell pkg-config --cflags htslib 2>/dev/null || echo "")
HTSLIB_LIBS=$(shell pkg-config --libs htslib 2>/dev/null || echo "-lhts -lz")

## new version of CRISP that calls genotypes using likelihood and MLE method
all:
	mkdir -p include bin;
	$(MAKE) -C . new
	#$(MAKE) -C . crisp
	#$(MAKE) -C . indel

new: include/readfasta.o include/bamsreader.o include/variant.o include/allelecounts.o include/chisquare.o include/contables.o include/newcrispcaller.o optionparser.c FET/lowcovFET.c include/bamread.o
	$(CC) -DPICALL=3 -Iparsebam $(HTSLIB_CFLAGS) -g -O2 include/readfasta.o include/variant.o include/allelecounts.o include/bamsreader.o include/newcrispcaller.o include/chisquare.o include/contables.o include/bamread.o -lm readmultiplebams.c -o bin/CRISP.binary $(HTSLIB_LIBS)


## this is the older version of CRISP that does not call genotypes
crisp: include/readfasta.o include/bamsreader.o include/variant.o include/allelecounts.o include/chisquare.o include/contables.o include/crispcaller.o include/bamread.o optionparser.c
	$(CC) -DPICALL=0 -Iparsebam $(HTSLIB_CFLAGS) -g -O2 include/readfasta.o include/variant.o include/allelecounts.o include/bamsreader.o include/crispcaller.o include/chisquare.o include/contables.o include/bamread.o -lm readmultiplebams.c -o bin/CRISP-AF $(HTSLIB_LIBS)

include/crispcaller.o:	crisp/crispcaller.c crisp/crispcaller.h FET/contables.h variant.h bamsreader.h FET/pooledFET.c crisp/crispprint.c
	$(CC) $(HTSLIB_CFLAGS) -c crisp/crispcaller.c -o include/crispcaller.o

## picall replacement to call indels only
indel: include/readfasta.o include/bamsreader.o include/variant.o include/allelecounts.o include/chisquare.o include/contables.o include/newcrispcaller.o optionparser.c FET/lowcovFET.c include/bamread.o
	$(CC) -DPICALL=3 -Iparsebam $(HTSLIB_CFLAGS) -g -O2 include/readfasta.o include/variant.o include/allelecounts.o include/bamsreader.o include/newcrispcaller.o include/chisquare.o include/contables.o include/bamread.o -lm readmultiplebams.c -o bin/CRISPindel $(HTSLIB_LIBS)

### shared objects ###

include/newcrispcaller.o:	crisp/newcrispcaller.c crisp/crispcaller.h FET/contables.h parsebam/variant.h parsebam/bamsreader.h parsebam/bamread.h crisp/crispEM.c FET/pooledFET.c crisp/newcrispprint.c
	$(CC) -Iparsebam $(HTSLIB_CFLAGS) -c crisp/newcrispcaller.c -o include/newcrispcaller.o

include/variant.o: parsebam/variant.c parsebam/variant.h parsebam/bamsreader.h parsebam/bamsreader.c parsebam/readfasta.h parsebam/readfasta.c parsebam/process_indel_variant.c
	$(CC) $(HTSLIB_CFLAGS) -c parsebam/variant.c -o include/variant.o

include/bamsreader.o: parsebam/readfasta.h parsebam/readfasta.c parsebam/bamsreader.h parsebam/bamsreader.c parsebam/bamread.h parsebam/bamread.c
	$(CC) $(HTSLIB_CFLAGS) -c parsebam/bamsreader.c -o include/bamsreader.o

include/bamread.o: parsebam/readfasta.h parsebam/readfasta.c parsebam/bamread.h parsebam/bamread.c
	$(CC) $(HTSLIB_CFLAGS) -c parsebam/bamread.c -o include/bamread.o

include/allelecounts.o: parsebam/readfasta.h parsebam/allelecounts.h parsebam/allelecounts.c parsebam/variant.h
	$(CC) $(HTSLIB_CFLAGS) -c parsebam/allelecounts.c -o include/allelecounts.o

include/readfasta.o: parsebam/readfasta.c parsebam/readfasta.h
	$(CC) -c parsebam/readfasta.c -o include/readfasta.o

include/chisquare.o: FET/chisquare.h FET/chisquare.c
	$(CC) -c FET/chisquare.c -o include/chisquare.o

include/contables.o: FET/contables.c FET/contables.h FET/chisquare.h FET/chisquare.c
	$(CC) -c FET/contables.c -o include/contables.o

clean:
	rm -rf include bin
