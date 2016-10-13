
#CC=gcc -Wall
CC=gcc -D_GNU_SOURCE
CFLAGS=-c -Wall -lz
SAMTOOLS=../parsebam/samtools-0.1.18
PARSE_BAM=../parsebam/

all:	
	$(MAKE) -C $(PARSE_BAM) hairs # compile code in extracthairs directory for hapCUT
	$(MAKE) -C . fosmid

fosmid: extracthairs-dev.c fosmidbam_hairs.c print_clusters.c
	$(CC) -I$(PARSE_BAM) -I$(SAMTOOLS) -g -O2 $(PARSE_BAM)/bamread.o $(PARSE_BAM)/hapfragments.o $(PARSE_BAM)/hashtable.o $(PARSE_BAM)/readfasta.o $(PARSE_BAM)/readvariant.o -o extractFOSMID extracthairs-dev.c  -L$(SAMTOOLS) -lbam -lm -lz

clean:
	rm -f extractFOSMID $(PARSE_BAM)/*.o
