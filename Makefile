
#CC=gcc -Wall
CC=gcc -D_GNU_SOURCE
CFLAGS=-c -Wall
SAMTOOLS=parsebam/samtools-0.1.18

all:	
	$(MAKE) -C parsebam/samtools-0.1.18 all
	$(MAKE) -C parsebam hairs
	$(MAKE) -C hapcut-src HAPCUT
	cp hapcut-src/HAPCUT parsebam/extractHAIRS .;


clean:
	$(MAKE) -C parsebam clean
	$(MAKE) -C hapcut-src clean
	rm -f extractHAIRS HAPCUT 
