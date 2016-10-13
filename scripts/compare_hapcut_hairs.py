#! /usr/bin/env python
import sys, os, glob, string, subprocess,time,math

# CODE created april 6 2012, last updated 08/27/2016
# compare phasing from hapcut output file to the original HAIRS file and output inconsistent Hairs and summary statistics 
# code is useful to assess what is the mode of errors in the input data, single sequencing errors, switch errors or a combination

"""
pipeline for indels:
1. include indels in original VCF file, extracthairs only SNPs, phase with hapcut 
2. extracthairs with indels as well.... 
3. for each indel, get the hairs that cover that indel and analyze consistency for two haplotypes using hapcut output file from SNPs alone 

in current implementation, unphased variants are absent in hapcut output file... need to use VCF file or something else...
weak SNPs and indels can be added onto strong haplotype scaffold, in some cases, this may not work 
two variant reads are sufficient to call SNP since consistency will 

"""

QVoffset = 33;

def read_hairs_file(input_file):
	hairlist = []; calls =0;
	## list of haplotype informative reads
	hairsfile = open(input_file,'r');
	for line in hairsfile:
		hair = line.strip().split(); qualitystr = hair[-1]; calls += len(qualitystr);
		varlist = []; q =0;
		for i in xrange(2,len(hair)-1,2):
			offset = int(hair[i]);
			for j in xrange(len(hair[i+1])): varlist.append([offset+j,hair[i+1][j],ord(qualitystr[q])-QVoffset]); q+=1; 
		varlist.append(line.strip()); ## last value is the original fragment string
		hairlist.append(varlist);
	hairsfile.close();
	print >>sys.stderr, "finished reading haplotype fragments file",len(hairlist),calls;
	return hairlist; 


def read_phased_vcf(input_file):
	phasedvariants = {}; varid = 1; blocks =0; phasedvars = 0;
	vcf = open(input_file,'r');
	for line in vcf:
		if line[0] == '#': continue; 
		else:
			var = line.strip().split('\t'); chrom = var[0]; position = int(var[1]); genotype = var[9].split(':'); 
			if genotype[0][1] != '/' and genotype[0][0] !='.' and genotype[0][2] !='.': 
				if var[8].split(':')[1] == 'BL': block = genotype[1]; 
				else: block = 1;
				phasedvariants[varid] = [position,genotype[0][0],genotype[0][2],var[3],var[4],var[9],block,[]];
				phasedvars +=1; 
			varid +=1; 
	vcf.close();
	print >>sys.stderr, "read",blocks,"blocks from VCF file",sys.argv[2],'and variants phased',phasedvars,varid;
	return phasedvariants;


def read_hapcut_output(input_file):
	phasedvariants = {}; 
	blocks=0; phasedvars=0;
	## read hapcut output file, block format
	hapcut = open(input_file,'r');
	for line in hapcut:
		if 'BLOCK' in line: blockinfo = line.strip(); var = line.split(); offset = int(var[2]); blocks +=1;
		elif '****' in line: pass;
		else:
			var = line.split(); 
			if var[1] != '-':
				phasedvariants[int(var[0])] = [var[4],var[1],var[2],var[5],var[6],var[7],var[8],blocks-1,offset,[]];
				phasedvars +=1; 
			offset +=1; 
	hapcut.close();
	print >>sys.stderr, "read",blocks,"blocks from file",sys.argv[2],'and variants phased',phasedvars;
	return phasedvariants;


def compare_blocks_hairs(hairlist,phasedvariants,pflag):
	total = 0; mismatches=0; STATS = [0.0,0.0,0,0.0,0.0,0.0,0]; 
	hairstats = []; 
	for t in xrange(len(hairlist)):
		hair = hairlist[t];
		h1 =0; h2 = 0;variants = 0; switchlist = []; m=0; switches=0;
		for i in xrange(len(hair)-1):
			try: 
				phasedvar = phasedvariants[hair[i][0]]; phasedvar[-1].append(t); ## list of hairs for each variant 
				if phasedvar[1] == hair[i][1]: 
					h2 +=1;
					if m==0: m = 1; 
					elif m == -1: 
						m = 1; switchlist.append(variants); # switch in phase
				else: 
					h1 +=1; 
					if m ==0: m = -1; 
					elif m ==1: m = -1; switchlist.append(variants);
				variants +=1;
			except KeyError: pass; 

		total += variants; switches = len(switchlist); ## position where phase between fragment and haplotype changes
		
		sw = 0; bp = 0; bp2=0; i=0; editstring = []; 
		while i < switches:
			if i < switches-1 and switchlist[i]+1 == switchlist[i+1]: bp +=1; editstring.append(i); i +=2;
			elif i < switches-1 and switchlist[i]+2 == switchlist[i+1]: bp2 +=1; editstring.append(i); i +=2;
			elif switchlist[i] ==1 or switchlist[i] == variants-1: bp +=1; editstring.append(i); i +=1; # count first or last vairant switch as bitflip
			else: sw +=1; editstring.append(-i); i +=1; 
		STATS[0] += sw; STATS[1] += bp; STATS[2] += variants-1; STATS[6] +=1; 

		hairstats.append([variants,sw,bp,bp2,h1,h2]); ## add scores to fragment 

		if sw == 0 and bp ==0 and bp2==0: continue; 
		elif sw == 0 and (bp ==1 or bp2 ==1) : STATS[3] +=1; continue; 
		elif sw ==1 and bp == 0 and bp2==0: STATS[4] +=1;
		else: STATS[5] +=1; 

		if pflag: ## print erroneous fragments
			print '\n','hair',hair;
			print 'STATS:',h1,h2,switches,'variants:',variants,'bitflips',bp,'bp-2',bp2,'switches',sw,switchlist,
			if h1 > 2 and h2 > 2 and switches < 2: print 'SWITCHERR',
			for i in xrange(len(hair)-1):
				try: 
					phasedvar = phasedvariants[hair[i][0]];
					print `hair[i][0]`+':'+phasedvar[1]+':'+hair[i][1] + ':' + `hair[i][2]`,
				except KeyError: pass;
			print;

	print >>sys.stderr,STATS[0:2],STATS[0]/STATS[2],STATS[1]/STATS[2];
	print >>sys.stderr,"fragment-errors","single-flips",STATS[3],"single-switch",STATS[4],"complex",STATS[5],"total",STATS[6];
	return hairstats;

################################################################################################


## for each variant, consider the hairs that cover that variant and assess consistency 
def phase_variant(hairlist,hairstats,phasedvariants):

	#phasedvariants[int(var[0])] = [var[0],var[1],var[2],blocks-1,offset];
	for variant in phasedvariants.iterkeys():
		info = phasedvariants[variant]; 
		hairs = info[-1]; 
		print info[0:5],len(hairs);
		for h in hairs: 	
			allele = '-'
			for i in xrange(len(hairlist[h])-1): 
				if variant == hairlist[h][i][0]: allele = hairlist[h][i][1];
			if hairstats[h][4] == 0 : hap = '0'; 
			elif hairstats[h][5] ==0: hap = '1'; 
			else: hap = '-';
			print variant,allele,hap,'|',hairstats[h],'|',hairlist[h][-1];
		print '\n';



#################################################################################################

if len(sys.argv) < 3: print 'python compare.py hairsfile genome.hapcut.phase'; sys.exit();

if 'vcf' in sys.argv[2] or 'VCF' in sys.argv[2]: 
	phasedvariants = read_phased_vcf(sys.argv[2]); hairlist = read_hairs_file(sys.argv[1]); 
	hairstats = compare_blocks_hairs(hairlist,phasedvariants,0);

else: 
	phasedvariants = read_hapcut_output(sys.argv[2]); hairlist = read_hairs_file(sys.argv[1]); 
	hairstats = compare_blocks_hairs(hairlist,phasedvariants,0);

## compare indels
#phase_variant(hairlist,hairstats,phasedvariants);
