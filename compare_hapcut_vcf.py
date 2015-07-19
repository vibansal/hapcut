#! /usr/bin/env python
import sys, os, glob, string, subprocess,time,math

# author Vikas Bansal, last updated July 14 2015 | This script performs three tasks: 


# 1. compare phasing from VCF file with phasing from hapcut output file 
# 2. output number of unphased variants that can now be linked to phased variants 
# 3. output new combined phased VCF file 

if len(sys.argv) < 3: print 'python compare.py phased.vcf hapcut.phase'; sys.exit();
snptable = {}; snps =0; indeltable = {}; phased =0; unphased=0; newphased=0;


########################## read VCF file ###############################

VCF_file = open(sys.argv[1],'r');
for line in VCF_file: 
	if line[0]== '#': continue;
	var = line.split(); 
	chrom = var[0]; position = var[1]; rsid = var[2]; allele1 = var[3]; alleles = var[4].split(','); allele2 = alleles[0]; 
	genotypes = var[9].split(':'); 

	if genotypes[0] == '0|1': snptable[(var[0],int(var[1]),'c')] = [0,1,genotypes,1]; phased +=1; 
	elif genotypes[0] == '1|0': snptable[(var[0],int(var[1]),'c')] = [1,0,genotypes,1]; phased +=1; 
	elif genotypes[0] == '0/1' or genotypes[0] == '1/0': snptable[(var[0],int(var[1]),'c')] = [-1,-1,genotypes,0]; unphased +=1; # in this case, genotype is unphased 
	else: pass; 

VCF_file.close();

print >>sys.stderr, 'comparing phase from VCF file to HAPCUT phase','unphased hets',unphased,'phased hets',phased;


################## read hapcut output file for phased genome ##################

FINAL_STATS = [0,0,0,0]; 
additionalphased=0; blocklist = [];
hapcut = open(sys.argv[2],'r');
for line in hapcut:
	if '****' in line: # end of previous block, print statistics 
		FINAL_STATS[0] += phased-1; FINAL_STATS[1] += switches
		FINAL_STATS[2] += phased-2; FINAL_STATS[3] += bitflips;

		print 'stats:',variants,'switcherrors',switches,'bitflips',bitflips;
		print blockinfo,'variants',phased+unphased,'unphasedhets',unphased,
		if newphased > 0 and phased > 0: 
			additionalphased += newphased; 
			print 'new-phased',newphased,'\n';
			#for snp in blocklist: snptable[snp][3] = 1;
		elif unphased > 0 and phased ==0: 
			print 'isolated block with',unphased,'hets\n';
		else: print '\n';
		#print blocklist;

	elif 'BLOCK' in line: # start of new block  
		blockinfo = line.strip();
		var = line.split(); blockid = var[2]; length = var[4]; 
		phased =0; unphased =0; variants = 0; newphased=0;
		blocklist = [];
		m = -1; bitflips = 0; switches = 0; runlength = 0;
	else:
		var = line.split(); #info = var[0].split('_'); 
		try: 
			snp = snptable[(var[3],int(var[4]),'c')];

			if snp[0] != -1 and var[1] != '-': ## phased in VCF file and hapcut output
				phased +=1; 
				blocklist.append([snp[0],snp[1],var]);
				flag = 0;
				if variants ==0 and snp[0] == int(var[1]): m = 0; runlength = 1;
				elif variants ==0 and snp[0] == int(var[2]): m = 1; runlength = 1; 
				elif m ==0 and snp[0] == int(var[1]): runlength +=1;
				elif m ==1 and snp[0] == int(var[2]): runlength +=1; 
				elif m ==0 and snp[0] == int(var[2]): ## change in phase
					if runlength > 0: switches +=1; m = 1; runlength = 0; flag = 1; 
					else: bitflips +=1; m = 1; runlength = 1; flag = 2; switches -=1; 
				elif m ==1 and snp[0] == int(var[1]): ## change in phase
					if runlength > 0: switches +=1; m = 0; runlength = 0; flag = 1; 
					else: bitflips +=1; m = 0; runlength = 1; flag = 2; switches -=1;
				variants +=1; 
				if flag > 0: print 'F' + `flag`,m,#runlength,
				else: print '--',m,#runlength,
				print line,

			else: # variant not phased in VCF file
				if var[1] != '-': newphased +=1; 
				unphased +=1; 
				#print 'SP',line,
				pass

		except KeyError: print line,
hapcut.close();
print >>sys.stderr, ' can phase',additionalphased,'additional variants in VCF file using sequence data';
print >>sys.stderr, 'FINAL_STATS_switches',FINAL_STATS[1],FINAL_STATS[0],float(FINAL_STATS[1])/FINAL_STATS[0],'bitflips:',FINAL_STATS[3],FINAL_STATS[2],float(FINAL_STATS[3])/FINAL_STATS[2];


