"""
this script converts heterozygous variants file previously used for HapCUT to a VCF file
input: hetsnps file for old version of hapCUT
output: VCF file used for new version of hapCUT 
author Vikas Bansal, last modified march 7 2012
"""
import sys, os, glob, string


print '##fileformat=VCFv4.0';
print '##source=createdforHapCUT';
print '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE';

File = open(sys.argv[1],'r');
for line in File:
		variant = line.strip().split();
		genotype = variant[5];
		if genotype[0] == variant[3] and genotype[2] == variant[3]: vcfgenotype = '0/0';
		elif genotype[0] == variant[4] and genotype[2] == variant[4]: vcfgenotype = '1/1';
		elif genotype[0] == variant[3] and genotype[2] == variant[4]: vcfgenotype = '0/1';
		elif genotype[0] == variant[4] and genotype[2] == variant[3]: vcfgenotype = '0/1';
		else: vcfgenotype = './.';
		
		print '%s\t%d\t.\t%s\t%s\t%d\tPASS\tNS=1\tGT\t%s' %(variant[1],int(variant[2]),variant[3],variant[4],int(variant[6]),vcfgenotype);
File.close();
		
