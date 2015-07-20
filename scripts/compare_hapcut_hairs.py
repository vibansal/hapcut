#! /usr/bin/env python
import sys, os, glob, string, subprocess,time,math

# CODE created april 6 2012 
# compared phasing from hapcut output file to the original HAIRS file and output for each HAIR whether it is consistent or inconsistent with phasing

if len(sys.argv) < 3: print 'python compare.py hairsfile genome.hapcut.phase'; sys.exit();

phasedvariants = {}; blocklist = [];
blocks=0; phasedvars=0;

hapcut = open(sys.argv[2],'r');

for line in hapcut:
	if 'BLOCK' in line: 
		blockinfo = line.strip(); var = line.split(); offset = int(var[2]);
		blocklist.append(var); blocks +=1;
	elif '****' in line: pass;
	else:
		var = line.split(); 
		if var[1] != '-':
			phasedvariants[offset] = [var[0],var[1],var[2],blocks-1];
			phasedvars +=1; 
		offset +=1; 

hapcut.close();
print >>sys.stderr, "read",blocks,"blocks from file",sys.argv[2],'and variants phased',phasedvars;

################################################################################################

hairlist = [];
hairs = open(sys.argv[1],'r');
for line in hairs:
	hair = line.strip().split();
	varlist = [int(hair[2])]; 
	for i in xrange(2,len(hair)-1,2):
		offset = int(hair[i]);
		for j in xrange(len(hair[i+1])): varlist.append([offset+j,hair[i+1][j]]); 
	varlist.append(line.strip());
	#print hair,varlist;
	hairlist.append(varlist);
hairs.close();

################################################################################################

hairlist.sort();

for hair in hairlist: 
	h1 =0; h2 = 0; bases = 0;
	#print hair[-1],
	for i in xrange(1,len(hair)-1,1):
		try: 
			phasedvar = phasedvariants[hair[i][0]];
			if phasedvar[1] == hair[i][1]: h2 +=1;
			else: h1 +=1; 
			bases +=1;
			#print phasedvar,
		except KeyError: pass; 
	if h1 == 0: pass; #print 'match';
	elif h2 ==0: pass; #print 'compmatch';
	else: 
		print hair[-1],'|',
		print h1,h2,
		for i in xrange(1,len(hair)-1,1):
			try: 
				phasedvar = phasedvariants[hair[i][0]];
				print phasedvar[0]+':'+phasedvar[1]+'|'+phasedvar[2],
			except KeyError: 'NotFound';
		print;



