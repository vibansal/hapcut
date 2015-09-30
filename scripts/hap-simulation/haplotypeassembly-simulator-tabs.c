// AUTHOR VIKAS BANSAL: last edited June 4 2009, Dec 8 2010
/* PROGRAM TO ESTIMATE HAPLOTYPE LENGTHS FROM SEQUENCING DATA */ 

#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<string.h>
#include<stdint.h>

struct SNP
{
	int* elist; int edges; int cc; int csize; int marked; int allele1,allele2; 
	int first,last;
};

struct hapblock { int first, last, size, span, adjustedspan; };

struct LIBRARY { int coverage; int readlength; int mean; int std; uint64_t reads; float fraction; /* fraction of total coverage */ };  // paired-end library with given mean, std and depth of coverage (x size of genome)


int compare_span_hapblock(const void *a,const void *b)
{
	const struct hapblock *h1 = (const struct hapblock*)a;
	const struct hapblock *h2 = (const struct hapblock*)b;
	if (h2->span > h1->span) return 1; // flipped this on Dec20/10 so that the blocks are sorted in descending order
	if (h1->span == h2->span) return 0;
	return -1; 
}
int compare_size_hapblock(const void *a,const void *b)
{
	const struct hapblock *h1 = (const struct hapblock*)a;
	const struct hapblock *h2 = (const struct hapblock*)b;
	if (h2->size > h1->size) return 1;
	if (h1->size == h2->size) return 0;
	return -1; 
}

int compare_adjustedspan_hapblock(const void *a,const void *b)
{
	const struct hapblock *h1 = (const struct hapblock*)a;
	const struct hapblock *h2 = (const struct hapblock*)b;
	if (h2->adjustedspan > h1->adjustedspan) return 1;
	if (h1->adjustedspan == h2->adjustedspan) return 0;
	return -1; 
}

float normal_variable(float mean, float std)
{
	float a = drand48(); float b= drand48(); float d = -2/log10(2.7182);
	float c = sqrt(log10(a)*d)*cos(2*3.14159*b); 
	c = c*std + mean;
	return c;
}

int DFS(int v, struct SNP* SNPgraph, int cc)
{
	if (SNPgraph[v].marked ==1) return 1;
	SNPgraph[v].marked = 1; SNPgraph[v].cc = cc; SNPgraph[cc].csize++; 
	if (v < SNPgraph[cc].first || SNPgraph[cc].first == -1) SNPgraph[cc].first = v; 								
	if (v > SNPgraph[cc].last || SNPgraph[cc].last == -1) SNPgraph[cc].last = v; 								
	int i=0;
	for (i=0;i<SNPgraph[v].edges;i++) 
	{
		if (SNPgraph[SNPgraph[v].elist[i]].marked  ==0) DFS(SNPgraph[v].elist[i],SNPgraph,cc); 
	} return 1; 
}


//////////////////////////////////////////// MAIN FUNCTION /////////////////////////////////////////////////////////////////


int read_snplocs(char* snpfile,int MATES,int coverage,struct LIBRARY* libraries, int ls)  
{
	time_t ts; time(&ts); srand48((long int)ts); 	
	FILE* sf;   char ch ='0';  int snps=0;

	sf = fopen(snpfile,"r"); if (sf == NULL) {fprintf(stdout,"error opening file \n"); exit(0); } 
	while (1) {  ch = fgetc(sf); if (ch == EOF) break; if (ch == '\n') snps++; } fclose(sf);

	int* snp_array = (int*)malloc(4*snps);
	struct SNP* SNPgraph = (struct SNP*)malloc(sizeof(struct SNP)*snps);
	char snpid[50];  char het[50]; int l1,l2; char c; char repeat[200];  char method[100];	int chr; int i=0,j=0,k=0,l=0;

	sf = fopen(snpfile,"r");
	for (i=0;i<snps;i++) 
	{
		fscanf(sf,"%d %s %s %d %d \n",&chr,snpid,het,&l1,&l2); //&c,&c,repeat,method);fprintf(stdout,"chr %d id %s loc %d \n",chr,snpid,l1); 
		snp_array[i] = l1; SNPgraph[i].cc = i; SNPgraph[i].csize =0; SNPgraph[i].edges =0;
		SNPgraph[i].allele1 =0; SNPgraph[i].allele2 = 0; 	SNPgraph[i].first = -1; SNPgraph[i].last = -1;
	} fclose(sf);

	int* snpmap; // map 1000bp of chromosome to first snp in that region indexed by snp_array 
	int blocks = (int)(snp_array[snps-1]/1000)+2; 
	snpmap = (int*)malloc(4*blocks); for (i=0;i<blocks;i++) snpmap[i] = -1;
	for (i=0;i<snps;i++)
	{
		j = (int)(snp_array[i]/1000); if (snpmap[j] == -1) snpmap[j] = i; 
	}
	for (i=1;i<blocks;i++) { if (snpmap[i] == -1) snpmap[i] = snpmap[i-1];} 

	//for (i=0;i<100;i++) fprintf(stdout,"snp locs %d %d \n",i,snp_array[i]);
	//for (i=0;i<100;i++) fprintf(stdout,"snp map %d %d \n",i,snpmap[i]);

	/************************************ CODE FOR SIMULATING PAIRED-READS ************************************/

	// simulate shotgun reads of length 1kb with mate pairs with certain span distribution
	//int readlength = 2000; int span1 = 3000, span2 = 1500, span3 = 36000,span;
	//int span1 = IS1, span2 = IS2, span3 = IS3, span; 
	//int no_reads = (int)(snp_array[snps-1]/readlength); no_reads *= coverage/2; 
	//if (MATES ==0) no_reads *=2; // for single reads no_reads should  be multiplied by 2 ?? 

	int wsnps=0; int good_reads=0;
	int start2,end2,ss,es, chromlength = snp_array[snps-1],start1,end1;
	int snpcoverage=0; int mincc; int MAX = 5000;
	int* snplist = (int*)malloc(sizeof(int)*MAX);
	int components =0,singletons=0,hap;
	srand(time(NULL));
	float random = rand(); 
	int pid = getpid();
	char fragfile[100]; sprintf(fragfile,"fragment.list.%d",pid);
	FILE* flist = fopen(fragfile,"w"); 
	
	uint64_t totalbases = (uint64_t)chromlength; totalbases *= coverage; 
	uint64_t lsreads =0,r=0;
	//fprintf(stderr,"snps %d %d\n",snps,chromlength); 

	for (l=0;l<ls;l++)
	{
		//lsreads = (uint64_t)(totalbases*libraries[l].fraction); lsreads /= libraries[l].readlength*2; libraries[l].reads = lsreads;
		lsreads = (uint64_t)(chromlength)*libraries[l].coverage; lsreads /= libraries[l].readlength*2; libraries[l].reads = lsreads;
		if (MATES ==0) libraries[l].reads *=2; lsreads *=2;
		fprintf(stderr,"reads for library %d %d mean %d std %d coverage %d #ofreads %ld\n",l+1,libraries[l].readlength,libraries[l].mean,libraries[l].std,libraries[l].coverage,libraries[l].reads);

		for (r=0;r<libraries[l].reads;r++) // iterate over the number of reads for this library 
		{
			start1 = lrand48()%chromlength; if (start1 > chromlength-50000) start1 = start1-50000;
			end1 = start1 + libraries[l].readlength; 
			j = (int)(start1/1000); ss = snpmap[j]; wsnps=0;
			if (ss > 0)
			{
				while(snp_array[ss] < start1) ss++; 
				while (snp_array[ss] <= end1 && wsnps < MAX )
				{
					snplist[wsnps] = ss;
					wsnps++; ss++; 
				}	
			}
			if (MATES > 0)
			{
				start2 = end1 + (int)(normal_variable(libraries[l].mean,libraries[l].std)) - libraries[l].readlength; end2 = start2 + libraries[l].readlength; 
				j = (int)(start2/1000); ss = snpmap[j];
				if (ss > 0)
				{
					while(snp_array[ss] < start2) ss++; 
					while (snp_array[ss] <= end2 && wsnps < MAX)
					{
						snplist[wsnps] = ss;
						wsnps++; ss++; 
					}
				}
			}
			if (drand48() <= 0.5) hap =0; else hap=1;
			if (wsnps > 0 && hap ==0) 
			{
				for (j=0;j<wsnps;j++) SNPgraph[snplist[j]].allele1++; 
			}
			if (wsnps > 0 && hap ==1) 
			{
				for (j=0;j<wsnps;j++) SNPgraph[snplist[j]].allele2++; 
			}
			if (wsnps > 1) good_reads++; 
			if (wsnps > MAX) fprintf(stderr,"wsnp %d is greater than MAX, need to increase value of MAX \n",wsnps);
			if (wsnps > 1)
			{
				snpcoverage += wsnps;
				fprintf(flist,"%d ",wsnps); 
				for (j=0;j<wsnps;j++) fprintf(flist,"%d ",snplist[j]); fprintf(flist,"\n"); 
				for (j=0;j<wsnps;j++) 
				{
					for (k=j+1;k<wsnps;k++)
					{	
						SNPgraph[snplist[j]].edges++; 	SNPgraph[snplist[k]].edges++; 	
					}
				}
			}
		} 
	}

	fclose(flist);
	flist = fopen(fragfile,"r");  
	for (i=0;i<snps;i++)
	{
		if (SNPgraph[i].edges > 0) SNPgraph[i].elist = (int*)malloc(4*SNPgraph[i].edges);
		SNPgraph[i].edges =0;
	}
	for (i=0;i<good_reads;i++)
	{
		fscanf(flist,"%d ",&wsnps); 
		for (j=0;j<wsnps;j++) snplist[j] = 0;
		for (j=0;j<wsnps;j++) fscanf(flist,"%d ",&snplist[j]); fscanf(flist,"\n");
		for (j=0;j<wsnps;j++) 
		{
			for (k=j+1;k<wsnps;k++)
			{	
				SNPgraph[snplist[j]].elist[SNPgraph[snplist[j]].edges] = snplist[k];
				SNPgraph[snplist[j]].edges++;
				SNPgraph[snplist[k]].elist[SNPgraph[snplist[k]].edges] = snplist[j];
				SNPgraph[snplist[k]].edges++;
			}
		}
	}
	fclose(flist);

	struct hapblock* hapblocklist = (struct hapblock*)malloc(sizeof(struct hapblock)*snps);
	int hapblocks =0; int distancespanned =0;

	for (i=0;i<snps;i++) SNPgraph[i].marked =0;
	for (i=0;i<snps;i++)
	{
		if (SNPgraph[i].marked ==0) DFS(i,SNPgraph,i); 
	}
	for (i=0;i<snps;i++)
	{
		if (SNPgraph[i].csize ==1) singletons++; 
		else if (SNPgraph[i].csize > 1)	
		{
			components++; 
		}
		if (SNPgraph[i].csize >=1) 
		{
			hapblocklist[hapblocks].first = SNPgraph[i].first;
			hapblocklist[hapblocks].last = SNPgraph[i].last;
			hapblocklist[hapblocks].size = SNPgraph[i].csize;
			distancespanned =snp_array[SNPgraph[i].last]-snp_array[SNPgraph[i].first];
			hapblocklist[hapblocks].span = distancespanned;
			hapblocklist[hapblocks].adjustedspan = distancespanned/(SNPgraph[i].last-SNPgraph[i].first+1);
			hapblocklist[hapblocks].adjustedspan *= SNPgraph[i].csize;

			//			if ( hapblocklist[hapblocks].adjustedspan != hapblocklist[hapblocks].span) fprintf(stdout,"badcomp SNPinfo %d loc %d edges %d cc %d size %d A1 %d A2 %d first %d last %d span %d snpspan %d LENGTH %d \n",i,snp_array[i],SNPgraph[i].edges,SNPgraph[i].cc,SNPgraph[i].csize,SNPgraph[i].allele1,SNPgraph[i].allele2,snp_array[SNPgraph[i].first],snp_array[SNPgraph[i].last],snp_array[SNPgraph[i].last]-snp_array[SNPgraph[i].first],SNPgraph[i].last-SNPgraph[i].first+1,(distancespanned*SNPgraph[i].csize)/(SNPgraph[i].last-SNPgraph[i].first+1)); 
			hapblocks +=1; 
		}
	}
	//fprintf(stdout,"no of hapblocks %d \n",hapblocks);
	// maybe we can the full analysis in C only 
	// declare an array of size SNPs with each entry having 'size-of-CC', span, reduced-span 
	int N50length_variants =0, N50length_adjustedspan=0,N50length_span=0; 
	double sum =0; double adjustedGS=0; double spanGS =0; int variantGS =0;
	qsort(hapblocklist,hapblocks,sizeof(struct hapblock),compare_adjustedspan_hapblock);
	for (i=0;i<hapblocks;i++)  { adjustedGS += hapblocklist[i].adjustedspan; }
	for (i=0;i<hapblocks;i++) 
	{
		sum += hapblocklist[i].adjustedspan;
		if (sum*2 > adjustedGS) 
		{
			N50length_adjustedspan  = hapblocklist[i].adjustedspan; i = hapblocks;
		}
	}
	sum =0;	qsort(hapblocklist,hapblocks,sizeof(struct hapblock),compare_span_hapblock);
	for (i=0;i<hapblocks;i++)  spanGS += hapblocklist[i].span;
	for (i=0;i<hapblocks;i++) 
	{
		sum += hapblocklist[i].span;
		if (sum*2 > spanGS) 
		{
			N50length_span  = hapblocklist[i].span; 	i = hapblocks; 
		}
	}

	sum =0; qsort(hapblocklist,hapblocks,sizeof(struct hapblock),compare_size_hapblock);
	for (i=0;i<hapblocks;i++)  variantGS += hapblocklist[i].size;
	for (i=0;i<hapblocks;i++) 
	{
		sum += hapblocklist[i].size;
		if (sum*2 > variantGS) 
		{
			N50length_variants  = hapblocklist[i].size; 	i = hapblocks; 
		}
	}
	//for (i=0;i<hapblocks;i++) fprintf(stdout,"block %d span %d snps %d adjusted %d last %d first %d\n",i,hapblocklist[i].span,hapblocklist[i].size,hapblocklist[i].adjustedspan,hapblocklist[i].last,hapblocklist[i].first); 
	//fprintf(stdout,"rlength PE cov IS1 IS2 %d IS3 %d IS1P %f IS2P %f N50-variants %d N50-span  %d N50-adjusted %d",readlength,MATES,coverage,IS1,IS2,IS3,clone1_fraction,clone2_fraction,N50length_variants,N50length_span,N50length_adjustedspan);
	//fprintf(stdout," components %d singletons %d \n",components,singletons); 
	coverage =0;
      char paramstring[100];

	// header
	//fprintf(stdout,"Input Params\t");
	//fprintf(stdout,"totalcoverage\tN50_vars\tN50_span\tN50_adjustedspan\t");
	//fprintf(stdout,"components\tphasedvariants\t");
	//fprintf(stdout,"adjustedGS\tspanGS\tactualsize\tvariantGS\n");

	// data
	for (l=0;l<ls;l++) 
	{
		fprintf(stdout,"RL %d IS %d cov %2d ",libraries[l].readlength,libraries[l].mean,libraries[l].coverage);
		coverage += libraries[l].coverage;
	}

      fprintf(stdout,"\t%3d\t%5d\t%6d\t%6d\t",coverage,N50length_variants,N50length_span,N50length_adjustedspan);
	fprintf(stdout,"%6d\t%6d\t",components,variantGS-singletons); 
	fprintf(stdout,"%f\t%f\t%d\t%d",adjustedGS,spanGS,chromlength,variantGS);
	

//fprintf(stderr,"%5d %1d %3d %5d %5d %5d %0.2f %0.2f %5d %6d %6d",readlength,MATES,coverage,IS1,IS2,IS3,clone1_fraction,clone2_fraction,N50length_variants,N50length_span,N50length_adjustedspan);
	//fprintf(stderr,"%6d %6d \n",components,singletons); 

	for (i=0;i<hapblocks;i++)
	{
		//	fprintf(stdout,"%d %d %d %d %d %d \n",i,hapblocklist[i].first,hapblocklist[i].last,hapblocklist[i].size,hapblocklist[i].span,hapblocklist[i].adjustedspan); 
	}
	free(snplist);

	/* CHECK WHY adjusted span != normal span for single reads....... May 5 2009 code modified to do qsort() and output N50 length directly from within C code itself. 
	*/

	char remove_cmd[100]; sprintf(remove_cmd,"rm -f %s",fragfile); system(remove_cmd);


	//fprintf(stderr,"reads %d good %d snpcoverage %d %f components %d singletons %d \n",no_reads,good_reads,snpcoverage,(float)snpcoverage/(float)snps,components,singletons); 
	//	fprintf(stderr,"snps %d chromlength %d reads %d good %d snpcoverage %d %f components %d singletons %d \n",snps,chromlength,no_reads,good_reads,snpcoverage,(float)snpcoverage/(float)snps,components,singletons); 
	//fprintf(stdout,"snps %d chromlength %d reads %d good %d snpcoverage %d %f components %d singletons %d \n",snps,chromlength,no_reads,good_reads,snpcoverage,(float)snpcoverage/(float)snps,components,singletons); 

	// HOW TO GET N50 haplotype length in kilobases from span of all connected components 
	// gcc -lm -o hapsim simulation_reads.c
	// ./hapsim chr1_HuRef_hetvariants_hetSNPs_24percent_homSNPs 150 1 20 3000 8000 20000 > a
	// ./hapsim chr1_HuRef_hetvariants_hetSNPs_24percent_homSNPs read_length mates_0/1 coverage IS1 IS2 IS3 IS1_pc IS2_pc  

	// sort the haplotype spans and add the number of SNPs in each connected component (cumulative dist) and take 50%
	// awk script: grep "SNPinfo" a | sort -g -k 20 | awk 'BEGIN { snps = 0; N50 = -1; } { snps += $10; if (snps > 153890/2 && N50 == -1 ) N50 = $20; } END { print N50; }'

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


int main (int argc, char** argv)
{
//	float m=10000; float s= 1000; float r= 0; int i=0; float sum =0;for (i=0;i<100000;i++) { r = normal_variable(m,s); printf("%f \n",r);  sum += r; } 
	// simulating reads of various lengths and paired ends for haplotype assembly 
	// fprintf(stdout,"random %f \n",normal_variable(10,0.5));  exit(0);
	// gcc -lm -o hapsim simulation_reads.c; ./hapsim chr1_HuRef_hetvariants_hetSNPs_24percent_homSNPs 100 1 30 3000 8000 20000 0 0
	int a =0; char snpfile[100]; int readlength =100; int MATES = 1; int IS1= 500, IS2=2000, IS3 =10000; int coverage=30;
	/*
	for (a=1;a<argc;a+=2)
	{
	  if (strcmp(argv[a],"--snpfile") ==0) strcpy(snpfile,argv[a+1]);
	  else if (strcmp(argv[a],"--cov") ==0) coverage = atoi(argv[a+1]);
	  else if (strcmp(argv[a],"--IS") ==0) coverage = atoi(argv[a+1]);
	}*/
	
	if (argc <= 5) 
	{
		for (a=0;a<100;a++) fprintf(stdout,"-");
		fprintf(stdout,"\n\nProgram for estimating N50 haplotype length using simulated reads\nProgram allows paired-end reads with an arbitrary number of insert libraries \n \n");
		fprintf(stdout,"input format is FILE.SNPlocations coverage RL1,IS1,CV1 RL2,IS2,CV2 .... readlength (bp) InsertSize (bp) depth of physical coverage using this library \n\n 1. snp file used to simulate heterozygous variant sites \n 2. total coverage, not required \n 3... three tuples corresponding to (readlength,insert_size,depth_of_coverage) \n"); 
		fprintf(stdout," ./new YRI.1000genomes.chrom1.snplocs 40 75 500 10 75 2500 70 75 9000 20 \n\n");
		for (a=0;a<100;a++) fprintf(stdout,"-"); 
		fprintf(stdout,"\n\n");
		exit(0); 
	}

	strcpy(snpfile,argv[1]); coverage = atoi(argv[2]); int ls = argc/3-1; 
	fprintf(stderr,"snpfile %s coverage %d libraries %d \n",snpfile,coverage,ls);
    struct LIBRARY* libraries = (struct LIBRARY*)malloc(sizeof(struct LIBRARY)*ls);
	for (a=0;a<ls;a++) 
	{
	  libraries[a].readlength = atoi(argv[3+3*a]); libraries[a].mean =  atoi(argv[3+3*a+1]);
		if (libraries[a].mean < 1) MATES =0;
	  libraries[a].coverage = (atoi(argv[3+3*a+2])); libraries[a].std = libraries[a].mean/10; 
	  //libraries[a].fraction = (float)(atoi(argv[3+3*a+2]))/100; libraries[a].std = libraries[a].mean/10; 
	}

	read_snplocs(argv[1],MATES,coverage,libraries,ls); return 1;

	// 3kb clones with light coverage of 8kb and 20kb clones improves it tremendously
	// mix of two clones is better especially at low coverage
	// even with 50bp reads and read lengths 1-10 kb we can get long haplotype lengths
	// plot 1: N50 vs coverage for 5 different mixes, 3kb, 8kb, 20kb, 3kb+8kb, 3kb+light sprinkling
	// N50 vs coverage for three different read lengths 
	// N50 as function of fraction of clones for 3kb + 8kb (20x 40x) same for 8kb +20kb 
	// talk about more DNA is required to sythesise 20kb clones 
	// for neighboring SNPs, use haplotypes to distinguish from error 

	//read_snplocs(argv[1],readlength,MATES,coverage,IS1,IS2,IS3,clone1_fraction,clone2_fraction); 

}


