#include<stdint.h>
int MIN_CLUSTER_DISTANCE = 10000;
int MAX_ISLAND_LENGTH = 200000;

int BLOCK_SIZE = 100; // 100 bp blocks
int BF = 10; // only for making bigger blocks for background density estimation..
#define BS1 20

// code for parsing fosmid pool bam files
// ./extractFOSMID --bam fosmid-data/SRR799544_GGACTCCTAAGGAGTA.rmdup.bam.chr1 --VCF fosmid-data/NA12878.1kg.hets.chr1.liftoverhg19.vcf --fosmids 1 --ref fosmid-data/chr1.fa > na
// ./extractFOSMID --bam /media/drive2/Haplotyping/NA12878-SOLID-fosmid/aligned-bams/pool_SPA1.novoalign.sorted.bam --VCF /media/drive2/Haplotyping/NA12878-SOLID-fosmid/NA12878.1KG.2010_03.vcf --fosmid 1 --out SPA1.frags.smooth.200bp --mask /media/drive2/Haplotyping/genome-mask/hg18-50bp-files/hg18.50bp.mask.fa --ref ~/Public/tools/reference-genomes/1000genomes-huref/human_b36_male.fa > SPA1.log

// block of reads in a window of size 'x' bp
struct BLOCK_READ
{
	int index; int cluster; 
	int firstread; int lastread; short counts; float reads; 
	int start, end; 
	float GC; // GC percentage of window // GC content of the full fragment is important |  both GC-rich fragments and AT-rich fragments are underrepresented 
	float mappability; // mappability of short reads in this window 
	short subcounts[BS1]; // count of reads in smaller blocks of 10 base pairs each 
	short adjusted_size;
	//int poslist[60];
	
	double bscore;  int previous; int previous_cluster; int reads_window; double score_window;
	double sumL,sumR;
	double mean,variance;
	char valid; 
	float GC_corr;
	//int nv,pv; // next valid and previous valid block.... 
	// bed file with 100bp windows -> read count for unique and non-uniquely mapping reads 
};

struct FOSMID // information for each fragment...
{
	int firstblock,lastblock; int firstread,lastread,ws; int hets,unique_variants; 
	double mean,reads; 
	int* varlist;
	//FRAGMENT fragment; fragment.variants =0; fragment.alist = (allele*)malloc(sizeof(allele)*4096);

};

int cmp_fosmid(const void *a,const void *b)
{
	const struct FOSMID* fa = (const struct FOSMID*)a;  const struct FOSMID* fb = (const struct FOSMID*)b; 
	return fa->firstblock- fb->firstblock; 
}

#include "print_clusters.c"

int calculate_block_stats(struct BLOCK_READ* RL, int blocks, REFLIST* reflist,REFLIST* genomemask)
{
	int i=0,j=0,k=0,c1,c2,c; int X0=0,X1=0,GC=0,valid=0; int subsize = BLOCK_SIZE/BS1;

	float counts[20]; for (i=0;i<20;i++) counts[i]=0; // readcounts for GC % blocks (0-5,5-10,10-15....95-100)
	int GCbins[20]; for (i=0;i<20;i++) GCbins[i] = 0; 
	float totalcounts =0; int totalbins =0;

	// use 0/1 counts for sub-bins to smooth count function | calculate mappability of block and GC content 
	for (i=0;i<blocks;i++) 
	{
		if (genomemask->ns > 0) 
		{
			X0 =0; X1 = 0; k=0;
			for (j=0;j<BLOCK_SIZE;j++)
			{
				c = genomemask->sequences[genomemask->current][i*BLOCK_SIZE+j]- 63; 
				c1 = c>>3; c2 = c & 7; 
				c1 = c1? 1<<(c1-1) : 0; c2 = c2? 1<<(c2-1) : 0;  
				if (c1 == 1 && c2 < 2) X0++; 
				if (c1 > 1 || c2 > 1 ) X1++; 
				if ((j+1)%subsize==0)  
				{
					if (X1 > 0 && RL[i].subcounts[k] > 0) { RL[i].subcounts[k] = 0; } 
					if (X1 > 0) RL[i].adjusted_size -= subsize; 
					k++; X1 = 0; 
				}
			}
			RL[i].mappability = (float)X0/BLOCK_SIZE; 
		}

		// smoothed read counts	
		RL[i].reads = 0; 
		for (k=0;k<BS1;k++) 
		{
			if (RL[i].subcounts[k] ==1) RL[i].reads +=1;
			else if (RL[i].subcounts[k] > 1) RL[i].reads +=1;
		}

		if (reflist->ns > 0) 
		{
			GC = 0; RL[i].GC = 0;
			for (j=0;j<BLOCK_SIZE;j++)
			{
				if (reflist->sequences[reflist->current][i*BLOCK_SIZE+j] == 'G' || reflist->sequences[reflist->current][i*BLOCK_SIZE+j] =='C') GC +=1;
			}
			RL[i].GC = (float)GC/BLOCK_SIZE;
			if (RL[i].mappability > 0.99) { counts[(GC/5)] += RL[i].reads; GCbins[GC/5] +=1; } 
			//fprintf(stderr,"%d %d mappability %f GC %d %d %d \n",i*BLOCK_SIZE,i*BLOCK_SIZE+j,RL[i].mappability,GC,RL[i].reads,RL[i].adjusted_size);	
		}
		//if (RL[i].mappability < 0.8) RL[i].valid = '0';
		if (RL[i].adjusted_size < subsize) RL[i].valid = '0'; 
		if (RL[i].GC <= 0.25 || RL[i].GC >= 0.7) RL[i].valid = '0'; 
		else if (reflist->ns > 0) 
		{
			totalcounts += RL[i].reads; totalbins++; 
		}

		if (RL[i].valid == '1') valid++;
	
		//for (k=0;k<BS1;k++) fprintf(stdout,"%d:",RL[i+j].subcounts[k]);
		//fprintf(stdout," block %d reads %d %d first %d last %d \n",i+j,RL[i+j].reads,empty,RL[i+j].firstread,RL[i+j].lastread);
	}
	double maxGC = 0; double meanRD = totalcounts/totalbins;
	for (i=1;i<19;i++) 
	{
		fprintf(stderr,"GC %2d-%2d reads %6.0f bins %6d average %0.3f GM %0.3f\n",i*5,i*5+5,counts[i],GCbins[i],(float)counts[i]/GCbins[i],meanRD); 
		counts[i] /= GCbins[i]; if ( counts[i] > maxGC) maxGC = counts[i]; 
	}
	for (i=0;i<blocks;i++) 
	{
		RL[i].GC_corr = counts[(int)(RL[i].GC*BLOCK_SIZE/5)]/maxGC; 
		RL[i].reads *= meanRD/counts[(int)(RL[i].GC*BLOCK_SIZE/5)]; // correction for GC content...
		//fprintf(stderr,"%d %d mappability %f GC %0.3f %d %d %0.3f\n",i*BLOCK_SIZE,i*BLOCK_SIZE+j,RL[i].mappability,RL[i].GC,RL[i].reads,RL[i].adjusted_size,RL[i].GC_corr);	
	}

	return valid;
}

// should be called on reads from single chromosome....
int cluster_reads(struct alignedread** readlist, int s,int e,FRAGMENT* flist,VARIANT* varlist,REFLIST* reflist,REFLIST* genomemask)
{
	int reads=0,i=0,j=0,k=0,first=0,last=0,empty=0,non_empty=0;
	double ws=0,ws_corrected=0;
	double blockscore=0,prior_size=0,mean=0,variance=0,W=0,mean_pb;

	double mean_block_size = 40000; double std_block_size = 8000; double pconst = log(std_block_size) + 0.5* log(2*3.14159); 
	double block_open_penalty = log(0.001);
	double exp_mean = 1.0/500000; double gap_length_penalty = 0;
	double score0,lambda,score_partial;
	
	int blocks = readlist[e-1]->position/BLOCK_SIZE+1;
	struct BLOCK_READ* RL = calloc(sizeof(struct BLOCK_READ),blocks); 
	for (i=0;i<blocks;i++) 
	{
		RL[i].reads=0; RL[i].counts = 0; 
		RL[i].firstread = -1; RL[i].lastread = -1; // init variables
		for (j=0;j<BS1;j++) RL[i].subcounts[j] =0;
	}
	int creads = 0; double empty_blocks=0,nonempty_blocks=0;
	int eb = 0; int block=0;
	double log_sizes[BLOCK_SIZE]; for (i=1;i<BLOCK_SIZE;i++) log_sizes[i] = log(i)-log(BLOCK_SIZE);

	// count # of reads in each block, starting read and last read for each block	
	for (i=s;i<e;i++)
	{
		if (readlist[i]->IS < 0 ||  ((readlist[i]->flag & 1024) ==1024)) continue;  
		if (readlist[i]->mquality < 20) continue;  // don't use these reads

		block = (readlist[i]->position/BLOCK_SIZE); 	
		//if (RL[block].counts < 60) RL[block].poslist[RL[block].counts] = readlist[i]->position; 
		RL[block].counts++; 
		RL[block].lastread = i; RL[block].end = readlist[i]->position; 
		if (RL[block].firstread < 0) { RL[block].firstread = i; RL[block].start = readlist[i]->position; } 
		readlist[i]->blockid = block;
		k = (readlist[i]->position%BLOCK_SIZE)/(BLOCK_SIZE/BS1); RL[block].subcounts[k]++;
		creads++;
	}
	
	for (i=0;i<blocks;i++) 
	{
		RL[i].reads = RL[i].counts;  RL[i].mappability = 1; RL[i].GC = 0.5; RL[i].valid = '1'; RL[i].adjusted_size = BLOCK_SIZE; 
		// default values 
	}
	int validblocks = calculate_block_stats(RL,blocks,reflist,genomemask); // GC and mappability of each window, smoothed read-depth
	//int validblocks = blocks; 


	int* BL = calloc(validblocks,sizeof(int)); j=0;
	for (i=0;i<blocks;i++) 
	{
		if (RL[i].valid == '0') continue; BL[j++] = i; 
	}

	// use bigger block size of 1000bp to calculate background density, if block of 1KB is empty -> due to background....
	// formula = lambda = -1*log(empty_blocks/total_blocks - heavy blocks)
	// some blocks could be empty due to mappability, exclude them or at least weight them appropriately
	for (i=0;i<blocks && i+BF < blocks;i+=BF) 
	{
		non_empty = 0; for (j=0;j<BF;j++) non_empty += RL[i+j].reads; 
		W = 0; for (j=0;j<BF;j++) 
		{ 
			if (RL[i+j].valid == '1') W += RL[i+j].mappability;
		}
		if (non_empty ==0) empty_blocks += W/BF; 
		else if (non_empty < 3*BF) nonempty_blocks += W/BF; // less than 3 reads per 100 bp BLOCK
	}
	double read_density = log((nonempty_blocks+empty_blocks)) - log(empty_blocks);  read_density /= BF;
	//read_density = 0.007; // average number of reads per 100 bp window 
	double logRD = log(read_density);
	fprintf(stderr,"1KB-blocks %0.1f %0.1f %0.8f valid blocks %d:%f\n",nonempty_blocks+empty_blocks,empty_blocks,read_density,validblocks,(float)validblocks/blocks);
	//////// block based analysis gives us probability that a block has 1 or more reads -> background poisson distribution for reads... ////////


	fprintf(stderr,"clustering reads using dynamic programming algorithm reads %d...%d %d %d %f \n",readlist[s]->position,readlist[e-1]->position,e-s+1,creads,read_density);
	fprintf(stdout,"\nclustering reads using dynamic programming algorithm reads %d %d %f \n\n",e-s+1,creads,read_density);

	// initialize all values	
	for (i=0;i<validblocks;i++) 
	{
		RL[BL[i]].previous = -1; RL[BL[i]].bscore = -10000000; RL[BL[i]].sumL = -10000000;

		/*
		//if (RL[BL[i]].reads > 1) 
		{
			fprintf(stdout,"BLOCK %d counts %d %0.1f %0.2f %0.2f %d ",BL[i]*BLOCK_SIZE,RL[BL[i]].counts,RL[BL[i]].reads,RL[BL[i]].GC,RL[BL[i]].mappability,RL[BL[i]].adjusted_size);
			//for (j=1;j<RL[BL[i]].counts && j < 60 ;j++) fprintf(stderr,"%d ",RL[BL[i]].poslist[j]-RL[BL[i]].poslist[j-1]); 
			fprintf(stdout,"\n");
		}
		*/
	}

	// dynamic programming loop for main code 	
	i = 0; score0 = RL[BL[i]].reads*(logRD - log_sizes[(int)RL[BL[i]].adjusted_size])-read_density*(float)RL[BL[i]].adjusted_size/BLOCK_SIZE;  
	RL[BL[1]].bscore = score0; RL[BL[i]].sumL = RL[BL[i]].bscore;

	for (i=1;i<validblocks;i++)
	{
		score0 = RL[BL[i]].reads*(logRD - log_sizes[(int)RL[BL[i]].adjusted_size])-read_density*(float)RL[BL[i]].adjusted_size/BLOCK_SIZE;  
		// score of empty segment... incorrect !! no penalty for non-mappability 
		RL[BL[i]].bscore = RL[BL[i-1]].bscore + score0;   RL[BL[i]].sumL = RL[BL[i]].bscore; 
		// does not matter if previous block was '1' or '0' (no penalty for switching from 1 -> 0), only penalty from 0 -> 1 for opening new block 
		RL[BL[i]].previous = i-1; 
		RL[BL[i]].previous_cluster = RL[BL[i-1]].previous_cluster; RL[BL[i]].score_window = score0;

		j = i-1; reads = RL[BL[i]].reads; eb = 0; ws =0; ws_corrected=0; 
		// mean is the average # of reads per 100 bp window = \lambda
		mean = RL[BL[i]].reads*(float)RL[BL[i]].adjusted_size/BLOCK_SIZE; W = (float)RL[BL[i]].adjusted_size/BLOCK_SIZE;
		variance = RL[BL[i]].reads*RL[BL[i]].reads*(float)RL[BL[i]].adjusted_size/BLOCK_SIZE; 
		//score_partial = RL[BL[i]].reads*log((float)RL[BL[i]].adjusted_size/BLOCK_SIZE);
		score_partial = RL[BL[i]].reads*log_sizes[(int)RL[BL[i]].adjusted_size];

		// consider segment that starts at 'j' and ends at 'i' 
		while ((BL[i]-BL[j]) < (int)(MAX_ISLAND_LENGTH/BLOCK_SIZE) && BL[j] > 0 && eb*BLOCK_SIZE < 10000) // segfault if we set BL[j] >= 0 
		{
			if (RL[BL[j]].reads > 0) eb = 0; else eb++; 
			reads += RL[BL[j]].reads; 
			
			mean += RL[BL[j]].reads*(float)RL[BL[j]].adjusted_size/BLOCK_SIZE; W += (float)RL[BL[j]].adjusted_size/BLOCK_SIZE;
			variance += RL[BL[j]].reads*RL[BL[j]].reads*(float)RL[BL[j]].adjusted_size/BLOCK_SIZE; 
			//score_partial += RL[BL[j]].reads*log((float)RL[BL[j]].adjusted_size/BLOCK_SIZE);
			score_partial += RL[BL[j]].reads*log_sizes[(int)RL[BL[j]].adjusted_size]; 

			ws_corrected += (float)RL[BL[j]].adjusted_size;//*RL[BL[j]].GC_corr;
			//ws_corrected += (float)RL[BL[j]].adjusted_size*RL[BL[j]].GC_corr;
			ws += BLOCK_SIZE;


			if (reads >= 5 && ws_corrected >= 1000 && mean >= 0.1*W)
			{
				//lambda = (double)reads/ws_corrected; 	blockscore = reads*(log(lambda)-1.0); 
				// 10/31/2014: //gap_length_penalty = log(1.0-exp((RL[j].previous_cluster-j-1)*exp_mean*BLOCK_SIZE)); 
				//fprintf(stderr,"penalty %d %f \n",RL[j].previous_cluster-j,gap_length_penalty);
				//if (mean/W > 5) fprintf(stderr,"window stats %d %d reads %d mean %f W %f \n",BL[i],BL[j],reads,mean/W,W*BLOCK_SIZE);
				blockscore = score_partial + reads*log(mean/W) - mean; 
				//mean_pb = mean/((float)BLOCK_SIZE*W); blockscore = reads*log(mean_pb) + (W*BLOCK_SIZE-reads)*log(1.0-mean_pb);

				if (j ==0) score0 = blockscore + block_open_penalty; // +prior_size + gap_length_penalty; 
				else score0 = RL[BL[j-1]].bscore + blockscore + block_open_penalty; // + prior_size + gap_length_penalty;


				if (score0 > RL[BL[i]].bscore)  // extra conditions to consider it a valid segment...
				{
					RL[BL[i]].bscore = score0; RL[BL[i]].previous = (j-1); RL[BL[i]].reads_window = reads; 
					RL[BL[i]].score_window = blockscore; 	RL[BL[i]].previous_cluster = i;
					RL[BL[i]].mean = mean/W; RL[BL[i]].variance = variance/W - RL[BL[i]].mean*RL[BL[i]].mean;
				}
				/*
				if (score0 < RL[BL[i]].sumL) RL[BL[i]].sumL += log(1.0 + exp(score0-RL[BL[i]].sumL)); 
				else 
				{
					//fprintf(stdout,"RL sum %f score0 %f %f ws %f mean_pb %f\n",RL[BL[i]].sumL,score0,blockscore,ws_corrected,mean_pb);
					RL[BL[i]].sumL = score0 + log(1.0 + exp(RL[BL[i]].sumL-score0));
				}*/
			}
			j--;
		}
		j= RL[BL[i]].previous; 
		//if (RL[BL[i]].reads >= 2) fprintf(stdout,"forwardscore1 for %d:%d %d:%d:%d %f %f \n",i,BL[i],RL[BL[i]].start,RL[BL[i]].end,(int)RL[BL[i]].reads,RL[BL[i]].bscore,RL[BL[i]].sumL);
		if (i%500000 ==0) fprintf(stderr,"done for block %d %d\n",i,blocks);
	}

	empty_blocks=0;
       	int reads_background = 0;


	struct FOSMID* fosmidlist = calloc(sizeof(struct FOSMID),100000); int fosmids = 0; 

	// DP backtracking to find optimal segmentation 
	int flag =0; int block0_start = -1; 
	i = validblocks-1; while (i >= 0)
	{
		j = RL[BL[i]].previous;   // both i and j are indexes in smaller list BL[i] and BL[j] are indexes into full list 
		ws = (BL[i]-BL[j])*BLOCK_SIZE;
		if (j < i-1) 
		{
			if (flag*BLOCK_SIZE < 3000) fprintf(stdout,"potential merge... ");
			if (reflist->current >= 0 && flag*BLOCK_SIZE <3000) 
			{
				for (k=0;k<flag;k++) fprintf(stdout,"%d:%0.2f:%0.2f | ",(k+block0_start)*BLOCK_SIZE,RL[block0_start-k].GC,RL[block0_start-k].mappability);
				fprintf(stdout," GC_mappability\n");
			}
			fprintf(stdout,"empty stretch %d %d-%d\n\n",flag*BLOCK_SIZE,(block0_start-flag)*BLOCK_SIZE,block0_start*BLOCK_SIZE);
			fprintf(stdout,"\n==block %d...%d length %0.0f score %0.2f %0.2f 1 | ",BL[j]*BLOCK_SIZE,BL[i]*BLOCK_SIZE,ws,RL[BL[i]].bscore,RL[BL[i]].score_window);
			//fprintf(stdout,"%d %0.2f | ",BLOCK_SIZE*(i-RL[i].previous1),RL[i].bscore1);
			fprintf(stdout,"%d....%d reads %d density %0.2f mean/variance %0.4f:%0.4f %0.1f ==\n",BL[j],BL[i],RL[BL[i]].reads_window,(float)ws/RL[BL[i]].reads_window,RL[BL[i]].mean,RL[BL[i]].variance,RL[BL[i]].mean*RL[BL[i]].mean/(RL[BL[i]].variance-RL[BL[i]].mean));

			// need to add code to add windows ignored due to GC content but adjacent to valid islands 
			// find the first non-empty window 
			first = RL[BL[j+1]].firstread; k=1; while(first < 0) { first = RL[BL[j+1]+k].firstread; k++; } 
			last = RL[BL[i]].lastread; k= 1;  while(last < 0) { last = RL[BL[i]-k].lastread; k++; } 
	
			fosmidlist[fosmids].firstblock = BL[j+1]; fosmidlist[fosmids].lastblock = BL[i]; 
			fosmidlist[fosmids].firstread = first; fosmidlist[fosmids].lastread = last;
			fosmidlist[fosmids].mean  = RL[BL[i]].mean; fosmidlist[fosmids].reads = RL[BL[i]].reads_window; fosmidlist[fosmids].ws = ws;
			fosmids++;

			flag =0;
			block0_start = -1;
		}
		else 
		{
			//fprintf(stderr,"empty block %d %d \n",BL[i],BL[j]); // this corresponds to the blocks not included in segmentation...
			if (block0_start < 0) block0_start = BL[j];
			// empty block... 
			flag++;
			empty_blocks++; reads_background += RL[BL[j]].reads; 
		}
		i = j;
	}
	fprintf(stderr,"empty blocks %0.1f reads %d density %f\n",empty_blocks,reads_background,(float)reads_background/(empty_blocks*BLOCK_SIZE));

	int distance_prev = 1000000; 
	qsort(fosmidlist,fosmids,sizeof(struct FOSMID),cmp_fosmid);
	for (i=0;i<fosmids;i++) 
	{
		if (i > 0) 
		{
			distance_prev = readlist[fosmidlist[i].firstread]->position - readlist[fosmidlist[i-1].lastread]->position; 
			//if (distance_prev > 1000) fprintf(fragment_file,"\n"); 
			//else fprintf(fragment_file,"%d ",distance_prev); 
		}
		generate_single_fragment(readlist,flist,varlist,&fosmidlist[i]);
		print_reads_window(readlist,fosmidlist[i].firstread,fosmidlist[i].lastread,flist,varlist,1);	
	}
	free(fosmidlist);

	free(RL); free(BL);

	return 1;

}


void process_chunk(struct alignedread** readlist,int s,int e,FRAGMENT* flist,VARIANT* varlist,REFLIST* reflist,REFLIST* genomemask )
{
	find_matepair(readlist,s,e);
	fprintf(stderr,"e %d s %d \n",e,s);
 
	int cl = init_clusters(readlist,s,e); // cluster using a maximum intra-cluster distance value

	//estimate_readdistance_distribution(readlist,s,e,cluster); // estimate distances between start positions of adjacent reads within same cluster 
	// cluster size distribution from data | probability of a read being a singleton read
	cluster_reads(readlist,s,e,flist,varlist,reflist,genomemask);

	fprintf(stdout,"\n\n");
	//print_clusters(readlist,s,e,flist,varlist);  // print clusters
}



int init_maskfile(char* maskfile,REFLIST* reflist)
{
        reflist->ns = 0; reflist->names = NULL; reflist->lengths = NULL; reflist->sequences = NULL; reflist->current = -1;
        if (strcmp(maskfile,"None") == 0) return -1;
	if (read_fastaheader(maskfile,reflist) > 0)
	{
		fprintf(stderr,"opening fasta file %s \n",maskfile);
		reflist->fp = fopen(maskfile,"r");
	}
	return 1;
}


// process all reads for a single chromosome since calculation of read-density requires data... 
// readlist is too big,not cleaned up after every chromosome... done 02/11/2014 
// extract haplotype informative reads from sorted bam file // 
int parse_bamfile_fosmid(char* bamfile,HASHTABLE* ht,CHROMVARS* chromvars,VARIANT* varlist,REFLIST* reflist,char* maskfile)
{
	REFLIST* genomemask = (REFLIST*)malloc(sizeof(REFLIST));  
	init_maskfile(maskfile,genomemask); 

	fprintf(stderr,"reading sorted bamfile %s for fosmid pool\n",bamfile);
	int reads=0;
	int MAX_READS = 5000000; // 10 million for now
	struct alignedread* read = (struct alignedread*)malloc(sizeof(struct alignedread));
	struct alignedread** readlist = calloc(MAX_READS,sizeof(struct alignedread*));
	for (reads=0;reads<MAX_READS;reads++) readlist[reads] = calloc(1,sizeof(struct alignedread)); 
	struct alignedread* read_pt;

	FRAGMENT fragment; fragment.variants =0; fragment.alist = (allele*)malloc(sizeof(allele)*10000);
	FRAGMENT* flist = (FRAGMENT*)malloc(sizeof(FRAGMENT)*MAX_READS/5); int fragments =0;
	
	int chrom=0; int r=0,i=0;
	int prevchrom=-1; int prevtid = -1; int prevposition = -1; // position of previous read in sorted bam file
	int lastread = 0;

	samfile_t *fp; 
	if ((fp = samopen(bamfile, "rb", 0)) == 0) { fprintf(stderr, "Fail to open BAM file %s\n", bamfile); return -1; }
	bam1_t *b = bam_init1();

	while (samread(fp, b) >= 0)
	{
		//readlist[r] = calloc(1,sizeof(struct alignedread));
		fetch_func(b, fp,readlist[r]); 
		if ((readlist[r]->flag & (BAM_FUNMAP|BAM_FSECONDARY|BAM_FQCFAIL|BAM_FDUP)))  // unmapped reads, PCR/optical dups are ignored
		{
			free_readmemory(readlist[r]); continue;
		}
		else if (readlist[r]->mquality < 20)  // ignore poorly mapped reads 
		{
			//free_readmemory(readlist[r]); continue;
		}
		// find the chromosome in reflist that matches read->chrom if the previous chromosome is different from current chromosome
		// if too many reads, break off when distance between adjacent reads is large 
		if (readlist[r]->tid != prevtid || r >= MAX_READS-1) 
		{

			if (r >= MAX_READS-1) fprintf(stderr,"limit on max reads %d exceeded.. need to clean up buffer \n",MAX_READS);
			if (prevtid >=0) 
			{
				fprintf(stderr,"reads in buffer %d process reads to identify islands for chromosome\n",r);
				process_chunk(readlist,lastread,r,flist,varlist,reflist,genomemask); 
				// free up reads in list and move up index 
				for (i=lastread;i<r;i++) free_readmemory(readlist[i]); 
				read_pt = readlist[0]; readlist[0] = readlist[r]; readlist[r] = read_pt; r = 0;
				for (i=0;i<fragments;i++) 
				{
					free(flist[i].alist); free(flist[i].id); 
				}
				fprintf(stderr,"free memory for reads from chrom %d cleaning up of fragment list %d\n",prevtid,fragments);
				fragments =0; 
				if (reflist->ns > 0 && prevtid < reflist->ns) free(reflist->sequences[prevtid]); 
				if (genomemask->ns > 0 && prevtid < reflist->ns) free(genomemask->sequences[prevtid]); 
				//return 1;  // need to remove this for all chromosomes
			}
			if (reflist->ns > 0 && readlist[r]->tid < reflist->ns) 
			{
				read_chromosome(reflist,readlist[r]->tid,reflist->fp); // read the next chromosome
				reflist->current = readlist[r]->tid; // assign current variable to tid
			}
			if (genomemask->ns > 0 && readlist[r]->tid < genomemask->ns)  // segfault here if the two genomes are not identical 
			{
				read_chromosome_mask(genomemask,readlist[r]->tid,genomemask->fp); // read the next chromosome
				genomemask->current = readlist[r]->tid; // assign current variable to tid
			}

			chrom = getindex(ht,readlist[r]->chrom); 
			get_chrom_name(readlist[r],ht,reflist); // reflist->current has chromosome index, -1 if no reflist 
			lastread = r;
		}
		else 	chrom = prevchrom;

		fragment.variants =0; //fragment.id = readlist[r]->readid;
                if (chrom >=0)  extract_variants_read(readlist[r],ht,chromvars,varlist,1,&fragment,chrom,reflist);
		if (fragment.variants > 0) 
		{ 
			add_fragment(flist,&fragment,readlist[r],fragments); readlist[r]->findex = fragments++; 
		}
		else readlist[r]->findex = -1;

		reads+=1; if (reads%2000000 ==0) fprintf(stderr,"processed %d reads, useful fragments \n",reads);
		prevchrom = chrom; prevtid = readlist[r]->tid; 
		if (readlist[r]->IS == 0) prevposition = readlist[r]->position; 
		else if (readlist[r]->IS > 0) prevposition = readlist[r]->position + readlist[r]->IS; // outer end of fragment r1....r2
		r++;
	}

	process_chunk(readlist,lastread,r,flist,varlist,reflist,genomemask); 
	for (i=lastread;i<r;i++) free_readmemory(readlist[i]); 
	for (reads=0;reads<MAX_READS;reads++) free(readlist[reads]); free(readlist);
	bam_destroy1(b);
	if (reflist->ns > 0 && prevtid < reflist->ns) { free(reflist->sequences[prevtid]); fclose(reflist->fp); } 
	if (genomemask->ns > 0 && prevtid < genomemask->ns) { free(genomemask->sequences[prevtid]); fclose(genomemask->fp); } 
	return 1;
}

