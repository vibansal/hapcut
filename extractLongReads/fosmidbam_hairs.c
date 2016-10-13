
#include<stdint.h>

#include "fragments.h"

#include "print_clusters.c"
#include "GC_repeat_analysis.c"
#include "analyze_clusters.c"
//#include "cluster_simple.c"
#include "sc-code.c"


// dynamic programming backtracking to find best segmentation 
int find_best_segmentation(struct BLOCK_READ* RL,int* BL,int validblocks,struct FOSMID* fosmidlist,int refflag)
{
	int i=validblocks-1,j=0,k=0,c=0; int fosmids=0; int flag =0; int block0_start = -1;  int first=0,last=0,l=0,t=0;
	double ws=0,mean=0;

	while (i >= 0) // start from last block in the list 
	{
		j = RL[BL[i]].previous[0];   mean = RL[BL[i]].mean;  
		// both i and j are indexes in smaller list BL[i] and BL[j] are indexes into full list 

		if (j < i-1) // (j-i) is a haplotype fragment
		{
			//if (VERBOSE ==2) print_block(RL,BL,i,j,flag,refflag,block0_start); 
			// need to add code to add windows ignored due to GC content but adjacent to valid islands 
			// find the first non-empty window 
			first = RL[BL[j+1]].firstread; k=1; while(first < 0) { first = RL[BL[j+1]+k].firstread; k++; } 
			last = RL[BL[i]].lastread; k= 1;  while(last < 0) { last = RL[BL[i]-k].lastread; k++; } 
	
			fosmidlist[fosmids].firstblock = BL[j+1]; fosmidlist[fosmids].lastblock = BL[i]; 
			fosmidlist[fosmids].start = BL[j+1]*BLOCK_SIZE; fosmidlist[fosmids].end = (BL[i]+1)*BLOCK_SIZE; 
			fosmidlist[fosmids].firstread = first; fosmidlist[fosmids].lastread = last;
			fosmidlist[fosmids].mean  = mean; fosmidlist[fosmids].reads_window = RL[BL[i]].reads_window; fosmidlist[fosmids].ws =  (BL[i]-BL[j])*BLOCK_SIZE;
			fosmidlist[fosmids].score_window = RL[BL[i]].score_window; 
			fosmids++;

			flag =0; block0_start = -1;
		}
		else // empty block background
		{
			if (block0_start < 0) block0_start = BL[j]; // empty block... 
			flag++; //empty_blocks++; reads_background += RL[BL[j]].reads; 
		}
		i = j;
	}
	return fosmids;
}

double calculate_bg_density(struct FOSMID* fosmidlist,int fosmids,struct BLOCK_READ* RL,int blocks)
{
        int i=0,k=0;
        double bg_reads=0,bg_bins=0; // read counts outside boundaries of fragments
        int genome_covered = 0; // genome spanned by fragments in basepairs 
        for (i=0;i<fosmids;i++)
        {
                genome_covered += fosmidlist[i].lastblock-fosmidlist[i].firstblock;
                if (i==fosmids-1) continue;
                //fprintf(stderr,"blocks %d .... %d \n",fosmidlist[i].lastblock,fosmidlist[i+1].firstblock);
                for (k=fosmidlist[i].lastblock+1;k<fosmidlist[i+1].firstblock;k++)
                {
                        if (RL[k].valid != '0') { bg_reads += RL[k].reads; bg_bins += 1; } // this is using reads while sc uses peaks 
                }
        }
        fprintf(stdout,"detected %d fosmid segments bgreads %0.0f %0.0f rd(perblock):%0.5f genome-cov %d frac %0.4f\n",fosmids,bg_reads,bg_bins,bg_reads/(bg_bins),genome_covered,(double)genome_covered/blocks);
        fprintf(stderr,"detected %d fosmid segments bgreads %0.0f %0.0f rd:%0.5f genome-cov %d frac %0.4f\n",fosmids,bg_reads,bg_bins,bg_reads/(bg_bins),genome_covered,(double)genome_covered/blocks);
        return bg_reads/bg_bins;

}


int print_fosmid_list(struct FOSMID* fosmidlist,int fosmids,struct alignedread** readlist,FRAGMENT* flist,VARIANT* varlist,struct BLOCK_READ* RL,int* BL,double frag_penalty)
{
	int i=0,k=0,l=0,j=0; double ws=0;
	int fragments_in_window= 0; int total_length =0; int distance_next=0;
	int clusters = 0; int curr_cluster = 0,prev=0;
	int cstart=-1,cend =-1; 
	double log_sizes[BLOCK_SIZE]; for (i=1;i<BLOCK_SIZE;i++) log_sizes[i] = log(i)-log(BLOCK_SIZE);
	double bg_ll,seg_ll,mean[3]; int vblocks =0;

	for (i=0;i<fosmids;i++) 
	{
		if (cstart <0) { cstart = fosmidlist[i].firstblock; cend = fosmidlist[i].lastblock; } 

		fosmidlist[i].cluster = curr_cluster; 
		k = fosmidlist[i].lastblock; 	ws = (k+1)*BLOCK_SIZE-fosmidlist[i].firstblock*BLOCK_SIZE;
		if (FOSMIDS != 2)
		{
			bg_ll = background_likelihood(RL,BL,fosmidlist[i].firstblock,fosmidlist[i].lastblock,log_sizes,read_density_global);
			//bg_ll =0; for (j=fosmidlist[i].firstblock;j<fosmidlist[i].lastblock;j++) { if (RL[j].valid =='1') bg_ll += RL[j].bgscore; }  
			seg_ll = segment_likelihood(RL,BL,fosmidlist[i].firstblock,fosmidlist[i].lastblock,log_sizes,mean);
			if (seg_ll+frag_penalty-bg_ll < 5) fprintf(stdout," small-delta "); fprintf(stdout,"comparison seg:%0.2f bg:%0.2f delta %0.2f \n",seg_ll,bg_ll,seg_ll+frag_penalty-bg_ll);
		}

		vblocks=0; for (j=fosmidlist[i].firstblock;j<fosmidlist[i].lastblock;j++) {  if (RL[j].valid == '1') vblocks++; } 
		print_block(RL,fosmidlist[i].firstblock,fosmidlist[i].lastblock,2);

		fprintf(stdout,"%d ==block %d...%d length %0.0f window:%0.2f 1 | ",fosmidlist[i].cluster,fosmidlist[i].firstblock*BLOCK_SIZE,k*BLOCK_SIZE,ws,fosmidlist[i].score_window);
		fprintf(stdout,"%d....%d reads %0.0f density(reads/bp) %0.3f mean(/kb) %0.3f val-blocks %0.3f\n",fosmidlist[i].firstblock,k,fosmidlist[i].reads_window,(float)fosmidlist[i].ws/fosmidlist[i].reads_window,fosmidlist[i].mean*1000/BLOCK_SIZE,(float)vblocks/(fosmidlist[i].lastblock-fosmidlist[i].firstblock));

		print_reads_window(readlist,fosmidlist[i].firstread,fosmidlist[i].lastread,flist,varlist,1);	
		generate_single_fragment(readlist,flist,varlist,&fosmidlist[i]);
		if (i < fosmids-1) distance_next = (fosmidlist[i+1].firstblock-fosmidlist[i].lastblock)*BLOCK_SIZE; else distance_next = 0;
		fosmidlist[i].distance_next = distance_next; 
		fragments_in_window++; total_length += ws;

		if (i < fosmids-1) print_block(RL,fosmidlist[i].lastblock+1,fosmidlist[i+1].firstblock-1,1);

		if (distance_next < MIN_SEPARATION) 
		{
			fprintf(stdout,"....MERGE.... %d\n\n",distance_next); total_length += distance_next; cend = fosmidlist[i+1].lastblock;  
		}
		else // analyze cluster and print output
		{
			if (fragments_in_window >=2 && fragments_in_window <=4) 
			{
				//fprintf(stdout,"calling function .... \n"); evaluate_triple_segments(RL,BL,fosmidlist[i-2].firstblock,fosmidlist[i-2].lastblock,fosmidlist[i-1].firstblock,fosmidlist[i-1].lastblock,fosmidlist[i].firstblock,fosmidlist[i].lastblock,log_sizes);
				//find_2seg(RL,BL,cstart,cend,log_sizes,frag_penalty);
			}
			if (fragments_in_window >=2 && fragments_in_window <= 4 && FIRSTPASS ==0)
			{
				find_3seg(RL,BL,cstart,cend,log_sizes,frag_penalty,fragments_in_window);
			}
			k = cend;
			prev = RL[k].previous[4]; 
			while (1)
			{
				l = (k-BL[prev+1]+1);
				if (l >=5) fprintf(stdout,"prev-4 %d-%d length:%d mean:%0.4f\n",BL[prev+1]*BLOCK_SIZE,k*BLOCK_SIZE,l*BLOCK_SIZE,RL[k].pmean[4]*1000/BLOCK_SIZE);
				if (prev <= 0 || BL[prev+1] <= cstart) break; 
				k = BL[prev]; prev = RL[k].previous[4];
			}

			fprintf(stdout,"\n-------- frags:%d ------ W:%d --------- %d ------------%d...%d\n\n\n\n\n",fragments_in_window,total_length,distance_next,cstart,cend);
			fragments_in_window = 0; total_length = 0; curr_cluster = i+1; clusters++; 
			cstart = fosmidlist[i+1].firstblock; cend = fosmidlist[i+1].lastblock;
		}
	}
	return clusters;
}

int dynamic_programming_loop(struct BLOCK_READ* RL,int blocks,int* BL,int validblocks,double read_density,double block_open_penalty,double length_penalty)
{
        int i = 0,j=0,eb=0,reads=0;
        double ws=0,ws_corrected=0,flag=0;
        double blockscore=0,prior_size=0,mean=0,variance=0,W=0,mean_pb;
        double score0,lambda,score_partial,sum_partial,score1,frag_density;
        double log_sizes[BLOCK_SIZE]; for (i=1;i<BLOCK_SIZE;i++) log_sizes[i] = log(i)-log(BLOCK_SIZE);
        double logRD = log(read_density);

        double block_open_penalty_1 = (block_open_penalty+length_penalty);
        double block_open_penalty_2 = (block_open_penalty+length_penalty)/3;
        double fragment_p = block_open_penalty; double fragment_p_1 = log(1.0-exp(fragment_p)); 

	double estimates[3]; estimate_parameters_sc(RL,blocks,estimates); fprintf(stderr,"estimating init parameters\n");

        score0 = RL[BL[i]].reads*(logRD + log_sizes[(int)RL[BL[i]].adjusted_size]-log(RL[BL[i]].GC_corr))-read_density*(float)RL[BL[i]].adjusted_size/(BLOCK_SIZE*RL[BL[i]].GC_corr);


	// shouldn't this be 0 instead of 1 ?? BUG 09/09/16 
        RL[BL[0]].bscore[0] = score0 + fragment_p_1;  RL[BL[0]].bscore[4] =  RL[BL[1]].bscore[0];
	RL[BL[i]].bgscore = score0; 
        RL[BL[i]].score_partial = RL[BL[i]].reads*(log_sizes[(int)RL[BL[i]].adjusted_size] -log(RL[BL[i]].GC_corr));
        RL[BL[i]].W = (float)RL[BL[i]].adjusted_size/(RL[BL[i]].GC_corr*BLOCK_SIZE);

        for (i=1;i<validblocks;i++)
        {
                score0 = RL[BL[i]].reads*(logRD + log_sizes[(int)RL[BL[i]].adjusted_size]-log(RL[BL[i]].GC_corr))-read_density*(float)RL[BL[i]].adjusted_size/(BLOCK_SIZE*RL[BL[i]].GC_corr);
		RL[BL[i]].bgscore = score0;
                // score of empty segment... incorrect !! no penalty for non-mappability 

                RL[BL[i]].bscore[0] = RL[BL[i-1]].bscore[0] + score0 + fragment_p_1; // fragment_p_1 is ln(1-p) where 'p' is probability of selecting bin as fragment start 
                RL[BL[i]].bscore[4] = RL[BL[i]].bscore[0];
                // does not matter if previous block was '1' or '0' (no penalty for switching from 1 -> 0), only penalty from 0 -> 1 for opening new block 
                RL[BL[i]].previous[0] = i-1; RL[BL[i]].previous_cluster = RL[BL[i-1]].previous_cluster; RL[BL[i]].score_window = score0;
                RL[BL[i]].previous[4] = i-1;

                // mean is the average # of reads per 100 bp window = \lambda
                RL[BL[i]].score_partial = RL[BL[i]].reads*(log_sizes[(int)RL[BL[i]].adjusted_size] -log(RL[BL[i]].GC_corr));
                RL[BL[i]].W = (float)RL[BL[i]].adjusted_size/(RL[BL[i]].GC_corr*BLOCK_SIZE);

                variance = RL[BL[i]].reads*RL[BL[i]].reads*(float)RL[BL[i]].adjusted_size/BLOCK_SIZE;
                score_partial = RL[BL[i]].score_partial; mean = RL[BL[i]].reads; W = RL[BL[i]].W; reads = RL[BL[i]].reads;
		//fprintf(stdout,"block:%d %d %f %f \n",i,BL[i],mean,score_partial);

                j = i-1; eb = 0; ws =0; ws_corrected=0;
                // consider segments that start at 'j' and ends at 'i'  (j decreasing) 
		//while ((BL[i]-BL[j]) < (int)(MAX_ISLAND_LENGTH/BLOCK_SIZE) && BL[j] > 0 && eb*BLOCK_SIZE < 10000) // segfault if we set BL[j] >= 0 
                while ((BL[i]-BL[j]) < (int)(MAX_FRAG_LENGTH/BLOCK_SIZE) && BL[j] > 0 && eb*BLOCK_SIZE < max_empty_stretch) // segfault if we set BL[j] >= 0 
                {
                        if (RL[BL[j]].reads > 0) eb = 0; else eb++;  reads += RL[BL[j]].reads;

                        mean += RL[BL[j]].reads; W += RL[BL[j]].W; score_partial += RL[BL[j]].score_partial;
                        variance += RL[BL[j]].reads*RL[BL[j]].reads*(float)RL[BL[j]].adjusted_size/BLOCK_SIZE;

                        ws_corrected += (float)RL[BL[j]].adjusted_size; // *RL[BL[j]].GC_corr;
                        ws += BLOCK_SIZE;

                        // default is 10 reads, 1 kb and 0.1= 1 read per 1kb | threshold should be based on background read 2-5 fold density...
                        frag_density = (double)ws/reads;
                        if (reads >= MIN_READS_PER_FRAGMENT && ws_corrected >= MIN_FRAG_LENGTH && mean*BLOCK_SIZE/W >= MIN_READ_DENSITY) // extra conditions to consider it a valid segment...
                        {
                                blockscore = score_partial + reads*log(mean/W) - mean;  // score of the full block 

                                if (j ==0) score0 = blockscore + block_open_penalty_1; // +prior_size + gap_length_penalty; 
                                else score0 = RL[BL[j-1]].bscore[0] + blockscore + block_open_penalty_1; // + prior_size + gap_length_penalty;
                                if (score0 > RL[BL[i]].bscore[0])
                                {
                                        RL[BL[i]].bscore[0] = score0; RL[BL[i]].previous[0] = (j-1); RL[BL[i]].reads_window = reads;
                                        RL[BL[i]].score_window = blockscore;    RL[BL[i]].previous_cluster = i;
                                        RL[BL[i]].mean = mean/W; RL[BL[i]].variance = variance/W - RL[BL[i]].mean*RL[BL[i]].mean;
                                }
                                RL[BL[j]].curr_score = RL[BL[j-1]].bscore[0] + blockscore;

                                if (j ==0) score0 = blockscore + block_open_penalty_2;
                                else score0 = RL[BL[j-1]].bscore[4] + blockscore + block_open_penalty_2; // + prior_size + gap_length_penalty;
                                if (score0 > RL[BL[i]].bscore[4])
                                {
                                        RL[BL[i]].bscore[4] = score0; RL[BL[i]].previous[4] = (j-1); RL[BL[i]].pmean[4] = mean/W;
                                }
                        }
                        j--;
                }
                if (i%500000 ==0) fprintf(stderr,"done for block %d %d\n",i,blocks);
        }
}

// should be called on reads from single chromosome....
int process_reads_chrom(struct alignedread** readlist,int s,int e,FRAGMENT* flist,VARIANT* varlist,REFLIST* reflist,REFLIST* genomemask )
{
	find_matepair(readlist,s,e,2000); fprintf(stderr,"start s:%d end read e:%d \n",s,e);	
	//int cl = init_clusters(readlist,s,e); // cluster using a maximum intra-cluster distance value
	//estimate_readdistance_distribution(readlist,s,e,cluster); // estimate distances between start positions of adjacent reads within same cluster 

	distance_nextread_dist(readlist,s, e); 

	int reads=0,i=0,j=0,k=0,empty=0,non_empty=0,eb=0,t=0;
	double blockscore=0,prior_size=0,mean=0,variance=0,W=0,mean_pb;

	int blocks = readlist[e-1]->position/BLOCK_SIZE+1; if (blocks < 100) return -1; // small number of blocks 
	struct BLOCK_READ* RL = calloc((blocks+1),sizeof(struct BLOCK_READ)); 

	init_blocklist(RL,blocks,readlist,s,e);

	int validblocks = calculate_block_stats(RL,blocks,reflist,genomemask); // calculate GC and mappability of each window, smoothed read-depth
	if (validblocks < 100) { free(RL); return -1; } // small number of blocks 
	int* BL = calloc(validblocks,sizeof(int)); j=0;
	for (i=0;i<blocks;i++) 
	{
		if (RL[i].valid == '0') continue; BL[j] = i; j++; 
	}

	//read_density = 0.011387; 
	double read_density = calculate_background_RD(RL,blocks,validblocks); double logRD = log(read_density); // per block
	if (read_density_global < 0) read_density_global = read_density;
	//////// block based analysis gives us probability that a block has 1 or more reads -> background poisson distribution for reads... ////////
	
	//double mean_block_size = 40000; double std_block_size = 8000; double pconst = log(std_block_size) + 0.5* log(2*3.14159); 
	// -7.72 = ln(0.000435) ?? penalty should be equal to #fragments/#blocks
	// frags = 1000, bins = 250,000,0, penalty = 0.0004
	// need better way to initialize block penalty 	
	double block_open_penalty = -7.72, length_penalty = -4.61; // -4.61 is penalty for segment length = 1/100 uniform distribution over 1-100 kb 
	//block_open_penalty = -6, length_penalty = -4; // -4.61 is penalty for segment length = 1/100 uniform distribution over 1-100 kb 
	block_open_penalty = block_penalty_global; length_penalty = length_penalty_global;


	fprintf(stderr,"clustering reads using dynamic programming algorithm reads %d...%d %d bg-read-density (per block) %f \n",readlist[s]->position,readlist[e-1]->position,e-s+1,read_density_global);
	fprintf(stdout,"\nclustering reads using dynamic programming algorithm reads %d bg-read-density (per block) %f \n\n",e-s+1,read_density_global);
	
	int refflag = 0;	if (reflist->current >=0) refflag = 1;
	struct FOSMID* fosmidlist = calloc(100000,sizeof(struct FOSMID)); 

	// dynamic programming loop for main code 	
	if (FOSMIDS ==2) dynamic_programming_sc(RL, blocks, BL,validblocks,read_density_global,block_open_penalty);
	else dynamic_programming_loop(RL,blocks,BL,validblocks,read_density_global,block_open_penalty,length_penalty);

	// call function to DP backtracking to find optimal segmentation 
	int fosmids =  find_best_segmentation(RL,BL,validblocks,fosmidlist,refflag);
	qsort(fosmidlist,fosmids,sizeof(struct FOSMID),cmp_fosmid);

	double frag_opening_penalty = log(fosmids+1)-log(validblocks); fprintf(stderr,"fragment penalty %0.2f \n",frag_opening_penalty);
	double new_read_density = calculate_bg_density(fosmidlist,fosmids, RL,blocks); 

	// print list of fragments identified, print fragments that can be merged and sum of lengths... print # of fragments with >=1 variants, others are useless
	int clusters = print_fosmid_list(fosmidlist, fosmids, readlist, flist,varlist,RL,BL,frag_opening_penalty);
	fprintf(stdout,"statistics for chrom %s fragments %d clusters %d blocks %d validblocks %d\n",readlist[s]->chrom,fosmids,clusters,blocks,validblocks);
	fprintf(stderr,"statistics for chrom %s fragments %d clusters %d blocks %d validblocks %d\n",readlist[s]->chrom,fosmids,clusters,blocks,validblocks);

	//calculate_GC_stats(RL,blocks,fosmidlist,fosmids); // this is not used to adjust read-depth but can be used in 2nd-pass
	free(fosmidlist); 
	free(RL); free(BL);

	return 1;
}
