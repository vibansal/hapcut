
// code implemented for analyzing single cell data but also useful in general 

int cmpfunc (const void * a, const void * b) {   return ( *(int*)a - *(int*)b );}

int USE_SL=0; // use sub-fragment lengths in fragment likelihood 

// function to estimate # of fragments, intra-fragment gap and background read density 
int estimate_parameters_sc1(struct BLOCK_READ* RL,int blocks,double estimates[])
{
	double seglen[2]={0,0}, gaplen[2] = {0,0}; int s=0,ps=0,empty=0,pempty=0;
	int j=0;
	// use peak variable to store block information on segment length/ gap length 
	for (j=0;j<blocks;j++) 
	{ 
		RL[j].peak = 0; 
		if (RL[j].valid == '0') continue; 

		if (RL[j].counts==0) 
		{ 
			if (s > 0) ps = s;  s=0; empty++; 
		}
		if (RL[j].counts > 0) 
		{
			if (empty > 0) 
			{
				if (ps > 0) fprintf(stdout,"seg-gap pair j:%d %d %d \n",j*BLOCK_SIZE,ps*BLOCK_SIZE,empty*BLOCK_SIZE); 
				pempty = empty; 
			}
			s+=1; empty = 0; 
		}
	}
}

// function to estimate # of fragments, intra-fragment gap and background read density 
int estimate_parameters_sc(struct BLOCK_READ* RL,int blocks,double estimates[])
{
	int i=0,j=0,k=0;
	int prev = 0, peaks = 0, sign = 0; // number of changes in coverage between adjacent bins, only positive, if two changes very close, disregard...
	// we can precalculate the locations of MDA clusters or peaks

	for (j=0;j<blocks;j++) 
	{ 
		RL[j].peak = 0; 
		if (RL[j].valid == '0') continue; 

		if (RL[j].counts > prev && sign ==0) { sign = 1; } 
		else if (RL[j].counts < prev && sign ==1) { sign = 0; RL[k].peak = 1; peaks++; } 
		else if (RL[j].counts < prev) sign = 0; 

		if (USE_SL ==1) 
		{
			// count segments instead of peaks 
			RL[j].peak = 0; 
			if (prev ==0 && RL[j].counts > 0) 
			{
				i = j+1; RL[j].peak =1; 
				while (i < blocks && (RL[i].valid == '0' || RL[i].counts > 0)) 
				{
					if (RL[i].valid == '1') RL[j].peak +=1; 		
					i++; 
				}
			}
			if (RL[j].peak >=1) peaks++;
		}

		prev = RL[j].counts; k = j; 
	}
	estimates[0] = 3000; estimates[1] = 75000;
	int* peakdist = calloc(peaks,sizeof(int)); peaks =0;
	prev = -1;  int prevpd = 0;
	for (j=0;j<blocks;j++) 
	{
		if (RL[j].peak >=1) 
		{ 
			if (prev >= 0) peakdist[peaks++] = j-prev; 
		//	fprintf(stdout,"peakpairs %d %d \n",prevpd*BLOCK_SIZE,(j-prev)*BLOCK_SIZE); 
			prevpd = j-prev;  
			prev = j; 
		} 
	}	
	qsort(peakdist,peaks,sizeof(int),cmpfunc); 

	// exponential distribution has tail prob = e^(-lambda x), use values of 75k and 150k to calculate ratio 
	double tail_counts[5] = {0,0,0,0,0}; double tail_points[5] = {50000,75000,100000,150000,200000}; // 50K,75K, 100K, 150K
	for (j=peaks-1;j>=0;j--) 
	{
		if (peakdist[j]*BLOCK_SIZE > tail_points[0]) tail_counts[0]+=1;
		if (peakdist[j]*BLOCK_SIZE >tail_points[1]) tail_counts[1]+=1;
		if (peakdist[j]*BLOCK_SIZE > tail_points[2]) tail_counts[2]+=1;
		if (peakdist[j]*BLOCK_SIZE >tail_points[3]) tail_counts[3]+=1;
		if (peakdist[j]*BLOCK_SIZE >tail_points[4]) tail_counts[4]+=1;
	}
	double exp_mean1  = (tail_points[2]-tail_points[0])/log(tail_counts[0]/tail_counts[2]); 
	double exp_mean2  = (tail_points[3]-tail_points[1])/log(tail_counts[1]/tail_counts[3]); 
	double exp_mean3  = (tail_points[4]-tail_points[2])/log(tail_counts[2]/tail_counts[4]); 
	
	double intragap_mean = peakdist[(peaks-(int)tail_counts[1])/2]*BLOCK_SIZE/log(2); // ln(2)/lambda 
	double intragap_mean1 = peakdist[(peaks)/2]*BLOCK_SIZE/log(2); // ln(2)/lambda 

	for (i=0;i<4;i++) fprintf(stderr,"%0.0f:%0.0f ",tail_points[i],tail_counts[i]);
	fprintf(stderr,"tail counts | mean of intra-gap %0.2f %0.2f %d mean_exp %f %f %f\n",intragap_mean,intragap_mean1,peaks,exp_mean1,exp_mean2,exp_mean3);
	for (i=0;i<4;i++) fprintf(stdout,"%0.0f:%0.0f ",tail_points[i],tail_counts[i]);
	fprintf(stdout,"tail counts mean1 %0.3f %d mean_exp %f %f %f\n",intragap_mean,peaks,exp_mean1,exp_mean2,exp_mean3);
	free(peakdist);

	estimates[0] = intragap_mean; estimates[1]= exp_mean1;  
	//if (tail_counts[3] >= 50) estimates[1] = exp_mean2; else estimates[1] = exp_mean1;
}


// mean = 200/ 3600 = 1/18 -> beta distribution... 
// B(x,y) = (x-1)! (y-1)! / (x+y-1)! 
// likelihood (k peaks | n, a,b) = C(n,k) B(k+a,n-k+b) / B(a,b) = (k+a-1)! x (n-k+b-1)! / (n+a+b-1)!  X term_computed / constant 
// mean = na/(a+b)
// problem is due to N(c,r) factor that has been added.... for background this factor is not there....  !!  
// use peak variable for all calculations..

double fact(int n) { int i=0; double ncr = 0.0; for (i=1;i<=n;i++) ncr += log(i); return ncr; } 

int dynamic_programming_sc(struct BLOCK_READ* RL,int blocks,int* BL,int validblocks,double read_density,double block_open_penalty)
{
        int i = 0,j=0,eb=0;
        double ws=0,ws_corrected=0,flag=0,reads=0;
        double blockscore=0,prior_size=0,mean=0,variance=0,W=0,mean_pb;
        double score0,lambda,score_partial,sum_partial,score1,frag_density;
        double log_sizes[BLOCK_SIZE]; for (i=1;i<BLOCK_SIZE;i++) log_sizes[i] = log(i)-log(BLOCK_SIZE);

	double estimates[3]; estimate_parameters_sc(RL,blocks,estimates);
	double intra_frag_gap = estimates[0], inter_frag_gap = estimates[1]; 

	//read_density = 2*(float)BLOCK_SIZE/inter_frag_gap;  // average frequency of non-empty blocks outside fragments 
	double logRD = log(read_density), logRD1 = log(1.0-read_density); fprintf(stderr,"bg read density %f %f \n",read_density,logRD);

	double ld_f = (double)BLOCK_SIZE/intra_frag_gap; // expected frequency of MDA clusters per BLOCK 
	double log_ld = log(ld_f), log_ld1 = log(1.0-ld_f); 

	MAX_FRAG_LENGTH = 1000000; MIN_READS_PER_FRAGMENT=4; MIN_FRAG_LENGTH=2000;  
	block_open_penalty = -8; 

        double block_open_penalty_1 = block_open_penalty;
        double fragment_p = block_open_penalty_1; double fragment_p_1 = log(1.0-exp(fragment_p)); 

	int a_b = 10,b_b = (int)((float)a_b*intra_frag_gap/(float)BLOCK_SIZE)-a_b; // parameters of beta binomial dist ratio = 1/15 
	double beta_const = fact(a_b+b_b-1)-fact(a_b-1)-fact(b_b-1); 
	fprintf(stderr,"inside DP loop for single cell MDA data, penalties Long-gaps %f perkb %f small-gaps %f %f %f beta priors: %d,%d\n",logRD,logRD+log(1000.0/BLOCK_SIZE),log_ld,intra_frag_gap,inter_frag_gap,a_b,b_b);
	//for (i=0;i<validblocks;i++) RL[BL[i]].peak = RL[BL[i]].reads; 

	// assume that mean segment length for background fragment is '1' -> likelihood term is -1.seglength + log(lambda) | lambda= 1 
	i = 0; 
	if (RL[BL[i]].peak >=1) score0 = logRD; else score0 = logRD1; if (USE_SL ==1) score0 -= RL[BL[i]].peak;
        RL[BL[i]].bscore[0] = score0;  // score if this block is outside fragments
        RL[BL[i]].W = 1;

	// reads is # of segments, lambda is sum of segment lengths, used to calculate mean segment length 
	// prior on # of segments per fragment and mean segment length is needed to avoid breaking
	// filter on # of segments is incorrect... 

        for (i=1;i<validblocks;i++)
        {
		// likelihood (k peaks | n, a,b) = n!/(n-k)!k!  X (k+a-1)! x (n-k+b-1)! / (n+a+b-1)!  X term_computed / constant 
		// n = 1, k=0/1
		//if (RL[BL[i]].peak ==1) RL[BL[i]].score_partial = log_ld; else RL[BL[i]].score_partial = log_ld1; score_partial = RL[BL[i]].score_partial; 

		if (RL[BL[i]].peak >=1) score0 = logRD; else score0 = logRD1; if (USE_SL ==1) score0 -= RL[BL[i]].peak;
                RL[BL[i]].bscore[0] = RL[BL[i-1]].bscore[0] + score0;
                RL[BL[i]].previous[0] = i-1; RL[BL[i]].score_window = score0; RL[BL[i]].mean = -1; 
               

                RL[BL[i]].W = 1; W = RL[BL[i]].W; 
		if (RL[BL[i]].peak >=1) { score_partial = fact(a_b) + fact(b_b-1) - fact(a_b+b_b); reads = 1; } 
		else { score_partial = fact(a_b-1) + fact(b_b) - fact(a_b+b_b); reads = 0;} 
		mean = reads; lambda = RL[BL[i]].peak;

                if (i%200000 ==0) fprintf(stderr,"done for block %d %d\n",i,blocks);
		// last+1 block of segment should have coverage 0 so that we don't break sub-fragments in middle....
		if (i+1<validblocks && RL[BL[i+1]].counts > 0) continue; // extra condition, seems to work

		// need to check for big low mappability gaps ...
                j = i-1; eb = 0; ws_corrected=0;
                while ((BL[i]-BL[j]) < (int)(MAX_FRAG_LENGTH/BLOCK_SIZE) && BL[j] > 0 && eb*BLOCK_SIZE < max_empty_stretch) // segfault if we set BL[j] >= 0 
                {
			if ( (BL[j+1]-BL[j])*BLOCK_SIZE > 20000) break; 
                        if (RL[BL[j]].counts > 0) eb = 0; else eb++;  
			W += RL[BL[j]].W; 

			if (RL[BL[j]].peak >=1) { reads +=1; mean +=1; score_partial += log((reads+a_b-1)/(W+a_b+b_b-1)); } // + log(W/reads);
			else score_partial += log((W-reads+b_b-1)/(W+a_b+b_b-1));// + log(W/(W-reads)); 
			lambda += RL[BL[i]].peak;
			
                        ws_corrected += (float)RL[BL[j]].adjusted_size; // *RL[BL[j]].GC_corr;
                        if (reads >= MIN_READS_PER_FRAGMENT && ws_corrected >= MIN_FRAG_LENGTH) // extra conditions to consider it a valid segment...
                        //if (reads >= 5 && ws_corrected >= 5000 && W/reads < intra_frag_gap*5 ) // extra conditions to consider it a valid segment...
                        {
                                blockscore = score_partial+beta_const;  
                                if (USE_SL ==1)  blockscore -= reads + reads*log(reads/lambda); // with segment length likelihood added 

                                if (j ==0) score0 = blockscore + block_open_penalty_1; // +prior_size + gap_length_penalty; 
                                else score0 = RL[BL[j-1]].bscore[0] + blockscore + block_open_penalty_1; // + prior_size + gap_length_penalty;
                                if (score0 > RL[BL[i]].bscore[0])
                                {
                                        RL[BL[i]].bscore[0] = score0; RL[BL[i]].previous[0] = (j-1); RL[BL[i]].reads_window = reads;
                                        RL[BL[i]].score_window = blockscore; RL[BL[i]].mean = mean/W; 
                                }
                        }
                        j--;
                }
		if (RL[BL[i]].mean < 0 && RL[BL[i]].mean < -2000) // empty block 
		{
			//fprintf(stdout,"empty block is preferred %d %d %f %d\n",i,RL[BL[i]].previous[0],RL[BL[i]].bscore[0],RL[BL[i]].peak);
			//fprintf(stdout,"fragment is preferred %d %d %d %f\n",i,BL[i],RL[BL[i]].previous[0],RL[BL[i]].bscore[0]);
			RL[BL[i]].mean = 0;eb = RL[BL[i]].previous[0];
			//j =i; while (j > eb) 
			{
			//	RL[BL[j]].previous[0] = j-1; j--; 
			}
		}
        }
}

// not used at the moment, instead of DP over fragments, do simple HMM 
int hmm_sc(struct BLOCK_READ* RL,int blocks,int* BL,int validblocks,double read_density,double block_open_penalty,double length_penalty,double intra_frag_gap,double inter_frag_gap)
{
        int i = 0,j=0,eb=0,reads=0;
        double ws=0,ws_corrected=0,blockscore=0,mean=0,score0,lambda,score_partial;
        double log_sizes[BLOCK_SIZE]; for (i=1;i<BLOCK_SIZE;i++) log_sizes[i] = log(i)-log(BLOCK_SIZE);
	
	double estimates[3]; estimate_parameters_sc(RL,blocks,estimates);

        read_density = BLOCK_SIZE/inter_frag_gap;
        double logRD = log(BLOCK_SIZE)-log(inter_frag_gap); // average frequency of non-empty blocks outside fragments
        double ld_f = (double)BLOCK_SIZE/intra_frag_gap; // expected frequency of MDA clusters per BLOCK
        double log_ld = log(ld_f), log_ld1 = log(1.0-ld_f);
        fprintf(stderr,"inside simple markov chain based loop, penalties Long-gaps %f small-gaps %f %f %f \n",logRD,log_ld,intra_frag_gap,inter_frag_gap);

        double block_open_penalty_1 = (block_open_penalty+length_penalty);
        double fragment_p = block_open_penalty_1; double fragment_p_1 = log(1.0-exp(fragment_p));
	

        i = 0;
        score0 = RL[BL[i]].reads*(logRD + log_sizes[(int)RL[BL[i]].adjusted_size])-read_density*(float)RL[BL[i]].adjusted_size/BLOCK_SIZE;
        RL[BL[i]].bscore[0] = score0;  // score if this block is outside fragments
        if (RL[BL[i]].reads > 0) RL[BL[i]].bscore[1] = fragment_p + log_ld; // score if this block is part of a fragment
        else RL[BL[i]].bscore[1] = -1000000;
        RL[BL[i]].W = (float)RL[BL[i]].adjusted_size/BLOCK_SIZE;
        for (i=1;i<validblocks;i++)
        {
                score0 = RL[BL[i]].reads*(logRD + log_sizes[(int)RL[BL[i]].adjusted_size])-read_density*(float)RL[BL[i]].adjusted_size/BLOCK_SIZE;
                RL[BL[i]].bscore[0] = RL[BL[i-1]].bscore[0] + score0; RL[BL[i]].previous[0] = 0;
                if (RL[BL[i-1]].bscore[1] + score0 > RL[BL[i]].bscore[0])
                {
                        RL[BL[i]].bscore[0] = RL[BL[i-1]].bscore[1] + score0; RL[BL[i]].previous[0] = 1;
                }
                RL[BL[i]].bscore[0] += fragment_p_1;


                if (RL[BL[i]].peak >=1) score0 = log_ld; else score0 = log_ld1;
                RL[BL[i]].bscore[1] = score0 + RL[BL[i-1]].bscore[1]; RL[BL[i]].previous[1] = 1;  // continue in fragment

                if (RL[BL[i-1]].bscore[0] + fragment_p + score0 > RL[BL[i]].bscore[1] && RL[BL[i]].reads > 0) // open new fragment
                {
                        RL[BL[i]].bscore[1] = RL[BL[i-1]].bscore[0] + fragment_p + score0; RL[BL[i]].previous[1] = 0;
                }

                if (i%500000 ==0) fprintf(stderr,"done for block %d %d\n",i,blocks);
        }

	int c=0,first=0,last=0,prevstart=0,fosmids=0;
	// quick estimate of # of fragments 
	i = validblocks-1; c = 0; // assume that we are outside fragment at last block
	while (i >= 0)
	{
		if (c == 0 && RL[BL[i]].previous[c] ==1)
		{
			 c = 1; last = BL[i];
		}
		else if (c ==1 && RL[BL[i]].previous[c] ==0)
		{
			c = 0; first = BL[i];
			fprintf(stdout,"new fragment chr?:%d-%d length %d gap %d\n",first*BLOCK_SIZE,last*BLOCK_SIZE,(last-first)*BLOCK_SIZE,(prevstart-last)*BLOCK_SIZE);
			prevstart = first; 
			fosmids++;
		}
		i--;
	}
	fprintf(stderr,"need to implement tracepack %d  fosmids\n",fosmids); 
}




// print intra-fragment gaps and coverage per block in compact format, useful for quick visualization 
void print_block(struct BLOCK_READ* RL,int start,int end,int flag)
{
	int j=0;int empty =0,empty1=0; float gaps=0; 
	int prev = 0, peaks = 0, sign = 0; // number of changes in coverage between adjacent bins, only positive, if two changes very close, disregard...
	if (flag ==0) 
	{
		for (j=start;j<=end;j++) 
		{ 
			if (RL[j].valid == '1' && RL[j].counts ==0) fprintf(stdout,"|."); 
			else if (RL[j].valid == '1') fprintf(stdout,"|%d",RL[j].counts); 
			else fprintf(stdout,"*"); 
		} fprintf(stdout," counts\n");
	}
	else
	{
		double seglen[3]={0,1.0e-6,0}, gaplen[3] = {0,1.0e-6,0}; int s=0;
		fprintf(stdout,"\ncounts-gap    ");
		for (j=start;j<=end;j++) 
		{ 
			if (RL[j].valid == '0') empty1++; 
			else 
			{
				// mean gap and segment lengths	
				if ((RL[j].counts==0 || j == end) && s > 0) { seglen[0] +=s; seglen[2] += s*s; seglen[1] += 1; s=0; } 
				if (RL[j].counts > 0) 
				{
					if (empty > 0) { gaplen[0] += empty; gaplen[2] += empty*empty; gaplen[1] +=1; } 
					s+=1;
				}

				if (RL[j].counts ==0) 
				{ 
					if (empty ==0) gaps += RL[j].mappability; empty++; 
				} 
				else
				{
					if (empty1+ empty < 2) fprintf(stdout,"|%d|",RL[j].counts);
					else if (empty1 >= 10) fprintf(stdout,"----%d+%d----|%d|",empty*BLOCK_SIZE,empty1*BLOCK_SIZE,RL[j].counts);  
					else fprintf(stdout,"----%d----|%d|",empty*BLOCK_SIZE,RL[j].counts);  
					empty =0; empty1=0; 
				}
				//if (RL[j].peak > 0) { fprintf(stdout,"+"); peaks++; } 
				if (RL[j].counts > prev && sign ==0) { sign = 1; } 
				else if (RL[j].counts < prev && sign ==1) { sign = 0; fprintf(stdout,"+"); peaks++; } 
				//if (RL[j].counts > prev && sign ==0) { fprintf(stdout,"+"); peaks++; sign = 1; } 
				//else if (RL[j].counts < prev) sign = 0; 
				prev = RL[j].counts; 
			}
		} 
		double std1 = sqrt(seglen[2]/seglen[1] - seglen[0]*seglen[0]/(seglen[1]*seglen[1]));
		double std2 = sqrt(gaplen[2]/gaplen[1] - gaplen[0]*gaplen[0]/(gaplen[1]*gaplen[1]));
		if (empty > 0) fprintf(stdout,"----%d----",empty*BLOCK_SIZE);  
		if (flag >= 2) fprintf(stdout," changes %d gaps %0.2f %d ld %0.2f seg:%0.3f:%d:%0.3f gap:%0.3f:%d:%0.2f \n\n",peaks,gaps,(end-start-empty1)*BLOCK_SIZE,(float)((end-start-empty1)*BLOCK_SIZE)/(peaks+1.0e-4),seglen[0]/seglen[1],(int)seglen[1],std1,BLOCK_SIZE*gaplen[0]/gaplen[1],(int)gaplen[1],std2);
		else fprintf(stdout,"inter-frag %0.2f %d %0.2f\n\n",gaps,(end-start-empty1)*BLOCK_SIZE,(float)((end-start-empty1)*BLOCK_SIZE)/(gaps+1.0e-4));
		
	}
}
