// RL[i].W incorporates GC correction 
double segment_likelihood(struct BLOCK_READ* RL, int* BL,int fb,int lb,double* log_sizes,double PARAMS[]) // first block .... last block
{
	int i=0; double mean=0,W=0,score_partial=0,reads=0,segment_score=0,sum_partial=0; 
	for (i=fb;i<lb;i++)
	{
		if (RL[i].valid != '0') { W += RL[i].W; reads += RL[i].reads;  segment_score += RL[i].score_partial; } 
	}
	mean = reads/W; segment_score += reads*log(mean)-mean*W; // same as below
	PARAMS[0] = reads; PARAMS[1] = W; PARAMS[2] = mean; 
	return segment_score; 
	//fprintf(stderr,"segscore %f mean %f reads %f W %f \n",segment_score,mean,reads,W);
}

// score of segment assuming no fragment overlaps
double background_likelihood(struct BLOCK_READ* RL, int* BL,int fb,int lb,double* log_sizes,double read_density)
{
	int i=0; double segment_score=0, reads=0,W=0;
	for (i=fb;i<lb;i++)
	{
		//segment_score += RL[i].reads*(logRD + log_sizes[(int)RL[i].adjusted_size]-log(RL[i].GC_corr))-read_density*(float)RL[i].adjusted_size/(BLOCK_SIZE*RL[i].GC_corr);
		if (RL[i].valid != '0') { W += RL[i].W; reads += RL[i].reads;  segment_score += RL[i].score_partial; } 
	}
	segment_score += reads*log(read_density)-read_density*W;
	return segment_score;
}

// expected  sum of two read-density | union instead of sum
double sum_coverages(double m1,double m2)
{
	int i=0,s=BLOCK_SIZE/BS1; double ne=0;
	for (i=0;i<1000;i++)
	{
		if (drand48() < m1/s || drand48() < m2/s) ne +=1;
	}
	ne /=1000; ne *=s;
	//if (m1+m2 < ne) return m1+m2; else 
	return ne;
}

// N=1, N=2,non-overlap/overlap
// best segmentation with 2 overlapping fragments  // this will also work for 1/2 non-overlapping fragments with minor changes 
double find_2seg(struct BLOCK_READ* RL, int* BL,int fb,int lb,double* log_sizes,double frag_penalty)
{
	double mean1[3]={0,0,0},mean2[3]={0,0,0},mean3[3]={0,0,0}; double logRD = log(read_density_global);
	int i=0,j=0,k=0;
	double l1,l2,l3,l2_1,l,l2_empty;
	double bestseg[7] = {-10000000,0,0,0,0,0,0};
	double best2seg[7] = {-10000000,0,0,0,0,0,0};
	double middlem = 0;
	int min_sl = 10; // minimum segment length units of 100 bp 
	for (i=fb+min_sl;i<lb-min_sl*2;i+=2)
	{
		if (RL[i].valid == '0') continue; 
		j = i+1; 
		l1 = 	segment_likelihood(RL,BL,fb,i,log_sizes,mean1); 
		l2_empty = background_likelihood(RL,BL,i,j,log_sizes,read_density_global);  // score of segment with no fragment 
		l2 = 	segment_likelihood(RL,BL,i,j,log_sizes,mean2);  l3 = 	segment_likelihood(RL,BL,j,lb,log_sizes,mean3);
		// further speedup by calculating updates

		while (j<lb-min_sl)
		{
			if (RL[j].valid != '0')
			{
				middlem= sum_coverages(mean1[2],mean3[2]);
				l2_1 = l2 - (mean2[0]*log(mean2[2]) - mean2[0]) + mean2[0]*log(middlem)-mean2[0]*(middlem)/mean2[2]; // using a + b as mean for middle segment
				l = l1+l3 + l2_1 + 2*frag_penalty;
				if (mean2[2] > mean1[2] && mean2[2] > mean3[2] && l > bestseg[0] && j-i >= min_sl) 
				{ 
					bestseg[0] = l; bestseg[1] = i; bestseg[2]= j; bestseg[3] = mean1[2]; bestseg[4] = mean3[2];bestseg[6] = middlem; bestseg[5] = mean2[2]; 
				} 
				//fprintf(stdout,"triple-seg %d...%d %f %f %f\n",i,j,l1,l2_1,l3);

				l = l1 + l2_empty + l3 + 2*frag_penalty; 
				if (l > best2seg[0]) // two non-overlapping fragments 
				{
					best2seg[0] = l; best2seg[1] = i; best2seg[2]= j; best2seg[3] = mean1[2]; best2seg[4] = mean3[2]; best2seg[5] = mean2[2]; 
				}

				// update l2,l3, mean2[x], mean3[x]  |  j-th bin is being added to segment2 and removed from segment3
				l2 -= mean2[0]*log(mean2[2])-mean2[0]; 
				mean2[1] += RL[j].W; mean2[0] += RL[j].reads;  mean2[2] = mean2[0]/mean2[1]; 
				l2 += log(RL[j].W)*RL[j].reads; l2 += mean2[0]*log(mean2[2])-mean2[0];

				l3 -= mean3[0]*log(mean3[2])-mean3[0];
				mean3[1] -= RL[j].W; mean3[0] -= RL[j].reads;  mean3[2] = mean3[0]/mean3[1]; 
				l3 -= log(RL[j].W)*RL[j].reads; l3 += mean3[0]*log(mean3[2])-mean3[0];
	
				l2_empty += RL[j].reads*(logRD + log_sizes[(int)RL[j].adjusted_size]-log(RL[j].GC_corr))-read_density_global*(float)RL[j].adjusted_size/(BLOCK_SIZE*RL[j].GC_corr);
			}
			j +=1; 
		}
			
	}

	// N=1
	double ls = segment_likelihood( RL,BL,fb,lb,log_sizes,mean2)+frag_penalty; fprintf(stdout,"single %0.2f %0.2f %d...%d\n",ls,mean2[2],fb,lb);

	// N=2, no overlap 
	fprintf(stdout,"best 2-seg %d....%d....%d....%d %0.2f | ",fb,(int)best2seg[1],(int)best2seg[2],lb,best2seg[0]);
	fprintf(stdout,"%0.0f:%0.2f %0.0f:%0.2f %0.0f:%0.2f \n",BLOCK_SIZE*(best2seg[1]-fb),best2seg[3],BLOCK_SIZE*(best2seg[2]-best2seg[1]),best2seg[5],BLOCK_SIZE*(lb-best2seg[2]),best2seg[4]);
	// N=2, overlapping by at least 1000 bases
	fprintf(stdout,"best overlap-seg %d....%d....%d....%d %0.2f | ",fb,(int)bestseg[1],(int)bestseg[2],lb,bestseg[0]);
	fprintf(stdout,"%0.0f:%0.2f %0.0f:%0.2f:middle:%0.2f %0.0f:%0.2f \n",BLOCK_SIZE*(bestseg[1]-fb),bestseg[3],BLOCK_SIZE*(bestseg[2]-bestseg[1]),bestseg[6],bestseg[5],BLOCK_SIZE*(lb-bestseg[2]),bestseg[4]);
	if (best2seg[0] > bestseg[0] && best2seg[0] > ls + 5) return best2seg[0]; 
	else if (bestseg[0] > best2seg[0] && bestseg[0] > ls + 5) return bestseg[0]; 
	return ls; // bogus
}


// 2 + 1 segmentation | N=3,non-overlap, N=3,1+2 should work...
double find_3seg(struct BLOCK_READ* RL, int* BL,int fb,int lb,double* log_sizes,double frag_penalty,int frags)
{
	// first find list of all potential breakpoints using weak penalty...  int* bkpts = calloc(sizeof(int),100); 
	int i=0; int k = lb; int prev = RL[lb].previous[0];
	int* bkpts = calloc(sizeof(int),frags*4); int b = 0; 

	while (1)
	{
		if (k-BL[prev] >=5 && BL[prev+1] > fb) bkpts[b++] = BL[prev]; 
		if (prev <= 0 || BL[prev+1] <= fb) break;
		k = BL[prev]; prev = RL[k].previous[0];
	}
	for (i=0;i<b;i++) fprintf(stdout,"bkpts %d:%d ",i,bkpts[i]); fprintf(stdout,"b %d\n",b);

	double mean[3]={0,0,0}; double ls1,ls2;

	if (frags <=4) 
	{
		ls2 = find_2seg(RL,BL,fb,lb,log_sizes,frag_penalty);
                //fprintf(stdout,"2way-seg ll %0.2f \n",ls2);
		// if the 2-way segmentation is much better, don't call the next loop
	}

	if (frags >=3) 
	{
		fprintf(stdout,"\n");
		ls2 = find_2seg(RL,BL,fb,bkpts[0],log_sizes,frag_penalty); 
		ls1 = segment_likelihood( RL,BL,bkpts[0],lb,log_sizes,mean)+frag_penalty;
		fprintf(stdout,"------ 2overlap+1 ll %0.2f %0.2f sum %0.2f-----------\n \n",ls1,ls2,ls1+ls2); 

		fprintf(stdout,"\n");
		ls1 = segment_likelihood( RL,BL,fb,bkpts[b-1],log_sizes,mean)+frag_penalty;
		ls2 = find_2seg(RL,BL,bkpts[b-1],lb,log_sizes,frag_penalty); 
		fprintf(stdout,"-----1+2overlap ll %0.2f %0.2f sum %0.2f---------------\n \n",ls1,ls2,ls1+ls2); 
	}
	
	free(bkpts);

}
