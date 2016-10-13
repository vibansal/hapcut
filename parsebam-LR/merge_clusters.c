

// functions for merging clusters or haplotype fragments identified from DP

// k! is missing -> we get both positive and negative likelihood values...
double segment_likelihood(struct BLOCK_READ* RL, int* BL,int fb,int lb,double* log_sizes,double PARAMS[]) // first block .... last block
{
	int i=0; double mean=0,W=0,score_partial=0,reads=0,segment_score=0,sum_partial=0; 
	for (i=fb;i<lb;i++)
	{
		if (RL[i].valid == '0') continue;
		W += RL[i].W; reads += RL[i].reads; 
		segment_score += log(RL[i].W)*RL[i].reads;
		//segment_score += log((float)RL[i].adjusted_size/(RL[i].GC_corr*BLOCK_SIZE))*RL[i].reads;
	}
	mean = reads/W; segment_score += reads*log(mean)-reads;
	//fprintf(stderr,"segscore %f mean %f reads %f W %f \n",segment_score,mean,reads,W);

	//segment_score = score_partial + reads*log(mean/W) - mean;
	//double mean_pb = mean/((float)BLOCK_SIZE*W); segment_score = reads*log(mean_pb) + (W*BLOCK_SIZE-reads)*log(1.0-mean_pb); 
	PARAMS[0] = reads; PARAMS[1] = W; PARAMS[2] = mean; 
	return segment_score; 
}

// score of segment assuming no fragment overlaps
double background_likelihood(struct BLOCK_READ* RL, int* BL,int fb,int lb,double* log_sizes,double read_density)
{
	int i=0; double segment_score=0; double logRD = log(read_density); 
	for (i=fb;i<lb;i++)
	{
		if (RL[i].valid == '0') continue;
		segment_score += RL[i].reads*(logRD + log_sizes[(int)RL[i].adjusted_size]-log(RL[i].GC_corr))-read_density*(float)RL[i].adjusted_size/(BLOCK_SIZE*RL[i].GC_corr);
	}
	return segment_score;
}

// check if a pair of segments/fragments can be merged into a single fragment
double evaluate_pair_segments(struct BLOCK_READ* RL, int* BL,int fb1,int lb1,int fb2,int lb2,double* log_sizes,double read_density,FILE* fp)
{
	double mean1[3]={0,0,0},mean2[3]={0,0,0},mean[3]={0,0,0};
	double segment1_score = segment_likelihood(RL,BL,fb1,lb1,log_sizes,mean1); 
	double segment2_score = segment_likelihood(RL,BL,fb2,lb2,log_sizes,mean2); 
	double middle_score = 0;
	if (lb1 < fb2) middle_score = background_likelihood(RL,BL,lb1,fb2,log_sizes,read_density); 
	double sum = segment1_score + segment2_score+ middle_score;
	double single_score = segment_likelihood(RL,BL,fb1,lb2,log_sizes,mean);  
	double delta = single_score - sum;
	if (delta > -13) fprintf(fp,"merge_fragments %d:%d:%f:%f %d:%d:%f:%f %f | %f %f delta %f\n",fb1,lb1,segment1_score,mean1[2],fb2,lb2,segment2_score,mean2[2],middle_score,single_score,mean[2],single_score-sum);
	if (delta + middle_score > -6) fprintf(fp,"merge_fragments_deletion %d:%d:%f:%f %d:%d:%f:%f %f | %f %f delta %f\n",fb1,lb1,segment1_score,mean1[2],fb2,lb2,segment2_score,mean2[2],middle_score,single_score,mean[2],single_score-sum+middle_score);
}

double evaluate_triple_segments(struct BLOCK_READ* RL, int* BL,int fb1,int lb1,int fb2,int lb2,int fb3,int lb3,double* log_sizes,double read_density,FILE* fp)
{
	// assess if triplet of consecutive segments correspond to two overlapping fragments... 
	// read density in middle fragment = a + b where 'a' and 'b' are read density in adjacent fragments... 
	//segment_score = reads*log(mean_pb) + (W*BLOCK_SIZE-reads)*log(1.0-mean_pb); 
	double mean1[3]={0,0,0},mean2[3]={0,0,0},mean3[3]={0,0,0};
	double segment1_score = segment_likelihood(RL,BL,fb1,lb1,log_sizes,mean1); 
	double segment2_score = segment_likelihood(RL,BL,fb2,lb2,log_sizes,mean2); 
	double segment3_score = segment_likelihood(RL,BL,fb3,lb3,log_sizes,mean3); 
	double sum = segment1_score + segment2_score + segment3_score; 
	double middle_score = 0;
	if (mean2[2] < mean1[2] || mean2[2] < mean3[2] || fb2-lb1 >=3 || fb3 -lb2 >=3) return 0;
	//	if (lb1 < fb2) middle_score = XXXX; if (lb2 < fb3) middle_score += XXXX

	// -lambda*w_i + r_i.log(lambda) + r_i (adjusted) 
	double segment2_score_new = segment2_score - (mean2[0]*log(mean2[2])-mean2[0]) + mean2[0]*log(mean1[2] + mean3[2]) - (mean1[2]+ mean3[2])*mean2[1]; 
	//mean2[0]*log(mean1[2]/BLOCK_SIZE + mean3[2]/BLOCK_SIZE) + (mean2[1]*BLOCK_SIZE-mean2[0])*log(1.0-(mean1[2]/BLOCK_SIZE + mean3[2]/BLOCK_SIZE)); 

	double delta = segment2_score_new-segment2_score;
	if (delta > -6) fprintf(fp,"overlapping_tripl_fragments segment2_score_orig %f:%f new %f:%f OP-DELTA %f\n",segment2_score,mean2[2],segment2_score_new,mean1[2]+mean3[2],segment2_score_new-segment2_score);
}

