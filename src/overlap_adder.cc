#include "overlap_adder.h"

// ============================================================================ //
//  < Pitch-Sync. Overlap-Add >
// ============================================================================ //
int OverlapAdder::PSOLA(float **frame, float *output, float *SentenceCoeff, float *windows,
					  int sample_index, float *Pitch, int WordStart, int WordMid, int WordEnd,
					  int frame_len, int frame_shift, int bias)
{
	int i,j,k;
	int syn_bias; // a shift before overlap-add for pitch-sync
	int total_syn_bias = 0;
	int syn_sample_no;
	float *frame_for_syn;
	float cross_corr, max_cross_corr;
	float T; // peroid in sample number

	for (i=0; i<WordMid-WordStart; i++)
	{
		for (j=0; j<frame_len && sample_index-bias+j>0; j++)
		{
			output[sample_index-bias+j] += frame[i][j]*windows[j];
			SentenceCoeff[sample_index-bias+j] += windows[j]*windows[j];
		}
		sample_index += frame_shift;
	}

	for (i=WordMid-WordStart; i<WordEnd-WordStart; i++)
	{
		T = floor((float)SamplingRate/Pitch[i-(WordMid-WordStart)]);
		syn_sample_no = bias;
		frame_for_syn = new float[syn_sample_no];

		if (sample_index-syn_sample_no >= 0)
		{
			for (j=0; j<syn_sample_no; j++)
				frame_for_syn[j] = output[sample_index-syn_sample_no+j];

			// find max cross-correlation
			max_cross_corr = 0;
			for (j=-1*(int)floor(T/2); j<(int)floor(T/2); j++)
			{
				cross_corr = 0;
				for (k=0; k<syn_sample_no; k++)
					cross_corr += frame_for_syn[k]*frame[i][frame_len-frame_shift-syn_sample_no-j+k];
			
				if (cross_corr > max_cross_corr)
				{
					max_cross_corr = cross_corr;
					syn_bias = j; // <0 if the latter frame shifts backward, and vice versa
				}
			}
		}
		else
			syn_bias = 0;

		for (j=0; j<frame_len && sample_index-bias+syn_bias+j>=0; j++)
		{
			output[sample_index-bias+syn_bias+j] += frame[i][j]*windows[j];
			SentenceCoeff[sample_index-bias+syn_bias+j] += windows[j]*windows[j];
		}
		sample_index += frame_shift+syn_bias;
		total_syn_bias += syn_bias;

		delete [] frame_for_syn;
	}

	return total_syn_bias;
}

// ============================================================================ //
//  < Non-Pitch-Sync. Overlap-Add >
// ============================================================================ //
int OverlapAdder::SOLA(float **frame, float *output, float *SentenceCoeff, float *windows,
							int sample_index, int WordStart, int WordMid, int WordEnd, 
							int frame_len, int frame_shift, int bias)
{
	int i,j;
	bool saturate;

	for (i=0; i<WordEnd-WordStart; i++)
	{
		// check if this frame saturates
		saturate = false;
		for (j=0; j<frame_len && sample_index-bias+j>0 && saturate==false; j++)
		{
			if (abs(frame[i][j]) >= 25000)
				saturate = true;
		}
		if (saturate == true)
		{
			for (j=0; j<frame_len && sample_index-bias+j>0; j++)
				frame[i][j] *= Atten_3dB;
		}

		for (j=0; j<frame_len && sample_index-bias+j>0; j++)
		{
			output[sample_index-bias+j] += frame[i][j]*windows[j];
			SentenceCoeff[sample_index-bias+j] += windows[j]*windows[j];
		}
		sample_index += frame_shift;
	}

	return 0;
}
