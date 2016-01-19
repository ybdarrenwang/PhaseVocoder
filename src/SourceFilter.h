#ifndef __SOURCEFILTER_H__
#define __SOURCEFILTER_H__

#include <math.h>
#include "VocoderFunctions.h"

#define order 16

class SourceFilter{
public:
	void clear_mem(int word_no);
	void SFDecomp_sync(double *data,double *excitation, int sample_no, 
						int word_no, int* WordStart, int* WordMid, 
						int* WordEnd, int frame_shift, int frame_len, double** Pitch);
	void SFDecomp_non_sync(double *data,double *excitation, int sample_no, 
							int word_no, int* WordStart, int* WordMid, 
							int* WordEnd, int frame_shift, int frame_len);
	void SFRecomb_sync(short* data, double* excitation, int sample_no, int word_no, 
					double* WordDurMul, int* WordStart, int* WordMid, int* WordEnd, 
					int frame_shift, int frame_len, int* syn_bias);
	void SFRecomb_FD(double** spec, int word_index, int frame_index, int WordMid);
	double LPC_From_Data(double *data, double *lpc, int n, int m);
	void LPC_Predict(double *coeff, double *prime, int m, double *data, long n);
	void LPC_SpecEnv(double *coeff, double *SpecEnv_mag, double *SpecEnv_pha);

private:
	int** vowel_mark; // vowel_mark[i][j]: the jth vowel f0 mark of ith word
	int** frame_no; // frame_no[i][j]: the number of frames of consonant(j=0) or vowel(j=1) of ith word
	double ****coeff; // coeff[i][j][k][m]: the mth LPC coefficient of kth frame in
					  // consonant(j=0) or vowel(j=1) of ith word
	VocoderFunc my_func;
};

#endif