#include "pitch_sync_source_filter.h"

// ============================================================================ //
//  < Source-Filter Decomposition >
//  < Pitch-Synchronous(vowel) LPA  >
// ============================================================================ //
void PitchSyncSourceFilter::Decomposition(float *data,float *excitation, int sample_no, int word_no, int* WordStart, int* WordMid, int* WordEnd, int frame_shift, int frame_len, float** Pitch)
{
	int i, j, k, m;
	int new_frame_len, vowel_mark_no;
	float error, *sf_frame;
	vowel_mark = new int *[word_no];
	frame_no = new int *[word_no];
	coeff = new float ***[word_no];
	VocoderFunc vocoder_function;
	float *windows;

	// Analysis: Frame Length = 20ms; Frame Shift = 10ms; Order = 16
	for (i=0; i<word_no; i++)
	{
		frame_no[i] = new int [2];
		coeff[i] = new float **[2];

		// Consonant part
		frame_no[i][0] = WordMid[i]-WordStart[i];
		coeff[i][0] = new float *[frame_no[i][0]];

		new_frame_len = frame_shift*2;
		sf_frame = new float[new_frame_len];

		windows = new float[new_frame_len];
		vocoder_function.Hamming(new_frame_len, windows);

		for (j=0; j<frame_no[i][0]; j++)
		{
			coeff[i][0][j] = new float[order];

			for (k=0; k<new_frame_len; k++) // read data
				sf_frame[k] = data[frame_shift*(WordStart[i]+j) + k]*windows[k];

			error = lpc.From_Data(sf_frame, coeff[i][0][j], new_frame_len, order);
		}
		delete [] sf_frame;
		delete [] windows;

		// Vowel part
		vowel_mark_no = lpc.CalPath(WordEnd[i]-WordMid[i], frame_shift, Pitch[i]);
		vowel_mark[i] = new int[vowel_mark_no];
		lpc.vowel_markPath(vowel_mark[i], vowel_mark_no, frame_shift, Pitch[i]);

		frame_no[i][1] = vowel_mark_no-2;
		coeff[i][1] = new float *[frame_no[i][1]];

		for (j=0; j<frame_no[i][1]; j++)
		{	
			coeff[i][1][j] = new float[order];
			new_frame_len = vowel_mark[i][j+2] - vowel_mark[i][j];
			sf_frame = new float[new_frame_len];

			windows = new float[new_frame_len];
			vocoder_function.Hamming(new_frame_len, windows);

			for (k=0; k<new_frame_len; k++) // read data
				sf_frame[k] = data[frame_shift*WordMid[i] + vowel_mark[i][j] + k]*windows[k];

			error = lpc.From_Data(sf_frame, coeff[i][1][j], new_frame_len, order);

			delete [] sf_frame;
			delete [] windows;
		}
	}

	// Filtering
	int sample_index;
	float *prime = new float[order];
	float sf_predict[1];
	float *int_coeff = new float[order]; // interpolated LPC coefficient
	float *curr_coeff = new float[order];
	float *prev_coeff = new float[order];
	float *next_coeff = new float[order];
	float b;

	for (i=0; i<word_no; i++)
	{
		for (j=0; j<frame_no[i][0]+frame_no[i][1]; j++)
		{
			if (j<frame_no[i][0]) // consonant
			{
				for (m=0; m<order; m++)
					curr_coeff[m] = coeff[i][0][j][m];

				if (j==0) // start
				{
					sample_index = WordStart[i]*frame_shift;
					new_frame_len = frame_shift*1.5;

					if (i==0 || WordStart[i]!=WordEnd[i-1])
					{
						for (m=0; m<order; m++)
							prime[m] = 0;
					}
					else
					{
						for (m=0; m<order; m++)
							prime[m] = data[sample_index-order+m];
					}
				}
				else if(j==frame_no[i][0]-1) // end
					new_frame_len = WordMid[i]*frame_shift - sample_index;
				else
					new_frame_len = frame_shift;

				// record previous LPC coefficient for interpolation
				if (j != 0)
				{
					for (m=0; m<order; m++)
						prev_coeff[m] = coeff[i][0][j-1][m];
				}
				else if (i != 0)
				{
					for (m=0; m<order; m++)
						prev_coeff[m] = coeff[i-1][1][frame_no[i-1][1]-1][m];
				}
				// record next LPC coefficient for interpolation
				if (j != frame_no[i][0]-1)
				{
					for (m=0; m<order; m++)
						next_coeff[m] = coeff[i][0][j+1][m];
				}
				else
				{
					for (m=0; m<order; m++)
						next_coeff[m] = coeff[i][1][0][m];
				}
			}
			else // vowel
			{
				for (m=0; m<order; m++)
					curr_coeff[m] = coeff[i][1][j-frame_no[i][0]][m];

				if (j==frame_no[i][0]) // start
				{
					sample_index = WordMid[i]*frame_shift;
					new_frame_len = (float)(vowel_mark[i][1]+vowel_mark[i][2])/2;			

					for (m=0; m<order; m++)
						prime[m] = data[sample_index-order+m];
				}
				else if (j==frame_no[i][0]+frame_no[i][1]-1) // end
					new_frame_len = (WordEnd[i]-1)*frame_shift + frame_len - sample_index;
				else
					new_frame_len = (float)(vowel_mark[i][j-frame_no[i][0]+1]+vowel_mark[i][j-frame_no[i][0]+2])/2 - (sample_index-WordMid[i]*frame_shift);

				// record previous LPC coefficient for interpolation
				if (j != frame_no[i][0])
				{
					for (m=0; m<order; m++)
						prev_coeff[m] = coeff[i][1][j-frame_no[i][0]-1][m];
				}
				else
				{
					for (m=0; m<order; m++)
						prev_coeff[m] = coeff[i][0][frame_no[i][0]-1][m];
				}
				// record next LPC coefficient for interpolation
				if (j != frame_no[i][0]+frame_no[i][1]-1)
				{
					for (m=0; m<order; m++)
						next_coeff[m] = coeff[i][1][j-frame_no[i][0]+1][m];
				}
				else if (i != word_no-1)
				{
					for (m=0; m<order; m++)
						next_coeff[m] = coeff[i+1][0][0][m];
				}
			}

			for (k=0; k<new_frame_len && sample_index<sample_no; k++)
			{
				// LPC coefficient interpolation
				if (k < new_frame_len-k)
				{
					b = (float)k/new_frame_len;
					if ((i==0||WordStart[i]!=WordEnd[i-1]) && j==0)
						b = 0;

					for (m=0; m<order; m++) 
						int_coeff[m] = (1.0-b)*curr_coeff[m] + b*prev_coeff[m];
				}
				else
				{
					b = (float)(new_frame_len-k)/new_frame_len;
					if ((i==word_no-1||WordStart[i+1]!=WordEnd[i]) && j==frame_no[i][0]+frame_no[i][1]-1)
						b = 0;

					for (m=0; m<order; m++) 
						int_coeff[m] = (1.0-b)*curr_coeff[m] + b*next_coeff[m];
				}

				// LPA
				lpc.Predict(int_coeff, prime, order, sf_predict, 1);
				excitation[sample_index] = data[sample_index] - sf_predict[0];

				// update prime value
				for (m=0; m<order-1; m++)
					prime[m] = prime[m+1];
				prime[m] = data[sample_index];
				
				sample_index++;
			} 
		}
	}
	delete [] prime;
	delete [] int_coeff;
	delete [] curr_coeff;
	delete [] prev_coeff;
	delete [] next_coeff;
}

// ============================================================================ //
//  < Source-Filter Recombination >
//  < Pitch-Synchronous(vowel) LPS  >
// ============================================================================ //
void PitchSyncSourceFilter::Recombination(short* data, float* excitation, int sample_no, int word_no, 
							float* WordDurMul, int* WordStart, int* WordMid, int* WordEnd, 
							int frame_shift, int frame_len, int* syn_bias)
{
	int i, j, k, m;
	int new_frame_len, sample_index;
	float sf_predict[1];
	float *int_coeff = new float[order]; // interpolated LPC coefficient
	float *curr_coeff = new float[order];
	float *prev_coeff = new float[order];
	float *next_coeff = new float[order];
	float b;
	float *prime = new float[order];
	float ave_syn_bias;
	bool saturate;

	for (i=0; i<word_no; i++)
	{
		ave_syn_bias = (float)syn_bias[i]/frame_no[i][1];

		for (j=0; j<frame_no[i][0]+frame_no[i][1]; j++)
		{
			saturate = false;

			if (j<frame_no[i][0]) // consonant
			{
				for (m=0; m<order; m++)
					curr_coeff[m] = coeff[i][0][j][m];

				if (j==0) // start
				{
					sample_index = WordStart[i]*frame_shift;
					new_frame_len = WordDurMul[i]*frame_shift*1.5;

					if (i==0 || WordStart[i]!=WordEnd[i-1])
					{
						for (m=0; m<order; m++)
							prime[m] = 0;
					}
					else
					{
						for (m=0; m<order; m++)
							prime[m] = data[sample_index-order+m];
					}
				}
				else if(j==frame_no[i][0]-1) // end
					new_frame_len = WordMid[i]*frame_shift - sample_index;
				else
					new_frame_len = WordDurMul[i]*frame_shift*(j+1.5) - (sample_index-WordStart[i]*frame_shift);

				// record previous LPC coefficient for interpolation
				if (j != 0)
				{
					for (m=0; m<order; m++)
						prev_coeff[m] = coeff[i][0][j-1][m];
				}
				else if (i != 0)
				{
					for (m=0; m<order; m++)
						prev_coeff[m] = coeff[i-1][1][frame_no[i-1][1]-1][m];
				}
				// record next LPC coefficient for interpolation
				if (j != frame_no[i][0]-1)
				{
					for (m=0; m<order; m++)
						next_coeff[m] = coeff[i][0][j+1][m];
				}
				else
				{
					for (m=0; m<order; m++)
						next_coeff[m] = coeff[i][1][0][m];
				}
			}
			else // vowel
			{
				for (m=0; m<order; m++)
					curr_coeff[m] = coeff[i][1][j-frame_no[i][0]][m];

				if (j==frame_no[i][0]) // start
				{
					sample_index = WordMid[i]*frame_shift;
					new_frame_len = WordDurMul[i]*(float)(vowel_mark[i][1]+vowel_mark[i][2])/2;

					for (m=0; m<order; m++)
						prime[m] = data[sample_index-order+m];
				}
				else if (j==frame_no[i][0]+frame_no[i][1]-1) // end
					new_frame_len = (WordEnd[i]-1)*frame_shift + frame_len - sample_index;
				else
					new_frame_len = WordDurMul[i]*(float)(vowel_mark[i][j-frame_no[i][0]+1]+vowel_mark[i][j-frame_no[i][0]+2])/2 - (sample_index-WordMid[i]*frame_shift);

				new_frame_len += floor(ave_syn_bias*(j-frame_no[i][0]+1)) - floor(ave_syn_bias*(j-frame_no[i][0]));

				// record previous LPC coefficient for interpolation
				if (j != frame_no[i][0])
				{
					for (m=0; m<order; m++)
						prev_coeff[m] = coeff[i][1][j-frame_no[i][0]-1][m];
				}
				else
				{
					for (m=0; m<order; m++)
						prev_coeff[m] = coeff[i][0][frame_no[i][0]-1][m];
				}
				// record next LPC coefficient for interpolation
				if (j != frame_no[i][0]+frame_no[i][1]-1)
				{
					for (m=0; m<order; m++)
						next_coeff[m] = coeff[i][1][j-frame_no[i][0]+1][m];
				}
				else if (i != word_no-1)
				{
					for (m=0; m<order; m++)
						next_coeff[m] = coeff[i+1][0][0][m];
				}
			}		

			for (k=0; k<new_frame_len && sample_index<sample_no; k++)
			{
				// LPC coefficient interpolation
				if (k < new_frame_len-k)
				{
					b = (float)k/new_frame_len;
					if ((i==0||WordStart[i]!=WordEnd[i-1]) && j==0)
						b = 0;

					for (m=0; m<order; m++) 
						int_coeff[m] = (1.0-b)*curr_coeff[m] + b*prev_coeff[m];
				}
				else
				{
					b = (float)(new_frame_len-k)/new_frame_len;
					if ((i==word_no-1||WordStart[i+1]!=WordEnd[i]) && j==frame_no[i][0]+frame_no[i][1]-1)
						b = 0;

					for (m=0; m<order; m++) 
						int_coeff[m] = (1.0-b)*curr_coeff[m] + b*next_coeff[m];
				}
			
				// LPS
				lpc.Predict(int_coeff, prime, order, sf_predict, 1);
				
				data[sample_index] = floor(sf_predict[0]+excitation[sample_index]);
				if (abs(data[sample_index]) >= 25000) // a-posteriori constrain of output signal
					saturate = true;

				// update prime value
				for (m=0; m<order-1; m++)
					prime[m] = prime[m+1];
				prime[m] = data[sample_index];

				sample_index++;
			}
			if (saturate == true)
			{
				for (k=0; k<new_frame_len; k++)
					data[sample_index-k] *= Atten_3dB;
			}
		}
	}

	// release memory
	for (i=0; i<word_no; i++)
		delete [] vowel_mark[i];
	delete [] vowel_mark;
	delete [] prime;
	delete [] int_coeff;
	delete [] curr_coeff;
	delete [] prev_coeff;
	delete [] next_coeff;
}
