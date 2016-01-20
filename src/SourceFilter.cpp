#include "SourceFilter.h"

void SourceFilter::clear_mem(int word_no)
{
	int i,j;

	for(i=0; i<word_no; i++)
	{
		for(j=0; j<frame_no[i][0]; j++)
			delete [] coeff[i][0][j];
		for(j=0; j<frame_no[i][1]; j++)
			delete [] coeff[i][1][j];

		delete [] coeff[i][0];
		delete [] coeff[i][1];
		delete [] coeff[i];
	}
	delete [] coeff;

	delete [] frame_no[0];
	delete [] frame_no[1];
	delete [] frame_no;
}

// ============================================================================ //
//  < Source-Filter Decomposition >
//  < Pitch-Synchronous(vowel) LPA  >
// ============================================================================ //
void SourceFilter::SFDecomp_sync(double *data,double *excitation, int sample_no,
							int word_no, int* WordStart, int* WordMid, int* WordEnd, 
							int frame_shift, int frame_len, double** Pitch)
{
	int i, j, k, m;
	int new_frame_len, vowel_mark_no;
	double error, *sf_frame;
	vowel_mark = new int *[word_no];
	frame_no = new int *[word_no];
	coeff = new double ***[word_no];
	VocoderFunc vocoder_function;
	double *windows;

	// Analysis: Frame Length = 20ms; Frame Shift = 10ms; Order = 16
	for (i=0; i<word_no; i++)
	{
		frame_no[i] = new int [2];
		coeff[i] = new double **[2];

		// Consonant part
		frame_no[i][0] = WordMid[i]-WordStart[i];
		coeff[i][0] = new double *[frame_no[i][0]];

		new_frame_len = frame_shift*2;
		sf_frame = new double[new_frame_len];

		windows = new double[new_frame_len];
		vocoder_function.Hamming(new_frame_len, windows);

		for (j=0; j<frame_no[i][0]; j++)
		{
			coeff[i][0][j] = new double[order];

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
		coeff[i][1] = new double *[frame_no[i][1]];

		for (j=0; j<frame_no[i][1]; j++)
		{	
			coeff[i][1][j] = new double[order];
			new_frame_len = vowel_mark[i][j+2] - vowel_mark[i][j];
			sf_frame = new double[new_frame_len];

			windows = new double[new_frame_len];
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
	double *prime = new double[order];
	double sf_predict[1];
	double *int_coeff = new double[order]; // interpolated LPC coefficient
	double *curr_coeff = new double[order];
	double *prev_coeff = new double[order];
	double *next_coeff = new double[order];
	double b;

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
					new_frame_len = (double)(vowel_mark[i][1]+vowel_mark[i][2])/2;			

					for (m=0; m<order; m++)
						prime[m] = data[sample_index-order+m];
				}
				else if (j==frame_no[i][0]+frame_no[i][1]-1) // end
					new_frame_len = (WordEnd[i]-1)*frame_shift + frame_len - sample_index;
				else
					new_frame_len = (double)(vowel_mark[i][j-frame_no[i][0]+1]+vowel_mark[i][j-frame_no[i][0]+2])/2 - (sample_index-WordMid[i]*frame_shift);

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
					b = (double)k/new_frame_len;
					if ((i==0||WordStart[i]!=WordEnd[i-1]) && j==0)
						b = 0;

					for (m=0; m<order; m++) 
						int_coeff[m] = (1.0-b)*curr_coeff[m] + b*prev_coeff[m];
				}
				else
				{
					b = (double)(new_frame_len-k)/new_frame_len;
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
//  < Source-Filter Decomposition >
//  < Fixed frame length LPA  >
// ============================================================================ //
void SourceFilter::SFDecomp_non_sync(double *data,double *excitation, int sample_no,
									int word_no, int* WordStart, int* WordMid, int* WordEnd, 
									int frame_shift, int frame_len)
{
	// setting
	int new_frame_len = frame_len;

	int i, j, k, m;
	double error, *sf_frame;
	vowel_mark = new int *[word_no];
	frame_no = new int *[word_no];
	coeff = new double ***[word_no];
	VocoderFunc vocoder_function;
	double *windows;

	// Analysis: Frame Length = 20ms; Frame Shift = 10ms; Order = 20
	for (i=0; i<word_no; i++)
	{
		frame_no[i] = new int [2];
		coeff[i] = new double **[2];
		sf_frame = new double[new_frame_len];
		windows = new double[new_frame_len];
		vocoder_function.Hamming(new_frame_len, windows);

		// Consonant part
		frame_no[i][0] = WordMid[i]-WordStart[i];
		coeff[i][0] = new double *[frame_no[i][0]];

		for (j=0; j<frame_no[i][0]; j++)
		{
			coeff[i][0][j] = new double[order];

			for (k=0; k<new_frame_len; k++) // read data
				sf_frame[k] = data[frame_shift*(WordStart[i]+j) + k]*windows[k];

			error = lpc.From_Data(sf_frame, coeff[i][0][j], new_frame_len, order);
		}

		// Vowel part
		frame_no[i][1] = WordEnd[i] - WordMid[i];
		coeff[i][1] = new double *[frame_no[i][1]];

		for (j=0; j<frame_no[i][1]; j++)
		{	
			coeff[i][1][j] = new double[order];

			for (k=0; k<new_frame_len; k++) // read data
				sf_frame[k] = data[frame_shift*(WordMid[i]+j) + k]*windows[k];

			error = lpc.From_Data(sf_frame, coeff[i][1][j], new_frame_len, order);
		}
		delete [] sf_frame;
		delete [] windows;
	}

	// Filtering
	int sample_index;
	double *prime = new double[order];
	double sf_predict[1];
	double *int_coeff = new double[order]; // interpolated LPC coefficient
	double *curr_coeff = new double[order];
	double *prev_coeff = new double[order];
	double *next_coeff = new double[order];
	double b;

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
					new_frame_len = frame_shift*1.5;			

					for (m=0; m<order; m++)
						prime[m] = data[sample_index-order+m];
				}
				else if (j==frame_no[i][0]+frame_no[i][1]-1) // end
					new_frame_len = (WordEnd[i]-1)*frame_shift + frame_len - sample_index;
				else
					new_frame_len = frame_shift;

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
					b = (double)k/new_frame_len;
					if ((i==0||WordStart[i]!=WordEnd[i-1]) && j==0)
						b = 0;

					for (m=0; m<order; m++) 
						int_coeff[m] = (1.0-b)*curr_coeff[m] + b*prev_coeff[m];
				}
				else
				{
					b = (double)(new_frame_len-k)/new_frame_len;
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
void SourceFilter::SFRecomb_sync(short* data, double* excitation, int sample_no, int word_no, 
							double* WordDurMul, int* WordStart, int* WordMid, int* WordEnd, 
							int frame_shift, int frame_len, int* syn_bias)
{
	int i, j, k, m;
	int new_frame_len, sample_index;
	double sf_predict[1];
	double *int_coeff = new double[order]; // interpolated LPC coefficient
	double *curr_coeff = new double[order];
	double *prev_coeff = new double[order];
	double *next_coeff = new double[order];
	double b;
	double *prime = new double[order];
	double ave_syn_bias;
	bool saturate;

	for (i=0; i<word_no; i++)
	{
		ave_syn_bias = (double)syn_bias[i]/frame_no[i][1];

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
					new_frame_len = WordDurMul[i]*(double)(vowel_mark[i][1]+vowel_mark[i][2])/2;

					for (m=0; m<order; m++)
						prime[m] = data[sample_index-order+m];
				}
				else if (j==frame_no[i][0]+frame_no[i][1]-1) // end
					new_frame_len = (WordEnd[i]-1)*frame_shift + frame_len - sample_index;
				else
					new_frame_len = WordDurMul[i]*(double)(vowel_mark[i][j-frame_no[i][0]+1]+vowel_mark[i][j-frame_no[i][0]+2])/2 - (sample_index-WordMid[i]*frame_shift);

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
					b = (double)k/new_frame_len;
					if ((i==0||WordStart[i]!=WordEnd[i-1]) && j==0)
						b = 0;

					for (m=0; m<order; m++) 
						int_coeff[m] = (1.0-b)*curr_coeff[m] + b*prev_coeff[m];
				}
				else
				{
					b = (double)(new_frame_len-k)/new_frame_len;
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

// ============================================================================ //
//  < Source-Filter Recombination at Frequency Domain >
// ============================================================================ //
void SourceFilter::SFRecomb_FD(double** spec, int word_index, int frame_index, int WordMid)
{
	double *SpecEnv_mag = new double[FFT_SIZE/2+1];
	double *SpecEnv_pha = new double[FFT_SIZE/2+1];
	double mag, pha;
	int i;

	if (frame_index >= WordMid)
		lpc.SpecEnv(coeff[word_index][1][frame_index-WordMid], SpecEnv_mag, SpecEnv_pha);
	else
		lpc.SpecEnv(coeff[word_index][0][frame_index], SpecEnv_mag, SpecEnv_pha);

	for (i=0; i<FFT_SIZE/2+1; i++)
	{
		mag = lpc.ABS2(spec[0][i], spec[1][i])*SpecEnv_mag[i];
		pha = atan2(spec[1][i], spec[0][i])+SpecEnv_pha[i];

		spec[0][i] = mag*cos(pha);
		spec[1][i] = mag*sin(pha);
	}

	delete [] SpecEnv_mag;
	delete [] SpecEnv_pha;
}
