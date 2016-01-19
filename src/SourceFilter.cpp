#include "stdafx.h"
#include "SourceFilter.h"

double MIN2(int a,int b)
{
	if (a<b)
		return (double)a;
	else
		return (double)b;
}

// calculate how many fundamental frequency vowel_marks needed
// return: the number of vowel_marks
int LPC_CalPath(int frame_no, int frame_shift, double *Pitch)
{
	int frame_ptr = 0, n = 0;
	double my_ptr = 0, T;

	while(my_ptr < frame_no*frame_shift)
	{
		if (Pitch[frame_ptr] == 0)
		{
			if (frame_ptr == 0)// if the first frame has no pitch
			{
				while(Pitch[frame_ptr] == 0)
					frame_ptr++;
				while(frame_ptr > 0)
				{
					Pitch[frame_ptr-1] = Pitch[frame_ptr];
					frame_ptr--;
				}
			}
			else
				Pitch[frame_ptr] = Pitch[frame_ptr-1];
		}

		T = (double)SamplingRate/Pitch[frame_ptr];
		my_ptr += T;
		frame_ptr = floor(my_ptr/frame_shift);
		n++;
	}
	
	return n;
}

// record the fundamental frequency vowel_marks
void LPC_vowel_markPath(int *vowel_mark, int n, int frame_shift, double *Pitch)
{
	int frame_ptr = 0;
	double my_ptr = 0, T;

	for (int i=0; i<n; i++)
	{
		vowel_mark[i] = floor(my_ptr);
		T = (double)SamplingRate/Pitch[frame_ptr];
		my_ptr += T;
		frame_ptr = floor(my_ptr/frame_shift);
	}
}

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

			error = LPC_From_Data(sf_frame, coeff[i][0][j], new_frame_len, order);
		}
		delete [] sf_frame;
		delete [] windows;

		// Vowel part
		vowel_mark_no = LPC_CalPath(WordEnd[i]-WordMid[i], frame_shift, Pitch[i]);
		vowel_mark[i] = new int[vowel_mark_no];
		LPC_vowel_markPath(vowel_mark[i], vowel_mark_no, frame_shift, Pitch[i]);

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

			error = LPC_From_Data(sf_frame, coeff[i][1][j], new_frame_len, order);

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
				LPC_Predict(int_coeff, prime, order, sf_predict, 1);
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

			error = LPC_From_Data(sf_frame, coeff[i][0][j], new_frame_len, order);
		}

		// Vowel part
		frame_no[i][1] = WordEnd[i] - WordMid[i];
		coeff[i][1] = new double *[frame_no[i][1]];

		for (j=0; j<frame_no[i][1]; j++)
		{	
			coeff[i][1][j] = new double[order];

			for (k=0; k<new_frame_len; k++) // read data
				sf_frame[k] = data[frame_shift*(WordMid[i]+j) + k]*windows[k];

			error = LPC_From_Data(sf_frame, coeff[i][1][j], new_frame_len, order);
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
				LPC_Predict(int_coeff, prime, order, sf_predict, 1);
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
				LPC_Predict(int_coeff, prime, order, sf_predict, 1);
				
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
		LPC_SpecEnv(coeff[word_index][1][frame_index-WordMid], SpecEnv_mag, SpecEnv_pha);
	else
		LPC_SpecEnv(coeff[word_index][0][frame_index], SpecEnv_mag, SpecEnv_pha);

	for (i=0; i<FFT_SIZE/2+1; i++)
	{
		mag = my_func.ABS2(spec[0][i], spec[1][i])*SpecEnv_mag[i];
		pha = atan2(spec[1][i], spec[0][i])+SpecEnv_pha[i];

		spec[0][i] = mag*cos(pha);
		spec[1][i] = mag*sin(pha);
	}

	delete [] SpecEnv_mag;
	delete [] SpecEnv_pha;
}

void SourceFilter::LPC_SpecEnv(double *coeff, double *SpecEnv_mag, double *SpecEnv_pha)
{
	double A_real, A_imag;
	int i,j;

	for (i=0; i<FFT_SIZE/2+1; i++)
	{
		A_real = 1;
		A_imag = 0;

		for (j=1; j<=order; j++)
		{
			A_real += coeff[j-1]*cos(i*j*2.0*PI/FFT_SIZE);
			A_imag -= coeff[j-1]*sin(i*j*2.0*PI/FFT_SIZE);
		}

		SpecEnv_mag[i] = 1.0/my_func.ABS2(A_real, A_imag);
		SpecEnv_pha[i] = -1*atan2(A_imag, A_real);
	}
}

/* Input : n elements of time domain data
   Output: m lpc coefficients, excitation energy */
double SourceFilter::LPC_From_Data(double *data,double *lpc,int n,int m)
{
  double *aut= new double[m+1];
  double error;
  int i,j;

  /* autocorrelation, p+1 lag coefficients */

  j=m+1;
  while(j--){
    double d=0;
    for(i=j;i<n;i++)d+=data[i]*data[i-j];
    aut[j]=d;
  }
  
  /* Generate lpc coefficients from autocorr values */

  error=aut[0];
  if(error==0){
    memset(lpc,0,m*sizeof(double));
    return 0;
  }
  
  for(i=0;i<m;i++){
    double r=-aut[i+1];

    /* Sum up this iteration's reflection coefficient; note that in
       Vorbis we don't save it.  If anyone wants to recycle this code
       and needs reflection coefficients, save the results of 'r' from
       each iteration. */

    for(j=0;j<i;j++)r-=lpc[j]*aut[i-j];
    r/=error; 

    /* Update LPC coefficients and total error */
    
    lpc[i]=r;
    for(j=0;j<i/2;j++){
      double tmp=lpc[j];
      lpc[j]+=r*lpc[i-1-j];
      lpc[i-1-j]+=r*tmp;
    }
    if(i%2)lpc[j]+=lpc[j]*r;
    
    error*=1.0-r*r;
  }
  
  /* we need the error value to know how big an impulse to hit the
     filter with later */
  delete [] aut;
  return error;
}

void SourceFilter::LPC_Predict(double *coeff,double *prime,int m,double *data,long n)
{
  /* in: coeff[0...m-1] LPC coefficients 
         prime[0...m-1] initial values (allocated size of n+m-1)
    out: data[0...n-1] data samples */

  long i,j,o,p;
  double y;
  double *work = new double[m+n];

  if(!prime)
    for(i=0;i<m;i++)
      work[i]=0.;
  else
    for(i=0;i<m;i++)
      work[i]=prime[i];

  for(i=0;i<n;i++){
    y=0;
    o=i;
    p=m;
    for(j=0;j<m;j++)
      y-=work[o++]*coeff[--p];
    data[i]=work[o]=y;
  }

  delete [] work;
}
