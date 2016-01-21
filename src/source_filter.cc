#include "source_filter.h"

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
//  < Fixed frame length LPA  >
// ============================================================================ //
void SourceFilter::Decomposition(float *data,float *excitation, int sample_no,
									int word_no, int* WordStart, int* WordMid, int* WordEnd, 
									int frame_shift, int frame_len)
{
	// setting
	int new_frame_len = frame_len;

	int i, j, k, m;
	float error, *sf_frame;
	vowel_mark = new int *[word_no];
	frame_no = new int *[word_no];
	coeff = new float ***[word_no];
    window_function = new Window(new_frame_len);
	float *windows = window_function->getHamming();

	// Analysis: Frame Length = 20ms; Frame Shift = 10ms; Order = 20
	for (i=0; i<word_no; i++)
	{
		frame_no[i] = new int [2];
		coeff[i] = new float **[2];
		sf_frame = new float[new_frame_len];

		// Consonant part
		frame_no[i][0] = WordMid[i]-WordStart[i];
		coeff[i][0] = new float *[frame_no[i][0]];

		for (j=0; j<frame_no[i][0]; j++)
		{
			coeff[i][0][j] = new float[order];

			for (k=0; k<new_frame_len; k++) // read data
				sf_frame[k] = data[frame_shift*(WordStart[i]+j) + k]*windows[k];

			error = lpc.From_Data(sf_frame, coeff[i][0][j], new_frame_len, order);
		}

		// Vowel part
		frame_no[i][1] = WordEnd[i] - WordMid[i];
		coeff[i][1] = new float *[frame_no[i][1]];

		for (j=0; j<frame_no[i][1]; j++)
		{	
			coeff[i][1][j] = new float[order];

			for (k=0; k<new_frame_len; k++) // read data
				sf_frame[k] = data[frame_shift*(WordMid[i]+j) + k]*windows[k];

			error = lpc.From_Data(sf_frame, coeff[i][1][j], new_frame_len, order);
		}
		delete [] sf_frame;
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
//  < Source-Filter Recombination at Frequency Domain >
// ============================================================================ //
void SourceFilter::Recombination(float** spec, int word_index, int frame_index, int WordMid)
{
	float *SpecEnv_mag = new float[FFT_SIZE/2+1];
	float *SpecEnv_pha = new float[FFT_SIZE/2+1];
	float mag, pha;
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
