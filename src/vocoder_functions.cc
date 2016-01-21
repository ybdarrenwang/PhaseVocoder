#include "vocoder_functions.h"

double VocoderFunctions::ABS2(double a,double b) {return sqrt(a*a+b*b);}

// calcuate the target pitch contour for pitch shifting
void VocoderFunctions::cal_target_contour(double **t_pitch, int *m, double **s_pitch, int *n, 
									 double **tar_pitch, int word_no)
{
	int i, j, PitchIndex;
	double PitchDurMul, b, int_pitch, pitch_quotient;
	double ave_t_pitch = 0, ave_s_pitch = 0;
	double tar_min_pitch, tar_max_pitch, s_min_pitch, s_max_pitch;
	double range_quot;
	int t_counter = 0, s_counter = 0;

	// calculate average log-f0
	for (i=0; i<word_no; i++)
	{
		for (j=0; j<m[i]; j++)
		{
			ave_t_pitch += log(t_pitch[i][j]);
			t_counter++;
		}
		for (j=0; j<n[i]; j++)
		{
			ave_s_pitch += log(s_pitch[i][j]);
			s_counter++;
		}
	}
	ave_t_pitch = exp(ave_t_pitch/t_counter);
	ave_s_pitch = exp(ave_s_pitch/s_counter);
	pitch_quotient = ave_s_pitch/ave_t_pitch;

	// calculate target pitch contour (minimize log-f0 difference)
	for (i=0; i<word_no; i++)
	{
		PitchDurMul = (double)m[i]/n[i];
		for (j=0; j<n[i]; j++)
		{			
			PitchIndex = floor(j*PitchDurMul);
			if(PitchIndex < m[i]-1)
				b = j*PitchDurMul - PitchIndex;
			else
				b = 0;
	
			int_pitch = (1-b)*t_pitch[i][PitchIndex] + b*t_pitch[i][PitchIndex+1];

			tar_pitch[i][j] = int_pitch*pitch_quotient;
		}
	}

	// find min. and max. pitch for normalization
/*	tar_min_pitch = tar_max_pitch = tar_pitch[0][0];
	s_min_pitch = s_max_pitch = s_pitch[0][0];
	for (i=0; i<word_no; i++)
	{
		for (j=0; j<n[i]; j++)
		{
			if (s_pitch[i][j] > s_max_pitch)
				s_max_pitch = s_pitch[i][j];

			if (s_pitch[i][j] < s_min_pitch)
				s_min_pitch = s_pitch[i][j];

			if (tar_pitch[i][j] > tar_max_pitch)
				tar_max_pitch = tar_pitch[i][j];

			if (tar_pitch[i][j] < tar_min_pitch)
				tar_min_pitch = tar_pitch[i][j];
		}
	}

	range_quot = MAX((tar_max_pitch-ave_s_pitch)/(s_max_pitch-ave_s_pitch), (tar_min_pitch-ave_s_pitch)/(s_min_pitch-ave_s_pitch));

	for (i=0; i<word_no; i++)
	{
		for (j=0; j<n[i]; j++)
			tar_pitch[i][j] = (tar_pitch[i][j]-ave_s_pitch)/range_quot + ave_s_pitch;
	}*/
}

// ============================================================================ //
// Source-filter Decomposition
//    First we represent the short-term spectrum in polar form,
//    then the spectral magnitude is decomposed into envelope and excitation.
//    The spectral magnitude envelope is obtained by low-passing the spectral magnitude, and
//    the excitation is obtained by directly dividing spectral magnitude by its envelope.
/*void VocoderFunctions::SFDecomp(double **Spectrum)
{
	int i,j,FreqHalfWinSize, FreqWinSize;

	for (i=0; i<FFT_SIZE/2+1; i++)
	{
		StudentMag[i] = ABS2(Spectrum[0][i],Spectrum[1][i]);
		StudentPha[i] = atan2(Spectrum[1][i], Spectrum[0][i]);
		StudentMagEnv[i] = 0;
	}

	for (i=0; i<FFT_SIZE/2+1; i++)
	{
		FreqHalfWinSize = MIN(4,abs(FFT_SIZE/4-abs(FFT_SIZE/4-i)));
		FreqWinSize = FreqHalfWinSize*2+1;

		for (j=0; j<FreqWinSize; j++)
			StudentMagEnv[i] += StudentMag[i+j-FreqHalfWinSize] / FreqWinSize;
	}
	
	for (i=0; i<FFT_SIZE/2+1; i++)
		StudentMagExc[i] = StudentMag[i] / StudentMagEnv[i];
}*/

// ============================================================================ //
// < Group the channels accoring to the local maxima and minima in each frame >
// ============================================================================ //
void VocoderFunctions::ChannelGrouping(double *Spec, int *ChannelGroupFlag)
{
	int i,j;
	int GroupBoundary, PrevGroupBoundary=0, PrevPeak=0;
	double MinTemp;

	for(i=2; i<FFT_SIZE/2-1; i++)
	{
		// find local mixima
		if (Spec[i] > Spec[i-1] && Spec[i] > Spec[i-2] && Spec[i] > Spec[i+1] && Spec[i] > Spec[i+2])
		{
			if (PrevPeak > 0)
			{
				// find minima between present and previous peak
				MinTemp = Spec[PrevPeak];
				for (j=PrevPeak+1; j<i; j++)
				{
					if (Spec[j] < MinTemp)
					{
						MinTemp = Spec[j];
						GroupBoundary = j;
					}
				}

				// annotate the previous group around the previous peak
				for (j = PrevGroupBoundary; j<GroupBoundary; j++)
					ChannelGroupFlag[j] = PrevPeak;

				PrevGroupBoundary = GroupBoundary;
			}
			else // annotate the first group
			{
				for (j=0; j<i; j++)
					ChannelGroupFlag[j] = i;
			}
			PrevPeak = i;
		}
	}

	// annotate the last group
	for (j = PrevGroupBoundary; j<FFT_SIZE/2+1; j++)
		ChannelGroupFlag[j] = PrevPeak;
}

void VocoderFunctions::new_ChannelGrouping(double *Spec, int *ChannelGroupFlag, double pitch)
{
	int i,j;
	int GroupBoundary, PrevGroupBoundary=0, PrevPeak=0;
	double MinTemp;

	// determine how much bins to look forward/backward
	double bin_period = FFT_SIZE/(SamplingRate/pitch); // f0 interval width counted by bin number
	int find_max_range = (int)floor(bin_period*0.5);
	bool find_max;

	// determine region of influence
	for(i=find_max_range; i<FFT_SIZE/2+1-find_max_range; i++)
	{
		// find local mixima
		find_max = true;
		for (j=i-find_max_range; j<i+find_max_range; j++)
		{
			if (Spec[i]<Spec[j])
				find_max = false;
		}

		// determine region of influence
		if (find_max == true)
		{
			if (PrevPeak > 0)
			{
				// find minima between present and previous peak
				MinTemp = Spec[PrevPeak];
				for (j=PrevPeak+1; j<i; j++)
				{
					if (Spec[j] < MinTemp)
					{
						MinTemp = Spec[j];
						GroupBoundary = j;
					}
				}

				// annotate the previous group around the previous peak
				for (j = PrevGroupBoundary; j<GroupBoundary; j++)
					ChannelGroupFlag[j] = PrevPeak;

				PrevGroupBoundary = GroupBoundary;
			}
			else // annotate the first group
			{
				for (j=0; j<i; j++)
					ChannelGroupFlag[j] = i;
			}
			PrevPeak = i;
		}
	}

	// annotate the last group
	for (j = PrevGroupBoundary; j<FFT_SIZE/2+1; j++)
		ChannelGroupFlag[j] = PrevPeak;
}

// ============================================================================ //
// < Voice cut-off frequency estimation >
// Reference: Applying the Harmonic Plus Noise Model in Concatenative Speech Synthesis, page3
// ============================================================================ //
int VocoderFunctions::VCO_estimation(double *Mag, int *harmonic_peak, double f0)
{
	// parameter setting
	double VCO_lb = 1000, VCO_ub = 8000; // the upper/lower bound of VCO
	int n = 3; // if n consecutive harmonics are unvoiced, then set VCO

	int i, j, harmonic_no=0, bin_index, VCO;
	int ChannelGroup[FFT_SIZE/2+1];
	double max_minor_mag, central_freq, residual;

	// record the bin index of harmonic bins
	for (i=0; i<FFT_SIZE/2+1; i++)
	{
		if (i == harmonic_peak[i])
			harmonic_no++;
	}
	
	int *harmonic_bin = new int[harmonic_no];
	bool *harmonic_voiced = new bool[harmonic_no];
	for (i=0; i<harmonic_no; i++)
		harmonic_voiced[i] = 0;

	j = 0;
	for (i=0; i<FFT_SIZE/2+1; i++)
	{
		if (i == harmonic_peak[i])
		{
			harmonic_bin[j] = i;
			j++;
		}
	}

	// determine how much bins to look forward/backward
	double bin_period = FFT_SIZE/(SamplingRate/f0); // f0 interval width counted by bin number
	int find_max_range = (int)floor(bin_period*0.5);

	// determine voiced/unvoiced harmonic peaks using peakiness
	ChannelGrouping(Mag, ChannelGroup); // locate minor peaks
	for (i=0; i<harmonic_no; i++)
	{
		// find max{A(omega-i)}
		max_minor_mag = 0;
		for (j=-1*find_max_range; j<=find_max_range; j++)
		{
			bin_index = harmonic_bin[i]+j;
			if (bin_index == ChannelGroup[bin_index] && j!=0)
			{
				if (Mag[bin_index] > max_minor_mag)
					max_minor_mag = Mag[bin_index];
			}
		}
		
		if (Mag[harmonic_bin[i]]/max_minor_mag > 4 // A(omega-c)-max(A(omega-i)) > 6dB
			|| max_minor_mag == 0) // there's no other peaks
			harmonic_voiced[i] = 1;
	}
/*	FILE* myfile=fopen("harmonic-1.txt","w");
	for (i=0; i<harmonic_no; i++)
		fprintf(myfile,"%d, %d\n",harmonic_bin[i], harmonic_voiced[i]);
	fclose(myfile);*/

	// determine voiced/unvoiced harmonic peaks using the location of harmonics
//	myfile=fopen("harmonic-2.txt","w");
	for (i=0; i<harmonic_no; i++)
	{
		central_freq = SamplingRate*((double)harmonic_bin[i]/FFT_SIZE);
		residual = central_freq/f0 - floor(central_freq/f0);
		if (residual > 0.5)
			residual = 1-residual;

		if (residual < 0.1)
		{
			harmonic_voiced[i] = 1;
//			fprintf(myfile,"%d, %d\n",harmonic_bin[i], 1);
		}
//		else
//			fprintf(myfile,"%d, %d\n",harmonic_bin[i], 0);
	}
//	fclose(myfile);

/*	myfile=fopen("harmonic.txt","w");
	for (i=0; i<harmonic_no; i++)
		fprintf(myfile,"%d, %d\n",harmonic_bin[i], harmonic_voiced[i]);
	fclose(myfile);*/

	// set the harmonics beneath lower bound be voiced
/*	i = 0;
	while(harmonic_bin[i] < (VCO_lb/SamplingRate)*FFT_SIZE)
	{
		harmonic_voiced[i] = 1;
		i++;
	}
	// set the harmonics exceed upper bound be unvoiced
	i = harmonic_no-1;
	while(harmonic_bin[i] > (VCO_ub/SamplingRate)*FFT_SIZE)
	{
		harmonic_voiced[i] = 0;
		i--;
	}*/

	// if n consecutive harmonics are unvoiced, then set VCO
	i = 0;
	while(harmonic_voiced[i] == 1 && i+n-1<harmonic_no)
	{
		i++;
		for (j=0; j<n; j++)
		{
			if (harmonic_voiced[i+j] == 1)
				harmonic_voiced[i] = 1;
		}
	}
//	VCO = SamplingRate*(double)harmonic_bin[i-1]/FFT_SIZE;
	VCO = harmonic_bin[i-1];
	if (VCO < (VCO_lb/SamplingRate)*FFT_SIZE || // lower bound
		VCO > (VCO_ub/SamplingRate)*FFT_SIZE) // upper bound
		VCO = 0;

//AfxMessageBox("VCO");

	delete [] harmonic_bin;
	delete [] harmonic_voiced;

	return VCO;
}

// ============================================================================ //
// < Pitch shifting, combining HNM >
// ============================================================================ //

void VocoderFunctions::VCO_contour(double ***Spectrum, int n, double *pitch)
{
	int i, j, counter;
	double a,average_bin_th = 0;
	int *subband_pitch;

	voiced_bin_th = new int[n];

	// calculate VCO
	for (i=0; i<n; i++)
	{
		for (j=0; j<FFT_SIZE/2+1; j++)
			StudentMag[j] = ABS2(Spectrum[i][0][j],Spectrum[i][1][j]);

		subband_pitch = new int[FFT_SIZE/2+1];
		new_ChannelGrouping(StudentMag, subband_pitch, pitch[i]);
		voiced_bin_th[i] = VCO_estimation(StudentMag, subband_pitch, pitch[i]);
		delete [] subband_pitch;
	}

	// VCO contour smoothing
/*	for (i=0; i<n; i++)
	{
		if (voiced_bin_th[i] > 0)
		{
			average_bin_th += voiced_bin_th[i];
			counter++;
		}
	}
	average_VCO/=counter;*/
	
	// interpolation
	for (i=0; i<n; i++)
	{
		if (voiced_bin_th[i] == 0)
		{
			if (i == 0) // if the first voiced_bin_th has no value
			{
				while(voiced_bin_th[i] == 0)
					i++;
				while(i > 0 && i < n)
				{
					voiced_bin_th[i-1] = voiced_bin_th[i];
					i--;
				}	
			}	
			else // interpolation
			{
				counter = 0;
				while(i+counter < n && voiced_bin_th[i+counter] == 0)
					counter++;

				if (i+counter < n)
					a = (double)(voiced_bin_th[i+counter]-voiced_bin_th[i-1])/(counter+1);
				else // if the last frame has no pitch
					a = 0;

				while(i < n && voiced_bin_th[i] == 0)
				{
					voiced_bin_th[i] = (int)floor(voiced_bin_th[i-1] + a);
					i++;
				}
			}
		}
	}
}

