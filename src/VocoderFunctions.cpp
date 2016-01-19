#include "stdafx.h"
#include "VocoderFunctions.h"

double my_PI = 3.14159265359;

int MIN(int a,int b)
{
	if (a<b)
		return a;
	else
		return b;
}

double MAX(double a, double b)
{
	return (a>b?a:b);
}

double VocoderFunc::ABS2(double a,double b)
{
	return sqrt(a*a+b*b);
}

// Pitch contour tracking and smoothing; reference: ªL°û©É master thesis
void VocoderFunc::pitch_contour_refine(double *pitch, int n)
{
	int i, counter;
	double a;
	double *pitch_temp = new double[n];
	for (i=0; i<n; i++)
		pitch_temp[i] = pitch[i];

	// tracking
	for (i=0; i<n-1; i++)
	{
		if (pitch_temp[i] != 0 && abs(1.0-pitch_temp[i+1]/pitch_temp[i]) > 0.12)
			pitch_temp[i+1] = 0;
	}

	// smoothing
	i = 0;
	while (i < n)
	{
		if (pitch_temp[i] != 0)
		{
			counter = 0;
			while (pitch_temp[i+counter]!=0 && i+counter<n)
				counter++;
	
			if (counter<=5 && i+counter<n)
			{	
				while (counter>0)
				{
					pitch_temp[i] = 0;
					i++;
					counter--;
				}
			}
			else
				i+=counter;
		}
		else
			i++;
	}

	// interpolation
	for (i=0; i<n; i++)
	{
		if (pitch_temp[i] == 0)
		{
			if (i == 0) // if the first frame has no pitch
			{
				while(pitch_temp[i] == 0)
					i++;
				while(i > 0 && i < n)
				{
					pitch_temp[i-1] = pitch_temp[i];
					i--;
				}	
			}	
			else // interpolation
			{
				counter = 0;
				while(i+counter < n && pitch_temp[i+counter] == 0)
					counter++;

				if (i+counter < n)
					a = (pitch_temp[i+counter]-pitch_temp[i-1])/(counter+1);
				else // if the last frame has no pitch
					a = 0;

				while(i < n && pitch_temp[i] == 0)
				{
					pitch_temp[i] = pitch_temp[i-1] + a;
					i++;
				}
			}
		}
	}

	for (i=0; i<n; i++)
		pitch[i] = pitch_temp[i];
	delete [] pitch_temp;
}

// calcuate the target pitch contour for pitch shifting
void VocoderFunc::cal_target_contour(double **t_pitch, int *m, double **s_pitch, int *n, 
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

// Generate Hamming window function
void VocoderFunc::Hamming(int len, double *window)
{
	int i, halflen = floor((double)len/2);
	if (len%2 == 1)
	{
		for (i=0; i<=halflen; i++)
			window[halflen+i] = window[halflen-i] = 0.53836 - 0.46164 * cos(2*my_PI*(halflen-i)/(len-1));
	}
	else
	{
		for (i=0; i<halflen; i++)
			window[halflen+i] = window[halflen-i-1] = 0.53836 - 0.46164 * cos(2*my_PI*(halflen-i)/len);
	}
}

// Calculate time path for time-scale mod.
void CalPath(double *path, double rate, int ElementCount)
{
	for (int i=0; i<ElementCount; i++)
	{
		if (i==0)
			path[i] = 0;
		else
			path[i] = path[i-1] + rate;
	}
}

// Find Max. Frequency Bin
int find_max(double* a, int n, double threshold)
{
	int i, max_bin;
	double max=0;
	for (i=0; i<n; i++)
	{
		if (a[i]>max)
		{
			max = a[i];
			max_bin = i;
		}
	}

	if (max > threshold)
		return max_bin;
	else
		return -1;
}

// ============================================================================ //
// Source-filter Decomposition
//    First we represent the short-term spectrum in polar form,
//    then the spectral magnitude is decomposed into envelope and excitation.
//    The spectral magnitude envelope is obtained by low-passing the spectral magnitude, and
//    the excitation is obtained by directly dividing spectral magnitude by its envelope.
/*void VocoderFunc::SFDecomp(double **Spectrum)
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
void ChannelGrouping(double *Spec, int *ChannelGroupFlag)
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

void new_ChannelGrouping(double *Spec, int *ChannelGroupFlag, double pitch)
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
int VCO_estimation(double *Mag, int *harmonic_peak, double f0)
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
// < Resampling on Frequency Dimension >
// ============================================================================ //
void VocoderFunc::PitchShifting(double **Spectrum, double PitchQuotient, int frame_shift, bool reset_ph)
{
	int i;
	int subband_pitch[FFT_SIZE/2+1];
	int DelFreqBin[FFT_SIZE/2+1];
	double PQResidual[FFT_SIZE/2+1]; // for interpolation
	double energy=0, new_energy=0; // for energy preservation

	double **SpecTemp = new double *[2];
	SpecTemp[0] = new double[FFT_SIZE/2+1];
	SpecTemp[1] = new double[FFT_SIZE/2+1];

	// initialization
	if (reset_ph == true)
	{
		for (i=0; i<FFT_SIZE/2+1; i++)
		{
			phasor_pitch[i] = 0;
			prev_subband_pitch[i] = i;
		}
	}

	for (i=0; i<FFT_SIZE/2+1; i++)
	{
		StudentMag[i] = ABS2(Spectrum[0][i],Spectrum[1][i]);
		StudentPha[i] = atan2(Spectrum[1][i], Spectrum[0][i]);
		energy += StudentMag[i]*StudentMag[i];
	}

	// allocate region-of-influence for sinusoid freq. tracking
	ChannelGrouping(StudentMag, subband_pitch);

	// calculate the phasor and bin number to shift for each peak
	for(i=0;i<FFT_SIZE/2+1;i++)
	{
		if (i==subband_pitch[i])
		{
			DelFreqBin[i] = floor((double)i*PitchQuotient) - i;
			PQResidual[i] = (double)i*PitchQuotient - floor((double)i*PitchQuotient);

			// cumulate the phasors
			phasor_pitch[i] = phasor_pitch[prev_subband_pitch[i]] + (i*(PitchQuotient-1))*(2.0*my_PI*frame_shift)/FFT_SIZE;

			while(phasor_pitch[i] >= my_PI)
				phasor_pitch[i] -= 2.0 * my_PI;
			while(phasor_pitch[i] < -1.0*my_PI)
				phasor_pitch[i] += 2.0 * my_PI;
		}
	}

	// perform phase shift for each region-of-influence
	for(i=0;i<FFT_SIZE/2+1;i++)
	{
		phasor_pitch[i] = phasor_pitch[subband_pitch[i]];
		StudentPha[i] += phasor_pitch[i];

		while(StudentPha[i] >= my_PI)
		   StudentPha[i] -= 2.0 * my_PI;
		while(StudentPha[i] < -1.0*my_PI)
		   StudentPha[i] += 2.0 * my_PI;

		SpecTemp[0][i] = StudentMag[i]*cos(StudentPha[i]);
		SpecTemp[1][i] = StudentMag[i]*sin(StudentPha[i]);

		Spectrum[0][i] = 0;
		Spectrum[1][i] = 0;
	}

	// Shift each group
	for(i=0;i<FFT_SIZE/2+1;i++)	
	{
		PQResidual[i] = PQResidual[subband_pitch[i]];
		DelFreqBin[i] = DelFreqBin[subband_pitch[i]];

		if (i+DelFreqBin[i] >= 0 && i+DelFreqBin[i] < FFT_SIZE/2)
		{
			Spectrum[0][i + DelFreqBin[i]] += (1-PQResidual[i])*SpecTemp[0][i];
			Spectrum[0][i + DelFreqBin[i] + 1] += PQResidual[i]*SpecTemp[0][i];
			Spectrum[1][i + DelFreqBin[i]] += (1-PQResidual[i])*SpecTemp[1][i];
			Spectrum[1][i + DelFreqBin[i] + 1] += PQResidual[i]*SpecTemp[1][i];
		}
	}

	// energy preservation
	for (i=0; i<FFT_SIZE/2+1; i++)
		new_energy += Spectrum[0][i]*Spectrum[0][i] + Spectrum[1][i]*Spectrum[1][i];
	double AmpQuotient = sqrt(energy/new_energy);
	
	for (i=0; i<FFT_SIZE/2+1; i++)
	{
		Spectrum[0][i] *= AmpQuotient;
		Spectrum[1][i] *= AmpQuotient;

		// record previous sub-bands
		prev_subband_pitch[i] = subband_pitch[i];
	}

	delete [] SpecTemp[0];
	delete [] SpecTemp[1];
	delete [] SpecTemp;
}

// ============================================================================ //
// < Pitch Shifting by KTH >
// ============================================================================ //
void VocoderFunc::PitchShifting_KTH(double **Spectrum, double PitchQuotient, int frame_shift, bool reset_ph, double *window, double pitch)
{
	int i, j;
//	int subband_pitch[FFT_SIZE/2+1];
	int DelFreqBin;
	double PQResidual; // for interpolation
	double energy=0, new_energy=0; // for energy preservation
	int max_bin;
	bool *mark_phasor = new bool[FFT_SIZE/2+1];
	double max_amp;
	int iteration;

	double **SpecTemp = new double *[2];
	SpecTemp[0] = new double[FFT_SIZE/2+1];
	SpecTemp[1] = new double[FFT_SIZE/2+1];
	double *MagTemp = new double[FFT_SIZE/2+1];

	// initialization
	for (i=0; i<FFT_SIZE/2+1; i++)
	{
		SpecTemp[0][i] = 0;
		SpecTemp[1][i] = 0;
		mark_phasor[i] = false;
	}

	if (reset_ph == true)
	{
		for (i=0; i<FFT_SIZE/2+1; i++)
		{
			phasor_pitch[i] = 0;
			prev_subband_pitch[i] = i;
		}
	}

	for (i=0; i<FFT_SIZE/2+1; i++)
	{
		StudentMag[i] = ABS2(Spectrum[0][i],Spectrum[1][i]);
		StudentPha[i] = atan2(Spectrum[1][i], Spectrum[0][i]);
		energy += StudentMag[i]*StudentMag[i];
	}

	// iteratively subtract window spectrum from speech spectrum
	iteration = floor(((double)SamplingRate/2.0)/pitch);
	for (i=0; i<iteration; i++)
	{
		max_bin = find_max(StudentMag, FFT_SIZE/2+1, 0);

		// calculate phasor
		DelFreqBin = floor((double)max_bin*PitchQuotient) - max_bin;
		if ((double)max_bin*PitchQuotient-floor((double)max_bin*PitchQuotient) > 0.5)
			DelFreqBin += 1;
		PQResidual = (double)max_bin*PitchQuotient - floor((double)max_bin*PitchQuotient);
		
		if (mark_phasor[max_bin] == false)
		{
			mark_phasor[max_bin] = true;
			phasor_pitch[max_bin] = phasor_pitch[prev_subband_pitch[max_bin]] + (i*(PitchQuotient-1))*(2.0*my_PI*frame_shift)/FFT_SIZE;

			while(phasor_pitch[max_bin] >= my_PI)
				phasor_pitch[max_bin] -= 2.0 * my_PI;
			while(phasor_pitch[max_bin] < -1.0*my_PI)
				phasor_pitch[max_bin] += 2.0 * my_PI;
		}

		// shift the entire window spectrum
		max_amp = StudentMag[max_bin];
		for(j=0; j<FFT_SIZE/2+1; j++)
		{
			if (j-max_bin+FFT_SIZE/2 >= 0 && j-max_bin+FFT_SIZE/2 < FFT_SIZE)
			{
				if (j+DelFreqBin >=0 && j+DelFreqBin < FFT_SIZE/2)
				{
/*					SpecTemp[0][j+DelFreqBin] += max_amp*window[j-max_bin+FFT_SIZE/2]*cos(StudentPha[j]+phasor_pitch[max_bin]);
					SpecTemp[1][j+DelFreqBin] += max_amp*window[j-max_bin+FFT_SIZE/2]*sin(StudentPha[j]+phasor_pitch[max_bin]);
*/
					// interpolation
					SpecTemp[0][j+DelFreqBin] += (1-PQResidual)*max_amp*window[j-max_bin+FFT_SIZE/2]*cos(StudentPha[j]+phasor_pitch[max_bin]);
					SpecTemp[0][j+DelFreqBin+1] += PQResidual*max_amp*window[j-max_bin+FFT_SIZE/2]*cos(StudentPha[j]+phasor_pitch[max_bin]);
					SpecTemp[1][j+DelFreqBin] += (1-PQResidual)*max_amp*window[j-max_bin+FFT_SIZE/2]*sin(StudentPha[j]+phasor_pitch[max_bin]);
					SpecTemp[1][j+DelFreqBin+1] += PQResidual*max_amp*window[j-max_bin+FFT_SIZE/2]*sin(StudentPha[j]+phasor_pitch[max_bin]);
				}

				Spectrum[0][j] -= max_amp*window[j-max_bin+FFT_SIZE/2]*cos(StudentPha[j]);
				Spectrum[1][j] -= max_amp*window[j-max_bin+FFT_SIZE/2]*sin(StudentPha[j]);

				StudentMag[j] -= max_amp*window[j-max_bin+FFT_SIZE/2];
			}
		}
	}

	// generate output spectrum
	for(i=0;i<FFT_SIZE/2+1;i++)
	{
//		Spectrum[0][i] += SpecTemp[0][i];
//		Spectrum[1][i] += SpecTemp[1][i];
		Spectrum[0][i] = SpecTemp[0][i];
		Spectrum[1][i] = SpecTemp[1][i];
	}

	// allocate region-of-influence for sinusoid freq. tracking
	for (i=0; i<FFT_SIZE/2+1; i++)
		MagTemp[i] = ABS2(Spectrum[0][i],Spectrum[1][i]);
	ChannelGrouping(MagTemp, prev_subband_pitch);

	// energy preservation; cumulate the phasors
	for (i=0; i<FFT_SIZE/2+1; i++)
		new_energy += Spectrum[0][i]*Spectrum[0][i] + Spectrum[1][i]*Spectrum[1][i];
	double AmpQuotient = sqrt(energy/new_energy);
	
	for (i=0; i<FFT_SIZE/2+1; i++)
	{
		Spectrum[0][i] *= AmpQuotient;
		Spectrum[1][i] *= AmpQuotient;

		phasor_pitch[i] = phasor_pitch[prev_subband_pitch[i]]; // cumulate the phasors
	}

	delete [] SpecTemp[0];
	delete [] SpecTemp[1];
	delete [] SpecTemp;
	delete [] mark_phasor;
}

// ============================================================================ //
// < Pitch shifting, combining HNM >
// ============================================================================ //

void VocoderFunc::VCO_contour(double ***Spectrum, int n, double *pitch)
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

void VocoderFunc::PitchShifting_KTH_HNM(double **Spectrum, double PitchQuotient, int frame_shift, bool reset_ph, 
										double *window_mag, double *window_pha, double pitch, 
										FILE* myfile1, FILE* myfile2, int vowel_frame_index)
{
	int i, j, window_index, spec_index;
	int DelFreqBin;
	double PQResidual; // for interpolation
	double max_amp;
	int iteration;
	int *peaks = new int[FFT_SIZE/2+1]; // record bin indices of spectral peaks
	int peak_no = 0; // the number of peaks

	int *subband_pitch = new int[FFT_SIZE/2+1];
	double **SpecTemp = new double *[2];
	SpecTemp[0] = new double[FFT_SIZE/2+1];
	SpecTemp[1] = new double[FFT_SIZE/2+1];

	// initialization
	for (i=0; i<FFT_SIZE/2+1; i++)
	{
		SpecTemp[0][i] = 0;
		SpecTemp[1][i] = 0;
	}

	if (reset_ph == true)
	{
		for (i=0; i<FFT_SIZE/2+1; i++)
		{
			phasor_pitch[i] = 0;
			prev_phasor_pitch[i] = 0;
			prev_subband_pitch[i] = i;
		}
	}

	for (i=0; i<FFT_SIZE/2+1; i++)
	{
		StudentMag[i] = ABS2(Spectrum[0][i],Spectrum[1][i]);
		StudentPha[i] = atan2(Spectrum[1][i], Spectrum[0][i]);
	}

	//ChannelGrouping(StudentMag, subband_pitch);
	new_ChannelGrouping(StudentMag, subband_pitch, pitch);
	
for (i=0; i<FFT_SIZE/2+1; i++)
{
	if (i == subband_pitch[i])
		fprintf(myfile1,"%d ",i);
}
fprintf(myfile2,"%d ",voiced_bin_th[vowel_frame_index]);

	// iteratively subtract window spectrum from speech spectrum
	for (i=0; i<voiced_bin_th[vowel_frame_index]; i++)
	{
		if (i == subband_pitch[i])
		{
			// calculate phasor
			DelFreqBin = floor((double)i*PitchQuotient) - i;
			PQResidual = (double)i*PitchQuotient - floor((double)i*PitchQuotient);
			if (PQResidual >= 0.5)
				DelFreqBin += 1;
			
			peaks[peak_no] = i+DelFreqBin;
			peak_no++;

			phasor_pitch[i] = prev_phasor_pitch[prev_subband_pitch[i]] + (i*(PitchQuotient-1))*(2.0*my_PI*frame_shift)/FFT_SIZE;
	
			while(phasor_pitch[i] >= my_PI)
				phasor_pitch[i] -= 2.0 * my_PI;
			while(phasor_pitch[i] < -1.0*my_PI)
				phasor_pitch[i] += 2.0 * my_PI;
			
			// shift the entire window spectrum
			max_amp = StudentMag[i];
			for(j=0; j<FFT_SIZE/2+1; j++)
			{
				window_index = j-i+FFT_SIZE/2;
				spec_index = j+DelFreqBin;
				if (window_index >= 0 && window_index < FFT_SIZE)
				{
					// shift to nearest bin
					if (spec_index >=0 && spec_index < FFT_SIZE/2+1)
					{
						SpecTemp[0][spec_index] += max_amp*window_mag[window_index]*cos(StudentPha[j]+phasor_pitch[i]);
						SpecTemp[1][spec_index] += max_amp*window_mag[window_index]*sin(StudentPha[j]+phasor_pitch[i]);
					}

					// interpolation
/*					if (j+DelFreqBin >=0 && j+DelFreqBin < FFT_SIZE/2)
					{
						SpecTemp[0][spec_index] += (1-PQResidual)*max_amp*window_mag[window_index]*cos(StudentPha[j]+phasor_pitch[i]);
						SpecTemp[0][spec_index+1] += PQResidual*max_amp*window_mag[window_index]*cos(StudentPha[j]+phasor_pitch[i]);
						SpecTemp[1][spec_index] += (1-PQResidual)*max_amp*window_mag[window_index]*sin(StudentPha[j]+phasor_pitch[i]);
						SpecTemp[1][spec_index+1] += PQResidual*max_amp*window_mag[window_index]*sin(StudentPha[j]+phasor_pitch[i]);
					}
*/
					Spectrum[0][j] -= max_amp*window_mag[window_index]*cos(StudentPha[j]);
					Spectrum[1][j] -= max_amp*window_mag[window_index]*sin(StudentPha[j]);
				}
			}
		}
	}

	// high-freq. regeneration
	i=peak_no;
	while (peaks[i-1] < voiced_bin_th[vowel_frame_index])
	{
		peaks[i] = peaks[i-peak_no] + peaks[peak_no-1];

		if (peaks[i] < voiced_bin_th[vowel_frame_index])
		{
			// shift the entire window spectrum
			window_index = j-peaks[i]+FFT_SIZE/2;
			for(j=0; j<FFT_SIZE/2+1; j++)
			{
				if (window_index >= 0 && window_index < FFT_SIZE)
				{
					// use the last available peak mag. (i.e.max_amp) as newly-generated peak mag.
					SpecTemp[0][j] += max_amp*window_mag[j-peaks[i]+FFT_SIZE/2]*cos(StudentPha[j]+phasor_pitch[i]);
					SpecTemp[1][j] += max_amp*window_mag[j-peaks[i]+FFT_SIZE/2]*sin(StudentPha[j]+phasor_pitch[i]);
				}
			}
		}
		i++;
	}

	// generate output spectrum
	for(i = 0; i < voiced_bin_th[vowel_frame_index]; i++)
	{
		Spectrum[0][i] = SpecTemp[0][i];
		Spectrum[1][i] = SpecTemp[1][i];
	}

	// record information about previous spectrum
	for (i=0; i<FFT_SIZE/2+1; i++)
	{
		prev_phasor_pitch[i] = phasor_pitch[subband_pitch[i]];
		prev_subband_pitch[i] = subband_pitch[i];
	}

	delete [] SpecTemp[0];
	delete [] SpecTemp[1];
	delete [] SpecTemp;
	delete [] subband_pitch;
	delete [] peaks;
}

// ============================================================================ //
// < Fixed KTH approach >
// ============================================================================ //

double *MUL2(double r1, double i1, double r2, double i2) //multiplation of complex numbers
{
	double *ans = new double[2];

	ans[0] = r1*r2-i1*i2;
	ans[1] = r1*i2+r2*i1;

	return ans;
}

void VocoderFunc::PS_KTH_new(double **Spectrum, double PitchQuotient, int frame_shift, bool reset_ph, 
							 double *window_real, double *window_imag)
{
	int i, j;
	int DelFreqBin;
	double PQResidual; // for interpolation
	double energy=0, new_energy=0; // for energy preservation
	int max_bin;
	bool *mark_phasor = new bool[FFT_SIZE/2+1];
	double max_amp_real, max_amp_imag, threshold;
//	double *ROI, *ROI_shift; // region of influence

	double **SpecTemp = new double *[2];
	SpecTemp[0] = new double[FFT_SIZE/2+1];
	SpecTemp[1] = new double[FFT_SIZE/2+1];
	double *MagTemp = new double[FFT_SIZE/2+1];
//AfxMessageBox("start");
	// initialization
	for (i=0; i<FFT_SIZE/2+1; i++)
	{
		SpecTemp[0][i] = 0;
		SpecTemp[1][i] = 0;
		mark_phasor[i] = false;
	}

	if (reset_ph == true)
	{
		for (i=0; i<FFT_SIZE/2+1; i++)
		{
			phasor_pitch[i] = 0;
			prev_subband_pitch[i] = i;
		}
	}

	for (i=0; i<FFT_SIZE/2+1; i++)
	{
		StudentMag[i] = ABS2(Spectrum[0][i],Spectrum[1][i]);
		StudentPha[i] = atan2(Spectrum[1][i], Spectrum[0][i]);
		energy += StudentMag[i]*StudentMag[i];
	}
//FILE* myfile = fopen("ROI.txt","w");
	// iteratively subtract window spectrum from speech spectrum
	max_bin = find_max(StudentMag, FFT_SIZE/2+1, 0);
	threshold = 0.007*StudentMag[max_bin];
	while (max_bin != -1)
	{
		// calculate phasor
		DelFreqBin = floor((double)max_bin*PitchQuotient) - max_bin;
//		if ((double)max_bin*PitchQuotient-floor((double)max_bin*PitchQuotient) > 0.5)
//			DelFreqBin += 1;
		PQResidual = (double)max_bin*PitchQuotient - floor((double)max_bin*PitchQuotient);
		
		if (mark_phasor[max_bin] == false)
		{
			mark_phasor[max_bin] = true;
			phasor_pitch[max_bin] = phasor_pitch[prev_subband_pitch[max_bin]] + (max_bin*(PitchQuotient-1))*(2.0*my_PI*frame_shift)/FFT_SIZE;

			while(phasor_pitch[max_bin] >= my_PI)
				phasor_pitch[max_bin] -= 2.0 * my_PI;
			while(phasor_pitch[max_bin] < -1.0*my_PI)
				phasor_pitch[max_bin] += 2.0 * my_PI;
		}

		// shift the entire window spectrum
		max_amp_real = Spectrum[0][max_bin];
		max_amp_imag = Spectrum[1][max_bin];
		for(j=0; j<FFT_SIZE/2+1; j++)
		{
			if (j-max_bin+FFT_SIZE/2 >= 0 && j-max_bin+FFT_SIZE/2 < FFT_SIZE) // subtract window spectrum from original spectrum
			{
				double *ROI = MUL2(max_amp_real, max_amp_imag, window_real[j-max_bin+FFT_SIZE/2], window_imag[j-max_bin+FFT_SIZE/2]);
//fprintf(myfile,"%f + %f*i\n",ROI[0],ROI[1]);
				if (j+DelFreqBin >=0 && j+DelFreqBin < FFT_SIZE/2) // shift the entire window spectrum
				{
					double *ROI_shift = MUL2(ROI[0], ROI[1], cos(phasor_pitch[max_bin]), sin(phasor_pitch[max_bin]));

/*					SpecTemp[0][j+DelFreqBin] += ROI_real_shift;
					SpecTemp[1][j+DelFreqBin] += ROI_imag_shift;*/

					// interpolation
					SpecTemp[0][j+DelFreqBin] += (1-PQResidual)*ROI_shift[0];
					SpecTemp[0][j+DelFreqBin+1] += PQResidual*ROI_shift[0];
					SpecTemp[1][j+DelFreqBin] += (1-PQResidual)*ROI_shift[1];
					SpecTemp[1][j+DelFreqBin+1] += PQResidual*ROI_shift[1];

					delete [] ROI_shift;
				}

				Spectrum[0][j] -= ROI[0];
				Spectrum[1][j] -= ROI[1];
				StudentMag[j] -= ABS2(ROI[0], ROI[1]);

				delete [] ROI;
			}
		}
		max_bin = find_max(StudentMag, FFT_SIZE/2+1, threshold);
	}

//FILE* myfile = fopen("output_spec.txt","w");
	// generate output spectrum
	for(i=0;i<FFT_SIZE/2+1;i++)
	{
		Spectrum[0][i] = SpecTemp[0][i];
		Spectrum[1][i] = SpecTemp[1][i];
//fprintf(myfile,"%f\n",ABS2(Spectrum[0][i],Spectrum[1][i]));
	}
//fclose(myfile);
//AfxMessageBox("end");
	// allocate region-of-influence for sinusoid freq. tracking
	for (i=0; i<FFT_SIZE/2+1; i++)
		MagTemp[i] = ABS2(Spectrum[0][i],Spectrum[1][i]);
	ChannelGrouping(MagTemp, prev_subband_pitch);

	// energy preservation; cumulate the phasors
	for (i=0; i<FFT_SIZE/2+1; i++)
		new_energy += Spectrum[0][i]*Spectrum[0][i] + Spectrum[1][i]*Spectrum[1][i];
	double AmpQuotient = sqrt(energy/new_energy);
	
	for (i=0; i<FFT_SIZE/2+1; i++)
	{
		Spectrum[0][i] *= AmpQuotient;
		Spectrum[1][i] *= AmpQuotient;

		phasor_pitch[i] = phasor_pitch[prev_subband_pitch[i]]; // cumulate the phasors
	}

	delete [] SpecTemp[0];
	delete [] SpecTemp[1];
	delete [] SpecTemp;
	delete [] mark_phasor;
}

// ============================================================================ //
// < Resampling on time axis >
// ============================================================================ //
void VocoderFunc::TimeScaleMod(double rate, int FrameNum, double ***spectrogram, 
							   int new_FrameNum, double ***new_spectrogram, bool reset_ph)
{
	int i,j,f;
	double a,b;
	double mag[FFT_SIZE/2+1], ph[FFT_SIZE/2+1];
	int subband_time[FFT_SIZE/2+1];
	double *time_path = new double[new_FrameNum];
	CalPath(time_path, rate, new_FrameNum);

	// initialize the first frame
	for(j=0;j<FFT_SIZE/2+1;j++)
		mag[j] = ABS2(spectrogram[0][0][j],spectrogram[0][1][j]);

	if (reset_ph == true)
	{
		for(j=0;j<FFT_SIZE/2+1;j++)
		{
			ph[j] = atan2(spectrogram[0][1][j],spectrogram[0][0][j]);
			new_spectrogram[0][0][j] = spectrogram[0][0][j];
			new_spectrogram[0][1][j] = spectrogram[0][1][j];
		}
		ChannelGrouping(mag, prev_subband_time);
	}
	else
	{
		ChannelGrouping(mag, subband_time);

		// calculate the phase of each peak
		for(j=0;j<FFT_SIZE/2+1;j++)
		{
			// for synchronizing the phases around the peaks
			// scaled PL
			if (j==subband_time[j])
				phasor_time[j] = atan2(spectrogram[0][1][j],spectrogram[0][0][j]) 
								- prev_phase[prev_subband_time[j]];
		}

		// calculate the phase of the other channels, and construct synthesis frame
		for(j=0;j<FFT_SIZE/2+1;j++)
		{
			ph[j] = prev_phase[j] + phasor_time[subband_time[j]];
			
			while(ph[j] >= my_PI)
				ph[j] -= 2.0 * my_PI;
			while(ph[j] < -1.0*my_PI)
				ph[j] += 2.0 * my_PI;
			
			new_spectrogram[0][0][j] = mag[j] * cos(ph[j]);
			new_spectrogram[0][1][j] = mag[j] * sin(ph[j]);

			// record previous sub-bands for scaled PL
			prev_subband_time[j] = subband_time[j];
		}
	}

	for (i=1; i<new_FrameNum; i++) 
	{
		f = (int)time_path[i];
		if (f >= FrameNum-1)
			f = FrameNum-2;
		b = time_path[i] - f;
		a = 1 - b;
		for(j=0;j<FFT_SIZE/2+1;j++)
			mag[j] = a*ABS2(spectrogram[f][0][j],spectrogram[f][1][j]) 
					+ b*ABS2(spectrogram[f+1][0][j],spectrogram[f+1][1][j]);

		ChannelGrouping(mag, subband_time);

		// calculate the phase of each peak
		for(j=0;j<FFT_SIZE/2+1;j++)
		{
			if (j==subband_time[j])
			{
				// for synchronizing the phases around the peaks
				// scaled PL
				phasor_time[j] = atan2(spectrogram[f+1][1][j],spectrogram[f+1][0][j]) 
							- atan2(spectrogram[f][1][prev_subband_time[j]],spectrogram[f][0][prev_subband_time[j]]);

				ph[j] += phasor_time[j];
				while(ph[j] >= my_PI)
					ph[j] -= 2 * my_PI;
				while(ph[j] < -1.0*my_PI)
					ph[j] += 2 * my_PI;
			}
		}

		// calculate the phase of the other channels, and construct synthesis frame
		for(j=0;j<FFT_SIZE/2+1;j++)
		{
			if (j!=subband_time[j])
			{
				ph[j] += phasor_time[subband_time[j]];
			
				while(ph[j] >= my_PI)
					ph[j] -= 2.0 * my_PI;
				while(ph[j] < -1.0*my_PI)
					ph[j] += 2.0 * my_PI;
			}

			new_spectrogram[i][0][j] = mag[j] * cos(ph[j]);
			new_spectrogram[i][1][j] = mag[j] * sin(ph[j]);

			// record previous sub-bands for scaled PL
			prev_subband_time[j] = subband_time[j];
		}
	}
	
	for(j=0;j<FFT_SIZE/2+1;j++)
		prev_phase[j] = ph[j];

	delete [] time_path;
}

// ============================================================================ //
//  < Pitch-Sync. Overlap-Add >
// ============================================================================ //
int VocoderFunc::PSOLA(double **frame, double *output, double *SentenceCoeff, double *windows,
					  int sample_index, double *Pitch, int WordStart, int WordMid, int WordEnd,
					  int frame_len, int frame_shift, int bias)
{
	int i,j,k;
	int syn_bias; // a shift before overlap-add for pitch-sync
	int total_syn_bias = 0;
	int syn_sample_no;
	double *frame_for_syn;
	double cross_corr, max_cross_corr;
	double T; // peroid in sample number

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
		T = floor((double)SamplingRate/Pitch[i-(WordMid-WordStart)]);
		syn_sample_no = bias;
		frame_for_syn = new double[syn_sample_no];

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
int VocoderFunc::non_PSOLA(double **frame, double *output, double *SentenceCoeff, double *windows,
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