#include "pitch_shifter.h"

float PitchShifter::ABS2(float a,float b) {return sqrt(a*a+b*b);}

// Find Max. Frequency Bin
int find_max(float* a, int n, float threshold)
{
	int i, max_bin;
	float max=0;
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
// < Resampling on Frequency Dimension >
// ============================================================================ //
void PitchShifter::PitchShifting(float **Spectrum, float PitchQuotient, int frame_shift, bool reset_ph)
{
	int i;
	int subband_pitch[FFT_SIZE/2+1];
	int DelFreqBin[FFT_SIZE/2+1];
	float PQResidual[FFT_SIZE/2+1]; // for interpolation
	float energy=0, new_energy=0; // for energy preservation

	float **SpecTemp = new float *[2];
	SpecTemp[0] = new float[FFT_SIZE/2+1];
	SpecTemp[1] = new float[FFT_SIZE/2+1];

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
			DelFreqBin[i] = floor((float)i*PitchQuotient) - i;
			PQResidual[i] = (float)i*PitchQuotient - floor((float)i*PitchQuotient);

			// cumulate the phasors
			phasor_pitch[i] = phasor_pitch[prev_subband_pitch[i]] + (i*(PitchQuotient-1))*(2.0*PI*frame_shift)/FFT_SIZE;

			while(phasor_pitch[i] >= PI)
				phasor_pitch[i] -= 2.0 * PI;
			while(phasor_pitch[i] < -1.0*PI)
				phasor_pitch[i] += 2.0 * PI;
		}
	}

	// perform phase shift for each region-of-influence
	for(i=0;i<FFT_SIZE/2+1;i++)
	{
		phasor_pitch[i] = phasor_pitch[subband_pitch[i]];
		StudentPha[i] += phasor_pitch[i];

		while(StudentPha[i] >= PI)
		   StudentPha[i] -= 2.0 * PI;
		while(StudentPha[i] < -1.0*PI)
		   StudentPha[i] += 2.0 * PI;

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
	float AmpQuotient = sqrt(energy/new_energy);
	
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
void PitchShifter::PitchShifting_KTH(float **Spectrum, float PitchQuotient, int frame_shift, bool reset_ph, float *window, float pitch)
{
	int i, j;
//	int subband_pitch[FFT_SIZE/2+1];
	int DelFreqBin;
	float PQResidual; // for interpolation
	float energy=0, new_energy=0; // for energy preservation
	int max_bin;
	bool *mark_phasor = new bool[FFT_SIZE/2+1];
	float max_amp;
	int iteration;

	float **SpecTemp = new float *[2];
	SpecTemp[0] = new float[FFT_SIZE/2+1];
	SpecTemp[1] = new float[FFT_SIZE/2+1];
	float *MagTemp = new float[FFT_SIZE/2+1];

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
	iteration = floor(((float)SamplingRate/2.0)/pitch);
	for (i=0; i<iteration; i++)
	{
		max_bin = find_max(StudentMag, FFT_SIZE/2+1, 0);

		// calculate phasor
		DelFreqBin = floor((float)max_bin*PitchQuotient) - max_bin;
		if ((float)max_bin*PitchQuotient-floor((float)max_bin*PitchQuotient) > 0.5)
			DelFreqBin += 1;
		PQResidual = (float)max_bin*PitchQuotient - floor((float)max_bin*PitchQuotient);
		
		if (mark_phasor[max_bin] == false)
		{
			mark_phasor[max_bin] = true;
			phasor_pitch[max_bin] = phasor_pitch[prev_subband_pitch[max_bin]] + (i*(PitchQuotient-1))*(2.0*PI*frame_shift)/FFT_SIZE;

			while(phasor_pitch[max_bin] >= PI)
				phasor_pitch[max_bin] -= 2.0 * PI;
			while(phasor_pitch[max_bin] < -1.0*PI)
				phasor_pitch[max_bin] += 2.0 * PI;
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
	float AmpQuotient = sqrt(energy/new_energy);
	
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

void PitchShifter::PitchShifting_KTH_HNM(float **Spectrum, float PitchQuotient, int frame_shift, bool reset_ph, 
										float *window_mag, float *window_pha, float pitch, int vowel_frame_index)
{
	int i, j, window_index, spec_index;
	int DelFreqBin;
	float PQResidual; // for interpolation
	float max_amp;
	int iteration;
	int *peaks = new int[FFT_SIZE/2+1]; // record bin indices of spectral peaks
	int peak_no = 0; // the number of peaks

	int *subband_pitch = new int[FFT_SIZE/2+1];
	float **SpecTemp = new float *[2];
	SpecTemp[0] = new float[FFT_SIZE/2+1];
	SpecTemp[1] = new float[FFT_SIZE/2+1];

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

	ChannelGrouping(StudentMag, subband_pitch);
	//new_ChannelGrouping(StudentMag, subband_pitch, pitch);
	
	// iteratively subtract window spectrum from speech spectrum
	for (i=0; i<voiced_bin_th[vowel_frame_index]; i++)
	{
		if (i == subband_pitch[i])
		{
			// calculate phasor
			DelFreqBin = floor((float)i*PitchQuotient) - i;
			PQResidual = (float)i*PitchQuotient - floor((float)i*PitchQuotient);
			if (PQResidual >= 0.5)
				DelFreqBin += 1;
			
			peaks[peak_no] = i+DelFreqBin;
			peak_no++;

			phasor_pitch[i] = prev_phasor_pitch[prev_subband_pitch[i]] + (i*(PitchQuotient-1))*(2.0*PI*frame_shift)/FFT_SIZE;
	
			while(phasor_pitch[i] >= PI)
				phasor_pitch[i] -= 2.0 * PI;
			while(phasor_pitch[i] < -1.0*PI)
				phasor_pitch[i] += 2.0 * PI;
			
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

float *MUL2(float r1, float i1, float r2, float i2) //multiplation of complex numbers
{
	float *ans = new float[2];

	ans[0] = r1*r2-i1*i2;
	ans[1] = r1*i2+r2*i1;

	return ans;
}

void PitchShifter::PS_KTH_new(float **Spectrum, float PitchQuotient, int frame_shift, bool reset_ph, 
							 float *window_real, float *window_imag)
{
	int i, j;
	int DelFreqBin;
	float PQResidual; // for interpolation
	float energy=0, new_energy=0; // for energy preservation
	int max_bin;
	bool *mark_phasor = new bool[FFT_SIZE/2+1];
	float max_amp_real, max_amp_imag, threshold;
//	float *ROI, *ROI_shift; // region of influence

	float **SpecTemp = new float *[2];
	SpecTemp[0] = new float[FFT_SIZE/2+1];
	SpecTemp[1] = new float[FFT_SIZE/2+1];
	float *MagTemp = new float[FFT_SIZE/2+1];
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
		DelFreqBin = floor((float)max_bin*PitchQuotient) - max_bin;
//		if ((float)max_bin*PitchQuotient-floor((float)max_bin*PitchQuotient) > 0.5)
//			DelFreqBin += 1;
		PQResidual = (float)max_bin*PitchQuotient - floor((float)max_bin*PitchQuotient);
		
		if (mark_phasor[max_bin] == false)
		{
			mark_phasor[max_bin] = true;
			phasor_pitch[max_bin] = phasor_pitch[prev_subband_pitch[max_bin]] + (max_bin*(PitchQuotient-1))*(2.0*PI*frame_shift)/FFT_SIZE;

			while(phasor_pitch[max_bin] >= PI)
				phasor_pitch[max_bin] -= 2.0 * PI;
			while(phasor_pitch[max_bin] < -1.0*PI)
				phasor_pitch[max_bin] += 2.0 * PI;
		}

		// shift the entire window spectrum
		max_amp_real = Spectrum[0][max_bin];
		max_amp_imag = Spectrum[1][max_bin];
		for(j=0; j<FFT_SIZE/2+1; j++)
		{
			if (j-max_bin+FFT_SIZE/2 >= 0 && j-max_bin+FFT_SIZE/2 < FFT_SIZE) // subtract window spectrum from original spectrum
			{
				float *ROI = MUL2(max_amp_real, max_amp_imag, window_real[j-max_bin+FFT_SIZE/2], window_imag[j-max_bin+FFT_SIZE/2]);
//fprintf(myfile,"%f + %f*i\n",ROI[0],ROI[1]);
				if (j+DelFreqBin >=0 && j+DelFreqBin < FFT_SIZE/2) // shift the entire window spectrum
				{
					float *ROI_shift = MUL2(ROI[0], ROI[1], cos(phasor_pitch[max_bin]), sin(phasor_pitch[max_bin]));

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
	float AmpQuotient = sqrt(energy/new_energy);
	
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
// < Voice cut-off frequency estimation >
// Reference: Applying the Harmonic Plus Noise Model in Concatenative Speech Synthesis, page3
// ============================================================================ //
int PitchShifter::VCO_estimation(float *Mag, int *harmonic_peak, float f0)
{
	// parameter setting
	float VCO_lb = 1000, VCO_ub = 8000; // the upper/lower bound of VCO
	int n = 3; // if n consecutive harmonics are unvoiced, then set VCO

	int i, j, harmonic_no=0, bin_index, VCO;
	int ChannelGroup[FFT_SIZE/2+1];
	float max_minor_mag, central_freq, residual;

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
	float bin_period = FFT_SIZE/(SamplingRate/f0); // f0 interval width counted by bin number
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
		central_freq = SamplingRate*((float)harmonic_bin[i]/FFT_SIZE);
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
//	VCO = SamplingRate*(float)harmonic_bin[i-1]/FFT_SIZE;
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

void PitchShifter::VCO_contour(float ***Spectrum, int n, float *pitch)
{
	int i, j, counter;
	float a,average_bin_th = 0;
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
					a = (float)(voiced_bin_th[i+counter]-voiced_bin_th[i-1])/(counter+1);
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

