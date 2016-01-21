#include "pitch_shifter.h"

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
// < Resampling on Frequency Dimension >
// ============================================================================ //
void PitchShifter::PitchShifting(double **Spectrum, double PitchQuotient, int frame_shift, bool reset_ph)
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
void PitchShifter::PitchShifting_KTH(double **Spectrum, double PitchQuotient, int frame_shift, bool reset_ph, double *window, double pitch)
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

void PitchShifter::PitchShifting_KTH_HNM(double **Spectrum, double PitchQuotient, int frame_shift, bool reset_ph, 
										double *window_mag, double *window_pha, double pitch, int vowel_frame_index)
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

double *MUL2(double r1, double i1, double r2, double i2) //multiplation of complex numbers
{
	double *ans = new double[2];

	ans[0] = r1*r2-i1*i2;
	ans[1] = r1*i2+r2*i1;

	return ans;
}

void PitchShifter::PS_KTH_new(double **Spectrum, double PitchQuotient, int frame_shift, bool reset_ph, 
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
