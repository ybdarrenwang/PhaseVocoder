#include "time_stretcher.h"

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


// ============================================================================ //
// < Resampling on time axis >
// ============================================================================ //
void TimeStretcher::Stretch(double rate, int FrameNum, double ***spectrogram, 
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
			
			while(ph[j] >= PI)
				ph[j] -= 2.0 * PI;
			while(ph[j] < -1.0*PI)
				ph[j] += 2.0 * PI;
			
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
				while(ph[j] >= PI)
					ph[j] -= 2 * PI;
				while(ph[j] < -1.0*PI)
					ph[j] += 2 * PI;
			}
		}

		// calculate the phase of the other channels, and construct synthesis frame
		for(j=0;j<FFT_SIZE/2+1;j++)
		{
			if (j!=subband_time[j])
			{
				ph[j] += phasor_time[subband_time[j]];
			
				while(ph[j] >= PI)
					ph[j] -= 2.0 * PI;
				while(ph[j] < -1.0*PI)
					ph[j] += 2.0 * PI;
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
