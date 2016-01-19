#ifndef __VOCODERFUNCTIONS_H__
#define __VOCODERFUNCTIONS_H__

#include "EngineDLL.h"
#include "WaveFrame.h"
#include <math.h>
#include "CPerformer.h"
#include "MyFFT.h"

#define FFT_SIZE 4096
#define SamplingRate 22050
#define Atten_3dB 0.707945784

class VocoderFunc{
public:
	double ABS2(double a,double b);
	void pitch_contour_refine(double *pitch, int n);
	void cal_target_contour(double **t_pitch, int *m, double **s_pitch, int *n, 
							double **tar_pitch, int word_no);
	void Hamming(int len, double *window);
	void VCO_contour(double ***Spectrum, int n, double *pitch);
//	void CalPath(double *path, double rate, int ElementCount);
//	void SFDecomp(double **Spectrum);
//	void ChannelGrouping(double *Spec, int *ChannelGroupFlag);
//	void new_ChannelGrouping(double *Spec, int *ChannelGroupFlag, double pitch);
	void PitchShifting(double **Spectrum, double PitchQuotient, int frame_shift, bool reset_ph);
	void PitchShifting_KTH(double **Spectrum, double PitchQuotient, int frame_shift, bool reset_ph, double *window, double pitch);
	void PitchShifting_KTH_HNM(double **Spectrum, double PitchQuotient, int frame_shift, bool reset_ph, 
								double *window_mag, double *window_pha, double pitch, 
								FILE* myfile1, FILE* myfile2, int vowel_frame_index);
	void PS_KTH_new(double **Spectrum, double PitchQuotient, int frame_shift, bool reset_ph, double *window_real, double *window_imag);
	void TimeScaleMod(double rate, int FrameNum, double ***spectrogram, int new_FrameNum, double ***new_spectrogram, bool reset_ph);
	int PSOLA(double **frame, double *output, double *SentenceCoeff, double *windows,
			int sample_index, double *Pitch, int WordStart, int WordMid, int WordEnd,
			int frame_len, int frame_shift, int bias);
	int non_PSOLA(double **frame, double *output, double *SentenceCoeff, double *windows,
				int sample_index, int WordStart, int WordMid, int WordEnd, 
				int frame_len, int frame_shift, int bias);

private:
	double StudentPha[FFT_SIZE/2+1];
	double StudentMag[FFT_SIZE/2+1];

	// for pitch shifting
	double phasor_pitch[FFT_SIZE/2+1];
	double prev_phasor_pitch[FFT_SIZE/2+1];
	int prev_subband_pitch[FFT_SIZE/2+1];
	int* voiced_bin_th;

	// for time scale modification
	double phasor_time[FFT_SIZE/2+1];
	int prev_subband_time[FFT_SIZE/2+1];
	double prev_phase[FFT_SIZE/2+1];
};

#endif