#include "stdafx.h"
#include "phase_vocoder.h"
#include <stdio.h>
#include <iostream>
using namespace std;

bool ini_teacher = false;

// ============================================================================ //
// Constructor/Destructor
// ============================================================================ //

Vocoder::Vocoder()
{
	frame_len = (int)(0.025*SamplingRate);
	frame_shift = (int)(0.010*SamplingRate);

	int i;

	// Hamming Window
	vocoder_function.Hamming(PV_frame_len, window);

	// expected phase advance in each bin
	dphi[0]=0;
	for(i=1;i<FFT_SIZE/2+1;i++)
		dphi[i] = (2.0*PI*frame_shift) / FFT_SIZE * i;
}

Vocoder::~Vocoder()
{
	if (ini_teacher == true) // Delete teacher information if needed
	{	
		delete [] TWordStart;
		delete [] TWordMid;
		delete [] TWordEnd;
		delete [] TWordFrameNum;
		for (int i=0; i<word_no; i++)
		{
			delete [] TPerceptedPitch[i];
			delete [] TPitch[i];
		}
		delete [] TPerceptedPitch;
		delete [] TPitch;
	}
}

// ============================================================================ //
// < Gather required information from teacher/student speech >
// ============================================================================ //
// 1. The start/end frame number of consonant/vowel of each syllable
// 2. The total frame count of each syllable
// 3. The normalized log-F0
// 4. The grade of pitch/timing
// 5. The time percentage of each syllable w.r.t. the whole utterance

void Vocoder::GetTeacherInfo(CUtterance& tmpUtterance, CUtteranceAlignment& UttAlign)
{	
	int counter = 0, i, j;
	TTotalLength = 0;
//	int pitch;

	if (ini_teacher == true) // Delete teacher information if needed
	{	
		delete [] TWordStart;
		delete [] TWordMid;
		delete [] TWordEnd;
		delete [] TWordFrameNum;
		for (i=0; i<word_no; i++)
		{
			delete [] TPerceptedPitch[i];
			delete [] TPitch[i];
		}
		delete [] TPerceptedPitch;
		delete [] TPitch;
	}

	ini_teacher = true;
	word_no = UttAlign.GetNumWordAlignments();

	TWordStart = new int[word_no];
	TWordMid = new int[word_no];
	TWordEnd = new int[word_no];
	TWordFrameNum = new int[word_no];
	TPerceptedPitch = new int *[word_no];
	TPitch = new double *[word_no];
//	TAveEng = new double *[word_no];

	tmpUtterance.FindMinMaxFramePitchForNorm(tmpUtterance.GetWaveFrameArray(), tmpUtterance.GetNumFrames(), TMinPitch, TMaxPitch);
	TPerceptedRange = log((double)TMaxPitch) - log((double)TMinPitch);

	for (i=0; i<word_no; i++)
	{	
		TWordStart[i] = (int)UttAlign[i].GetWordStartFrame();
		TWordMid[i] = (int)UttAlign[i].GetPhonemeStartFrame(1);
		TWordEnd[i] = (int)UttAlign[i].GetWordEndFrame();
		TWordFrameNum[i] = TWordEnd[i] - TWordStart[i];
		TPerceptedPitch[i] = new int [TWordEnd[i]-TWordMid[i]];
		TPitch[i] = new double [TWordEnd[i]-TWordMid[i]];
//FILE* myfile=fopen("pitch_contour.txt","w");
		for (j=0; j<TWordEnd[i]-TWordMid[i]; j++)
		{
			CWaveFrame& Frame = tmpUtterance.GetWaveFrame(TWordMid[i]+j);
			TPerceptedPitch[i][j] = (int)Frame.GetPerceptedPitch();
			TPitch[i][j] = exp(log((double)TMinPitch)+TPerceptedRange*((double)Frame.GetPerceptedPitch()-20.0)/100.0);
//fprintf(myfile,"%f/n",TPitch[i][j]);
//			pitch = (int)Frame.GetPitch();
//			TPitch[i][j] = (double)pitch;
		}
		TTotalLength += TWordFrameNum[i];
//fclose(myfile);
		// pitch contour smoothing
		vocoder_function.pitch_contour_refine(TPitch[i], TWordEnd[i]-TWordMid[i]);
	}
}

void Vocoder::GetStudentInfo(CUtterance& tmpUtterance, CUtteranceAlignment& UttAlign)
{
	int counter = 0, i, j;
	STotalLength = 0;

	SWordStart = new int[word_no];
	SWordMid = new int[word_no];
	SWordEnd = new int[word_no];
	SWordFrameNum = new int[word_no];
	SPerceptedPitch = new int *[word_no];
	SPitch = new double *[word_no];
//	SAveEng = new double *[word_no];

	tmpUtterance.FindMinMaxFramePitchForNorm(tmpUtterance.GetWaveFrameArray(), tmpUtterance.GetNumFrames(), SMinPitch, SMaxPitch);
	SPerceptedRange = log((double)SMaxPitch) - log((double)SMinPitch);

	for (i=0; i<word_no; i++)
	{	
		SWordStart[i] = (int)UttAlign[i].GetWordStartFrame();
		SWordMid[i] = (int)UttAlign[i].GetPhonemeStartFrame(1);
		SWordEnd[i] = (int)UttAlign[i].GetWordEndFrame();
		SWordFrameNum[i] = SWordEnd[i] - SWordStart[i];
		SPerceptedPitch[i] = new int [SWordEnd[i] - SWordMid[i]];
		SPitch[i] = new double [SWordEnd[i] - SWordMid[i]];

		for (j=0; j<SWordEnd[i] - SWordMid[i]; j++)
		{
			CWaveFrame& Frame = tmpUtterance.GetWaveFrame(SWordMid[i]+j);
			SPerceptedPitch[i][j] = (int)Frame.GetPerceptedPitch();
			SPitch[i][j] = exp(log((double)SMinPitch)+SPerceptedRange*((double)Frame.GetPerceptedPitch()-20.0)/100.0);
//			SPitch[i][j] = (double)Frame.GetPitch();
		}
		STotalLength += SWordFrameNum[i];

		// pitch contour smoothing
		vocoder_function.pitch_contour_refine(SPitch[i], SWordEnd[i]-SWordMid[i]);
	}
}

void Vocoder::GetGrades(CUtteranceAlignment& UttAlign)
{
//	PitchGrade = new int[word_no];
//	TimingGrade = new int[word_no];
	TPercent = new int[word_no];
	SPercent = new int[word_no];
//	TEnergy = new int[word_no];
//	SEnergy = new int[word_no];
//FILE* myfile = fopen("energy.txt","w");
	for (int i=0; i<word_no; i++)
	{	
		const CWordGrade& WordGrade = UttAlign[i].GetWordGrade();		
//		PitchGrade[i] = WordGrade.GetPitchGrade(0);
//		TimingGrade[i] = WordGrade.GetTimingGrade(0);

		TPercent[i] = WordGrade.GetTeacherTimingPercent(0);
        SPercent[i] = WordGrade.GetStudentTimingPercent(0);

//		TEnergy[i] = WordGrade.GetTeacherIntensityPercent();
//		SEnergy[i] = WordGrade.GetStudentIntensityPercent();
//fprintf(myfile,"%d	%d\n",TEnergy[i],SEnergy[i]);
	}
//fclose(myfile);
}

// ============================================================================ //
//  < Synthesis >
// ============================================================================ //
//	1. Do STFT on each frame. 
//	2. Performe pitch-scaling, which is in fact simply resampling on the frequency axle, w.r.t. each frame. 
//	3. Performe time-scaling on the whole syllable spectrogram.
//	4. Finally, by combining all the words into one sentence, we get the pitch-scaled, time-scaled utterance.

void Vocoder::Synthesis(const short* pnSamples)
{	
	/* Initialization */
	int i, j, k;
	int in_sample_index, pitch_index,new_sample_no = 0;
	int	bias = (int)((PV_frame_len-frame_len)/2);
	int sample_no = (SWordEnd[word_no-1]-SWordStart[0])*frame_shift+(frame_len-frame_shift)+2*frame_shift; // add 20ms blank signals at start
	MyFFT myfft;
	double b;

	// data
	double ***StudentSpec, ***new_StudentSpec;
	double *StudentWave, **NewStudentWave;
	double *WaveImag; // just a temp for imaginary part of IFFT wave output
	double ***Spec_vowel;

	// for pitch-sync. overlap add of output frames
	double **new_SPitch = new double *[word_no]; // pitch contour after pitch-shifting
	double **new2_SPitch = new double *[word_no]; // pitch contour after time-scaling
	int *syn_bias = new int[word_no];

	// for phase reset
	bool *reset_ph_time = new bool[word_no];
	bool reset_ph_pitch;

	// for gain control
	double threshold_GC = 25000;
	bool saturate;
	double max_amp, gain_ctrl;

	// for transient detection
	int find_tran_range = 4; // determine how many frames around syllible boundary to check
	double threshold_TD = 10;

	/* Memory allocation */
	double *Excitation = new double[sample_no];
	double *NewExcitation = new double[2*sample_no];
	memset(NewExcitation,0.0,2*sample_no*sizeof(double));
	double *NewSentenceCoeff = new double[2*sample_no];
	memset(NewSentenceCoeff,0.0,2*sample_no*sizeof(double));
	double *WordDurMul = new double[word_no];
	double *spectral_real = new double[FFT_SIZE];
	double *spectral_imag = new double[FFT_SIZE];
	new_SWordFrameNum = new int[word_no];
	new_SWordStart = new int[word_no];
	new_SWordMid = new int[word_no];
	new_SWordEnd = new int[word_no];

	/* Read data */
	double *my_samples = new double[sample_no];
	saturate = false;
	for (i=0; i<sample_no-2*frame_shift; i++)
	{
		my_samples[i+2*frame_shift] = (double)pnSamples[SWordStart[0]*frame_shift+i];
		
		if (abs(my_samples[i+2*frame_shift]) > threshold_GC) // gain control
		{
			saturate = true;
			max_amp = abs(my_samples[i+2*frame_shift]);
		}
	}

	for (i=0; i<2*frame_shift; i++) // add 20ms blank signals at start
		my_samples[i] = 0;
	SWordStart[0] -= 2;
	SWordFrameNum[0] += 2;

	for (i=word_no-1; i>=0; i--) // set the 1st frame of student utt. as beginning of samples
	{
		SWordEnd[i] -= SWordStart[0];
		SWordMid[i] -= SWordStart[0];
		SWordStart[i] -= SWordStart[0];
	}

	/* Transient Detection */
	double *frame_eng = new double[find_tran_range];
	for (i=0; i<find_tran_range; i++)
		frame_eng[i] = 0;

	reset_ph_time[0] = true;
	for (i=1; i<word_no; i++)
	{
		if (SWordEnd[i-1] != SWordStart[i])
			reset_ph_time[i] = true;
		else
		{
			// calculate frame energy
			for (j=0; j<find_tran_range; j++)
			{
				for (k=0; k<frame_len; k++)
				{
					in_sample_index = (SWordStart[i]-find_tran_range/2.0+j)*frame_shift-bias+k;
					frame_eng[j] += (double)my_samples[in_sample_index]*my_samples[in_sample_index]/frame_len;
				}
			}
			
			// transient detection
			reset_ph_time[i] = false;
			for (j=0; j<find_tran_range-1; j++)
			{
				if (frame_eng[j+1]/frame_eng[j] > threshold_TD)
					reset_ph_time[i] = true;
			}
		}
	}
	delete [] frame_eng;

	/* Gain control */
	if (saturate == true) // if max amplitude exceed 25000, then attenuate signal such that max amplitude be 25000.
	{
		gain_ctrl = threshold_GC/max_amp;
		for (i=0; i<sample_no; i++)
			my_samples[i] *= gain_ctrl;
	}

	/* Source-Filter Decomposition */
//	s_f_model.SFDecomp_sync(my_samples, Excitation, sample_no, word_no, SWordStart, SWordMid, SWordEnd, frame_shift, frame_len, SPitch);
	s_f_model.SFDecomp_non_sync(my_samples, Excitation, sample_no, word_no, SWordStart, SWordMid, SWordEnd, frame_shift, frame_len);
	delete [] my_samples;

	/* Calcuate the target pitch contour for pitch shifting */
	int *vowel_frame_no = new int[word_no];
	int *new_vowel_frame_no = new int[word_no];

	for (i=0; i<word_no; i++)
	{
		new_SPitch[i] = new double[SWordEnd[i]-SWordMid[i]];
		vowel_frame_no[i] = TWordEnd[i]-TWordMid[i];
		new_vowel_frame_no[i] = SWordEnd[i]-SWordMid[i];
	}

	vocoder_function.cal_target_contour(TPitch, vowel_frame_no, SPitch, new_vowel_frame_no, new_SPitch, word_no);
	
	delete [] vowel_frame_no;
	delete [] new_vowel_frame_no;

	/* Calculate the spectrum of hamming window */
	double time_hamming[FFT_SIZE], spec_hamming_real[FFT_SIZE], spec_hamming_imag[FFT_SIZE];
	for (i=0; i<PV_frame_len; i++)
		time_hamming[i] = window[i];

	if (PV_frame_len < FFT_SIZE) // zero-crossing
	{
		for (i=PV_frame_len; i<FFT_SIZE; i++)
			time_hamming[i] = 0;
	}

	myfft.fft_double(FFT_SIZE,0,window,NULL,spectral_real,spectral_imag);
	double normalizer = sqrt(spectral_real[0]*spectral_real[0]+spectral_imag[0]*spectral_imag[0]);
	
	for (i=0; i<FFT_SIZE/2; i++)
	{
		spec_hamming_real[i+FFT_SIZE/2] = spectral_real[i]/normalizer;
		spec_hamming_imag[i+FFT_SIZE/2] = spectral_imag[i]/normalizer;
	}
	for (i=FFT_SIZE/2; i<FFT_SIZE; i++)
	{
		spec_hamming_real[i-FFT_SIZE/2] = spectral_real[i]/normalizer;
		spec_hamming_imag[i-FFT_SIZE/2] = spectral_imag[i]/normalizer;
	}

	/* Phase Vocoder */
	for (i=0; i<word_no; i++)
	{	
		// Initialize Word Spectrogram
		StudentSpec = new double **[SWordFrameNum[i]];
		for (j=0; j<SWordFrameNum[i]; j++)
			StudentSpec[j] = new double *[2];
		for (j=0; j<SWordFrameNum[i]; j++)
		{
			StudentSpec[j][0] = new double [FFT_SIZE/2+1];
			StudentSpec[j][1] = new double [FFT_SIZE/2+1];
		}

		/* Windowed STFT */
		for (j=0; j<SWordFrameNum[i]; j++)
		{	
			StudentWave = new double[FFT_SIZE];
			for (k=0; k<PV_frame_len; k++)
			{
				in_sample_index = (SWordStart[i]+j)*frame_shift-bias+k;
				if (in_sample_index<sample_no && in_sample_index>=0)
					StudentWave[k] = window[k]*Excitation[in_sample_index];
				else
					StudentWave[k] = 0;
			}
			// zero-crossing
			if (PV_frame_len < FFT_SIZE)
			{
				for (k=PV_frame_len; k<FFT_SIZE; k++)
					StudentWave[k] = 0;
			}

			myfft.fft_double(FFT_SIZE,0,StudentWave,NULL,spectral_real,spectral_imag);
			delete [] StudentWave;

			for (k=0; k<FFT_SIZE/2+1; k++)
			{
				StudentSpec[j][0][k] = spectral_real[k];
				StudentSpec[j][1][k] = spectral_imag[k];
			}
		}

		/* Calculate VCO contour */
/*		double*** Spec_vowel = new double** [SWordEnd[i]-SWordMid[i]];
		for (j=0; j<SWordEnd[i]-SWordMid[i]; j++) // record vowel spectrum
		{
			Spec_vowel[j] = new double* [2];
			Spec_vowel[j][0] = new double [FFT_SIZE/2+1];
			Spec_vowel[j][1] = new double [FFT_SIZE/2+1];

			for (k=0; k<FFT_SIZE/2+1; k++)
			{
				Spec_vowel[j][0][k] = StudentSpec[j+(SWordMid[i]-SWordStart[i])][0][k];
				Spec_vowel[j][1][k] = StudentSpec[j+(SWordMid[i]-SWordStart[i])][1][k];
			}
		}
		
		vocoder_function.VCO_contour(Spec_vowel, SWordEnd[i]-SWordMid[i], SPitch[i]);

		for (j=0; j<SWordEnd[i]-SWordMid[i]; j++) // clear vowel spectrum
		{	
			delete [] Spec_vowel[j][0];
			delete [] Spec_vowel[j][1];
			delete [] Spec_vowel[j];
		}
		delete [] Spec_vowel;*/

		/* Pitch Shifting on Vowel */
		reset_ph_pitch = true; // reset phase only at the beginning of vowel
		for (j=0; j<SWordEnd[i]-SWordMid[i]; j++)
		{
//			vocoder_function.PitchShifting(StudentSpec[j+(SWordMid[i]-SWordStart[i])], new_SPitch[i][j]/SPitch[i][j], frame_shift, reset_ph_pitch);
//			vocoder_function.PitchShifting_KTH(StudentSpec[j+(SWordMid[i]-SWordStart[i])], new_SPitch[i][j]/SPitch[i][j], frame_shift, reset_ph_pitch, spec_hamming, SPitch[i][j]);
//			vocoder_function.PitchShifting_KTH_HNM(StudentSpec[j+(SWordMid[i]-SWordStart[i])], new_SPitch[i][j]/SPitch[i][j], 
//												frame_shift, reset_ph_pitch, spec_hamming_mag, spec_hamming_pha, SPitch[i][j],
//												myfile1, myfile2, j);
			vocoder_function.PS_KTH_new(StudentSpec[j+(SWordMid[i]-SWordStart[i])], new_SPitch[i][j]/SPitch[i][j], 
										frame_shift, reset_ph_pitch, spec_hamming_real, spec_hamming_imag);

			reset_ph_pitch = false;
		}

		/* Source-Filter Recombination (freq. domain) */
		for (j=0; j<SWordEnd[i]-SWordStart[i]; j++)
			s_f_model.SFRecomb_FD(StudentSpec[j], i, j, SWordMid[i]-SWordStart[i]);

		/* Time-scale modification */
		new_SWordFrameNum[i] = floor(((double)TWordFrameNum[i]/TTotalLength) * STotalLength);
		WordDurMul[i] = (double)new_SWordFrameNum[i]/SWordFrameNum[i];
//		new_SWordFrameNum[i] = SWordFrameNum[i];
//		WordDurMul[i] = 1;

		new_StudentSpec = new double **[new_SWordFrameNum[i]];	
		for (j=0; j<new_SWordFrameNum[i]; j++)
		{
			new_StudentSpec[j] = new double *[2];
			new_StudentSpec[j][0] = new double [FFT_SIZE/2+1];
			new_StudentSpec[j][1] = new double [FFT_SIZE/2+1];
		}

		vocoder_function.TimeScaleMod(1.0/WordDurMul[i], SWordFrameNum[i], StudentSpec,
									new_SWordFrameNum[i], new_StudentSpec, reset_ph_time[i]);

/*		for (j=0; j<new_SWordFrameNum[i]; j++)
		{
			for (k=0; k<FFT_SIZE/2+1; k++)
			{
				new_StudentSpec[j][0][k] = StudentSpec[j][0][k];
				new_StudentSpec[j][1][k] = StudentSpec[j][1][k];
			}
		}*/

		for (j=0; j<SWordFrameNum[i]; j++)
		{
			delete [] StudentSpec[j][0];
			delete [] StudentSpec[j][1];
			delete [] StudentSpec[j];
		}	
		delete [] StudentSpec;

		/* Generate Synthesized Sentence */

		// calculate new pitch contour (for pitch-synchronous overlap-add)
		new_SWordStart[i] = floor((double)new_sample_no/frame_shift);
		new_SWordMid[i] = new_SWordStart[i] + floor(WordDurMul[i]*(double)(SWordMid[i]-SWordStart[i]));
		if (WordDurMul[i]*(double)(SWordMid[i]-SWordStart[i])-floor(WordDurMul[i]*(double)(SWordMid[i]-SWordStart[i]))>=0.5)
			new_SWordMid[i] += 1;
		new_SWordEnd[i] = new_SWordStart[i] + new_SWordFrameNum[i];
				
		new2_SPitch[i] = new double[new_SWordEnd[i]-new_SWordMid[i]];
		for (j=0; j<new_SWordEnd[i]-new_SWordMid[i]; j++)
		{
			pitch_index = floor((double)j/WordDurMul[i]);
			if(pitch_index < SWordEnd[i]-SWordMid[i]-1)
				b = (double)j/WordDurMul[i] - pitch_index;
			else
				b = 0;

			new2_SPitch[i][j] = (1-b)*new_SPitch[i][pitch_index] + b*new_SPitch[i][pitch_index+1];
		}

		// Combine words into sentence
		NewStudentWave = new double *[new_SWordFrameNum[i]];
		for (j=0; j<new_SWordFrameNum[i]; j++) // j: Frame count
		{
			NewStudentWave[j] = new double[FFT_SIZE];
			WaveImag = new double[FFT_SIZE]; // just a temp for imaginary part of IFFT wave output

			for(k=0; k<FFT_SIZE/2+1; k++)
			{
				spectral_real[k] = new_StudentSpec[j][0][k];
				spectral_imag[k] = new_StudentSpec[j][1][k];
			}

			for(k=FFT_SIZE-1; k>FFT_SIZE/2; k--) // 1023 ~ 513
			{  
				spectral_real[k] = spectral_real[FFT_SIZE-k];
				spectral_imag[k] = -1.0 * spectral_imag[FFT_SIZE-k];
			}

			myfft.fft_double(FFT_SIZE,1,spectral_real,spectral_imag,NewStudentWave[j],WaveImag);

			delete [] WaveImag;
		}

		// overlap-add frames
/*		syn_bias[i] = vocoder_function.PSOLA(NewStudentWave, NewExcitation, NewSentenceCoeff, window, new_sample_no, 
											new_SWordStart[i], new_SWordMid[i], new_SWordEnd[i], 
											PV_frame_len, frame_shift, bias);*/
		syn_bias[i] = vocoder_function.non_PSOLA(NewStudentWave, NewExcitation, NewSentenceCoeff, window, new_sample_no, 
												new_SWordStart[i], new_SWordMid[i], new_SWordEnd[i], 
												PV_frame_len, frame_shift, bias);

		new_sample_no += new_SWordFrameNum[i]*frame_shift - syn_bias[i];

		for (j=0; j<new_SWordFrameNum[i]; j++) // j: Frame count
			delete [] NewStudentWave[j];
		delete [] NewStudentWave;

		// insert silence between words
		if (i<word_no-1 && TWordStart[i+1] != TWordEnd[i])
		{	
			int SilenceFrameNum = TWordStart[i+1]-TWordEnd[i];

			for (j=0; j<SilenceFrameNum; j++)
			{
				for(k=0; k<PV_frame_len; k++)
					NewSentenceCoeff[new_sample_no-bias+k + j*frame_shift] += window[k]*window[k];
			}
			for (j=0; j<SilenceFrameNum*frame_shift; j++)
				NewExcitation[new_sample_no+j]=0;

			new_sample_no += SilenceFrameNum * frame_shift;
		}

		// release memory
		for (j=0; j<new_SWordFrameNum[i]; j++)
		{
			delete [] new_StudentSpec[j][1];
			delete [] new_StudentSpec[j][0];
			delete [] new_StudentSpec[j];
		}
		delete [] new_StudentSpec;
	}
	new_sample_no += (frame_len-frame_shift);

	// divide by the sum of square of window coefficients
	for (i=0; i<new_sample_no; i++)
	{
		if (NewSentenceCoeff[i] != 0)
			NewExcitation[i] /= NewSentenceCoeff[i];
	}
	for(i=0; i<2*frame_shift-1; i++)
		NewExcitation[i] = 0;

	short *NewSentence = new short[new_sample_no];
	memset(NewSentence,0,new_sample_no*sizeof(short));

	/* Source-Filter Recombination (time domain) */
/*	s_f_model.SFRecomb(NewSentence, NewExcitation, new_sample_no, word_no, WordDurMul, 
						new_SWordStart, new_SWordMid, new_SWordEnd, frame_shift, frame_len, syn_bias);
	
	for(i=0; i<2*frame_shift-1; i++)
		NewSentence[i] = 0;*/

	/* Gain control */
	max_amp = 0;
	for (i=0; i<new_sample_no; i++)
	{
		if (abs(NewExcitation[i]) > max_amp)
			max_amp = abs(NewExcitation[i]);
	}

	// output data
	if (max_amp < threshold_GC) // amplify the signal so that the max amplitude be 25000.
		gain_ctrl = threshold_GC/max_amp;
	else
		gain_ctrl = 1;
	
	for (i=0; i<new_sample_no; i++)
		NewSentence[i] =  (short)floor(NewExcitation[i]*gain_ctrl);	
		
	/* Generate output wave file */
	COWave SynthesisWave;
	WORD wNumChannels = 1;   // mono
	DWORD nNumSamplesPerSec = 22050;  // sampling rate
	WORD wBitsPerSample = 16;   // sample resolution
	bool bCopyBuffer = true;
	SynthesisWave.BuildFormat(wNumChannels, nNumSamplesPerSec, wBitsPerSample);
	
	SynthesisWave.SetBuffer((BYTE *)NewSentence, new_sample_no, bCopyBuffer);
	SynthesisWave.Save(_T("test.wav"));

	delete [] NewSentence;

	// output excitation
/*	short *short_excitation = new short[sample_no];
	for (i=0; i<sample_no; i++)
		short_excitation[i] = (short)floor(Excitation[i]);

	SynthesisWave.SetBuffer((BYTE *)short_excitation, sample_no, bCopyBuffer);
	SynthesisWave.Save(_T("excitation.wav"));

	delete [] short_excitation;*/

	// output new_excitation
/*	short_excitation = new short[new_sample_no];
	for (i=0; i<new_sample_no; i++)
		short_excitation[i] = (short)floor(NewExcitation[i]);

	SynthesisWave.SetBuffer((BYTE *)short_excitation, new_sample_no, bCopyBuffer);
	SynthesisWave.Save(_T("new_excitation.wav"));

	delete [] short_excitation;*/

	// Release memory
	s_f_model.clear_mem(word_no);
	delete [] Excitation;
	delete [] NewExcitation;
	delete [] reset_ph_time;
	delete [] NewSentenceCoeff;
	delete [] spectral_real;
	delete [] spectral_imag;

	// Delete student information
	delete [] WordDurMul;
	delete [] SWordStart;
	delete [] SWordMid;
	delete [] SWordEnd;
	delete [] SWordFrameNum;

	for (i=0; i<word_no; i++)
	{
		delete [] SPerceptedPitch[i];
		delete [] SPitch[i];
		delete [] new_SPitch[i];
		delete [] new2_SPitch[i];
	}
	delete [] SPerceptedPitch;
	delete [] SPitch;
	delete [] new_SPitch;
	delete [] new2_SPitch;
	delete [] syn_bias;

	delete [] TPercent;
	delete [] SPercent;
//	delete [] TEnergy;
//	delete [] SEnergy;
	delete [] new_SWordFrameNum;
	delete [] new_SWordStart;
	delete [] new_SWordMid;
	delete [] new_SWordEnd;
}
