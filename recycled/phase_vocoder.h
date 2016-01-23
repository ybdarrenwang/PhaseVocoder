#ifndef __PHASEVOCODER_H__
#define __PHASEVOCODER_H__

#include "SourceFilter.h"

#define PV_frame_len 1024

class Vocoder{
public:
	Vocoder();
	virtual ~Vocoder();
	void GetTeacherInfo(CUtterance& tmpUtterance, CUtteranceAlignment& UttAlign);
	void GetStudentInfo(CUtterance& tmpUtterance, CUtteranceAlignment& UttAlign);
	void GetGrades(CUtteranceAlignment& UttAlign);
	void Synthesis(const short *pnSamples);

private:
	double window[FFT_SIZE];
	double dphi[FFT_SIZE/2+1];
	int frame_len;
	int frame_shift;
	int word_no;

	int *TWordFrameNum;
	int *TWordStart, *TWordMid, *TWordEnd;
	int TTotalLength;
//	double *TAveEng;

	int *SWordFrameNum;
	int *SWordStart, *SWordMid, *SWordEnd;
	int STotalLength;
//	double *SAveEng;

	double **TPitch, **SPitch;
	int **TPerceptedPitch, **SPerceptedPitch;
	double TPerceptedRange, SPerceptedRange;
	unsigned int TMinPitch, TMaxPitch, SMinPitch, SMaxPitch;

	int *PitchGrade, *TimingGrade;
	int *TPercent, *SPercent;
//	int *TEnergy, *SEnergy;

	int *new_SWordFrameNum;
	int *new_SWordStart, *new_SWordMid, *new_SWordEnd;
	VocoderFunc vocoder_function;
	SourceFilter s_f_model;
};

#endif