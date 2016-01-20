#ifndef ENHANCE_H
#define ENHANCE_H

class CEnhance  
{
public:
	CEnhance();
	virtual ~CEnhance();

	// Enable/disable the performing of spectral subtraction on the signal.
	// This reduces the amount of background noise in the signal by subtracting
	// a fraction of the quietest frames from all frames (in the frequency domain).
	// The degree (ie fSubDegree) is a percentage factor on the energy found in the
	// quietest frames. A degree of 100 gives modest background noise removal.
	void EnableSpectralSubtraction(double fSubDegree);
	void EnableSpectralSubtraction();		// use default sub-degree
	void DisableSpectralSubtraction();
	double GetSubtractionDegree() const { return m_fSubDegree; };

	double median (const double* frame, int n);
//	void IniNoiseSpec(const short *npInSampleInputBuffer, int nNumSamples);
	int VAD(const short *npInSampleInputBuffer,int nNumSamples,double fSecondPerSample,short *npSampleOutputBuffer);
	// Do signal enhancement, return 0 if success, < 0 if failed;
	// fSecondsPerSample is the inverse of sampling rate
	int DoEnhance(const short *npInSampleInputBuffer,int nNumSamples,double fSecondPerSample,short *npSampleOutputBuffer,bool *mark);

private:
	int m_nDospectsub;
//	int docompression;
//	int dofixscale;
	double m_fSubDegree;		// subtraction strength
//	short noisemag[1024];
};

#endif
