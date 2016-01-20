#include <math.h>
#include <stdlib.h>
#include "Enhance.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif
/* manifest constants */
#define WINTIME	0.030		/* process in 30ms windows */
#define DBRANGE	40.0		/* dynamic range */

#define buffsize 1024
#define SPEECH 1
#define	NONSPEECH 0

/* buffering */
short	*buff;
double	*window;
double	maxenergy;

/* operations selected */
const double	mu=10.0;				/* mu value for compression */

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CEnhance::CEnhance()
{
	EnableSpectralSubtraction();	// use default subdegree
}

CEnhance::~CEnhance()
{

}

void CEnhance::EnableSpectralSubtraction()
{
	const double fDefaultSubDegree = 100.;
	EnableSpectralSubtraction(fDefaultSubDegree);
}

void CEnhance::EnableSpectralSubtraction(double fSubDegree)
{
	m_nDospectsub = 1;
	m_fSubDegree = fSubDegree;
}

void CEnhance::DisableSpectralSubtraction()
{
	m_nDospectsub = 0;
}

int compare (const void* a, const void* b)
{
  return ( *(int*)a - *(int*)b );
}

double CEnhance::median (const double* frame, int n)
{
	double* nframe = new double[n];
	for (int i=0; i<n; i++)
		nframe[i] = frame[i];

	qsort(nframe, n, sizeof(double), compare);
	if (n % 2 == 1)
		return (nframe[(int)floor((double)n/2)]);
	else
		return (0.5*(nframe[(int)((double)n/2)]+nframe[(int)((double)n/2-1)]));
}

//////////////////////////////////////////////////////////////////////
// Voice Active Detection
//////////////////////////////////////////////////////////////////////
int CEnhance::VAD(const short *pInSamples,int nNumSamples,double frameduration,short *pOutSamples)
{
	/* Initialize */
	int i,j,counter = 0;
	double quantile = 0.9, p = 0.1;
	double total_SNR = 0;
	int frame_shift = buffsize/4, 
		noise_frame = 18,
		threshold = 3,    // unit : dB
		L = 10,
		r = 2*L*quantile,
		k = floor(r),
		f = r-k;
//FILE* MyFile = fopen("VAD.txt","w");

	int frame_no = floor(((int)nNumSamples - buffsize)/frame_shift) + 1;
	double* energy = new double[frame_no]; //log energy
	double* noise = new double[frame_no];
	double* signal = new double[frame_no];
	bool* mark = new bool[frame_no];
	double* energy_sort = new double[2*L+1];
	double* SNR = new double[frame_no];

	/* read in speech samples */
	short* speech = new short[nNumSamples];
	for (i = 0 ; i < (int)nNumSamples; i++)
		speech[i] = pInSamples[i];	

	/* framing & computing log energy */
	for (i = 0 ; i < frame_no ; i++)
	{
		energy[i] = 0;

		for (j = i*frame_shift ; j < (i*frame_shift + buffsize) ; j++)
			energy[i] += speech[j]*speech[j];
	      
		energy[i] = log10(energy[i]);
	}

	/* OS filter */
	// initialize noise energy as the median of the beginning 18 frames
	for (i=0 ; i<=noise_frame ; i++)
		noise[i]=median(energy, noise_frame); 

	// VAD for the beginning 18 frames
	int rr, kk, ff;
	for (i=0 ; i<noise_frame ; i++)
	{
	    for (j=0 ; j<=i ; j++)
			energy_sort[j] =  energy[j];

		qsort(energy_sort, i+1, sizeof(double), compare);
	    rr = quantile*i;
	    kk = floor(rr);
	    ff = rr-kk;

	    if (kk!=0)
			signal[i]=(1-ff)*energy_sort[kk]+ff*energy_sort[kk+1];
		else
			signal[i]=ff*energy_sort[kk+1];

		mark[i]=NONSPEECH;
	}

	// VAD for the rest of the frames
	for (i = noise_frame ; i < (frame_no-L) ; i++)
	{
	    for (j = 0; j < (2*L+1); j++)
			energy_sort[j] = energy[i-L+j];

		qsort(energy_sort, 2*L+1, sizeof(double), compare);
	    signal[i] = (1-f)*energy_sort[k] + f*energy_sort[k+1]; //quantile filter
	    SNR[i] = 10*(signal[i]-noise[i]);

	    if (SNR[i] >= threshold)
		{
			mark[i] = SPEECH;
			noise[i+1] = noise[i];
			total_SNR += SNR[i];
			counter++;
		}
		else
		{
			mark[i] = NONSPEECH;
			noise[i+1] = (1-p)*noise[i] + p*energy_sort[L];
		}
	}
	
	for (i = frame_no-L ; i < frame_no ; i++)
		mark[i]=NONSPEECH;

	// check if nonspeech segments are at least 5 frames
	i = noise_frame;
	while (i < frame_no-L)
	{
		j=0;

		while (mark[i+j] == mark[i] && i+j < frame_no-L)
			j++;
	
		if (j<=5 && i+j < frame_no-L)
		{
			while(j>0)
			{
				mark[i]++;
				i++;
				j--;
			}
		}
		else
			i+=j;
	}

	if (total_SNR/counter < 15)
		DoEnhance(pInSamples,nNumSamples,frameduration,pOutSamples,mark);
	else
	{
		for (i = 0; i<nNumSamples ;i++)
			pOutSamples[i] = pInSamples[i];
	}

	delete [] energy;
	delete [] noise;
	delete [] signal;
	delete [] mark;
	delete [] energy_sort;
	delete [] SNR;

	return 0;
}

//////////////////////////////////////////////////////////////////////
// Spectral Subtraction
//////////////////////////////////////////////////////////////////////
int CEnhance::DoEnhance(const short *pInSamples,int nNumSamples,double frameduration,short *pOutSamples,bool *mark)
{
	extern int	optind;		/* option index */
	extern char	*optarg;	/* option argument */
	REAL		*fsp,*wsp,*mag2,*minmag2,*lstmag2;
	COMPLEX		*spect;
	double		smooth_mag, smooth_noise;
	double		factor;
	int frame_shift = buffsize/4;
	int frame_num = floor((double)nNumSamples/frame_shift) + 1;

	/* get buffers */
	short* buff = new short[nNumSamples];
	int* hold = new int[buffsize];
	short* noisemag2 = new short[buffsize];

	for (int i = 0 ; i < (int)nNumSamples; i++)
		buff[i] = pInSamples[i];

//	memset(hold,0,buffsize*sizeof(short));
	window = (double *)calloc(buffsize,sizeof(double));

	/* get buffers */
	fsp = (REAL *)calloc(buffsize,sizeof(REAL));
	wsp = (REAL *)calloc(buffsize,sizeof(REAL));
	mag2 = (REAL *)calloc(buffsize,sizeof(REAL));
	lstmag2 = (REAL *)calloc(buffsize,sizeof(REAL));
	minmag2 = (REAL *)calloc(buffsize,sizeof(REAL));
	spect = (COMPLEX *)calloc(buffsize,sizeof(COMPLEX));

	/* calculate raised cosine window */
	for (i=0;i<buffsize;i++)
		window[i] = cos(PI * (i-buffsize/2) / buffsize);

	/* set smoothing factor */
	smooth_mag = 0; // no smoothing
	smooth_noise = 0.9;

	/* process file */
	for (i=0; i<frame_num && i*frame_shift+buffsize<(int)nNumSamples; i++)
	{
		// framing, windowing
		for (int j=0;j<buffsize;j++)
			fsp[j] = buff[i*frame_shift+j] * window[j];
		rfft(fsp,spect,buffsize);

		for (j=0;j<buffsize;j++)
			mag2[j] = sqrt(spect[j].re*spect[j].re + spect[j].im*spect[j].im);
	
		if (i==0) // initialize
		{
			for (j=0;j<buffsize;j++)
			{
				lstmag2[j] = mag2[j];
				noisemag2[j] = mag2[j]; 
				hold[j] = 0;
			}
		}
		else
		{
			for (j=0;j<buffsize;j++)
				lstmag2[j] = smooth_mag*lstmag2[j] + (1-smooth_mag)*mag2[j]; // smoothing of magnitude spectrum
		}

		/* Spectral Subtraction */
		for (j=0;j<buffsize;j++)
		{
			if (mark[i] == SPEECH)
			{
				factor = (double)(lstmag2[j]-0.5*noisemag2[j])/(double)lstmag2[j];
				if (factor < 0)
					factor = 0;
			}
			else
			{
				factor = 0.5;
				noisemag2[j] = smooth_noise*noisemag2[j] + (1-smooth_noise)*lstmag2[j]; // update noise spectrum
			}

			spect[j].re *= factor;
			spect[j].im *= factor;
		}
		irfft(fsp,spect,buffsize);

		/* Use a hold buffer to sum up the frames and then output*/
		for (j=0; j<buffsize; j++)
			hold[j] = (int)(hold[j] + fsp[j]*window[j]);

		/* write out result of processing */
		for (j=0; j<frame_shift; j++)
			buff[i*frame_shift+j] = (short)((double)hold[j]/2);

		/* shift hold buffer and reset */
		for (j=0; j<buffsize-frame_shift; j++)
			hold[j] = hold[frame_shift + j];
		for (j=buffsize-frame_shift; j<buffsize; j++)
			hold[j] = 0;
	}
	// process the last frame
	for (int j=0;j < (int)buffsize-frame_shift;j++)
	{
		if ( (i*frame_shift+j) < (int)nNumSamples)
			buff[i*frame_shift+j] = (short)((double)hold[j]/2);
	}

	/* write out result */
	for (int l = 0; l < (int) nNumSamples ;l++)
		pOutSamples[l] = buff[l];

	// release memory
	if (buff)
		delete [] buff;
	if (hold)
		delete [] hold;
	free(fsp);
	free(wsp);
	free(mag2);
	free(minmag2);
	free(spect);
	// Chialin: more to free
	free(window);
	free(lstmag2);
	return 0;
}
