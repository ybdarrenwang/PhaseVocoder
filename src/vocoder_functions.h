#ifndef __VOCODERFUNCTIONS_H__
#define __VOCODERFUNCTIONS_H__

#include <math.h>
#include "window.h"

class VocoderFunctions {
    public:
        VocoderFunctions(){}
        virtual ~VocoderFunctions(){}

        void cal_target_contour(double **t_pitch, int *m, double **s_pitch, int *n, 
                double **tar_pitch, int word_no);
        void VCO_contour(double ***Spectrum, int n, double *pitch);
        //	void CalPath(double *path, double rate, int ElementCount);
        //	void SFDecomp(double **Spectrum);
        void ChannelGrouping(double *Spec, int *ChannelGroupFlag);
        void new_ChannelGrouping(double *Spec, int *ChannelGroupFlag, double pitch);
        int VCO_estimation(double *Mag, int *harmonic_peak, double f0);
        double ABS2(double a,double b);

    protected:
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
