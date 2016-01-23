#ifndef __VOCODERFUNCTIONS_H__
#define __VOCODERFUNCTIONS_H__

#include <math.h>
#include "window.h"

class VocoderFunctions {
    public:
        VocoderFunctions(){}
        virtual ~VocoderFunctions(){}
        void ChannelGrouping(double *Spec, int *ChannelGroupFlag);
        void new_ChannelGrouping(double *Spec, int *ChannelGroupFlag, double pitch);
        void VCO_contour(double ***Spectrum, int n, double *pitch);
        int VCO_estimation(double *Mag, int *harmonic_peak, double f0);
        double ABS2(double a,double b);

    protected:
        double StudentPha[FFT_SIZE/2+1];
        double StudentMag[FFT_SIZE/2+1];
        int* voiced_bin_th;
};

#endif
