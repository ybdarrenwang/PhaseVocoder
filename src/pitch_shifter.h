#ifndef __PITCHSHIFTER_H__
#define __PITCHSHIFTER_H__

#include "vocoder_functions.h"

class PitchShifter : public VocoderFunctions
{
    public:
        void PitchShifting(double **Spectrum, double PitchQuotient, int frame_shift, bool reset_ph);
        void PitchShifting_KTH(double **Spectrum, double PitchQuotient, int frame_shift, bool reset_ph, double *window, double pitch);
        void PitchShifting_KTH_HNM(double **Spectrum, double PitchQuotient, int frame_shift, bool reset_ph, 
                double *window_mag, double *window_pha, double pitch, int vowel_frame_index);
        void PS_KTH_new(double **Spectrum, double PitchQuotient, int frame_shift, bool reset_ph, double *window_real, double *window_imag);

    private:
        double phasor_pitch[FFT_SIZE/2+1];
        double prev_phasor_pitch[FFT_SIZE/2+1];
        int prev_subband_pitch[FFT_SIZE/2+1];
};

#endif
