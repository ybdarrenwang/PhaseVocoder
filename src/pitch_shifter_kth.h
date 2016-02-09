#ifndef __PITCHSHIFTER_H__
#define __PITCHSHIFTER_H__

#include "vocoder_functions.h"

class PitchShifter : public VocoderFunctions
{
    public:
        void PitchShifting(float **Spectrum, float PitchQuotient, int frame_shift, bool reset_ph);
        void PitchShifting_KTH(float **Spectrum, float PitchQuotient, int frame_shift, bool reset_ph, float *window, float pitch);
        void PitchShifting_KTH_HNM(float **Spectrum, float PitchQuotient, int frame_shift, bool reset_ph, float *window_mag, float *window_pha, float pitch, int vowel_frame_index);
        void PS_KTH_new(float **Spectrum, float PitchQuotient, int frame_shift, bool reset_ph, float *window_real, float *window_imag);

    private:
        float phasor[FFT_SIZE/2+1];
        float prev_phasor[FFT_SIZE/2+1];
        int prev_subband[FFT_SIZE/2+1];
};

#endif
