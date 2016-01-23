#ifndef __TIMESTRETCHER_H__
#define __TIMESTRETCHER_H__

#include "vocoder_functions.h"

class TimeStretcher : public VocoderFunctions
{
    public:
        void Stretch(double rate, int FrameNum, double ***spectrogram, int new_FrameNum, double ***new_spectrogram, bool reset_ph);

    private:
        double phasor_time[FFT_SIZE/2+1];
        int prev_subband_time[FFT_SIZE/2+1];
        double prev_phase[FFT_SIZE/2+1];
};

#endif
