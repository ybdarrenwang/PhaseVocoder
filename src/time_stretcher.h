#ifndef __TIMESTRETCHER_H__
#define __TIMESTRETCHER_H__

#include "vocoder_functions.h"

class TimeStretcher : public VocoderFunctions
{
    public:
        void Stretch(double rate, int FrameNum, double ***spectrogram, int new_FrameNum, double ***new_spectrogram, bool reset_ph);
};

#endif
