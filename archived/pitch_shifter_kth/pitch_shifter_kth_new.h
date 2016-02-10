#ifndef __PITCHSHIFTERKTH2_H__
#define __PITCHSHIFTERKTH2_H__

#include "pitch_shifter_kth.h"

using namespace std;

class PitchShifterKTH2 : public PitchShifterKTH
{
    public:
        void PS_KTH_new(float **Spectrum, float PitchQuotient, int frame_shift, bool reset_ph, float *window_real, float *window_imag);
};

#endif
