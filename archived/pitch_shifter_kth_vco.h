#ifndef __PITCHSHIFTERVCO_H__
#define __PITCHSHIFTERVCO_H__

#include "pitch_shifter_kth.h"

using namespace std;

class PitchShifterVCO : public PitchShifterKTH
{
    public:
        void PitchShifting_KTH_HNM(float **Spectrum, float PitchQuotient, int frame_shift, bool reset_ph, float *window_mag, float *window_pha, float pitch, int vowel_frame_index);
};

#endif
