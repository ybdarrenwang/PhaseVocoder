#ifndef __PITCHSHIFTERKTH_H__
#define __PITCHSHIFTERKTH_H__

#include "pitch_shifter.h"

using namespace std;

/**
 * Pitch Shifting proposed by KTH
 * 1. assume known pitch
 * 2. iteratively subtract window spectrum from speech spectrum
 */
class PitchShifterKTH : public PitchShifter
{
    public:
        virtual void UpdatePhase(vector<float>& mag, vector<float>& synth_ph, float factor);
        virtual void SynthesizeFrame(vector<float>& mag, vector<float>& ph, Frame *f);

    protected:
        vector< pair<int, int> > sorted_mag; // (magnitude, bin_index)
};

#endif
