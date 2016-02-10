#ifndef __PITCHSHIFTERPL_H__
#define __PITCHSHIFTERPL_H__

#include "pitch_shifter.h"

using namespace std;

class PitchShifterPL : public PitchShifter
{
    public:
        PitchShifterPL(int n, int s) : PitchShifter(n, s) {
            for (int i=0; i<n/2+1; ++i)
                prev_subband.push_back(i);
        }

        void UpdatePhase(vector<float>& mag, vector<float>& synth_ph, float factor);

    protected:
        vector<int> prev_subband;
};

#endif
