#ifndef __PITCHSHIFTERPL_H__
#define __PITCHSHIFTERPL_H__

#include "pitch_shifter.h"

using namespace std;

/**
* Phase-locked phase vocoder:
* 1. Gather channels around each peak value into 1 group
* 2. calculate the phase of each peak
* 3. calculate the phase of the other channels with the peak's phasor
* 
* Note: so far sounds far worse than ordinary pitch shifting; to find out why
*/
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
