#ifndef __TIMESTRETCHERPL_H__
#define __TIMESTRETCHERPL_H__

#include "time_stretcher.h"

using namespace std;

/**
* Phase-locked phase vocoder:
* 1. Gather channels around each peak value into 1 group
* 2. calculate the phase of each peak
* 3. calculate the phase of the other channels with the peak's phasor
*/
class TimeStretcherPL : public TimeStretcher
{
    public:
        TimeStretcherPL(int n) : TimeStretcher(n) {
            for (int i=0; i<FFT_SIZE/2+1; ++i)
                prev_subband.push_back(i);
        }
        virtual ~TimeStretcherPL() {}

        void SynthesizePhase(vector<float> mag, vector<float> prev_phase, vector<float> phase, vector<float>& synth_ph);

    protected:
        vector<int> prev_subband; // the sub-band information from previous frame
};

#endif
