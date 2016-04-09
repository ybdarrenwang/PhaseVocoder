#ifndef TIMESTRETCHERPL_H
#define TIMESTRETCHERPL_H

#include "time_stretcher.h"

using namespace std;

/**
* Phase-locked phase vocoder:
* 1. Gather channels around each peak value into 1 group
* 2. calculate the phase of each peak
* 3. calculate the phase of the other channels by maintaining the difference between their phases to peak's phasor
* 
* Note: sounds worse than ordinary time stretching; to find out why
*/
class TimeStretcherPL : public TimeStretcher
{
    public:
        TimeStretcherPL(int n, int s) : TimeStretcher(n, s) {
            for (int i=0; i<FFT_SIZE/2+1; ++i)
                prev_local_peaks.push_back(i);
        }
        virtual ~TimeStretcherPL() {}
        virtual void UpdatePhase(const double& ts_factor, vector<double> mag, vector<double> prev_phase, vector<double> next_phase, vector<double>& synth_ph);

    protected:
        vector<int> prev_local_peaks; // the sub-band information from previous frame
};

#endif
