#ifndef TIMESTRETCHER_FD_PL_H
#define TIMESTRETCHER_FD_PL_H

#include "time_stretcher_fd.h"

using namespace std;

/**
 * Interpolate/extrapolate on spectrogram, and perform phase locking as well
*/
class TimeStretcherFDPL : public TimeStretcherFD
{
    public:
        TimeStretcherFDPL(int n, int s) : TimeStretcherFD(n, s) {
            for (int i=0; i<FFT_SIZE/2+1; ++i)
                prev_local_peaks.push_back(i);
        }
        virtual ~TimeStretcherFDPL(){}
        virtual void UpdatePhase(const double& ts_factor, vector<double> mag, vector<double> prev_phase, vector<double> next_phase, vector<double>& synth_ph);

    protected:
        vector<int> prev_local_peaks; // the sub-band information from previous frame
};

#endif
