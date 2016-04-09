#ifndef TIMESTRETCHER_H
#define TIMESTRETCHER_H

#include "util.h"
#include "frame.h"

using namespace std;

/**
 * Common approach for time stretching:
 * - number of frames remain the same
 * - the output synthesis-rate/frame-shift/hop-size is modified correspondingly
*/
class TimeStretcher
{
    public:
        TimeStretcher(int n, int s) : FFT_SIZE(n), FRAME_SHIFT(s) {
            cached_phase = vector<double>(n/2+1, 0.0);
        }
        virtual ~TimeStretcher() {}

        virtual void UpdatePhase(const double& ts_factor, vector<double> mag, vector<double> prev_phase, vector<double> phase, vector<double>& synth_ph);
        virtual void SynthesizeFrame(vector<double>&, vector<double>&, Frame*);
        virtual void Stretch(const double& ts_factor, vector<Frame*>& input_spec, vector<Frame*>& output_spec, int &synthesis_frame_shift, bool reset_phase);

    protected:
        int FFT_SIZE;
        int FRAME_SHIFT;
        vector<double> cached_phase; // the last phase spectrum from previous Stretch execution
};

#endif
