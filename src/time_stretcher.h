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
        TimeStretcher(double f, int n, int s) : ts_factor(f), FFT_SIZE(n), FRAME_SHIFT(s) {
            cached_phase = vector<double>(n/2+1, 0.0);
        }
        virtual ~TimeStretcher() {}

        virtual void UpdatePhase(vector<double> mag, vector<double> prev_phase, vector<double> phase, vector<double>& synth_ph);
        virtual void SynthesizeFrame(vector<double>&, vector<double>&, Frame*);
        virtual void Stretch(vector<Frame*>& input_spec, vector<Frame*>& output_spec, int &synthesis_frame_shift, bool reset_phase);

    protected:
        double ts_factor;
        int FFT_SIZE;
        int FRAME_SHIFT;
        vector<double> cached_phase; // the last phase spectrum from previous Stretch execution
};

#endif
