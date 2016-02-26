#ifndef __TIMESTRETCHER_H__
#define __TIMESTRETCHER_H__

#include "vocoder_functions.h"
#include "frame.h"

using namespace std;

class TimeStretcher
{
    public:
        TimeStretcher(int n, int s) : FFT_SIZE(n), FRAME_SHIFT(s) {
            cached_phase = vector<double>(n/2+1, 0.0);
            vocoder_func = new VocoderFunctions(n, s);
        }
        virtual ~TimeStretcher() { delete vocoder_func; }

        virtual void UpdatePhase(vector<double> mag, vector<double> prev_phase, vector<double> phase, vector<double>& synth_ph);
        void SynthesizeFrame(vector<double>&, vector<double>&, Frame*);
        void Stretch(double rate, vector<Frame*>& input_spec, vector<Frame*>& output_spec, bool reset_phase);

    protected:
        int FFT_SIZE;
        int FRAME_SHIFT;
        vector<double> cached_phase; // the last phase spectrum from previous Stretch execution
        VocoderFunctions* vocoder_func;
};

#endif
