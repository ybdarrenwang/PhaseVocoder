#ifndef __TIMESTRETCHER_H__
#define __TIMESTRETCHER_H__

#include "vocoder_functions.h"
#include "frame.h"

using namespace std;

class TimeStretcher
{
    public:
        TimeStretcher(int n, int s) : FFT_SIZE(n), FRAME_SHIFT(s) {
            cached_phase = vector<float>(n/2+1, 0.0);
            vocoder_func = new VocoderFunctions(n, s);
        }
        virtual ~TimeStretcher() { delete vocoder_func; }

        virtual void UpdatePhase(vector<float> mag, vector<float> prev_phase, vector<float> phase, vector<float>& synth_ph);
        void SynthesizeFrame(vector<float>&, vector<float>&, Frame*);
        void Stretch(float rate, vector<Frame*>& input_spec, vector<Frame*>& output_spec, bool reset_phase);

    protected:
        int FFT_SIZE;
        int FRAME_SHIFT;
        vector<float> cached_phase; // the last phase spectrum from previous Stretch execution
        VocoderFunctions* vocoder_func;
};

#endif
