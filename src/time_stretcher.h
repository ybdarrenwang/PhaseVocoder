#ifndef __TIMESTRETCHER_H__
#define __TIMESTRETCHER_H__

#include "vocoder_functions.h"
#include "frame.h"

using namespace std;

class TimeStretcher : public VocoderFunctions
{
    public:
        TimeStretcher(int n) : VocoderFunctions(n) {
            prev_subband = vector<int>(n/2+1,0);
            cached_phase = vector<float>(n/2+1,0);
            phasor = vector<float>(n/2+1,0);
        }
        virtual ~TimeStretcher() {}
        void SynthesizePhase(vector<float> mag, vector<float> prev_phase, vector<float> phase, vector<float>& synth_ph);
        void SynthesizeFrame(vector<float>&, vector<float>&, Frame*);
        void Stretch(float rate, vector<Frame*>& input_spec, vector<Frame*>& output_spec, bool reset_phase);

    private:
        vector<int> prev_subband; // the sub-band information from previous frame
        vector<float> cached_phase; // the last phase spectrum from previous Stretch execution
        vector<float> phasor; // the delta term for the phase of every frequency bin
};

#endif
