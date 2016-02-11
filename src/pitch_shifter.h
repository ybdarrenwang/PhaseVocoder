#ifndef __PITCHSHIFTER_H__
#define __PITCHSHIFTER_H__

#include "vocoder_functions.h"
#include "frame.h"

using namespace std;

/**
 * Resample complex spectrum along frequency axis while taking care of phases
*/
class PitchShifter
{
    public:
        PitchShifter(int n, int s) : FFT_SIZE(n), FRAME_SHIFT(s) {
            cached_phase = vector<float>(n/2+1, 0.0);
            bin_shift_residual = vector<float>(n/2+1, 0.0);
            synth_freq_bin = vector<int>(n/2+1, 0);
            vocoder_func = new VocoderFunctions(n, s);
        }

        virtual void UpdatePhase(vector<float>& mag, vector<float> prev_phase, vector<float> next_phase, vector<float>& synth_ph, float factor);
        virtual void SynthesizeFrame(vector<float>& mag, vector<float>& ph, Frame *f);
        virtual void Shift(float factor, vector<Frame*>& input_spec, vector<Frame*>& output_spec, bool reset_phase);

    protected:
        int FFT_SIZE;
        int FRAME_SHIFT;
        vector<float> cached_phase; // the last phase spectrum from previous Shift execution
        vector<int> synth_freq_bin;
        vector<float> bin_shift_residual; // for interpolation
        VocoderFunctions* vocoder_func;
};

#endif
