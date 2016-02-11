#ifndef __VOCODERFUNCTIONS_H__
#define __VOCODERFUNCTIONS_H__

#include "global.h"

using namespace std;

class VocoderFunctions {
    public:
        VocoderFunctions(int n, int s) : FFT_SIZE(n), FRAME_SHIFT(s) {}
        float ABS2(float a,float b);
        float phaseUnwrapping(float delta_phase, int freq_bin);
        vector<float> vectorWeightedSum(vector<float> v1, vector<float> v2, float w1, float w2);
        vector<int> groupChannel(vector<float>& spec);
        //void new_ChannelGrouping(float *Spec, int *ChannelGroupFlag, float pitch);

    private:
        int FFT_SIZE;
        int FRAME_SHIFT;
};

#endif
