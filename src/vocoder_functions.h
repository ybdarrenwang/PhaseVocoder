#ifndef __VOCODERFUNCTIONS_H__
#define __VOCODERFUNCTIONS_H__

#include "global.h"

using namespace std;

class VocoderFunctions {
    public:
        float ABS2(float a,float b);
        vector<float> vectorWeightedSum(vector<float> v1, vector<float> v2, float w1, float w2);
        vector<int> groupChannel(vector<float>& spec);
        //void new_ChannelGrouping(float *Spec, int *ChannelGroupFlag, float pitch);
};

#endif
