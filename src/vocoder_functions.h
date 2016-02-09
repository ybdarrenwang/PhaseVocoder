#ifndef __VOCODERFUNCTIONS_H__
#define __VOCODERFUNCTIONS_H__

#include "global.h"

using namespace std;

class VocoderFunctions {
    public:
        VocoderFunctions(int n) : FFT_SIZE(n) {}
        virtual ~VocoderFunctions(){}
        vector<float> vectorWeightedSum(vector<float> v1, vector<float> v2, float w1, float w2);
        vector<int> groupChannel(vector<float>& spec);
        //void new_ChannelGrouping(float *Spec, int *ChannelGroupFlag, float pitch);
        float ABS2(float a,float b);

    protected:
        int FFT_SIZE;
        //float StudentPha[FFT_SIZE/2+1];
        //float StudentMag[FFT_SIZE/2+1];
        //int* voiced_bin_th;
};

#endif
