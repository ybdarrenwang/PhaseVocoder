#ifndef __VOCODERFUNCTIONS_H__
#define __VOCODERFUNCTIONS_H__

#include "global.h"

using namespace std;

class VocoderFunctions {
    public:
        VocoderFunctions(int n, int s) : FFT_SIZE(n), FRAME_SHIFT(s) {}
        double ABS2(double a,double b);
        double unwrapPhase(double delta_phase, int freq_bin);
        vector<double> vectorWeightedSum(const vector<double> &, const vector<double> &, double, double);
        vector<int> getLocalPeaks(vector<double> &spec);
        //void new_ChannelGrouping(double *Spec, int *ChannelGroupFlag, double pitch);

    private:
        int FFT_SIZE;
        int FRAME_SHIFT;
};

#endif
