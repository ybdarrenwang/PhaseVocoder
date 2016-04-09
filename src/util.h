#ifndef VOCODERFUNCTIONS_H
#define VOCODERFUNCTIONS_H

#include "global.h"

using namespace std;

void resample(short* data, const int& len, const double& rate, vector<short>& new_data);
double unwrapPhase(double delta_phase, const int& freq_bin, const int& FRAME_SHIFT, const int& FFT_SIZE);
vector<double> vectorWeightedSum(const vector<double> &, const vector<double> &, double, double);
vector<int> getLocalPeaks(vector<double> &spec);
//void new_ChannelGrouping(double *Spec, int *ChannelGroupFlag, double pitch);

#endif
