#include "util.h"

void resample(short* data, const int& len, const double& rate, vector<short>& new_data)
{
    double idx = 0;
    while (idx<len)
    {
        new_data.push_back(data[(int)idx]);
        idx += rate;
    }
}

double unwrapPhase(double delta_phase, const int& freq_bin, const int& FRAME_SHIFT, const int& FFT_SIZE) {
    delta_phase -= 2*PI*freq_bin*FRAME_SHIFT/FFT_SIZE;
    delta_phase = fmod(delta_phase, 2.0*PI);
    return delta_phase + 2*PI*freq_bin*FRAME_SHIFT/FFT_SIZE;
}

vector<double> vectorWeightedSum(const vector<double> &v1, const vector<double> &v2, double w1, double w2) {
    if (w2==0) return v1;
    if (w1==0) return v2;
    vector<double> ans;
    if (v1.size()!=v2.size()) return ans;
    for (unsigned int i=0; i<v1.size(); ++i)
        ans.push_back(v1[i]*w1+v2[i]*w2);
    return ans;
}

/**
 * Input:  a double vector
 *         e.g. [0,1,2,1,0,3,6,5,0]
 * Output: an int vector indication the corresponding peak index of each bin;
 *         e.g. [2,2,2,2,7,7,7,7,7]
 * Definition:
 * - peak: larger then preceding 2 and following 2 elements
 * - valley: the minimum between 2 peaks
 */
vector<int> getLocalPeaks(vector<double>& spec) {
    vector<int> peak_idx;
    int peak_ptr = 0;
    int valley_ptr = 0;
    peak_idx.push_back(0);
    peak_idx.push_back(0);
    for (unsigned int i=2; i<spec.size()-2; ++i) {
        if (spec[i] > spec[i-1] && spec[i] > spec[i-2] && spec[i] > spec[i+1] && spec[i] > spec[i+2]) {
            // peak: overwrite from valley to current position
            for (unsigned int j=valley_ptr; j<i; ++j)
                peak_idx[j] = i;
            peak_ptr = i;
            valley_ptr = i;
        }
        else // otherwise: update valley if necessary
            if (spec[i]<spec[valley_ptr])
                valley_ptr = i;
        peak_idx.push_back(peak_ptr);
    }
    peak_idx.push_back(peak_ptr);
    peak_idx.push_back(peak_ptr);
    return peak_idx;
}
/*
void new_ChannelGrouping(double *Spec, int *ChannelGroupFlag, double pitch) {
    int i,j;
    int GroupBoundary, PrevGroupBoundary=0, PrevPeak=0;
    double MinTemp;

    // determine how much bins to look forward/backward
    double bin_period = FFT_SIZE/(SamplingRate/pitch); // f0 interval width counted by bin number
    int find_max_range = (int)floor(bin_period*0.5);
    bool find_max;

    // determine region of influence
    for(i=find_max_range; i<FFT_SIZE/2+1-find_max_range; i++)
    {
        // find local mixima
        find_max = true;
        for (j=i-find_max_range; j<i+find_max_range; j++)
        {
            if (Spec[i]<Spec[j])
                find_max = false;
        }

        // determine region of influence
        if (find_max == true)
        {
            if (PrevPeak > 0)
            {
                // find minima between present and previous peak
                MinTemp = Spec[PrevPeak];
                for (j=PrevPeak+1; j<i; j++)
                {
                    if (Spec[j] < MinTemp)
                    {
                        MinTemp = Spec[j];
                        GroupBoundary = j;
                    }
                }

                // annotate the previous group around the previous peak
                for (j = PrevGroupBoundary; j<GroupBoundary; j++)
                    ChannelGroupFlag[j] = PrevPeak;

                PrevGroupBoundary = GroupBoundary;
            }
            else // annotate the first group
            {
                for (j=0; j<i; j++)
                    ChannelGroupFlag[j] = i;
            }
            PrevPeak = i;
        }
    }

    // annotate the last group
    for (j = PrevGroupBoundary; j<FFT_SIZE/2+1; j++)
        ChannelGroupFlag[j] = PrevPeak;
}*/
