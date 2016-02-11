#include "time_stretcher_pl.h"

void TimeStretcherPL::UpdatePhase(vector<float> mag, vector<float> prev_phase, vector<float> next_phase, vector<float>& synth_ph) {
    vector<int> subband = vocoder_func->groupChannel(mag);

    for(int freq_bin=0; freq_bin<FFT_SIZE/2+1; ++freq_bin)
        if (freq_bin==subband[freq_bin]) {
            phasor[freq_bin] = next_phase[freq_bin]-prev_phase[freq_bin];
            //phasor[freq_bin] = (next_phase[freq_bin]-prev_phase[prev_subband[freq_bin]]) + (freq_bin-prev_subband[freq_bin])*2.0*PI/FFT_SIZE;
    }

    for(int freq_bin=0; freq_bin<FFT_SIZE/2+1; ++freq_bin)
    {
        synth_ph[freq_bin] += phasor[subband[freq_bin]];
        while(synth_ph[freq_bin] >= PI) synth_ph[freq_bin] -= 2.0 * PI;
        while(synth_ph[freq_bin] < -1.0*PI) synth_ph[freq_bin] += 2.0 * PI;
        prev_subband[freq_bin] = subband[freq_bin];
    }
}
