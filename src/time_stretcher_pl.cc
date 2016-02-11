#include "time_stretcher_pl.h"

void TimeStretcherPL::UpdatePhase(vector<float> mag, vector<float> prev_phase, vector<float> next_phase, vector<float>& synth_ph) {
    vector<int> subband = vocoder_func->groupChannel(mag);

    for(int freq_bin=0; freq_bin<FFT_SIZE/2+1; ++freq_bin)
        if (freq_bin==subband[freq_bin]) {
            //float delta_phase = next_phase[freq_bin]-prev_phase[prev_subband[freq_bin]];
            float delta_phase = next_phase[freq_bin]-prev_phase[freq_bin];
            synth_ph[freq_bin] = fmod(synth_ph[freq_bin]+vocoder_func->phaseUnwrapping(delta_phase, freq_bin), 2.0*PI);
    }

    for(int freq_bin=0; freq_bin<FFT_SIZE/2+1; ++freq_bin)
        if (freq_bin!=subband[freq_bin])
            synth_ph[freq_bin] = fmod(synth_ph[subband[freq_bin]]+next_phase[freq_bin]-next_phase[subband[freq_bin]], 2.0*PI);

    prev_subband = subband;
}
