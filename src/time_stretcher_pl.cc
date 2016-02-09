#include "time_stretcher_pl.h"

/**
* Phase-locked phase vocoder:
* 1. Gather channels around each peak value into 1 group
* 2. calculate the phase of each peak
* 3. calculate the phase of the other channels with the peak's phasor
*/
void TimeStretcherPL::UpdatePhase(vector<float> mag, vector<float> prev_phase, vector<float> phase, vector<float>& synth_ph) {
    vector<int> subband = vocoder_func->groupChannel(mag);

    for(int freq=0; freq<FFT_SIZE/2+1; ++freq)
        if (freq==subband[freq])
            phasor[freq] = phase[freq] - prev_phase[prev_subband[freq]];

    for(int freq=0; freq<FFT_SIZE/2+1; ++freq)
    {
        synth_ph[freq] += phasor[subband[freq]];
        while(synth_ph[freq] >= PI) synth_ph[freq] -= 2.0 * PI;
        while(synth_ph[freq] < -1.0*PI) synth_ph[freq] += 2.0 * PI;
        prev_subband[freq] = subband[freq];
    }
}
