#include "pitch_shifter_pl.h"

/**
* Phase-locked phase vocoder:
* 1. Gather channels around each peak value into 1 group
* 2. calculate the phase of each peak
* 3. calculate the phase of the other channels with the peak's phasor
* 
* Note: so far sounds far worse than ordinary pitch shifting; to find out why
*/
void PitchShifterPL::UpdatePhase(vector<float>& mag, vector<float>& synth_ph, float factor) {
    // allocate region-of-influence for sinusoid freq. tracking
    vector<int> subband = vocoder_func->groupChannel(mag);

    // calculate the phasor and bin number to shift for each peak
    for(int i=0;i<FFT_SIZE/2+1;i++)
        if (i==subband[i]) {
            freq_bin_shift[i] = floor(i*factor) - i;
            phase_shift_residual[i] = i*factor-floor(i*factor);
            phasor[i] = phasor[prev_subband[i]] + (i*(factor-1))*(2.0*PI*FRAME_SHIFT)/FFT_SIZE;
            while(phasor[i] >= PI) phasor[i] -= 2.0 * PI;
            while(phasor[i] < -1.0*PI) phasor[i] += 2.0 * PI;
        }

    // update phase and other cached info
    for(int i=0;i<FFT_SIZE/2+1;i++)
    {
        phasor[i] = phasor[subband[i]];
        phase_shift_residual[i] = phase_shift_residual[subband[i]];
        freq_bin_shift[i] = freq_bin_shift[subband[i]];

        synth_ph[i] += phasor[i];
        while (synth_ph[i] >= PI) synth_ph[i] -= 2.0 * PI;
        while (synth_ph[i] < -1.0*PI) synth_ph[i] += 2.0 * PI;

        prev_subband[i] = subband[i];
    }
}
