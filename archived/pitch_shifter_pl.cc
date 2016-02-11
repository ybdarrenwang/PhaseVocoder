#include "pitch_shifter_pl.h"

void PitchShifterPL::UpdatePhase(vector<float>& mag, vector<float>& synth_ph, float factor) {
    // allocate region-of-influence for sinusoid freq. tracking
    vector<int> subband = vocoder_func->groupChannel(mag);

    // calculate the phasor and bin number to shift for each peak
    for(int i=0;i<FFT_SIZE/2+1;i++)
        if (i==subband[i]) {
            freq_bin_shift[i] = floor(i*factor) - i;
            phase_shift_residual[i] = i*factor-floor(i*factor);
            phasor[i] = fmod(phasor[prev_subband[i]]+(i*(factor-1))*(2.0*PI*FRAME_SHIFT)/FFT_SIZE, 2.0*PI);
        }

    // update phase and other cached info
    for(int i=0;i<FFT_SIZE/2+1;i++)
    {
        phasor[i] = phasor[subband[i]];
        phase_shift_residual[i] = phase_shift_residual[subband[i]];
        freq_bin_shift[i] = freq_bin_shift[subband[i]];
        synth_ph[i] = fmod(synth_ph[i]+phasor[i], 2.0*PI);
        prev_subband[i] = subband[i];
    }
}
