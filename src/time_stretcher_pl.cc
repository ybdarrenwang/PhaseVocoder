#include "time_stretcher_pl.h"

void TimeStretcherPL::UpdatePhase(const double& ts_factor, vector<double> mag, vector<double> prev_phase, vector<double> next_phase, vector<double>& synth_ph) {
    vector<int> local_peaks = getLocalPeaks(mag);

    for(int freq_bin=0; freq_bin<FFT_SIZE/2+1; ++freq_bin)
        if (freq_bin==local_peaks[freq_bin])
        {
            double delta_phase = (next_phase[freq_bin]-prev_phase[prev_local_peaks[freq_bin]])*ts_factor;
            synth_ph[freq_bin] = fmod(synth_ph[freq_bin]+unwrapPhase(delta_phase, freq_bin, FRAME_SHIFT, FFT_SIZE), 2.0*PI);
        }

    for(int freq_bin=0; freq_bin<FFT_SIZE/2+1; ++freq_bin)
        if (freq_bin!=local_peaks[freq_bin])
        {
            double delta_phase = next_phase[freq_bin]-next_phase[local_peaks[freq_bin]];
            synth_ph[freq_bin] = fmod(synth_ph[local_peaks[freq_bin]]+delta_phase, 2.0*PI);
        }

    prev_local_peaks = local_peaks;
}
