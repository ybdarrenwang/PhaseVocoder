#include "pitch_shifter_kth.h"

// Find Max. Frequency Bin
bool sort(pair<int, int> i, pair<int, int> j){
    return (i.first>j.first);
}

void PitchShifterKTH::UpdatePhase(vector<float>& mag, vector<float>& synth_ph, float factor) {
    for (int i=0; i<mag.size(); ++i)
        sorted_mag.push_back(pair<int,int>(mag[i], i));
    sort(sorted_mag.begin(), sorted_mag.end(), compare_mag_bin); // in descending order

    for (int order=0; order<SamplingRate/2.0/pitch; ++order) { // until reach theoretical number of pitch peaks in spectrum
        int i = sorted_mag[order].second;

        freq_bin_shift[i] = floor(i*factor) - i;
        phase_shift_residual[i] = i*factor-floor(i*factor);
        phasor[i] += (i*(factor-1))*(2.0*PI*FRAME_SHIFT)/FFT_SIZE;
        while(phasor[i] >= PI) phasor[i] -= 2.0 * PI;
        while(phasor[i] < -1.0*PI) phasor[i] += 2.0 * PI;

        synth_ph[i] += phasor[i];
        while (synth_ph[i] >= PI) synth_ph[i] -= 2.0 * PI;
        while (synth_ph[i] < -1.0*PI) synth_ph[i] += 2.0 * PI;
    }
}

void PitchShifterKTH::SynthesizeFrame(vector<float>& mag, vector<float>& ph, Frame *f) {
    vector<Complex> synth_spec(FFT_SIZE/2+1, Complex(0.0, 0.0));
    float energy=0, new_energy=0; // for energy preservation
    vector<bool> mark_phasor(FFT_SIZE/2+1, false);
    float max_amp;

    for (int order=0; order<SamplingRate/2.0/pitch; ++order) { // until reach theoretical number of pitch peaks in spectrum
        int max_bin = sorted_mag[order].second;

        // shift the entire window spectrum
        max_amp = mag[max_bin];
        for(j=0; j<FFT_SIZE/2+1; j++) {
            if (j-max_bin+FFT_SIZE/2 >= 0 && j-max_bin+FFT_SIZE/2 < FFT_SIZE) {
                if (j+freq_bin_shift[max_bin] >=0 && j+freq_bin_shift[max_bin] < FFT_SIZE/2) {
                    // shift to nearest bin
                    // synth_spec[0][j+freq_bin_shift[max_bin]] += max_amp*window[j-max_bin+FFT_SIZE/2]*cos(ph[j]+phasor[max_bin]);
                    // synth_spec[1][j+freq_bin_shift[max_bin]] += max_amp*window[j-max_bin+FFT_SIZE/2]*sin(ph[j]+phasor[max_bin]);
                    // interpolation
                    synth_spec[j+freq_bin_shift[max_bin]].real += (1-phase_shift_residual[max_bin])*max_amp*window[j-max_bin+FFT_SIZE/2]*cos(ph[j]+phasor[max_bin]);
                    synth_spec[j+freq_bin_shift[max_bin]+1].real += phase_shift_residual[max_bin]*max_amp*window[j-max_bin+FFT_SIZE/2]*cos(ph[j]+phasor[max_bin]);
                    synth_spec[j+freq_bin_shift[max_bin]].imag += (1-phase_shift_residual[max_bin])*max_amp*window[j-max_bin+FFT_SIZE/2]*sin(ph[j]+phasor[max_bin]);
                    synth_spec[j+freq_bin_shift[max_bin]+1].imag += phase_shift_residual[max_bin]*max_amp*window[j-max_bin+FFT_SIZE/2]*sin(ph[j]+phasor[max_bin]);
                }
                mag[j] -= max_amp*window[j-max_bin+FFT_SIZE/2];
            }
        }
    }

    // allocate region-of-influence for sinusoid freq. tracking
    // for (i=0; i<FFT_SIZE/2+1; i++)
    //     MagTemp[i] = ABS2(Spectrum[0][i],Spectrum[1][i]);
    // ChannelGrouping(MagTemp, prev_subband);

    // energy preservation; cumulate the phasors
    for (int i=0; i<FFT_SIZE/2+1; i++)
        new_energy += synth_spec[i].getEnergy();
    float amplitude_norm = sqrt(energy/new_energy);
    for (int i=0; i<FFT_SIZE/2+1; i++) {
        synth_spec[i].real *= amplitude_norm;
        synth_spec[i].imag *= amplitude_norm;
    }
}
