#include "pitch_shifter.h"

void PitchShifter::UpdatePhase(vector<float>& mag, vector<float>& synth_ph, float factor) {
    // update phase and other cached info
    for(int i=0;i<FFT_SIZE/2+1;i++)
    {
        freq_bin_shift[i] = floor(i*factor) - i;
        phase_shift_residual[i] = i*factor-floor(i*factor);
        phasor[i] = fmod(phasor[i]+(i*(factor-1))*(2.0*PI*FRAME_SHIFT)/FFT_SIZE, 2.0*PI);
        synth_ph[i] = fmod(synth_ph[i]+phasor[i], 2.0*PI);
    }
}

void PitchShifter::SynthesizeFrame(vector<float>& mag, vector<float>& ph, Frame *f) {
    vector<Complex> synth_spec(FFT_SIZE/2+1, Complex(0.0, 0.0));
    float energy=0, new_energy=0; // for energy preservation

    // complex number interpolation (seems problematic?)
    for(int i=0;i<FFT_SIZE/2+1;i++) {
        energy += mag[i]*mag[i];
        if (i+freq_bin_shift[i]>=0 && i+freq_bin_shift[i]<FFT_SIZE/2+1) {
            synth_spec[i+freq_bin_shift[i]].real += (1-phase_shift_residual[i])*mag[i]*cos(ph[i]);
            synth_spec[i+freq_bin_shift[i]].imag += (1-phase_shift_residual[i])*mag[i]*sin(ph[i]);
        }
        if (i+freq_bin_shift[i]+1>=0 && i+freq_bin_shift[i]+1<FFT_SIZE/2+1) {
            synth_spec[i+freq_bin_shift[i]+1].real += phase_shift_residual[i]*mag[i]*cos(ph[i]);
            synth_spec[i+freq_bin_shift[i]+1].imag += phase_shift_residual[i]*mag[i]*sin(ph[i]);
        }
    }

    // energy preservation
    for (int i=0; i<FFT_SIZE/2+1; i++)
        new_energy += synth_spec[i].getEnergy();
    float amplitude_norm = sqrt(energy/new_energy);
    for (int i=0; i<FFT_SIZE/2+1; i++) {
        synth_spec[i].real *= amplitude_norm;
        synth_spec[i].imag *= amplitude_norm;
    }

    // output
    f->setSpectrum(&synth_spec[0]);
}

// ============================================================================ //
// < Resampling on Frequency Dimension >
// ============================================================================ //
void PitchShifter::Shift(float factor, vector<Frame*>& input_spec, vector<Frame*>& output_spec, bool reset_phase) {
    vector<int> subband;
    vector<float> mag, ph;
    Frame *f;

    if (reset_phase)
        for (int i=0; i<FFT_SIZE/2+1; i++) {
            phasor[i] = 0;
            //prev_subband[i] = i;
        }

    for (int sample_idx=0; sample_idx<input_spec.size(); ++sample_idx) {
        mag = input_spec[sample_idx]->getMagnitude();
        ph = input_spec[sample_idx]->getPhase();
        UpdatePhase(mag, ph, factor);
        f = new Frame(FFT_SIZE);
        SynthesizeFrame(mag, ph, f);
        output_spec.push_back(f);
    }
}
