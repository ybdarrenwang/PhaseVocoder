#include "pitch_shifter.h"

void PitchShifter::UpdatePhase(vector<double>& mag, vector<double> prev_phase, vector<double> next_phase, vector<double>& synth_ph, double factor) {
    for(int i=0;i<FFT_SIZE/2+1;i++) {
        synth_freq_bin[i] = floor(i*factor);
        bin_shift_residual[i] = i*factor-synth_freq_bin[i];
        synth_ph[i] = fmod(synth_ph[i]+factor*vocoder_func->unwrapPhase(next_phase[i]-prev_phase[i], i), 2.0*PI);
    }
}

void PitchShifter::SynthesizeFrame(vector<double>& mag, vector<double>& ph, Frame *f) {
    vector<Complex> synth_spec(FFT_SIZE/2+1, Complex(0.0, 0.0));
    double energy=0, new_energy=0; // for energy preservation

    // interpolate complex spectrum along frequency axis
    for(int i=0;i<FFT_SIZE/2+1;i++) {
        energy += mag[i]*mag[i];
        if (synth_freq_bin[i]>=0 && synth_freq_bin[i]<FFT_SIZE/2+1) {
            synth_spec[synth_freq_bin[i]].real += (1-bin_shift_residual[i])*mag[i]*cos(ph[i]);
            synth_spec[synth_freq_bin[i]].imag += (1-bin_shift_residual[i])*mag[i]*sin(ph[i]);
        }
        if (synth_freq_bin[i]+1>=0 && synth_freq_bin[i]+1<FFT_SIZE/2+1) {
            synth_spec[synth_freq_bin[i]+1].real += bin_shift_residual[i]*mag[i]*cos(ph[i]);
            synth_spec[synth_freq_bin[i]+1].imag += bin_shift_residual[i]*mag[i]*sin(ph[i]);
        }
    }

    // energy preservation
    for (int i=0; i<FFT_SIZE/2+1; i++)
        new_energy += synth_spec[i].getEnergy();
    double amplitude_norm = sqrt(energy/new_energy);
    for (int i=0; i<FFT_SIZE/2+1; i++) {
        synth_spec[i].real *= amplitude_norm;
        synth_spec[i].imag *= amplitude_norm;
    }

    // output
    f->setSpectrum(&synth_spec[0]);
}

void PitchShifter::Shift(double factor, vector<Frame*>& input_spec, vector<Frame*>& output_spec, bool reset_phase) {
    vector<int> subband;
    vector<double> mag, ph;
    Frame *f;

    for (int sample_idx=0; sample_idx<input_spec.size(); ++sample_idx) {
        mag = input_spec[sample_idx]->getMagnitude();
        if (sample_idx)
            UpdatePhase(mag, input_spec[sample_idx-1]->getPhase(), input_spec[sample_idx]->getPhase(), ph, factor);
        else {
            if (reset_phase) {
                ph = input_spec[0]->getPhase();
                for (int i=0; i<FFT_SIZE/2+1; i++)
                    synth_freq_bin[i] = (int)(i*factor);
            }
            else {
                ph = cached_phase;
                UpdatePhase(mag, cached_phase, input_spec[0]->getPhase(), ph, factor);
            }
        }
        f = new Frame(FFT_SIZE);
        SynthesizeFrame(mag, ph, f);
        output_spec.push_back(f);
    }
    cached_phase = ph;
}
