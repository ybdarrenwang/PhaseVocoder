#include "time_stretcher.h"

void TimeStretcher::UpdatePhase(vector<float> mag, vector<float> prev_phase, vector<float> next_phase, vector<float>& synth_ph) {
    for(int freq_bin=0; freq_bin<FFT_SIZE/2+1; ++freq_bin) {
        phasor[freq_bin] = next_phase[freq_bin] - prev_phase[freq_bin];
        synth_ph[freq_bin] = fmod(synth_ph[freq_bin]+phasor[freq_bin], 2.0*PI);
    }
}

void TimeStretcher::SynthesizeFrame(vector<float>& mag, vector<float>& ph, Frame* f){
    for(int freq_bin=0; freq_bin<FFT_SIZE/2+1; ++freq_bin)
        f->setSpectrum(freq_bin, Complex(mag[freq_bin]*cos(ph[freq_bin]), mag[freq_bin]*sin(ph[freq_bin])));
    for(int freq_bin=FFT_SIZE/2+1; freq_bin<FFT_SIZE; ++freq_bin)
        f->setSpectrum(freq_bin, f->getSpectrum(FFT_SIZE-freq_bin).getConjugate());
}

void TimeStretcher::Stretch(float rate, vector<Frame*>& input_spec, vector<Frame*>& output_spec, bool reset_phase) {
    vector<float> mag, ph;
    Frame *f;
    float sample_ptr = 0.0; // the pointer to the old spectrum, where the new magnitude/phase should be synthesized from.
    while ((int)sample_ptr+1<input_spec.size()) {
        int prev_frame_idx = (int)sample_ptr;
        float prev_frame_weight = 1-(sample_ptr-(int)sample_ptr);
        int next_frame_idx = prev_frame_idx+1;
        float next_frame_weight = 1-prev_frame_weight;

        mag = vocoder_func->vectorWeightedSum(input_spec[prev_frame_idx]->getMagnitude(), input_spec[next_frame_idx]->getMagnitude(), prev_frame_weight, next_frame_weight);

        f = new Frame(FFT_SIZE);
        if (output_spec.size()>0)
            UpdatePhase(mag, input_spec[prev_frame_idx]->getPhase(), input_spec[next_frame_idx]->getPhase(), ph);
        else // initialize the first frame
            if (reset_phase)
                ph = input_spec[0]->getPhase();
            else {
                ph = cached_phase;
                UpdatePhase(mag, cached_phase, input_spec[0]->getPhase(), ph);
            }
        SynthesizeFrame(mag, ph, f);
        output_spec.push_back(f);

        sample_ptr += 1.0/rate;
    }

    for(int freq_bin=0; freq_bin<FFT_SIZE/2+1; ++freq_bin)
        cached_phase[freq_bin] = ph[freq_bin];
}
