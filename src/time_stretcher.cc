#include "time_stretcher.h"

void TimeStretcher::SynthesizePhase(vector<float> mag, vector<float> prev_phase, vector<float> phase, vector<float>& synth_ph) {
    for(int freq=0; freq<FFT_SIZE/2+1; ++freq) {
        phasor[freq] = phase[freq] - prev_phase[freq];
        synth_ph[freq] += phasor[freq];
        while(synth_ph[freq] >= PI) synth_ph[freq] -= 2.0 * PI;
        while(synth_ph[freq] < -1.0*PI) synth_ph[freq] += 2.0 * PI;
    }
}

void TimeStretcher::SynthesizeFrame(vector<float>& mag, vector<float>& ph, Frame* f){
    for(int freq=0; freq<FFT_SIZE/2+1; ++freq)
        f->setSpectrum(freq, Complex(mag[freq]*cos(ph[freq]), mag[freq]*sin(ph[freq])));
    for(int freq=FFT_SIZE/2+1; freq<FFT_SIZE; ++freq)
        f->setSpectrum(freq, f->getSpectrum(FFT_SIZE-freq));
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
            SynthesizePhase(mag, input_spec[next_frame_idx]->getPhase(), input_spec[prev_frame_idx]->getPhase(), ph);
        else // initialize the first frame
            if (reset_phase)
                ph = input_spec[0]->getPhase();
            else {
                ph = cached_phase;
                SynthesizePhase(mag, input_spec[0]->getPhase(), cached_phase, ph);
            }
        SynthesizeFrame(mag, ph, f);
        output_spec.push_back(f);
        sample_ptr += 1.0/rate;
    }

    for(int freq=0; freq<FFT_SIZE/2+1; ++freq)
        cached_phase[freq] = ph[freq];
}
