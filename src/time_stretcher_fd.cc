#include "time_stretcher_fd.h"

void TimeStretcherFD::UpdatePhase(const double& ts_factor, vector<double> mag, vector<double> prev_phase, vector<double> next_phase, vector<double>& synth_ph) {
    for(int freq_bin=0; freq_bin<FFT_SIZE/2+1; ++freq_bin)
        synth_ph[freq_bin] = fmod(synth_ph[freq_bin]+next_phase[freq_bin]-prev_phase[freq_bin], 2.0*PI);
}

void TimeStretcherFD::Stretch(const double& ts_factor, vector<Frame*>& input_spec, vector<Frame*>& output_spec, int &synthesis_frame_shift, bool reset_phase) {
    vector<double> mag, ph;
    Frame *f;
    double sample_ptr = 0.0; // the pointer to the old spectrum, where the new magnitude/phase should be synthesized from.
    while ((unsigned int)sample_ptr+1<input_spec.size())
    {
        int prev_frame_idx = (int)sample_ptr;
        double prev_frame_weight = 1-(sample_ptr-(int)sample_ptr);
        int next_frame_idx = prev_frame_idx+1;
        double next_frame_weight = 1-prev_frame_weight;

        mag = vectorWeightedSum(input_spec[prev_frame_idx]->getMagnitude(), input_spec[next_frame_idx]->getMagnitude(), prev_frame_weight, next_frame_weight);

        f = new Frame(FFT_SIZE);
        if (output_spec.size()>0)
            UpdatePhase(ts_factor, mag, input_spec[prev_frame_idx]->getPhase(), input_spec[next_frame_idx]->getPhase(), ph);
        else // initialize the first frame
        {
            if (reset_phase)
                ph = input_spec[0]->getPhase();
            else
            {
                ph = cached_phase;
                UpdatePhase(ts_factor, mag, cached_phase, input_spec[0]->getPhase(), ph);
            }
        }
        SynthesizeFrame(mag, ph, f);
        output_spec.push_back(f);

        sample_ptr += 1.0/ts_factor;
    }

    for(int freq_bin=0; freq_bin<FFT_SIZE/2+1; ++freq_bin)
        cached_phase[freq_bin] = ph[freq_bin];
}
