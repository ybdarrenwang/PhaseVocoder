#include "time_stretcher.h"

void TimeStretcher::UpdatePhase(const double& ts_factor, vector<double> mag, vector<double> prev_phase, vector<double> next_phase, vector<double>& synth_ph) {
    for(int freq_bin=0; freq_bin<FFT_SIZE/2+1; ++freq_bin)
        synth_ph[freq_bin] = fmod(synth_ph[freq_bin]+(next_phase[freq_bin]-prev_phase[freq_bin])*ts_factor, 2.0*PI);
}

void TimeStretcher::SynthesizeFrame(vector<double>& mag, vector<double>& ph, Frame* f){
    for(int freq_bin=0; freq_bin<FFT_SIZE/2+1; ++freq_bin)
        f->setSpectrum(freq_bin, Complex(mag[freq_bin]*cos(ph[freq_bin]), mag[freq_bin]*sin(ph[freq_bin])));
    for(int freq_bin=FFT_SIZE/2+1; freq_bin<FFT_SIZE; ++freq_bin)
        f->setSpectrum(freq_bin, f->getSpectrum(FFT_SIZE-freq_bin).getConjugate());
}

void TimeStretcher::Stretch(const double& ts_factor, vector<Frame*>& input_spec, vector<Frame*>& output_spec, int &synthesis_frame_shift, bool reset_phase) {
    vector<double> mag, ph;
    Frame *f;
    for (unsigned frame_idx=0; frame_idx<input_spec.size(); ++frame_idx) {
        mag = input_spec[frame_idx]->getMagnitude();
        f = new Frame(FFT_SIZE);
        if (frame_idx > 0)
            UpdatePhase(ts_factor, mag, input_spec[frame_idx-1]->getPhase(), input_spec[frame_idx]->getPhase(), ph);
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
    }

    for(int freq_bin=0; freq_bin<FFT_SIZE/2+1; ++freq_bin)
        cached_phase[freq_bin] = ph[freq_bin];

    synthesis_frame_shift*=ts_factor;
}
