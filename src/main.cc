#include "time_stretcher.h"
#include "time_stretcher_pl.h"
#include "pitch_shifter.h"
#include "pitch_shifter_pl.h"
#include "complex.h"
#include "window.h"
#include "hamming_window.h"
#include "frame.h"
#include "wav_io.h"
#include <iostream>

using namespace std;

int main(int argc, char **argv) {
    // Initialize parameters and functions
    cout<<"Initialize parameters and functions"<<endl;
    string input_file = "test/HungarianDanceNo5.wav";
    string output_file = "test/tmp_pl.wav";
    int FRAME_LENGTH = 4096;
    int FRAME_SHIFT = 1024;
    int FFT_SIZE = FRAME_LENGTH;
    float ts_rate = 2;
    //float ps_rate = 1.125;
    float ps_rate = 1;

    Window *window = new HammingWindow(FRAME_LENGTH);
    MyFFT *fft = new MyFFT(FFT_SIZE);

    // Read wave file
    cout<<"Read wave file"<<endl;
    WavFileIO* wav = new WavFileIO(input_file.c_str());

    // Framing
    cout<<"Framing"<<endl;
    vector<Frame*> recording;
    int num_frame = (wav->myDataSize/2-FRAME_LENGTH)/FRAME_SHIFT+1;
    for (int frame_idx=0; frame_idx<num_frame; ++frame_idx) {
        Frame *f = new Frame(FRAME_LENGTH);
        f->loadSample(wav->myData_short, frame_idx*FRAME_SHIFT);
        f->applyWindow(window);
        f->runFFT(fft);
        recording.push_back(f);
    }

    // Time stretching
    cout<<"Time stretching"<<endl;
    TimeStretcher *ts = new TimeStretcherPL(FFT_SIZE, FRAME_SHIFT);
    vector<Frame*> tmp_recording;
    ts->Stretch(ts_rate, recording, tmp_recording, true);

    // Pitch shifting
    cout<<"Pitch shifting"<<endl;
    PitchShifter *ps = new PitchShifter(FFT_SIZE, FRAME_SHIFT);
    vector<Frame*> synth_recording;
    ps->Shift(ps_rate, tmp_recording, synth_recording, true);

    // Synthesis
    cout<<"Synthesis"<<endl;

    vector<float> square_window(FRAME_LENGTH, 1.0); // for the denominator in synthesis
    window->applyWindow(&square_window[0]);
    window->applyWindow(&square_window[0]);

    int synth_size = wav->myDataSize/2*ts_rate;
    vector<short> synth_signal(synth_size, 0);
    vector<float> synth_normalize_coeff(synth_size, 0.0);
    for (int frame_idx=0; frame_idx<synth_recording.size(); ++frame_idx) {
        synth_recording[frame_idx]->runIFFT(fft);
        synth_recording[frame_idx]->applyWindow(window);
        float *frame = synth_recording[frame_idx]->getFrame();
        for (int sample_idx=0; sample_idx<FRAME_LENGTH && frame_idx*FRAME_SHIFT+sample_idx<synth_size; ++sample_idx) {
            synth_signal[frame_idx*FRAME_SHIFT+sample_idx]+=frame[sample_idx];
            synth_normalize_coeff[frame_idx*FRAME_SHIFT+sample_idx]+=square_window[sample_idx];
        }
    }
    for (int i=0; i<synth_size; ++i)
        synth_signal[i]/=synth_normalize_coeff[i];

    wav->setPath(output_file.c_str());
    wav->myDataSize = synth_size*2;
    wav->myData_short = &synth_signal[0];
    wav->save();

    cout<<"end"<<endl;
    return 0;
}
