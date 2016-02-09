#include "time_stretcher.h"
#include "overlap_adder.h"
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
    int FRAME_LENGTH = 4096;
    int FRAME_SHIFT = 1024;
    int FFT_SIZE = FRAME_LENGTH;

    Window *window = new HammingWindow(FRAME_LENGTH);
    MyFFT *fft = new MyFFT(FFT_SIZE);

    // Read wave file
    cout<<"Read wave file"<<endl;
    string test_file = "test/HungarianDanceNo5.wav";
    WavFileIO* wav = new WavFileIO(test_file.c_str());

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
    TimeStretcher ts(FFT_SIZE);
    vector<Frame*> new_recording;
    ts.Stretch(2, recording, new_recording, true);

    // Synthesis
    cout<<"Synthesis"<<endl;
    for (int frame_idx=0; frame_idx<new_recording.size(); ++frame_idx) {
        new_recording[frame_idx]->runIFFT(fft);
    }

    cout<<"end"<<endl;
    return 0;
}
