#include "time_stretcher.h"
#include "pitch_shifter.h"
#include "overlap_adder.h"
#include "complex.h"
#include "window.h"
#include "frame.h"
#include "wav_io.h"
#include <iostream>

using namespace std;

int main(int argc, char **argv)
{
    string test_file = "test/HungarianDanceNo5.wav";

    WavFileIO* wav = new WavFileIO(test_file.c_str());

    // Framing
    vector<Frame*> recording;
    int num_frame = (wav->myDataSize/2-FRAME_LENGTH)/FRAME_SHIFT+1;
    for (int frame_idx=0; frame_idx<num_frame; ++frame_idx) {
        Frame* f = new Frame(FRAME_LENGTH);
        f->loadSample(wav->myData_short, frame_idx*FRAME_SHIFT);
        f->applyWindow();
        f->runFFT();
        recording.push_back(f);
    }

    cout<<"end"<<endl;
    return 0;
}
