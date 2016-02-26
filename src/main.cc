#include "time_stretcher.h"
#include "time_stretcher_pl.h"
#include "pitch_shifter.h"
#include "complex.h"
#include "window.h"
#include "hamming_window.h"
#include "frame.h"
#include "wav_io.h"
#include <iostream>

using namespace std;

void usage( const char *prog ) {
    cout << endl
        << "Usage: "<<prog<<" [options] ..."<<endl<<endl
        << "options:"<<endl
        << " -t TIME_SCALE_FACTOR  (>1 for slowdown, <1 for speedup) [default:1]"<<endl
        << " -p PITCH_SCALE_FACTOR (>1 for higher, <1 for lower)     [default:1]"<<endl
        << " -i INPUT_WAVE_FILE_PATH"<<endl
        << " -o SYNTHESIZED_WAVE_FILE_PATH"<<endl<<endl
        << " --phaselock           (enable phase-locking)            [default:off]"<<endl<<endl;
    exit(1);
}

void readConfig(vector<string> &Args, double &ts_rate, double &ps_rate, string &input_file, string &output_file, bool &phase_lock) {
    for(int i=0; i<Args.size(); ++i) {
        if( Args[i] == "-t" && Args.size() > i ) {
            ++i;
            stringstream ss(Args[i]);
            ss >> ts_rate;
        }
        else if( Args[i] == "-p" && Args.size() > i )
        {
            ++i;
            stringstream ss(Args[i]);
            ss >> ps_rate;
        }
        else if( Args[i] == "-i" && Args.size() > i ) {
            ++i;
            input_file = Args[i];
        }
        else if( Args[i] == "-o" && Args.size() > i ) {
            ++i;
            output_file = Args[i];
        }
        else if( Args[i] == "--phaselock" && Args.size() > i ) {
            phase_lock = true;
        }
    }
}

int main(int argc, char **argv) {
    if( argc < 2 ){
        usage( argv[0] );
        exit(0);
    }

    // Read arguments
    string input_file = "";
    string output_file = "";
    double ts_rate = 1;
    double ps_rate = 1;
    bool phase_lock = false;

    vector<string> Args;
    for(int i=1; i<argc; ++i)
    {
        string tmpstr(argv[i]);
        Args.push_back( tmpstr );
    }
    readConfig(Args, ts_rate, ps_rate, input_file, output_file, phase_lock);
    if (input_file=="" || output_file=="") {
        cerr<<"ERROR: input or output file not specified"<<endl;
        exit(1);
    }
    cout<<endl<<"[ "<<input_file<<" -> "<<output_file<<" ]"<<endl;

    // Initialize parameters and functions
    cout<<"Initialize parameters and functions"<<endl;
    int FRAME_LENGTH = 4096;
    int FRAME_SHIFT = 1024;
    int FFT_SIZE = FRAME_LENGTH;

    Window *window = new HammingWindow(FRAME_LENGTH);
    MyFFT *fft = new MyFFT(FFT_SIZE);

    // Read wave file
    cout<<"Read wave file"<<endl;
    WavFileIO* input_wav = new WavFileIO(input_file);
    int sampling_rate = input_wav->mySampleRate;

    // Framing
    cout<<"Framing"<<endl;
    vector<Frame*> recording;
    int num_frame = (input_wav->myDataSize/2-FRAME_LENGTH)/FRAME_SHIFT+1;
    for (int frame_idx=0; frame_idx<num_frame; ++frame_idx) {
        Frame *f = new Frame(FRAME_LENGTH);
        f->loadSample(input_wav->myData_short, frame_idx*FRAME_SHIFT);
        f->applyWindow(window);
        f->runFFT(fft);
        recording.push_back(f);
    }

    // Time stretching
    cout<<"Time stretching"<<endl;
    TimeStretcher *ts;
    if (phase_lock)
        ts = new TimeStretcherPL(FFT_SIZE, FRAME_SHIFT);
    else
        ts = new TimeStretcher(FFT_SIZE, FRAME_SHIFT);
    vector<Frame*> tmp_recording;
    ts->Stretch(ts_rate, recording, tmp_recording, true);

    // Pitch shifting
    cout<<"Pitch shifting"<<endl;
    PitchShifter *ps = new PitchShifter(FFT_SIZE, FRAME_SHIFT);
    vector<Frame*> synth_recording;
    ps->Shift(ps_rate, tmp_recording, synth_recording, true);

    // Synthesis
    cout<<"Synthesis"<<endl;

    vector<double> square_window(FRAME_LENGTH, 1.0); // for the denominator in synthesis
    window->applyWindow(&square_window[0]);
    window->applyWindow(&square_window[0]);

    int synth_size = input_wav->myDataSize/2*ts_rate;
    vector<short> synth_signal(synth_size, 0);
    vector<double> synth_normalize_coeff(synth_size, 0.0);
    for (int frame_idx=0; frame_idx<synth_recording.size(); ++frame_idx) {
        synth_recording[frame_idx]->runIFFT(fft);
        synth_recording[frame_idx]->applyWindow(window);
        double *frame = synth_recording[frame_idx]->getFrame();
        for (int sample_idx=0; sample_idx<FRAME_LENGTH && frame_idx*FRAME_SHIFT+sample_idx<synth_size; ++sample_idx) {
            synth_signal[frame_idx*FRAME_SHIFT+sample_idx]+=frame[sample_idx];
            synth_normalize_coeff[frame_idx*FRAME_SHIFT+sample_idx]+=square_window[sample_idx];
        }
    }
    for (int i=0; i<synth_size; ++i)
        synth_signal[i]/=synth_normalize_coeff[i];

    // Write wav file
    WavFileIO* output_wav = new WavFileIO(*input_wav);
    output_wav->setPath(output_file);
    output_wav->myDataSize = synth_size*2;
    output_wav->myData_short = &synth_signal[0]; // note: delete output_wav cause double deletion!!!
    output_wav->save();

    // Release memory
    delete window;
    delete fft;
    delete input_wav;
    delete ts;
    delete ps;
    for (vector<Frame*>::iterator f=recording.begin(); f!=recording.end(); ++f) delete (*f);
    for (vector<Frame*>::iterator f=tmp_recording.begin(); f!=tmp_recording.end(); ++f) delete (*f);
    for (vector<Frame*>::iterator f=synth_recording.begin(); f!=synth_recording.end(); ++f) delete (*f);
    cout<<"Complete!"<<endl;
    return 0;
}
