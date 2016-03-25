#include "phasevocoder.h"

using namespace std;

PhaseVocoder::PhaseVocoder(int frame_length, int frame_shift, bool pl)
{
    #ifdef DEBUG
    cout<<"Initialize parameters and functions"<<endl;
    #endif
    analysis_frame_shift = synthesis_frame_shift = frame_shift;
    FFT_SIZE = frame_length;
    phase_lock = pl;
    window = new HammingWindow(FFT_SIZE);
    fft = new MyFFT(FFT_SIZE);
    ts_rate = 1.0;
    ps_rate = 1.0;
    if (phase_lock)
        ts = new TimeStretcherPL(FFT_SIZE, analysis_frame_shift);
    else
        ts = new TimeStretcherFD(FFT_SIZE, analysis_frame_shift);
    ps = new PitchShifter(FFT_SIZE, analysis_frame_shift);
}
    
PhaseVocoder::~PhaseVocoder()
{
    delete ts;
    delete ps;
    delete window;
    delete fft;
    if (input_wav) delete input_wav;
    if (output_wav) delete output_wav;
    for (vector<Frame*>::iterator f=spectrogram.begin(); f!=spectrogram.end(); ++f) delete (*f);
}

void PhaseVocoder::ReadWave(string input_file)
{
    #ifdef DEBUG
    cout<<"Read "<<input_file<<endl;
    #endif
    input_wav = new WavFileIO(input_file);
    sampling_rate = input_wav->mySampleRate;
    #ifdef DEBUG
    cout<<input_wav->getSummary()<<endl;
    #endif
}

void PhaseVocoder::Analysis()
{
    int num_frame = (input_wav->myDataSize/2-FFT_SIZE)/analysis_frame_shift+1;
    for (int frame_idx=0; frame_idx<num_frame; ++frame_idx) {
        Frame *f = new Frame(FFT_SIZE);
        f->loadSample(input_wav->myData_short, frame_idx*analysis_frame_shift);
        f->applyWindow(window);
        f->runFFT(fft);
        spectrogram.push_back(f);
    }
}

void PhaseVocoder::TimeStretching(float rate)
{
    ts_rate = rate;
    if (ts_rate!=1) {
        #ifdef DEBUG
        cout<<"Time stretching"<<endl;
        #endif
        vector<Frame*> new_spectrogram;
        ts->Stretch(ts_rate, spectrogram, new_spectrogram, true);
        for (vector<Frame*>::iterator f=spectrogram.begin(); f!=spectrogram.end(); ++f) delete (*f);
        spectrogram = new_spectrogram;
    }
}

void PhaseVocoder::PitchShifting(float rate)
{
    ps_rate = rate;
    if (ps_rate!=1) {
        #ifdef DEBUG
        cout<<"Pitch shifting"<<endl;
        #endif
        PitchShifter *ps = new PitchShifter(FFT_SIZE, analysis_frame_shift);
        vector<Frame*> new_spectrogram;
        ps->Shift(ps_rate, spectrogram, new_spectrogram, true);
        for (vector<Frame*>::iterator f=spectrogram.begin(); f!=spectrogram.end(); ++f) delete (*f);
        spectrogram = new_spectrogram;
    }
}

void PhaseVocoder::Synthesis()
{
    #ifdef DEBUG
    cout<<"Synthesis"<<endl;
    #endif
    vector<double> square_window(FFT_SIZE, 1.0); // for the denominator in synthesis
    window->applyWindow(&square_window[0]);
    window->applyWindow(&square_window[0]);

    synth_size = input_wav->myDataSize/2*ts_rate;
    synth_signal = new short[synth_size];
    for (unsigned int i=0; i<synth_size; ++i) synth_signal[i]=0;

    vector<double> synth_normalize_coeff(synth_size, 0.0);
    for (unsigned int frame_idx=0; frame_idx<spectrogram.size(); ++frame_idx) {
        spectrogram[frame_idx]->runIFFT(fft);
        spectrogram[frame_idx]->applyWindow(window);
        double *frame = spectrogram[frame_idx]->getFrame();
        for (unsigned int sample_idx=0; sample_idx<FFT_SIZE && frame_idx*synthesis_frame_shift+sample_idx<synth_size; ++sample_idx) {
            synth_signal[frame_idx*synthesis_frame_shift+sample_idx]+=frame[sample_idx];
            synth_normalize_coeff[frame_idx*synthesis_frame_shift+sample_idx]+=square_window[sample_idx];
        }
    }
    for (unsigned int i=0; i<synth_size; ++i)
        synth_signal[i]/=synth_normalize_coeff[i];
}

void PhaseVocoder::WriteWave(string output_file)
{
    #ifdef DEBUG
    cout<<"Write "<<output_file<<endl;
    #endif
    output_wav = new WavFileIO(*input_wav);
    output_wav->setPath(output_file);
    output_wav->myDataSize = synth_size*2;
    output_wav->myData_short = synth_signal;
    output_wav->save();
}
