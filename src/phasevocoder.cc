#include "phasevocoder.h"

using namespace std;

PhaseVocoder::PhaseVocoder(int frame_length, int frame_shift, bool phase_lock, bool fd_time_stretch, bool fd_pitch_shift, double _ts_factor, double _ps_factor)
{
    #ifdef DEBUG
    cout<<"Initialize parameters and functions"<<endl;
    #endif
    analysis_frame_shift = synthesis_frame_shift = frame_shift;
    FFT_SIZE = frame_length;
    window = new HammingWindow(FFT_SIZE);
    fft = new MyFFT(FFT_SIZE);
    ts_factor = _ts_factor;
    ps_factor = _ps_factor;
    if (phase_lock)
    {
        if (fd_time_stretch)
            ts = new TimeStretcherFDPL(FFT_SIZE, analysis_frame_shift);
        else
            ts = new TimeStretcherPL(FFT_SIZE, analysis_frame_shift);
    }
    else
    {
        if (fd_time_stretch)
            ts = new TimeStretcherFD(FFT_SIZE, analysis_frame_shift);
        else
            ts = new TimeStretcher(FFT_SIZE, analysis_frame_shift);
    }
    if (fd_pitch_shift)
        ps = new PitchShifter(FFT_SIZE, analysis_frame_shift);
    else
        ps = NULL;
}
    
PhaseVocoder::~PhaseVocoder()
{
    delete ts;
    if (ps) delete ps;
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
    synth_size = input_wav->myDataSize/2*ts_factor;
    #ifdef DEBUG
    cout<<input_wav->getSummary()<<endl;
    #endif
}

void PhaseVocoder::Analysis()
{
    vector<short> data;
    if (ps==NULL && ps_factor!=1) {
        resample(input_wav->myData_short, input_wav->myDataSize/2, ps_factor, data); // if ps_factor==1, this is equal to copying data
        ts_factor*=ps_factor;
        ps_factor = 1;
    }
    else
        resample(input_wav->myData_short, input_wav->myDataSize/2, 1, data);

    int num_frame = (data.size()-FFT_SIZE)/analysis_frame_shift+1;
    for (int frame_idx=0; frame_idx<num_frame; ++frame_idx) {
        Frame *f = new Frame(FFT_SIZE);
        f->loadSample(&(data[0]), frame_idx*analysis_frame_shift);
        f->applyWindow(window);
        f->runFFT(fft);
        spectrogram.push_back(f);
    }
}

void PhaseVocoder::TimeStretching()
{
    if (ts!=NULL && ts_factor!=1) {
        #ifdef DEBUG
        cout<<"Time stretching"<<endl;
        #endif
        vector<Frame*> new_spectrogram;
        ts->Stretch(ts_factor, spectrogram, new_spectrogram, synthesis_frame_shift, true);
        for (vector<Frame*>::iterator f=spectrogram.begin(); f!=spectrogram.end(); ++f) delete (*f);
        spectrogram = new_spectrogram;
    }
}

void PhaseVocoder::PitchShifting()
{
    if (ps!=NULL && ps_factor!=1) {
        #ifdef DEBUG
        cout<<"Pitch shifting"<<endl;
        #endif
        vector<Frame*> new_spectrogram;
        ps->Shift(ps_factor, spectrogram, new_spectrogram, true);
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
