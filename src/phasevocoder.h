#ifndef PHASEVOCODER_H
#define PHASEVOCODER_H

#include "time_stretcher.h"
#include "time_stretcher_fd.h"
#include "time_stretcher_pl.h"
#include "pitch_shifter.h"
#include "complex.h"
#include "window.h"
#include "hamming_window.h"
#include "frame.h"
#include "wav_io.h"

using namespace std;

class PhaseVocoder {
    public:
        PhaseVocoder(int frame_length, int frame_shift, bool pl);
        virtual ~PhaseVocoder();

        virtual void ReadWave(string input_file);
        virtual void Analysis();
        virtual void TimeStretching(float ts_rate);
        virtual void PitchShifting(float ps_rate);
        virtual void Synthesis();
        virtual void WriteWave(string output_file);

    protected:
        TimeStretcher *ts;
        PitchShifter *ps;
        int analysis_frame_shift;
        int synthesis_frame_shift;
        unsigned int FFT_SIZE;
        bool phase_lock;
        Window *window;
        MyFFT *fft;
        WavFileIO *input_wav, *output_wav;
        int sampling_rate;
        double ts_rate;
        double ps_rate;
        vector<Frame*> spectrogram;
        short* synth_signal;
        unsigned int synth_size;
};

#endif
