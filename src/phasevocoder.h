#ifndef PHASEVOCODER_H
#define PHASEVOCODER_H

#include "time_stretcher.h"
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

        void ReadWave(string input_file);
        void Analysis();
        void TimeStretching(float ts_rate);
        void PitchShifting(float ps_rate);
        void Synthesis();
        void WriteWave(string output_file);

    private:
        int analysis_frame_shift;
        int synthesis_frame_shift;
        int FFT_SIZE;
        bool phase_lock;
        Window *window;
        MyFFT *fft;
        WavFileIO *input_wav, *output_wav;
        int sampling_rate;
        double ts_rate;
        double ps_rate;
        vector<Frame*> spectrogram;
        short* synth_signal;
        int synth_size;
};

#endif
