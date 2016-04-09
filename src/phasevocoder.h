#ifndef PHASEVOCODER_H
#define PHASEVOCODER_H

#include "time_stretcher.h"
#include "time_stretcher_fd.h"
#include "time_stretcher_pl.h"
#include "time_stretcher_fd_pl.h"
#include "pitch_shifter.h"
#include "complex.h"
#include "window.h"
#include "hamming_window.h"
#include "frame.h"
#include "wav_io.h"

using namespace std;

class PhaseVocoder {
    public:
        PhaseVocoder(int frame_length, int frame_shift, bool phase_lock, bool fd_time_stretch, bool fd_pitch_shift, double _ts_factor, double _ps_factor);
        virtual ~PhaseVocoder();

        virtual void ReadWave(string input_file);
        virtual void Analysis();
        virtual void TimeStretching();
        virtual void PitchShifting();
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
        double ts_factor;
        double ps_factor;
        vector<Frame*> spectrogram;
        short* synth_signal;
        unsigned int synth_size;
};

#endif
