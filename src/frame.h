#ifndef __FRAME_H__
#define __FRAME_H__

#include "my_fft.h"
#include "complex.h"
#include "window.h"

class Frame
{
    public:
        Frame(int l) : length(l) {
            frame = new float[length];
            spectrum = new Complex[length];
        }
        virtual ~Frame() {
            delete frame;
            delete spectrum;
        }
        void loadSample(short* samples, int begin);
        void applyWindow() {window->applyWindow(frame);}
        void runFFT();
        void runIFFT();

    private:
        int length;
        float* frame;
        Complex* spectrum;
        MyFFT fft;
        static Window *window;
};

#endif
