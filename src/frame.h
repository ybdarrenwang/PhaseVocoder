#ifndef __FRAME_H__
#define __FRAME_H__

#include "my_fft.h"
#include "complex.h"

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

        void RunFFT();
        void RunIFFT();

    private:
        int length;
        float* frame;
        Complex* spectrum;
        MyFFT fft;
};

#endif
