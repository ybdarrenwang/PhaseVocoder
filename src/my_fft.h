#ifndef __MYFFT_H__
#define __MYFFT_H__

#include "global.h"

class MyFFT  
{
    public:
        MyFFT(unsigned n) : NumSamples(n) {}
        virtual ~MyFFT(){}
        void fft_float(bool InverseTransform, float *RealIn, float *ImagIn, float *RealOut, float *ImagOut);

    private:
        unsigned NumSamples;
        unsigned NumberOfBitsNeeded ( unsigned PowerOfTwo );
        unsigned ReverseBits (unsigned index, unsigned NumBits);
};

#endif
