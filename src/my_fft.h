#ifndef __MYFFT_H__
#define __MYFFT_H__

#include "global.h"

class MyFFT  
{
    public:
        MyFFT(){}
        virtual ~MyFFT(){}
        void fft_float(unsigned  NumSamples, bool InverseTransform, float *RealIn, float *ImagIn, float *RealOut, float *ImagOut);
        unsigned NumberOfBitsNeeded ( unsigned PowerOfTwo );
        unsigned ReverseBits (unsigned index, unsigned NumBits);
};

#endif
