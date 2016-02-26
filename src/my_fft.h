#ifndef __MYFFT_H__
#define __MYFFT_H__

#include "global.h"

class MyFFT  
{
    public:
        MyFFT(unsigned n) : NumSamples(n) {}
        virtual ~MyFFT(){}
        void fft_double(bool InverseTransform, double *RealIn, double *ImagIn, double *RealOut, double *ImagOut);

    private:
        unsigned NumSamples;
        unsigned NumberOfBitsNeeded ( unsigned PowerOfTwo );
        unsigned ReverseBits (unsigned index, unsigned NumBits);
};

#endif
