#include "frame.h"
        
void Frame::RunFFT()
{
    float *real;
    float *imag;
    fft.fft_float(length, false, frame, NULL, real, imag);
    delete spectrum;
    spectrum = new Complex[length];
    for (int i=0; i<length; ++i) {
        spectrum[i].real = real[i];
        spectrum[i].imag = imag[i];
    }
}

void Frame::RunIFFT()
{
    float *real = new float[length];
    float *imag = new float[length];
    float *imag_out = new float[length];
    delete frame;
    frame = new float[length];
    for (int i=0; i<length; ++i) {
        real[i] = spectrum[i].real;
        imag[i] = spectrum[i].imag;
    }
    fft.fft_float(length, true, real, imag, frame, imag_out);
    delete real;
    delete imag;
    delete imag_out;
}
