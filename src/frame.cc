#include "frame.h"

using namespace std;
        
void Frame::loadSample(short* samples, int begin) {
    for (int i=0; i<length; ++i)
        frame[i] = (float)((int)samples[begin+i]);
}

void Frame::runFFT(MyFFT* fft) {
    float *real = new float[length];
    float *imag = new float[length];
    fft->fft_float(false, frame, NULL, real, imag);
    for (int i=0; i<length; ++i) {
        spectrum[i].real = real[i];
        spectrum[i].imag = imag[i];
    }
    delete real;
    delete imag;
}

void Frame::runIFFT(MyFFT* fft) {
    float *real = new float[length];
    float *imag = new float[length];
    float *imag_out = new float[length];
    for (int i=0; i<length; ++i) {
        real[i] = spectrum[i].real;
        imag[i] = spectrum[i].imag;
    }
    fft->fft_float(true, real, imag, frame, imag_out);
    delete real;
    delete imag;
    delete imag_out;
}

vector<float> Frame::getMagnitude() {
    vector<float> ans;
    for (int i=0; i<length/2+1; ++i)
        ans.push_back(spectrum[i].getMagnitude());
    return ans;
}

vector<float> Frame::getPhase() {
    vector<float> ans;
    for (int i=0; i<length/2+1; ++i)
        ans.push_back(spectrum[i].getPhase());
    return ans;
}
