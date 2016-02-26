#include "frame.h"

using namespace std;
        
void Frame::loadSample(short* samples, int begin) {
    for (int i=0; i<length; ++i)
        frame[i] = (double)((int)samples[begin+i]);
}

void Frame::runFFT(MyFFT* fft) {
    double *real = new double[length];
    double *imag = new double[length];
    fft->fft_double(false, frame, NULL, real, imag);
    for (int i=0; i<length; ++i) {
        spectrum[i].real = real[i];
        spectrum[i].imag = imag[i];
    }
    delete real;
    delete imag;
}

void Frame::runIFFT(MyFFT* fft) {
    double *real = new double[length];
    double *imag = new double[length];
    double *imag_out = new double[length];
    for (int i=0; i<length; ++i) {
        real[i] = spectrum[i].real;
        imag[i] = spectrum[i].imag;
    }
    fft->fft_double(true, real, imag, frame, imag_out);
    delete real;
    delete imag;
    delete imag_out;
}

vector<double> Frame::getMagnitude() {
    vector<double> ans;
    for (int i=0; i<length/2+1; ++i)
        ans.push_back(spectrum[i].getMagnitude());
    return ans;
}

vector<double> Frame::getPhase() {
    vector<double> ans;
    for (int i=0; i<length/2+1; ++i)
        ans.push_back(spectrum[i].getPhase());
    return ans;
}
