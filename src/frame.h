#ifndef __FRAME_H__
#define __FRAME_H__

#include "my_fft.h"
#include "complex.h"
#include "window.h"

using namespace std;

class Frame
{
    public:
        Frame(int l) : length(l) {
            frame = new float[length];
            spectrum = new Complex[length];
        }

        virtual ~Frame() {
            cout<<"delete frame"<<endl;
            if (frame) delete frame;
            if (spectrum) delete spectrum;
        }

        void loadSample(short* samples, int begin);
        void applyWindow(Window* window) {window->applyWindow(frame);}
        void runFFT(MyFFT* fft);
        void runIFFT(MyFFT* fft);

        vector<float> getMagnitude();
        vector<float> getPhase();
        Complex getSpectrum(int freq) {return spectrum[freq];}
        void setSpectrum(Complex* _spectrum) {
            for (int freq=0; freq<length/2+1; ++freq)
                spectrum[freq] = _spectrum[freq];
            for (int freq=length/2+1; freq<length; ++freq)
                spectrum[freq] = spectrum[length-freq].getConjugate();
        }
        void setSpectrum(int freq, Complex spec) {spectrum[freq] = spec;}

        float* getFrame() {return frame;}

    private:
        int length;
        float* frame;
        Complex* spectrum;
};

#endif
