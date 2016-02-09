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
            if (frame) delete frame;
            if (spectrum) delete spectrum;
        }

        void loadSample(short* samples, int begin);
        void applyWindow(Window* window) {window->applyWindow(frame);}
        void runFFT(MyFFT* fft);
        void runIFFT(MyFFT* fft);

        vector<float> getMagnitude() {
            vector<float> ans;
            for (int i=0; i<length/2+1; ++i)
                ans.push_back(spectrum[i].getMagnitude());
            return ans;
        }

        vector<float> getPhase() {
            vector<float> ans;
            for (int i=0; i<length/2+1; ++i)
                ans.push_back(spectrum[i].getPhase());
            return ans;
        }

    private:
        int length;
        float* frame;
        Complex* spectrum;
};

#endif
