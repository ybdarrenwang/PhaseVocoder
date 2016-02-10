#ifndef __OVERLAPADDER_H__
#define __OVERLAPADDER_H__

#include "vocoder_functions.h"

class OverlapAdder : public VocoderFunctions
{
    public:
        int PSOLA(float **frame, float *output, float *SentenceCoeff, float *windows,
                int sample_index, float *Pitch, int WordStart, int WordMid, int WordEnd,
                int frame_len, int frame_shift, int bias);
        int SOLA(float **frame, float *output, float *SentenceCoeff, float *windows,
                int sample_index, int WordStart, int WordMid, int WordEnd, 
                int frame_len, int frame_shift, int bias);
};

#endif
