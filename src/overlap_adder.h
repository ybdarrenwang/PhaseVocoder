#ifndef __OVERLAPADDER_H__
#define __OVERLAPADDER_H__

#include "vocoder_functions.h"

class OverlapAdder : public VocoderFunctions
{
    public:
        int PSOLA(double **frame, double *output, double *SentenceCoeff, double *windows,
                int sample_index, double *Pitch, int WordStart, int WordMid, int WordEnd,
                int frame_len, int frame_shift, int bias);
        int SOLA(double **frame, double *output, double *SentenceCoeff, double *windows,
                int sample_index, int WordStart, int WordMid, int WordEnd, 
                int frame_len, int frame_shift, int bias);
};

#endif
