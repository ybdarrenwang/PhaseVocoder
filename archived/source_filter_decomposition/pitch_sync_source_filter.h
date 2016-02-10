#ifndef __PITCHSYNCSOURCEFILTER_H__
#define __PITCHSYNCSOURCEFILTER_H__

#include "SourceFilter.h"

class PitchSyncSourceFilter : public SourceFilter
{
    public:
        virtual void Decomposition(float *data, float *excitation,
                int sample_no, int word_no, int* WordStart, int* WordMid,
                int* WordEnd, int frame_shift, int frame_len, float** Pitch);
        virtual void Recombination(short* data, float* excitation, int sample_no, int word_no, 
                        float* WordDurMul, int* WordStart, int* WordMid, int* WordEnd, 
                        int frame_shift, int frame_len, int* syn_bias);
};

#endif
