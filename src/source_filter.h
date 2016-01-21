#ifndef __SOURCEFILTER_H__
#define __SOURCEFILTER_H__

#include <math.h>
#include "vocoder_functions.h"
#include "lpc.h"

/**
* For derived classes "PitchSync" and "FreqDomain" in the future, Decomposition
* and recombination need additiona information as input; need to figure out how
* to implement.
*/
class SourceFilter
{
    public:
        SourceFilter(){}
        virtual ~SourceFilter(){}
        virtual void clear_mem(int word_no);
        virtual void Decomposition(float *data, float *excitation,
                int sample_no, int word_no, int* WordStart, int* WordMid,
                int* WordEnd, int frame_shift, int frame_len);
        /*virtual void Recombination(short* data, float* excitation, int sample_no, int word_no, 
                float* WordDurMul, int* WordStart, int* WordMid, int* WordEnd, 
                int frame_shift, int frame_len, int* syn_bias) = 0;*/
        virtual void Recombination(float** spec, int word_index, int frame_index, int WordMid);

    protected:
        int** vowel_mark; // vowel_mark[i][j]: the jth vowel f0 mark of ith word
        int** frame_no; // frame_no[i][j]: the number of frames of consonant(j=0) or vowel(j=1) of ith word
        float ****coeff; // coeff[i][j][k][m]: the mth LPC coefficient of kth frame in
        // consonant(j=0) or vowel(j=1) of ith word
        Window* window_function;
        LPC lpc;
};

#endif
