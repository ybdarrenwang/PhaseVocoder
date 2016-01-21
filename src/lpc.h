#ifndef __LPC_H__
#define __LPC_H__

#include "global.h"

class LPC{
public:
    float ABS2(float a,float b) {return sqrt(a*a+b*b);}
    int CalPath(int frame_no, int frame_shift, float *Pitch);
    void vowel_markPath(int *vowel_mark, int n, int frame_shift, float *Pitch);
    void SpecEnv(float *coeff, float *SpecEnv_mag, float *SpecEnv_pha);
	float From_Data(float *data, float *lpc, int n, int m);
	void Predict(float *coeff, float *prime, int m, float *data, long n);
};

#endif
