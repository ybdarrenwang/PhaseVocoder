#ifndef __LPC_H__
#define __LPC_H__

#include "global.h"

class LPC{
public:
    double ABS2(double a,double b) {return sqrt(a*a+b*b);}
    int CalPath(int frame_no, int frame_shift, double *Pitch);
    void vowel_markPath(int *vowel_mark, int n, int frame_shift, double *Pitch);
    void SpecEnv(double *coeff, double *SpecEnv_mag, double *SpecEnv_pha);
	double From_Data(double *data, double *lpc, int n, int m);
	void Predict(double *coeff, double *prime, int m, double *data, long n);
};

#endif
