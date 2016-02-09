#ifndef __COMPLEX_H__
#define __COMPLEX_H__

#include <math.h>

class Complex
{
    public:
        Complex() : real(0.0), imag(0.0) {}
        Complex(float r, float i) : real(r), imag(i) {}
        float getEnergy(){return real*real+imag*imag;}
        float getMagnitude(){return sqrt(getEnergy());}
        float getPhase(){return atan2(imag, real);}
        float real;
        float imag;
};

#endif
