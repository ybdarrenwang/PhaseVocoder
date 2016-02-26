#ifndef __COMPLEX_H__
#define __COMPLEX_H__

#include <math.h>

class Complex
{
    public:
        Complex() : real(0.0), imag(0.0) {}
        Complex(double r, double i) : real(r), imag(i) {}
        double getEnergy(){return real*real+imag*imag;}
        double getMagnitude(){return sqrt(getEnergy());}
        double getPhase(){return atan2(imag, real);}
        Complex getConjugate(){return Complex(real, -1*imag);}
        double real;
        double imag;
};

#endif
