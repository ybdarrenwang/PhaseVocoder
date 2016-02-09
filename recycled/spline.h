#ifndef SPLINE_H
#define SPLINE_H

#include "matrix.h"

using namespace std;

// Cubic spline piecewise polynomial
// pp(x) = a_i*x^3 + b_i*x^2 + c_i*x + d_i, for knots[i] < x < knots[i+1]
class Spline
{
	public:
		Spline();
		virtual ~Spline();
		void CubicSpline(vector<double>&, vector<double>&);
		void GetPPValue(int, int, vector<double>&);
		void Reset();

	private:
		vector<double> a;
		vector<double> b;
		vector<double> c;
		vector<double> d;
		vector<double> knots;
};

#endif
