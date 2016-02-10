#include "spline.h"

using namespace std;

Spline::Spline()
{
}

Spline::~Spline()
{
}

// Performe cubic spline interpolation (natural spline)
// Ex. if pitch contour stylization is performed......
//       x: time index of each available pitch value
//       y: corresponding pitch values
//       return: interpolated pitch contour
void Spline::CubicSpline(vector<double>& x, vector<double>& y)
{
//	cout<<"Performing cubic spline interpolation"<<endl;

	int i;

	// construct the matrix equation for the second deritives M:
	// coeff * M = sol -> M = coeff^-1 * sol
	vector<double> dx, dy;
	for (i=0; i<(int)x.size()-1; i++) dx.push_back(x[i+1]-x[i]);
	for (i=0; i<(int)y.size()-1; i++) dy.push_back(y[i+1]-y[i]);

	Matrix coeff(x.size()-2, x.size()-2);
	for (i=0; i<(int)x.size()-2; i++)	
	{
		if (i>0) coeff.SetEntry(i, i-1, dx[i]);
		coeff.SetEntry(i, i, 2*(dx[i]+dx[i+1]));
		if (i<(int)x.size()-3) coeff.SetEntry(i, i+1, dx[i+1]);
	}

	Matrix sol(x.size()-2, 1);
	for (i=0; i<(int)x.size()-2; i++)
		sol.SetEntry(i, 0, 6.0*(dy[i+1]/dx[i+1] - dy[i]/dx[i]));

	// calculate second deritives
	vector<double> M;

	Matrix_GaussJordan(coeff, sol);
	M.push_back(0);
	for (i=0; i<sol.GetDim(0); i++)
		M.push_back(sol.GetEntry(i, 0));
	M.push_back(0);

	// calculate piecewise polynomial coefficients
	for (i=0; i<(int)M.size()-1; i++)
	{
		a.push_back((M[i+1]-M[i])/6/dx[i]);
		b.push_back(M[i]/2);
		c.push_back(dy[i]/dx[i] - dx[i]*(M[i+1]+2*M[i])/6);
		d.push_back(y[i]);
	}

	knots = x;
}

void Spline::GetPPValue(int start, int end, vector<double>& y)
{
//	cout<<"Get the interpolated value of interval ["<<start<<", "<<end<<")"<<endl;

	int x, i = 0;
	y.clear();

	for (x=start; x<end; x++)
	{
		while (i < (int)knots.back() && x >= knots[i])
			i++;
		y.push_back(a[i-1]*pow(x-knots[i-1],3) + b[i-1]*pow(x-knots[i-1],2) + c[i-1]*(x-knots[i-1]) + d[i-1]);
	}
	y.push_back(a[i-1]*pow(end-knots[i-1],3) + b[i-1]*pow(end-knots[i-1],2) + c[i-1]*(end-knots[i-1]) + d[i-1]);
}

void Spline::Reset()
{
	a.clear();
	b.clear();
	c.clear();
	d.clear();
	knots.clear();
}
