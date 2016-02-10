#ifndef	MATRIX_H
#define MATRIX_H

#include "mymath.h"

using namespace std;

class Matrix
{
	public:
		Matrix();
		Matrix(int, int);
		virtual ~Matrix();

		void Init(int, int);
		void Init(vector< vector<double> >&);

		double GetEntry(int, int);
		void SetEntry(int, int, double);
		void GetRow(int, vector<double>&);
		void SetRow(int, vector<double>&);
		void GetColumn(int, vector<double>&);
		void SetColumn(int, vector<double>&);
		int GetDim(int);

		void RowAddition(int, int, double);
		void RowMul(int, double);
		void RowSwitch(int, int);

		bool IsIdentity();

	private:
		double** entry;
		int dim[2];
		static const double precision = 1e-6; // precision interval for checking entries
};

Matrix Mat_Mul(Matrix, Matrix);
void Matrix_GaussJordan(Matrix&, Matrix&);

#endif
