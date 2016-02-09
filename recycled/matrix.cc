#include "matrix.h"

using namespace std;

Matrix::Matrix()
{
	dim[0] = dim[1] = 0;
}

Matrix::Matrix(int dim1, int dim2)
{
	Init(dim1, dim2);
}

void Matrix::Init(int dim1, int dim2)
{
	entry = new double*[dim1];
	for (int i=0; i<dim1; i++)
	{
		entry[i] = new double[dim2];
		for (int j=0; j<dim2; j++)
		{
			entry[i][j] = 0.0;
		}
	}
	dim[0] = dim1;
	dim[1] = dim2;
}

void Matrix::Init(vector< vector<double> >& data)
{
	dim[0] = (int)data.size();
	dim[1] = (int)data[0].size();

	entry = new double*[dim[0]];
	for (int i=0; i<dim[0]; i++)
	{
		if ((int)data[i].size() != dim[1])
		{
			cerr<<"[Error] instance sizes for matrix initialization unequal"<<endl;
			exit(1);
		}

		entry[i] = new double[dim[1]];
		for (int j=0; j<dim[1]; j++)
		{
			entry[i][j] = data[i][j];
		}
	}
}

Matrix::~Matrix()
{
	if (dim[0] > 0)
	{
		for (int i=0; i<dim[0]; i++)
		{
			delete [] entry[i];
		}
		delete [] entry;
	}
}
		
double Matrix::GetEntry(int i, int j)
{
	if (i>=dim[0] || i<0 || j>=dim[1] || j<0)
	{
		cerr<<"[Error] request exceed matrix dimension"<<endl;
		exit(1);
	}
	
	return entry[i][j];
}

void Matrix::SetEntry(int i, int j, double value)
{
	if (i>=dim[0] || i<0 || j>=dim[1] || j<0)
	{
		cerr<<"[Error] request exceed matrix dimension"<<endl;
		exit(1);
	}
	
	entry[i][j] = value;
}

void Matrix::GetRow(int i, vector<double>& tmp_row)
{
	if (i>=dim[0] || i<0)
	{
		cerr<<"[Error] request exceed matrix dimension"<<endl;
		exit(1);
	}

	tmp_row.clear();
	for (int j=0; j<dim[1]; j++)
		tmp_row.push_back(entry[i][j]);
}

void Matrix::SetRow(int i, vector<double>& row)
{
	if (i>=dim[0] || i<0)
	{
		cerr<<"[Error] request exceed matrix dimension"<<endl;
		exit(1);
	}
	if ((int)row.size() != dim[1])
	{
		cerr<<"[Error] vector and matrix dimension does not match"<<endl;
		exit(1);
	}

	for (int j=0; j<dim[1]; j++)
		entry[i][j] = row[j];
}

void Matrix::GetColumn(int j, vector<double>& col)
{
	if (j>=dim[1] || j<0)
	{
		cerr<<"[Error] request exceed matrix dimension"<<endl;
		exit(1);
	}

	col.clear();
	for (int i=0; i<dim[0]; i++)
		col.push_back(entry[i][j]);
}

void Matrix::SetColumn(int j, vector<double>& col)
{
	if (j>=dim[1] || j<0)
	{
		cerr<<"[Error] request exceed matrix dimension"<<endl;
		exit(1);
	}
	if ((int)col.size() != dim[0])
	{
		cerr<<"[Error] vector and matrix dimension does not match"<<endl;
		exit(1);
	}

	for (int i=0; i<dim[0]; i++)
		entry[i][j] = col[i];
}

int Matrix::GetDim(int i)
{
	if (i != 0 && i != 1)
	{
		cerr<<"Wrong request for matrix dimension ("<<i<<")"<<endl;
		exit(1);
	}

	return dim[i];
}

// Perform matrix row addition
// Ri <- Ri + mul*Rj
void Matrix::RowAddition(int i, int j, double mul)
{
	if (mul != 0)
	{
		if (mul == 1)
			for (int p=0; p<dim[1]; p++)
				entry[i][p] += entry[j][p];
		else
		{
			for (int p=0; p<dim[1]; p++)
				entry[i][p] += entry[j][p]*mul;
		}
	}
}

// Perform matrix row multiplication
// Ri <- mul*Ri
void Matrix::RowMul(int i, double mul)
{
	if (mul != 1)
	{
		if (mul == 0)
		{
			for (int j=0; j<dim[1]; j++)
				entry[i][j] = 0;
		}
		else
		{
			for (int j=0; j<dim[1]; j++)
				entry[i][j]*=mul;
		}
	}
}

// Perform matrix row switching
void Matrix::RowSwitch(int i, int j)
{
	double tmp;
	for (int p=0; p<dim[1]; p++)
	{
		tmp = entry[i][p];
		entry[i][p] = entry[j][p];
		entry[j][p] = tmp;
	}
}

// Check whether this is an identity matrix
bool Matrix::IsIdentity()
{
	if (dim[0] != dim[1])
	{
		cerr<<"[Warning] non-square matrix asked to check identity"<<endl;
		return false;
	}

	for (int i=0; i<dim[0]; i++)
	{
		for (int j=0; j<dim[1]; j++)
		{
			if (i != j)
			{
				if (entry[i][j] < 0-precision || entry[i][j] > 0+precision) return false;
			}
			else
			{
				if (entry[i][j] < 1-precision || entry[i][j] > 1+precision) return false;
			}
		}
	}

	return true;
}

// ==================================================================
// Matrix manipulation
// ==================================================================
// Perform matrix multiplication
Matrix Matrix_Mul(Matrix A, Matrix B)
{
	Matrix answer(A.GetDim(0),B.GetDim(1));
	vector<double> tmp_row, tmp_col;

	if (A.GetDim(1) != B.GetDim(0))
	{
		cerr<<"[Error] matrix dimension does not match for multiplication"<<endl;
		exit(1);
	}

	for (int i=0; i<answer.GetDim(0); i++)
	{
		A.GetRow(i, tmp_row);
		for (int j=0; j<answer.GetDim(1); j++)
		{
			B.GetColumn(j, tmp_col);
			answer.SetEntry(i, j, Vector_InnerProd(tmp_row, tmp_col));
		}
	}

	return answer;
}

// Perform Gauss-Jordan elimination on matrix A,
// and do same operation on matrix B.
void Matrix_GaussJordan(Matrix& A, Matrix& B)
{
	if (A.GetDim(0) != B.GetDim(0))
	{
		cerr<<"[Error] matrix dimension does not match for Gauss-Jordan elimination"<<endl;
		exit(1);
	}

	int i,j;
	int row_num = A.GetDim(0);
	vector<double> tmp_row1, tmp_row2, tmp_row_sum, tmp_col;

	// horizontally concatenate matrix A and B
	Matrix tmp_mat(row_num, A.GetDim(1)+B.GetDim(1));
	for (i=0; i<A.GetDim(1); i++)
	{
		A.GetColumn(i, tmp_col);
		tmp_mat.SetColumn(i, tmp_col);
	}
	for (i=0; i<B.GetDim(1); i++)
	{
		B.GetColumn(i, tmp_col);
		tmp_mat.SetColumn(A.GetDim(1)+i, tmp_col);
	}

	// top-down eliminate the i-th variable
	for(i=0; i<row_num; i++)
	{
		// ensure the i-th element on i-th row is not 0
		j = i;
		if(tmp_mat.GetEntry(i, i) == 0)
		{
			j = i+1;
			while(j<row_num && tmp_mat.GetEntry(j, i) == 0) j++;
			if (j == row_num)
			{
				cerr<<"[Error] not full rank square matrix, cannot perform Gauss-Jordan elimination"<<endl;
				exit(1);
			}
			else
				tmp_mat.RowSwitch(i, j);
		}

		tmp_mat.RowMul(i, 1.0/tmp_mat.GetEntry(i,i));
		for (j=i+1; j<row_num; j++)
			tmp_mat.RowAddition(j, i, -1.0*tmp_mat.GetEntry(j,i));
	}

	// buttom-up eliminate the i-th variable
	for(i=row_num-1; i>=0; i--)
	{
		for (j=0; j<i; j++)
			tmp_mat.RowAddition(j, i, -1.0*tmp_mat.GetEntry(j,i));
	}

	// check if succeeded
	Matrix tmp_A(A.GetDim(0), A.GetDim(1));
	for (i=0; i<A.GetDim(1); i++)
	{
		tmp_mat.GetColumn(i, tmp_col);
		tmp_A.SetColumn(i, tmp_col);
	}

	if (!tmp_A.IsIdentity())
	{
		cerr<<"[Error] input matrix for Gauss-Jordan elimination is not invertible"<<endl;
		exit(1);
	}
	else
	{
		for (i=0; i<B.GetDim(1); i++)
		{
			tmp_mat.GetColumn(A.GetDim(1)+i, tmp_col);
			B.SetColumn(i, tmp_col);
		}
	}
}
