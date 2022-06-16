#pragma once
#include "sparsematrix.h"

using namespace sparsematrix;
typedef double real;
namespace slae
{
	struct SLAE
	{
		SparseMatrixAsym A;
		SparseMatrixAsym LU;
		double* b;

		int count_LOS(real* x, int maxiter, real eps);
		int count_LOS_Simple(real* x, int maxiter, real eps);
		void init(SparseMatrixAsym A);
		void init(SparseMatrixAsym A, double* b);

		void setOneVariableSolve(int iVar, double varMean);
	private:

		real* r;;
		real* z; 
		real* p; 
		real* f; 
		real* buf_v; 
		real* buf_v1;

		int calc_Lx(SparseMatrixAsym& LU, real* x, real* f, bool is_diag_1);

		int calc_Ux(SparseMatrixAsym& LU, real* x, real* f, bool is_diag_1);	
	};

	int solveSLAU3(double A[3][3], double b[3], double x[3]);
}
