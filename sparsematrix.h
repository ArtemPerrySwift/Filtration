#pragma once
//#include "mesh.h"
#include "gridbuilder.h"

namespace sparsematrix
{
	struct SparseMatrixSym;

	struct SparseMatrixPortrait
	{
		int* ig;
		int* jg;
		double* di;
		int n;
		SparseMatrixPortrait();
		//SparseMatrixPortrait(SparseMatrixPortrait& portarait);

		double* operator * (double* v);
		void mult(double* v, double* res);
		void buildPortrait(CoordStorage coordStorage, FinitElementStore finitElementStore);
		virtual void fillMatrix(double mean) = 0;
	protected:
		
		void copy(SparseMatrixPortrait& M);
		virtual void multMatElemAndVecElem(int i, int j, int ij, double* v, double*res) = 0;
		//virtual void allocateMemoryForElems() = 0;	
		virtual void allmemory() = 0;

		
		//virtual void printInFullFormat() = 0;
	};

	struct SparseMatrixAsym :  SparseMatrixPortrait
	{
		double* ggl;
		double* ggu;

		SparseMatrixAsym();

		void copy(SparseMatrixAsym& M);

		void copy(SparseMatrixSym& M);

		int decomp_mat_LU(SparseMatrixAsym& LU);

		void printFullMatrix();

		void fillMatrix(double mean) override;
	private:

		void multMatElemAndVecElem(int i, int j, int ij, double* v, double* res) override;
		void allocateMemoryForElems();
		virtual void allmemory() override;
		
	};

	struct SparseMatrixSym : SparseMatrixPortrait
	{
		double* gg;

		SparseMatrixSym();

		void copy(SparseMatrixSym& M);

		//operator SparseMatrixAsym() const;

		void fillMatrix(double mean) override;

	private:
		void multMatElemAndVecElem(int i, int j, int ij, double* v, double* res) override;
		void allocateMemoryForElems();
		virtual void allmemory() override;
	};

	
}

