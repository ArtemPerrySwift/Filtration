#include "fem.h"
#include "DifferentEquParams.h"
#include "array.h"
#include "GlobalMatAssembler.h"
#include "gridbuilder.h"
#include "facebasic.h"
#include <iostream>
#include <iomanip>
#include <fstream>

using namespace arrayspace;

void FEM::init(CalculationArea& calculationArea)
{
	this->calculationArea = calculationArea;
	SparseMatrixAsym A;
	std::cout << "Построение портрета матрицы" << std::endl;
	A.buildPortrait(calculationArea.coordsStore, calculationArea.finitElementStore);
	slae.init(A);
	q = new double[A.n];
	//qPrev = new double[A.n];
	//bForCheck = new double[A.n];
	
	//descr = bForCheck;
	arrayspace::fill_vec(q, A.n, 0.0);
	LinearTask();
}

void FEM::addFirstConditions()
{
	std::cout << "Учёт первых краевых условий" << std::endl;
	int* Iknots = calculationArea.knots1Cond.IKnots;
	int kt1 = calculationArea.knots1Cond.kt1;

	Coord* knots = calculationArea.coordsStore.coords;
	double uMean;
	int ind;

	for (int i = 0; i < kt1; i++)
	{
		ind = Iknots[i];
		uMean = DifferentEquParams::u1(knots[ind].x, knots[ind].y, knots[ind].z);
		slae.setOneVariableSolve(ind, uMean);
	}
}
/*
void FEM::addEmptyConditions()
{
	std::cout << "Учёт фиктивных узлов" << std::endl;
	int* Iknots = calculationArea.knotsEmptyCond.IKnots;
	int kt1 = calculationArea.knotsEmptyCond.kt1;

	Coord* knots = calculationArea.coordsStore.coords;
	int ind;

	for (int i = 0; i < kt1; i++)
	{
		ind = Iknots[i];
		slae.setOneVariableSolve(ind, 0);
	}
}
*/
void FEM::addSecondConditions()
{
	std::cout << "Учёт вторых краевых условий" << std::endl;
	Face* faces = calculationArea.Faces2CondStore.faces;
	int n = calculationArea.Faces2CondStore.count;
	int n_basic = FaceBasic2Cond::N_BASIC;
	FaceBasic2Cond faceBasic2Cond;
	double LocalB[FaceBasic2Cond::N_BASIC];

	int* iFaceCoord;
	int i, j;

	double* b = slae.b;
	for (i = 0; i < n; i++)
	{
		faceBasic2Cond.init(faces[i], calculationArea.coordsStore);
		faceBasic2Cond.calcAddLocalB(LocalB);
		iFaceCoord = faces[i].knots;
		for (j = 0; j < n_basic; j++)
		{
			/*
			if (iFaceCoord[j] == 7670)
			{
				for (int k = 0; k < n_basic; k++)
				{
					std::cout << "nan " << LocalB[k] << " ";
				}
				std::cout << std::endl;
			}
			*/
			//if (!i) std::cout << LocalB[j] << std::endl;
			b[iFaceCoord[j]] += LocalB[j];
		}
			
	}
}


void FEM::LinearTask()
{
	//if (binfile_in.gcount() != sizeof(double)) return false;
	slae.A.fillMatrix(0.0);
	fill_vec(slae.b, slae.A.n, 0);
	/*
	cout << "di Before" << endl;
	for (int i = 0; i < 50; i++)
		cout << slae.A.di[i] << " ";
	cout << endl;
	*/
	std::cout << "Сборка глобальной матрицы" << std::endl;
	simpleAssembler.assembleGlobalMatrix(slae, calculationArea, q);
	
	
	
	/*
	cout << "ggl " << endl;
	int n_gg = slae.A.ig[slae.A.n];
	for (int i = 0; i < n_gg; i++)
		cout << slae.A.di[i] << " ";
	cout << endl;
	*/
	//slae.A.printFullMatrix();
	addSecondConditions();
	//addEmptyConditions();
	addFirstConditions();
	//slae.A.printFullMatrix();
	/*
	cout << "di" << endl;
	for (int i = 0; i < slae.A.n; i++)
		cout << slae.A.di[i] << std::endl;
	cout << endl;
	*/
	/*
	std::ofstream fileCheck;
	fileCheck.open("CheckFile.txt");
	double sum = 0;

	
	fileCheck << "b" << std::endl;
	for (int i = 0; i < slae.A.n; i++)
	{
		sum += slae.b[i];
		fileCheck << "b[" << i << "] = " << slae.b[i] << "; sum = " << sum << std::endl;
		//std::cout << slae.b[i] << std::endl;
	}
	fileCheck.close();
	std::cout << endl;
	*/
	//double sum = 0;
	//for (int i = 0; i < slae.A.n; i++)
	//	sum += slae.b[i];
	//std::cout << "Sum " <<  sum << endl;
	//double* b = new double[slae.A.n];
	//double* bPrev = new double[slae.A.n];

	/*
	ifstream binfile_in;
	binfile_in.open("v2.dat", ios::binary);
	for (int i = 0; i < mesh.knotsStorage.kuzlov; i++)
		binfile_in.read((char*)&(q[i]), sizeof(double));
	arrayspace::copy(bPrev, q, slae.A.n);
	*/
	//slae.A.mult(q, bPrev);
//	slae.A.mult(q, b);

//	for (int i = 17; i < mesh.knotsStorage.kuzlov; i+=35)
//	{
//		slae.b[i] = b[i];
//	}
	
	std::cout << "Решение СЛАУ" << std::endl;
	//slae.count_LOS(q, 100000, 1e-14);
	slae.count_LOS_Simple(q, 100000, 1e-14);
	/*
	slae.A.fillMatrix(0.0);
	fill_vec(slae.b, slae.A.n, 0);
	simpleAssembler.assembleGlobalMatrix(slae, calculationArea, q);
	addFirstConditions();
	slae.A.mult(q, b);
	/*
	cout << "B check" << endl;
	double sum1, sum2;
	for (int i = 0; i < slae.A.n; i++)
	{
		cout << i << " " << mesh.knotsStorage.knots[i].x << " " << mesh.knotsStorage.knots[i].y << " " << bPrev[i] << " " << q[i] << " " << (bPrev[i] - q[i]) / bPrev[i] * 100 << "%" << endl;
	}
		
	cout << endl;
	*/
}

double* FEM::getSolutWeights()
{
	double* qRes = new double[slae.A.n];
	arrayspace::copy(qRes, q, slae.A.n);

	return qRes;	
}

void FEM::getSolutWeights(double* &qRes)
{
	arrayspace::copy(qRes, q, slae.A.n);
}

#if NON_LINEAR
void FEM::NonLinearTask()
{
	double nonLinearErrMax = 1e-12;
	double nonLinearErr = nonLinearErrMax * 10;
	int maxiter = 1000;
	double minSolutionChange = 1e-10;
	double solutionChange = minSolutionChange * 10;
	for (int i = 0; i < maxiter && nonLinearErr > nonLinearErrMax && solutionChange > minSolutionChange; i++)
	{
		arrayspace::copy(qPrev, q, slae.A.n);
		slae.A.fillMatrix(0.0);
		newtAssembler.assembleGlobalMatrix(slae, mesh, q);
		addFirstConditions();
		slae.count_LOS(q, 1000, 1e-12);
		nonLinearErr = countDescr();
		solutionChange = countChange();
		cout << "Iteration " << i << " Descr = " << nonLinearErr << " Solution Change = " << solutionChange << endl;
	}
}
#endif

double FEM::countDescr()
{
	simpleAssembler.assembleGlobalMatrix(slae, calculationArea, q);
	
	addFirstConditions();
	
	slae.A.mult(q, bForCheck);
	arrayspace::minus(slae.b, bForCheck, descr, slae.A.n);
	return abs(scal(descr, descr, slae.A.n) / scal(slae.b, slae.b, slae.A.n));
}

/*
double FEM::solution(double x, double y, double&A, double&Bx, double&By)
{
	FinalElemBase finalElemBase;
	int materialInd;
	Material material;
	for (int i = 0; i < mesh.finalElementStorage.ktr; i++)
	{
		materialInd = mesh.finalElementMaterialStorage.finalElemMaterials[i];
		material = mesh.materialStorage.findMaterialByNum(materialInd);
		finalElemBase.setFinalElement(mesh.finalElementStorage.endElements[i], mesh.knotsStorage, material, q);
		if (finalElemBase.isIn(x, y))
		{
			A = finalElemBase.countSolut(x, y);
			Bx = -finalElemBase.du_dy(x, y);
			By = finalElemBase.du_dx(x, y);
			return 0;
		}
			
	}
	return 0;
}
*/

double FEM::countChange()
{
	double* dif = qPrev;
	arrayspace::minus(q, qPrev, dif, slae.A.n);

	return abs(scal(dif, dif, slae.A.n) / scal(q, q, slae.A.n));
}

void FEM::printSolution(std::ofstream &out)
{
	int nWeigts = calculationArea.coordsStore.count;
	Coord* coords = calculationArea.coordsStore.coords;
	out << "_________________________________________________________________" << std::endl;
	out << "|-------x-------|-------y-------|-------z-------|---Solution----|" << std::endl;
	out << "_________________________________________________________________" << std::endl;
	//out.width(15);
	for (int i = 0; i < nWeigts; i++)
	{
		out << "|";
		out.width(15);
		out << coords[i].x;
		out << "|";
		out.width(15);
		out << coords[i].y;
		out << "|";
		out.width(15);
		out << coords[i].z; 
		out << "|";
		out.width(15);
		out << q[i];
		out << "|" << std::endl;
	}
		
	out << "_________________________________________________________________" << std::endl;
}



