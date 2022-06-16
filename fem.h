#pragma once
#include "slae.h"
//#include "mesh.h"
#include "FinalElemetsBasic.h"
#include "GlobalMatAssembler.h"
#define NON_LINEAR false

//using namespace meshspace;
using namespace slae;

class FEM
{
	SLAE slae;
	//Mesh mesh;
	//NewtAssembler newtAssembler;
	CalculationArea calculationArea;
	SimpleAssembler simpleAssembler;


	double* descr;
	double* bForCheck;
	void addFirstConditions();
	void addSecondConditions();
	//void addEmptyConditions();
	void LinearTask();
#if NON_LINEAR
	void NonLinearTask();
#endif

	/// <summary>
	/// Расчёт невязки
	/// </summary>
	double countDescr();
	double countChange();
	double* qPrev;
public:
	double* q;
	//double solution(double x, double y, double& A, double& Bx, double& By);
	void init(CalculationArea& calculationArea);
	void printSolution(std::ofstream& out);
	double* getSolutWeights();
	void getSolutWeights(double*& qRes);
};

