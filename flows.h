#pragma once
#include "gridbuilder.h"
#include "FinalElemetsBasic.h"

/// <summary>
/// ����� ������� ������� ����� �����
/// </summary>
class FlowCulcer
{
	CalculationArea calculationArea;
	double* q;
	int qN;
	double* flows;

public:
	FlowCulcer();
	bool isInnerN(int iLocalFace);
	void init(CalculationArea calculationArea, double* q);
	void calcFlows();
	double* getflows();
	void printFlows(std::ofstream &out);

};
