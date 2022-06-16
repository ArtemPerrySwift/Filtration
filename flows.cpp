#include "flows.h"
#include "array.h"
#include "facebasic.h"
#include "DifferentEquParams.h"
#include "flowelemculcer.h"
#include <iostream>

FlowCulcer::FlowCulcer()
{
	qN = 0;
	q = new double[qN];
	flows = NULL;
}

void FlowCulcer::init(CalculationArea calculationArea, double* q)
{
	this->calculationArea = calculationArea;
	int n = calculationArea.coordsStore.count;
	flows = calculationArea.flowStore.flows;

	if (qN != n)
	{
		delete[] this->q;
		this->q = new double[n];
		qN = n;
	}

	arrayspace::copy(this->q, q, n);

	arrayspace::fill_vec(flows, calculationArea.flowStore.count, 0);

}

bool FlowCulcer::isInnerN(int iLocalFace)
{
	switch (iLocalFace)
	{
	case 0:
		return true;
	case 1:
		return false;
	case 2:
		return true;
	case 3:
		return true;
	case 4:
		return false;
	case 5:
		return false;
	}
}

void FlowCulcer::calcFlows()
{
	FinitElement* finitElements = calculationArea.finitElementStore.finitElements;
	int nElems = calculationArea.finitElementStore.nFinitElement;
	Face* faces = calculationArea.faceStore.faces;
	int i, j;
	int iFace;
	FaceBasicFlow faceBasicFlow;
	FinitElement finitElement;
	double flow;

	//double flowSigns[FACES_NUM];
	double flowLoc[FACES_NUM];
	FlowElemCulcer flowElemCulcer;
	for (i = 0; i < nElems; i++)
	{
		finitElement = finitElements[i];
		flowElemCulcer.init(finitElement, calculationArea.coordsStore, q);
		flowElemCulcer.getFlows(flowLoc, finitElements[i].flowSign);
		/*if (!i)
		{
			std::cout << "Flow of 0 element" << std::endl;
			for (j = 0; j < FACES_NUM; j++)
			{
				std::cout << flowLoc[j] << " ";
			}
			std::cout << std::endl;
		}*/
		
		for (j = 0; j < FACES_NUM; j++)
		{
			iFace = finitElement.faces[j];
			flows[iFace] += flowLoc[j] / 2;
			//faceBasicFlow.init(faces[iFace], calculationArea.coordsStore, q, isInnerN(i), DifferentEquParams::lambda(finitElement.K, finitElement.phases.phases, finitElement.phases.count));
			//low = faceBasicFlow.getFlow();
			//
			//flowSigns[j] = flow < 0 ? -1: 1;
			
		}

		//arrayspace::copy(finitElements[i].flowSign, flowSigns, FACES_NUM);
		
	}
}

double* FlowCulcer::getflows() 
{
	double* flowsRes;
	flowsRes = new double[calculationArea.faceStore.count];
	return flowsRes;
}

void FlowCulcer::printFlows(std::ofstream &out)
{
	out << "____________________________________" << std::endl;
	out << "|----N of Face----|------Flow-------|" << std::endl;
	out << "____________________________________" << std::endl;
	out.width(17);
	for (int i = 0; i < qN; i++)
		out << "|" << i << "|" << flows[i] << "|" << std::endl;
	out << "___________________________________" << std::endl;
}