#include "satculcer.h"
#include "DifferentEquParams.h"

const double SatCulcer::S_MIN = 0.01;
const double SatCulcer::S_MAX = 0.05;

void SatCulcer::init(CalculationArea calculationArea, double dtMax)
{
	this->calculationArea = calculationArea;
	flows = calculationArea.flowStore.flows;
	this->dtMax = dtMax;
}

double SatCulcer::choose_dt()
{
	FinitElement* finitElems = calculationArea.finitElementStore.finitElements;
	int nFinElems = calculationArea.finitElementStore.nFinitElement;

	int i , j;
	FinitElement finitElementBuf;
	Phase phase;
	int nPhase = finitElems[0].phaseStorage.count;

	double dtMin, dt;
	double coeffSum, phaseCoeff;
	double outQSum;

	PhaseOut phaseOut;
	dtMin = dtMax;
	for (i = 0; i < nFinElems; i++)
	{
		finitElementBuf = finitElems[i];
		coeffSum = calcCoeffSum(finitElementBuf.phaseStorage);
		outQSum = calcOutQSum(finitElementBuf);

		for (j = 0; j < nPhase; j++)
		{
			phase = finitElementBuf.phaseStorage.phases[j];

			if (phase.saturation < S_MIN)
			{
				phaseOut.iFinElem = i;
				phaseOut.iPhase = j;
				phaseSatSmin.push_back(phaseOut);
				continue;
			}

			if (phase.saturation < S_MAX)
			{
				phaseOut.iFinElem = i;
				phaseOut.iPhase = j;
				phaseSatSmax.push_back(phaseOut);
				continue;
			}

			phaseCoeff = calcPhaseCoeff(finitElementBuf.phaseStorage.phases[j], coeffSum);
			dt = (finitElementBuf.sqare * finitElementBuf.FI * phase.saturation) / (phaseCoeff * outQSum);

			if (dt < dtMin) dtMin = dt;
		}
	}

	return dtMin;
}

double SatCulcer::calcOutQSum(FinitElement finitElement)
{
	double sum = 0;
	for (int i = 0; i < FACES_NUM; i++)
	{
		if (finitElement.flowSign[i] > 0)
			sum += flows[finitElement.faces[i]];
	}

	return sum;
}

double SatCulcer::calcPhaseCoeff(Phase phase, double coeffSum)
{
	return (DifferentEquParams::km(phase.saturation) / phase.dynamicVisc) / coeffSum;
}

double SatCulcer::calcCoeffSum(PhaseStorage phaseStorage)
{
	int nPhase = phaseStorage.count;
	double sum = 0;
	Phase phase;
	for (int i = 0; i < nPhase; i++)
	{
		phase = phaseStorage.phases[i];
		sum += DifferentEquParams::km(phase.saturation) / phase.dynamicVisc;
	}
	
	return sum;
}

double  SatCulcer::calcAllVol(double dt)
{
	FinitElement* finitElems = calculationArea.finitElementStore.finitElements;
	int nFinElems = calculationArea.finitElementStore.nFinitElement;
	calcBorderVol(calculationArea.borderFacesStore, dt);
	ContainerPhaseVol containerPhaseVol = calculationArea.containerPhaseVol;
	int i, j, k;
	FinitElement finitElementBuf;
	Phase phase;
	int nPhase = finitElems[0].phaseStorage.count;
	int nOutPhase;
	int volInd[FACES_NUM];
	Phase phasesL[FACES_NUM];

	double coeffSum, phaseCoeff;

	PhaseVolStore* phaseVolStore = calculationArea.containerPhaseVol.phaseVolStore;
	PhaseVolStore phaseVolStoreMix = calculationArea.containerPhaseVol.phaseVolStoreMix;
	int nVolStore = calculationArea.containerPhaseVol.count;

	for (i = 0; i < nFinElems; i++)
	{
		finitElementBuf = finitElems[i];
		coeffSum = calcCoeffSum(finitElementBuf.phaseStorage);
		nOutPhase = 0;
		for (j = 0; j < FACES_NUM; j++)
		{
			if (finitElementBuf.flowSign[j] > 0)
			{
				volInd[nOutPhase] = finitElementBuf.faces[j];
				phasesL[nOutPhase] = finitElementBuf.phaseStorage.phases[j];
				nOutPhase++;
			}
		}

		for (j = 0; j < nOutPhase; j++)
		{
			for (k = 0; k < nPhase; k++)
			{
				phaseCoeff = calcPhaseCoeff(finitElementBuf.phaseStorage.phases[k], coeffSum);
				phaseVolStore[k].PhaseVol[volInd[j]] = phaseCoeff * flows[volInd[j]] * dt;
			}
			phaseVolStoreMix.PhaseVol[volInd[j]] = flows[volInd[j]] * dt;
		}
	}
	return 0;
}

void SatCulcer::pushOutPhases()
{
	int nSmax = phaseSatSmax.size();
	int i, j;
	PhaseOut phaseOutBuf;
	FinitElement* finitElems = calculationArea.finitElementStore.finitElements;
	FinitElement finitElemBuf;
	PhaseVolStore* phaseVolStore = calculationArea.containerPhaseVol.phaseVolStore;
	double sumV;
	double Vm;

	for (i = 0; i < nSmax; i++)
	{
		phaseOutBuf = phaseSatSmax[i];
		finitElemBuf = finitElems[phaseOutBuf.iFinElem];
		sumV = 0;
		for (j = 0; j < FACES_NUM; j++)
			if (finitElemBuf.flowSign[j] > 0)
				sumV += phaseVolStore[phaseOutBuf.iPhase].PhaseVol[finitElemBuf.faces[j]];
		Vm = finitElemBuf.sqare * finitElemBuf.FI * finitElemBuf.phaseStorage.phases[phaseOutBuf.iPhase].saturation;
		if(sumV > Vm)
			phaseOut.insert(phaseSatSmax[i]);
	}
		
	phaseSatSmax.clear();

	int nSmin = phaseSatSmin.size();

	for (i = 0; i < nSmin; i++)
		phaseOut.insert(phaseSatSmin[i]);

	phaseSatSmin.clear();

	PhaseVolStore phaseVolStoreMix = calculationArea.containerPhaseVol.phaseVolStoreMix;
	int nPhases = finitElems[0].phaseStorage.count;
	bool* outPhasesInd = new bool[nPhases]; // ������ ����� ��� ������� ������������� ��� � �������� ��������
	int nOutPhase = 0; // ���������� ������������� ��� � �������� ��������
	int iElem; // ����� �������� � ������� ����� ������������� ����
	

	int k, l;
	int iPhase; // ����� ������������� ����
	double sumCoeffPhase, // ����� ������������ ������������� ��� ���, ������� �� �������������
		del, // �������� ��� ������� ������� ������������� � �������� ��������
		outQSum; // ����� ��������� �������
	
	int* faceInd; //������ ���������� ������� ������ ��������� ��������
	int iFace; // ���������� ����� ����� ��������� ��������
	double faceVolBuf; // ������ ��� �������� ������ ������������� ��������

	double poreVol; // ����� ��� � �������� ��������

	Phase* phasesEl; // ������ ��� � ��������

	std::set<PhaseOut>::iterator elem = phaseOut.begin();
	nOutPhase = 0;
	iElem = elem->iFinElem;
	l = 0;
	for (; elem != phaseOut.end(); elem++)
	{
		if (iElem == elem->iFinElem)
		{
			iPhase = elem->iPhase;
			for (; l < iPhase; l++)
				outPhasesInd[l] = false;
			outPhasesInd[l] = true;
			l++;
		}
		else // ����� ����������� ���������� � ������������� ����� � ������� �������� ��������
		{
			finitElemBuf = finitElems[iElem];
			poreVol = finitElemBuf.FI * finitElemBuf.sqare;
			del = calcCoeffSum(finitElemBuf.phaseStorage);
			sumCoeffPhase = 0;
			for (i = 0; i < nPhases; i++)
				sumCoeffPhase += outPhasesInd[i] ? 0: calcPhaseCoeff(finitElemBuf.phaseStorage.phases[i], del);

			outQSum = calcOutQSum(finitElemBuf);
			
			faceInd = finitElemBuf.faces;
			phasesEl = finitElemBuf.phaseStorage.phases;
			for (i = 0; i < FACES_NUM; i++)
			{
				iFace = faceInd[i];
				for (j = 0; j < nPhases; j++)
				{
					if (outPhasesInd[j])
					{
						faceVolBuf = phaseVolStore[j].PhaseVol[iFace] = flows[iFace] / outQSum * poreVol * phasesEl[j].saturation;
						phaseVolStoreMix.PhaseVol[iFace] -= faceVolBuf;
					}
				}
			}

			for (i = 0; i < FACES_NUM; i++)
			{
				iFace = faceInd[i];
				for (j = 0; j < nPhases; j++)
				{
					if (!outPhasesInd[j])
						phaseVolStore[j].PhaseVol[iFace] = calcPhaseCoeff(finitElemBuf.phaseStorage.phases[i], del) * sumCoeffPhase;
				}
			}
		}
	}
}

void SatCulcer::calcNewSat()
{
	FinitElement* finitElems = calculationArea.finitElementStore.finitElements;
	PhaseVolStore* phaseVolStore = calculationArea.containerPhaseVol.phaseVolStore;

	int nFinEl = calculationArea.finitElementStore.nFinitElement; 
	int i, j, k;

	int nPhases = finitElems[0].phaseStorage.count; // ���������� ���

	FinitElement finitElemBuf;
	double poreVol; // ����� ��� � �������� ��������
	int* elFaces; // ������ ���������� ������� ������ ��������� ��������
	Phase* elPhases; // ������ ��� ��������� ��������
	Phase jElPhase; // J-�� ���� ��������� ��������
	signed char* elFlowSign;

	double phaseElVol; // ����� ����� ���� � ������
	for (i = 0; i < nFinEl; i++)
	{
		finitElemBuf = finitElems[i];
		poreVol = finitElemBuf.FI * finitElemBuf.sqare;
		elFaces = finitElemBuf.faces;
		elPhases = finitElemBuf.phaseStorage.phases;
		elFlowSign = finitElemBuf.flowSign;

		for (j = 0; j < nPhases; j++)
		{
			jElPhase = elPhases[j];
			phaseElVol = poreVol * jElPhase.saturation;

			for (k = 0; k < FACES_NUM; k++)
				phaseElVol -= elFlowSign[k] * phaseVolStore[j].PhaseVol[elFaces[k]];

			elPhases[j].saturation = phaseElVol / poreVol;
		}
	}
}

void SatCulcer::calcBorderVol(BorderFacesStore borderFacesStore, double dt)
{
	int* iFaces = borderFacesStore.iFaces;
	int nFaces = borderFacesStore.nFaces;
	PhaseStorage phaseStorage = borderFacesStore.phaseStorage;
	Phase* phases = phaseStorage.phases;
	int nPhases = phaseStorage.count;
	double coeffSum = calcCoeffSum(phaseStorage);
	double phaseCoeff;
	int i, j;
	PhaseVolStore phaseVolStoreMix = calculationArea.containerPhaseVol.phaseVolStoreMix;
	PhaseVolStore* phaseVolStore = calculationArea.containerPhaseVol.phaseVolStore;
	for (i = 0; i < nPhases; i++)
	{
		phaseCoeff = calcPhaseCoeff(phases[i], coeffSum);
		for (j = 0; j < nFaces; j++)
		{
			phaseVolStore[i].PhaseVol[iFaces[j]] = phaseCoeff * 2 * flows[iFaces[j]] * dt;
		}
	}

	for (j = 0; j < nFaces; j++)
	{
		phaseVolStoreMix.PhaseVol[iFaces[i]] = 2 * flows[iFaces[i]] * dt;
	}
}

double SatCulcer::reculcSat()
{
	double dt = choose_dt();
	calcAllVol(dt);
	pushOutPhases();
	calcNewSat();
	return dt;
}

void SatCulcer::printSat(std::ofstream &out)
{
	int nPhases = calculationArea.nPhases;
	int i;
	int nFinElems = calculationArea.finitElementStore.nFinitElement;
	FinitElement* finitElements = calculationArea.finitElementStore.finitElements;
	out << "___________________________________" << std::endl;
	for (i = 0; i < nPhases; i++)
		out << "____________________";

	out << "|----N of Face----|-----Flow------|" << std::endl;
	out.width(16);
	out.setf(out.left);
	out.fill('-');
	for (i = 0; i < nPhases; i++)
		out << "----Phase " << i << "|";
	out << "___________________________________" << std::endl;
	for (i = 0; i < nPhases; i++)
		out << "____________________";

	out.width(15);
	int j;
	Phase* phasesFinElem;
	for (int i = 0; i < nFinElems; i++)
	{
		phasesFinElem = finitElements[i].phaseStorage.phases;
		out << "|" << i << "|";
		for(j = 0; j < nPhases; j++)
			out << phasesFinElem[j].saturation << "|" << std::endl;
		out << std::endl;
	}
	out << "___________________________________" << std::endl;
	for (i = 0; i < nPhases; i++)
		out << "____________________";
}

bool operator <(const PhaseOut& left, const PhaseOut& right)
{
	if (left.iFinElem == right.iFinElem)
		return left.iPhase < right.iPhase;
	return left.iFinElem < right.iFinElem;
}