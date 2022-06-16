#include "finitelemculcer.h"
#include "array.h"
double FinitElemCulcer::culcPhi(ECoord p)
{
	double phi = 0;
	for (int i = 0; i < N; i++)
		phi += phi_i(p, i)*q[i];

	return phi;
}

ECoord FinitElemCulcer::culcGrad(ECoord p)
{
	ECoord grad(0, 0, 0);
	for (int i = 0; i < N; i++)
		grad += countGrad_i(p, i) * q[i];
	return grad;
}

void FinitElemCulcer::init(FinitElement finitElement, CoordStorage coordsStore, double* q)
{
	this->finitElement = finitElement;

	coordsStore.copyElemsByInd(finitElement.ver, coords, VER_NUM);
	arrayspace::copyElemsByInd(finitElement.ver, this->q, q, VER_NUM);
}

double FinitElemCulcer::equ1(Coord p)
{

}
double FinitElemCulcer::equ2(Coord p)
{

}
double FinitElemCulcer::equ3(Coord p)
{

}