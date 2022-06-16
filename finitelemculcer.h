#pragma once
#include "basicfunct3d.h"
#include "gridbuilder.h"
#include "sne.h"

class FinitElemCulcer: public BasicFunct3D, private SNE3D
{
public:
	double culcPhi(ECoord p);
	ECoord culcGrad(ECoord p);
	virtual void init(FinitElement finitElem, CoordStorage coordsStore, double* q);

protected:
	double q[N];
	FinitElement finitElement;

private:
	double equ1(Coord p) override;
	double equ2(Coord p) override;
	double equ3(Coord p) override;
};
