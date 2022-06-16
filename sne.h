#pragma once
#include "coord.h"
class SNE3D
{
public:
	SNE3D();
	bool countSolut(Coord initAppr, Coord& ans);
	bool setMaxiter(int maxiter);
	bool setErr(double err);
	bool setBettaEps(double bettaEps);
	bool setEpsDif(double epsDif);
	bool setMinEpsDif(double minEpsDif);

	bool isNumericalDif;
private:
	static const int DIM = 3;
	double J[DIM][DIM];
	double F[DIM];
	double ans[DIM];
	int maxiter;
	double err;
	double bettaEps;
	double epsDif; // Шаг чисенного расчёта производной
	double minEpsDif;

	void fillJacobMatrix(double J[DIM][DIM], Coord p);
	void countF(double F[DIM], Coord p);
	Coord getIterSolutUsingBetta(Coord p, double& betta, double normF);

	virtual double equ1(Coord p) = 0;
	virtual double equ2(Coord p) = 0;
	virtual double equ3(Coord p) = 0;

	virtual double dEqu1dx(Coord p);
	virtual double dEqu2dx(Coord p);
	virtual double dEqu3dx(Coord p);

	virtual double dEqu1dy(Coord p);
	virtual double dEqu2dy(Coord p);
	virtual double dEqu3dy(Coord p);

	virtual double dEqu1dz(Coord p);
	virtual double dEqu2dz(Coord p);
	virtual double dEqu3dz(Coord p);
};
