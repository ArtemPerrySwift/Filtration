#include "sne.h"
#include "programlog.h"
#include "array.h"
#include "slae.h"

SNE3D::SNE3D()
{
	maxiter = 1000;
	err = 1e-14;
	bettaEps = 1e-14;
	epsDif = 1e-14;
	isNumericalDif = true;
}
void SNE3D::fillJacobMatrix(double J[DIM][DIM], Coord p)
{
	J[0][0] = dEqu1dx(p);
	J[0][1] = dEqu1dy(p);
	J[0][2] = dEqu1dz(p);

	J[1][0] = dEqu2dx(p);
	J[1][1] = dEqu2dy(p);
	J[1][2] = dEqu2dz(p);

	J[2][0] = dEqu3dx(p);
	J[2][1] = dEqu3dy(p);
	J[2][2] = dEqu3dz(p);
}

void SNE3D::countF(double F[DIM], Coord p)
{
	F[0] = equ1(p);
	F[1] = equ2(p);
	F[2] = equ3(p);
}

bool SNE3D::setMaxiter(int maxiter)
{
	if (maxiter < 1)
	{
		programlog::writeErr("SNE3D - Attempting to set the maximum number of iterations parameter to less than 1");
		return false;
	}
	this->maxiter = maxiter;
	return true;
}
bool SNE3D::setErr(double err)
{
	if (err <= 0)
	{
		programlog::writeErr("SNE3D - An attempt to assign a value equal to 0 or less than 0 to the parameter of the permissible error in solving a system of nonlinear equations");
		return false;
	}
	this->err = err;
	return true;
}
bool SNE3D::setBettaEps(double bettaEps)
{
	if (bettaEps <= 0)
	{
		programlog::writeErr("SNE3D - An attempt to assign to the parameter of the minimum allowable change in the solution of a system of nonlinear equations at each iteration a value equal to 0 or less than 0");
		return false;
	}
	this->bettaEps = bettaEps;
	return true;
}
bool SNE3D::setEpsDif(double epsDif)
{
	if (epsDif <= 0)
	{
		programlog::writeErr("SNE3D - An attempt to assign a value equal to 0 or less than 0 to the step parameter of the numerical calculation of derivatives");
		return false;
	}
	this->epsDif = epsDif;
	return true;
}

double SNE3D::dEqu1dx(Coord p)
{
	Coord p1 = p;
	p1.x += epsDif;
	return (equ1(p1) - equ1(p)) / epsDif;
}

double SNE3D::dEqu1dy(Coord p)
{
	Coord p1 = p;
	p1.y += epsDif;
	return (equ1(p1) - equ1(p)) / epsDif;
}

double SNE3D::dEqu1dz(Coord p)
{
	Coord p1 = p;
	p1.z += epsDif;
	return (equ1(p1) - equ1(p)) / epsDif;
}

double SNE3D::dEqu2dx(Coord p)
{
	Coord p1 = p;
	p1.x += epsDif;
	return (equ2(p1) - equ2(p)) / epsDif;
}

double SNE3D::dEqu2dy(Coord p)
{
	Coord p1 = p;
	p1.y += epsDif;
	return (equ2(p1) - equ2(p)) / epsDif;
}

double SNE3D::dEqu2dz(Coord p)
{
	Coord p1 = p;
	p1.z += epsDif;
	return (equ2(p1) - equ2(p)) / epsDif;
}

double SNE3D::dEqu3dx(Coord p)
{
	Coord p1 = p;
	p1.x += epsDif;
	return (equ3(p1) - equ3(p)) / epsDif;
}

double SNE3D::dEqu3dy(Coord p)
{
	Coord p1 = p;
	p1.y += epsDif;
	return (equ3(p1) - equ3(p)) / epsDif;
}

double SNE3D::dEqu3dz(Coord p)
{
	Coord p1 = p;
	p1.z += epsDif;
	return (equ3(p1) - equ3(p)) / epsDif;
}


bool SNE3D::countSolut(Coord initAppr, Coord& answer)
{
	countF(F, initAppr);
	double curErr = arrayspace::scal(F, F, DIM);
	double betta;
	Coord p = initAppr;
	int i;
	if(isNumericalDif)
	{
		double begEpsDif = epsDif;
		for (i = 0; i < maxiter && curErr > err && betta > bettaEps && epsDif > minEpsDif; i++)
		{
			betta = 1;
			fillJacobMatrix(J, p);
			if (slae::solveSLAU3(J, F, ans) == -1)
			{	
				while (epsDif > minEpsDif && slae::solveSLAU3(J, F, ans) == -1)
				{
					epsDif /= 2.0;
					fillJacobMatrix(J, p);
				}

				if (epsDif > minEpsDif)
				{
					epsDif = begEpsDif;
					p = getIterSolutUsingBetta(p, betta, curErr);
				}
				else
				{
					programlog::writeErr("SLAE made in SNU cannot be solved");
					return false;
				}
					
					
			}
			else
				p = getIterSolutUsingBetta(p, betta, curErr);
		}
	}
	else
	{
		for (i = 0; i < maxiter && curErr > err && betta > bettaEps; i++)
		{
			betta = 1;
			fillJacobMatrix(J, p);
			if (slae::solveSLAU3(J, F, ans) == -1)
			{
				programlog::writeErr("SLAE made in SNU cannot be solved");
				return false;
			}
				
			p = getIterSolutUsingBetta(p, betta, curErr);
		}
	}

	if (i == maxiter)
	{
		programlog::writeErr("Solving of SNE reached the limit of iterations");
		return false;
	}

	if (betta <= bettaEps)
	{
		programlog::writeErr("Solving of SNE stocked");
		return false;
	}

	answer = p;
	return true;
}

Coord SNE3D::getIterSolutUsingBetta(Coord p,double& betta, double normF)
{
	Coord p1;
	Coord dp;
	dp = ans;
	double bettaBuf = betta;
	p1 = p + dp * bettaBuf;
	countF(F, p1);
	double normCurF = normCurF = arrayspace::scal(F, F, DIM);
	while (normCurF > normF && bettaBuf > bettaEps);
	{
		bettaBuf /= 2;
		p1 = p + dp * bettaBuf;
		countF(F, p1);
		normCurF = arrayspace::scal(F, F, DIM);
	}

	if (normCurF > normF)
	{
		betta = bettaBuf;
		return p1;
	}
	else
	{
		betta = 0;
		return p;
	}
}

