#include "DifferentEquParams.h"

double DifferentEquParams::f(double x, double y, double z)
{
	return 0;
}

double DifferentEquParams::u1(double x, double y, double z)
{
	return 100;
}

double DifferentEquParams::du_dn(double x, double y, double z)
{
	return -100;
}

double DifferentEquParams::lambda(double K, Phase* phases, int nPhases)
{
	double sum = 0;
	for (int i = 0; i < nPhases; i++)
	{
		sum += km(phases[i].saturation) / phases[i].dynamicVisc;
	}

	return K * sum;
}

double DifferentEquParams::km(double str)
{
	return str;
}
