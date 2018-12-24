#include "pch.h"

vec f1(double e, vec x)
{
	double sigma = sig0_a * pi / 180 ;
	double r = x(0);
	double gamma = x(1);
	double V = sqrt(2 * (1 / r - e));
	vec aero = calcu_aero(r, V);
	double L = aero(0);
	double D = aero(1);
	vec xdot(2);
	xdot(0) = sin(gamma) / D;
	xdot(1) = 1 / (D*V*V)*(L*cos(sigma) + (V*V / r - 1 / r / r)*cos(gamma));
	return xdot;
}

