#include "pch.h"

vec dxde1(double e, vec x)
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

vec dxde2(double e, vec x)
{
	double sigma = sig0;

	double r = x(0);
	double theta = x(1);
	double phi = x(2);
	double gamma = x(3);
	double psi = x(4);
	double V = sqrt(2 * (1 / r - e));
	vec aero = calcu_aero(r, V);
	double L = aero(0);
	double D = aero(1);
	vec xdot(6);
	xdot(0) = sin(gamma) / D;
	xdot(1) = cos(gamma)*sin(psi) / (r*D*cos(phi));
	xdot(2) = cos(gamma)*cos(psi) / (r*D);
	xdot(3) = 1 / (D*V*V)*(L*cos(sigma) + (V*V / r - 1 / r / r)*cos(gamma));
	xdot(4) = 1 / (D*V)*(L*sin(sigma) / V / cos(gamma) + V / r * cos(gamma)*sin(psi)*tan(phi));
	xdot(5) = 1 / (D*V);
	return xdot;
}