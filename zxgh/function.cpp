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

vec dxdt1(double t, vec x)//³õÊ¼¶Î·ÂÕæ
{
	double sigma = sig0;

	double r = x(0);
	double theta = x(1);
	double phi = x(2);
	double V = x(3);
	double gamma = x(4);
	double psi = x(5);
	vec aero = calcu_aero(r, V);
	double L = aero(0);
	double D = aero(1);
	vec xdot(6);
	xdot(0) = V * sin(gamma);
	xdot(1) = V * cos(gamma)*sin(psi) / (r*cos(phi));
	xdot(2) = V * cos(gamma)*cos(psi) / r;
	xdot(3) = -D - (sin(gamma) / r / r);
	xdot(4) = 1 / V * (L * cos(sigma) + (V * V / r - 1 / r / r)*cos(gamma));
	xdot(5) = (L*sin(sigma) / V / cos(gamma) + V / r * cos(gamma)*sin(psi)*tan(phi));
	return xdot;
}