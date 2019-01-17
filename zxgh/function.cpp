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

double dvds1(double s, double v)
{
	
	double sigma0=sig0_a/180*pi;
	double sigma;
	double vm = 5600 / Vc;
	double vf = Vf / Vc;
	if (v > vm)
	{
		sigma = limit(sigma0 + (sigma_mid - sigma0) / (vm - vv0)*(v - vv0), v);
	}
	else
	{
		sigma = limit(sigmaf + (sigma_mid - sigmaf) / (vm - vf)*(v - vf), v);
	}


	double r = 1.01;
	vec aero = calcu_aero(r, v);
	double Cl = aero(2);
	double Cd = aero(3);
	double vdot = (1 / r - v * v)*(Cd / Cl) / (v*cos(delta_psi)*cos(sigma));
	return vdot;
}

vec dvds2(double s, vec x,double sigma)
{
	double r = x(0);
	double v = x(1);
	double gamma = x(2);
	

	vec aero = calcu_aero(r, v);
	double Cl = aero(2);
	double Cd = aero(3);
	double D = aero(1);
	double L = aero(0);

	double rdot = -r * tan(gamma);
	//double vdot = (-D - (sin(gamma) / r / r)) / (-v * cos(gamma)/ r);
	double vdot = D * r / (v * cos(gamma)) + tan(gamma) / v / r;
	//double gammadot = (1 / v * (L * cos(sigma) + (v * v / r - 1 / r / r)*cos(gamma))) / (-v * cos(gamma) *cos(delta_psi) / r);
	double gammadot = -1 / v / v * (L*r*cos(sigma) / cos(gamma) + (v*v - 1 / r));

	vec xdot = { rdot,vdot,gammadot };
	return xdot;

}
