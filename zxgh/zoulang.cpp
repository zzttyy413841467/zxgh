#include "pch.h"

void zoulang()
{
	vec V = linspace(850, 7800, 140);
	vec aero(4);
	mat h_bianjie(140, 4);
	for (uword i = 0; i < 140; i++)
	{
		aero = calcu_aero(1, V(i) / Vc);
		h_bianjie(i, 0) = 2 * hs*log(sqrt(1.225)*pow(V(i), 3.15)*9.4369e-5 / q_max);
		h_bianjie(i, 1) = hs * log(1.225*pow(V(i), 2) / 2 / 1.5e4);
		//ЗЈЯђЙ§ди
		h_bianjie(i, 2) = hs * log(1.225*S*V(i)*V(i)*sqrt(aero(2)*aero(2) + aero(3)*aero(3)) / 2 / 2.5 / g0 / m);
		h_bianjie(i, 3) = newton_h(1, V(i) / Vc, aero(2));
	}
	h_bianjie.save("E:/zaixianguihua/pre_corre/h_bianjie.data", raw_ascii);
}

double newton_h(double x0, double v, double Cl)
{
	double f = 0;
	double fdot = 0;
	double x1 = 1.01;
	const double eps = 1e-10;
	while (fabs(x1 - x0) / fabs(x1) >= eps)
	{
		x0 = x1;
		f = 1 / x0 / x0 - v * v / x0 - 1.0 / 2 * rho0*exp(-(x0*Re - Re) / hs)*v*v*Cl*S / m * Re *cos(10 * pi / 180);
		fdot = -2 / pow(x0, 3) + v * v / pow(x0, 2) + 1.0 / 2 * rho0*exp(-(x0*Re - Re) / hs)*v*v*Cl*S / m * Re * Re / hs * cos(10 * pi / 180);
		x1 = x0 - f / fdot;
	}
	return x1*Re-Re;
}