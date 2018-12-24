#include "pch.h"

//全局变量
double sig0_a = 0;

//主函数
int main()
{
	double sigma0;

	sigma0=calcu_sigma0();

	cout << sigma0 * 57.3 << endl;
	return 0;
}

//计算初始段倾侧角
double calcu_sigma0()
{
	double r0, rf, v0, vf, e0, ef;
	double sigma0 = 0;
	r0 = (Re + H0) / Re;
	rf = (Re + Hf) / Re;
	v0 = V0 / Vc;
	vf = Vf / Vc;
	e0 = 1 / r0 - v0 * v0 / 2;
	ef = 1 / rf - vf * vf / 2;
	vec espan = linspace(e0, ef, 5000);

	double gamma0 = -0.5*pi / 180;//初始弹道倾角
	vec x0 = { r0,gamma0 };
	vec sigma = zeros(2000, 1);
	vec qdot = zeros(2000, 1);
	sigma(0) = 30 * pi / 180;
	sigma(1) = 35 * pi / 180;
	uword i;
	mat ex;
	int num_r = 0;
	double rho = 0;
	double v = 0;
	double qdot_delta = 0;
	for (i = 0; i < 1999; i++)
	{
		sig0_a = sigma(i) / pi * 180 ;
		ex = rk1(espan, x0);
		num_r = ex.n_rows;
		ef = ex(num_r - 1, 0);
		rf = ex(num_r - 1, 1);
		rho = rho0 * exp(-(rf*Re - Re) / hs);
		v = sqrt(2 * (1 / rf - ef));
		qdot(i) = k_q * sqrt(rho)*pow(v, 3.15);
		qdot_delta = qdot(i) - q_max + 20000;
		if (abs(qdot_delta) < 0.001)
		{
			break;
		}
		if (i > 0)
		{
			sigma(i + 1) = magnitude(sigma(i) - (sigma(i) - sigma(i - 1)) / (qdot(i) - qdot(i - 1))*qdot_delta);
		}
	}
	sigma0 = sigma(i);
	//cout.precision(12);
	//cout.setf(ios::fixed);
	//sigma.rows(0,5).raw_print(cout,"");
	return sigma0;
}



vec calcu_aero(double r, double V)
{
	double ma = V * Vc / 340;
	double alpha;
	if (ma >= 15)
	{
		alpha = 45 * pi / 180;
	}
	else
	{
		alpha = (45 - 0.21*(ma - 15)*(ma - 15)) * pi / 180;
	}
	double Cl = cl0 + cl1 * alpha + cl2 * alpha*alpha;
	double Cd = cd0 + cd1 * Cl + cd2 * Cl*Cl;
	double h = r * Re - Re;
	double rho = rho0 * exp(-h / hs);

	double q = 0;
	q = 1.0 / 2 * rho * pow(V,2);
	vec aero(4);
	aero(2) = Cl;
	aero(3) = Cd;
	aero(0) = q * Cl * S / m * Re;
	aero(1) = q * Cd * S / m * Re;

	return aero;
}


int trans(double e, vec x)
{
	int result = 0;
	double r = x(0);
	double gamma = x(1);
	double v = sqrt(2 * (1 / r - e));
	vec aero=calcu_aero(r, v);
	double L = aero(0);
	double D = aero(1);
	double drdv = v * sin(gamma) / (-D - sin(gamma) / r / r);
	double r0 = 1;
	double r1 = newton_r(r0, v, aero(3));
	double rho1 = rho0 * exp(-(r1*Re - Re) / hs);
	double drhodr = -Re / hs * rho1;
	double drdvq = -6.3*rho1 / (drhodr*v);
	const double delta = 0.01;
	if (abs(drdv - drdvq) < delta)
	{
		result = 1;
	}
	else
	{
		result = 0;
	}
	return result;
}

double newton_r(double x0, double v, double Cl)
{
	double f = 0;
	double fdot = 0;
	double x1 = 0.5;
	const double eps = 1e-10;
	while (fabs(x1 - x0) / fabs(x1) >= eps)
	{
		x0 = x1;
		f = 1 / x0 / x0 - v * v / x0 - 1.0 / 2 * rho0*exp(-(x0*Re - Re) / hs)*v*v*Cl*S / m * Re;
		fdot = -2 / pow(x0, 3) + v * v / pow(x0, 2) + 1.0 / 2 * rho0*exp(-(x0*Re - Re) / hs)*v*v*Cl*S / m * Re*Re / hs;
		x1 = x0 - f / fdot;
	}
	return x1;
}


double magnitude(double x)
{
	if (x > 85 * pi / 180)
	{
		x = 85 * pi / 180;
	}
	if (x < 20 * pi / 180)
	{
		x = 20 * pi / 180;
	}
	return x;
}






