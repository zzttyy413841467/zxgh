#include "pch.h"

//全局变量
double sig0_a = 0;
Noflyzone nfz1, nfz2;
double sig0 = 0;
double delta_psi;
double sigma_mid;
double vv0;
int item = 0;
vec sigma_qeg(501);
double sigmaf;


//主函数
int main()
{
	double sigma0;

	nfz1.theta = 25.0 / 180 * pi;
	nfz1.phi = 8.5 / 180 * pi;
	nfz1.radius = 2.3 / 180 * pi;
	nfz2.theta = 45.0 / 180 * pi;
	nfz2.phi = 20.0 / 180 * pi;
	nfz2.radius = 3.8 / 180 * pi;

	double ee0 = 1.0 / ((Re + H0) / Re) - pow(V0 / Vc, 2) / 2;
	double eef = 1.0 / ((Re + Hf) / Re) - pow(Vf / Vc, 2) / 2;

	double r0 = (Re + H0) / Re;
	double v0 = V0 / Vc;
	double theta0 = 0;
	double phi0 = 0;
	double gamma0 = -0.5*pi / 180;
	double psi0 = 65 * pi / 180;
	double tau0 = 0;
	double psi_t = 0;

	double rf = (Re + Hf) / Re;
	double vf = Vf / Vc;
	vec rv = calcu_aero(rf, vf);
	sigmaf = acos((1 / rf / rf - vf * vf / rf) / rv(0));

	if (phif - phi0 > 0)
	{
		psi_t= asin(cos(phif)*sin(thetaf - theta0) / sin(acos(cos(phi0)*cos(phif)*cos(theta0 - thetaf) + sin(phi0)*sin(phif))));
	}
	else
	{
		psi_t = pi - asin(cos(phif)*sin(thetaf - theta0) / sin(acos(cos(phi0)*cos(phif)*cos(theta0 - thetaf) + sin(phi0)*sin(phif))));
	}

	sigma0 = calcu_sigma0();

	sig0 = -sigma0 * sign(psi0 - psi_t);
	
	vec xxx0 = { r0,theta0,phi0,v0,gamma0,psi0};
	mat ex1 = rk2(xxx0);

	int ex1_n_r = ex1.n_rows;
	cout << ex1(ex1_n_r - 1, 1)*Re - Re << endl;

	ex1.save("E:/zaixianguihua/pre_corre/ex1.data", raw_ascii);

	zoulang();//再入走廊
	sigma_limit();
	vec  x1 = ex1.row(ex1_n_r - 1).st();
	vv0 = x1(4);

	calcu_qeg(x1);
	vec x10 = { x1(1),x1(4),x1(5) };

	double theta1 = x1(2);
	double phi1 = x1(3);
	double s_togo = acos(cos(phi1)*cos(phif)*cos(theta1 - thetaf) + sin(phi1)*sin(phif));

	mat ex2 = rk5(s_togo, vv0);
	ex2.save("E:/zaixianguihua/pre_corre/ex2.data", raw_ascii);

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
	vec espan = linspace(e0, ef, 2000);

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
		qdot_delta = qdot(i) - q_max ;
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


//计算气动
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
	double Cl = aero(2);
	double drdv = v * sin(gamma) / (-D - sin(gamma) / r / r);
	double r0 = 1;
	double r1 = newton_r(r0, v, Cl);
	double rho1 = rho0 * exp(-(r1*Re - Re) / hs);
	double drhodr = -Re / hs * rho1;
	double K = Re * S * Cl / 2 / m;
	double drdvq = -6.3*rho1 / (drhodr*v);
	//double drdvq = -2 * v * (K* rho1 * pow(r1, 3) + pow(r1, 2)) / (K*pow(r1, 3)*drhodr*pow(v, 2) - pow(v, 2)*r1 + 2);
	//double drdvq = -2 * v * pow((K* rho1 * pow(r1, 2) + r1), 2) / (K * pow(r1, 2)*drhodr + 2 * K*rho1*r1 + 1);

	const double delta = 0.001;
	if (fabs(drdv - drdvq) < delta)
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
	double x1 = 1.01;
	const double eps = 1e-12;
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
	if (x < 15 * pi / 180)
	{
		x = 15 * pi / 180;
	}
	return x;
}

void calcu_qeg(vec x1)
{

	double theta = x1(2);
	double phi = x1(3);
	double psi = x1(6);
	double s_togo= acos(cos(phi)*cos(phif)*cos(theta - thetaf) + sin(phi)*sin(phif));
	double psi_t = 0;
	if (phif - phi > 0)
	{
		psi_t = asin(cos(phif)*sin(thetaf - theta) / sin(acos(cos(phi)*cos(phif)*cos(theta - thetaf) + sin(phi)*sin(phif))));
	}
	else
	{
		psi_t = pi - asin(cos(phif)*sin(thetaf - theta) / sin(acos(cos(phi)*cos(phif)*cos(theta - thetaf) + sin(phi)*sin(phif))));
	}
	delta_psi = psi - psi_t;

	double vf = Vf / Vc;
	uword n = 1000;
	vec sigma_m(n);
	vec v_f(n);
	sigma_m(0) = 40 * pi / 180;
	sigma_m(1) = 45 * pi / 180;
	uword i;
	double delta_v;
	for (i = 0; i < 999; i++)
	{
		sigma_mid = sigma_m(i);
		v_f(i) = rk3(s_togo, vv0);

		delta_v = v_f(i) - vf;
		if (fabs(delta_v) < 1e-10)
		{
			break;
		}

		if (i > 0)
		{
			sigma_m(i + 1) = sigma_m(i) - (sigma_m(i) - sigma_m(i - 1)) / (v_f(i) - v_f(i - 1))*delta_v;
		}
	}
}


double limit(double sigma, double v)
{
	double sigma0 = 15 * pi / 180;
	double r = 1.01;
	vec aero = calcu_aero(r, v);
	double Cl = aero(2);
	double Cd = aero(3);
	double K1 = Re * Cl * S / 2 / m;

	double sigmaqdot = acos(k_q * k_q * (1 - v * v) * pow(v, 4.3) / (K1 * q_max * q_max));
	double sigmaq = acos(Vc * Vc * (1 - v * v) / (2 * K1 * 1.5e4));
	double sigman = acos((1 - v * v) / 2.5*(sqrt(Cl * Cl + Cd * Cd) / Cl));

	double min = 0;

	if (sigmaqdot < sigmaq)
		min = sigmaqdot;
	else
		min = sigmaq;

	if (min > sigman)
		min = sigman;
	if (sigma > min)
		sigma = min;
	if (sigma < sigma0)
		sigma = sigma0;

	return sigma;
}
