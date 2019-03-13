#include "pch.h"

//初始段规划微分方程
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

//初始段仿真微分方程
vec dxdt1(double t, vec x)
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

//滑翔段规划(单)微分方程
double dvds1(double s, double v)
{
	
	double sigma0=sig0_a/180*pi;
	double sigma;
	double vm = 4200 / Vc;
	double vf = vvf;
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
	double vdot = (1 / r - v * v)*(Cd / Cl) / (v*cos(1 * delta_psi)*cos(sigma));
	return vdot;
}

//滑翔段规划(多)微分方程
vec dvds2(double s, vec x, double sigma)
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


vec dxde2(double e, vec x, double sigma)
{
	double r = x(0);
	double gamma = x(1);
	double s_togo = x(2);

	double v = sqrt(2 * (1 / r - e));

	vec aero = calcu_aero(r, v);
	double Cl = aero(2);
	double Cd = aero(3);
	double D = aero(1);
	double L = aero(0);

	
	vec xdot(3);
	xdot(0) = sin(gamma) / D;
	xdot(1) = 1 / (D*v*v)*(L*cos(sigma) + (v*v / r - 1 / r / r)*cos(gamma));
	xdot(2) = -cos(gamma) / D / r;
	return xdot;
}


//运动微分方程
vec dxdt2(double t, vec x,double sigma,double alpha)
{
	
	double r = x(0);
	double theta = x(1);
	double phi = x(2);
	double V = x(3);
	double gamma = x(4);
	double psi = x(5);

	double Cl = cl0 + cl1 * alpha + cl2 * alpha*alpha;
	double Cd = cd0 + cd1 * Cl + cd2 * Cl*Cl;
	double h = r * Re - Re;
	double rho = rho0 * exp(-h / hs);

	double q = 1.0 / 2 * rho * pow(V, 2);
	
	double L = q * Cl * S / m * Re;
	double D = q * Cd * S / m * Re;

	vec xdot(6);
	xdot(0) = V * sin(gamma);
	xdot(1) = V * cos(gamma)*sin(psi) / (r*cos(phi));
	xdot(2) = V * cos(gamma)*cos(psi) / r;
	xdot(3) = -D - (sin(gamma) / r / r);
	xdot(4) = 1 / V * (L * cos(sigma) + (V * V / r - 1 / r / r)*cos(gamma));
	xdot(5) = (L*sin(sigma) / V / cos(gamma) + V / r * cos(gamma)*sin(psi)*tan(phi));
	return xdot;
}

//分段二次插值
vec interp2(double s_togo, mat Ref)
{
	uword i = abs(Ref.col(0) - s_togo).index_min();
	uword n_r = Ref.n_rows;
	uword n_c = Ref.n_cols;
	vec out(n_c - 1);
	if (i == 0)
	{
		i = 1;
	}
	if (i == n_r - 1)
	{
		i = n_r - 2;
	}
	for (uword j = 1; j < n_c; j++)
	{
		out(j - 1) = (s_togo - Ref(i, 0))*(s_togo - Ref(i + 1, 0)) / (Ref(i - 1, 0) - Ref(i, 0)) / (Ref(i - 1, 0) - Ref(i + 1, 0))*Ref(i - 1, j) +
			(s_togo - Ref(i - 1, 0))*(s_togo - Ref(i + 1, 0)) / (Ref(i, 0) - Ref(i - 1, 0)) / (Ref(i, 0) - Ref(i + 1, 0))*Ref(i, j) +
			(s_togo - Ref(i - 1, 0))*(s_togo - Ref(i, 0)) / (Ref(i + 1, 0) - Ref(i - 1, 0)) / (Ref(i + 1, 0) - Ref(i, 0))*Ref(i + 1, j);
	}

	return out;

}


//控制量变化幅度
double d_mag(double theta)
{
	if (theta>5 * pi / 180)
	{
		theta = 5 * pi / 180;
	}
	if (theta < -5 * pi / 180)
	{
		theta = -5 * pi / 180;
	}
	return theta;
}
