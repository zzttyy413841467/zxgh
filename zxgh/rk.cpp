#include "pch.h"

mat rk1(vec espan, vec x0)//��ʼ�ι滮
{
	uword n = espan.n_elem;
	double ef = espan(n-1);
	double e0 = espan(0);
	double h = (ef - e0) / (n - 1);
	double e = 0;
	int result;
	vec x = x0;
	mat x1(n, x0.n_elem);
	vec k1(2), k2(2), k3(2), k4(2);
	uword i;
	for (i = 0; i < n; i++)
	{
		e = espan(i);
		x1.row(i) = x.st();
		result = trans(e, x);
		if (result == 1)
		{
			break;
		}
		k1 = dxde1(e, x);
		k2 = dxde1(e + h / 2, x + h * k1 / 2);
		k3 = dxde1(e + h / 2, x + h * k2 / 2);
		k4 = dxde1(e + h, x + h * k3);
		x = x + (h / 6)*(k1 + 2 * k2 + 2 * k3 + k4);
	}
	mat xx(n, 1 + x0.n_elem);
	xx = join_rows(espan, x1);
	return xx.rows(0, i-1);
}
/*
mat rk2(double e0,double ef, vec x0)
{
	uword n = 2000;
	vec espan = linspace(e0, ef, n);
	double h = (ef - e0) / (n - 1);
	double e = 0;

	int result;
	vec x = x0;
	mat x1(n, x0.n_elem);
	vec k1(6), k2(6), k3(6), k4(6);
	uword i;
	vec xxx(2);
	for (i = 0; i < n; i++)
	{
		e = espan(i);
		x1.row(i) = x.st();
		xxx(0) = x(0);
		xxx(1) = x(3);
		result = trans(e, xxx);
		if (result == 1)
		{
			break;
		}
		k1 = dxde2(e, x);
		k2 = dxde2(e + h / 2, x + h * k1 / 2);
		k3 = dxde2(e + h / 2, x + h * k2 / 2);
		k4 = dxde2(e + h, x + h * k3);
		x = x + (h / 6)*(k1 + 2 * k2 + 2 * k3 + k4);
	}
	mat xx(n, 1 + x0.n_elem);
	xx = join_rows(espan, x1);
	return xx.rows(0, i - 1);
}*/

mat rk2(vec x0)
{
	
	uword n = 10000;
	double h = 0.1;
	double hh = h / 2;
	double t0 = 0;
	int result;
	vec x = x0;
	mat x1(n, x0.n_elem);
	vec s_togo(n);
	vec k1(6), k2(6), k3(6), k4(6);
	uword i;
	vec xxx(2);
	vec x_1(6), x_2(6);
	const double eps_y = 4e-5;
	double yan=0;
	double e = 0;
	for (i = 0; i < n; i++)
	{
		x1.row(i) = x.st();
		s_togo(i) = acos(cos(x(2))*cos(phif)*cos(x(1) - thetaf) + sin(x(2))*sin(phif));;
		xxx(0) = x(0);
		xxx(1) = x(4);
		e = 1 / x(0) - x(3)*x(3) / 2;
		result = trans(e, xxx);
		if (result == 1)
		{
			break;
		}
		s1:
		k1 = dxdt1(t0, x);
		k2 = dxdt1(t0 + h / 2, x + h * k1 / 2);
		k3 = dxdt1(t0 + h / 2, x + h * k2 / 2);
		k4 = dxdt1(t0 + h, x + h * k3);
		x_1 = x + (h / 6)*(k1 + 2 * k2 + 2 * k3 + k4);

		k1 = dxdt1(t0, x);
		k2 = dxdt1(t0 + hh / 2, x + hh * k1 / 2);
		k3 = dxdt1(t0 + hh / 2, x + hh * k2 / 2);
		k4 = dxdt1(t0 + hh, x + hh * k3);
		x_2 = x + (hh / 6)*(k1 + 2 * k2 + 2 * k3 + k4);

		yan = norm_2(x_2-x_1);
		
		if (yan > eps_y)
		{
			h = h / 2;
			hh = h / 2;
			goto s1;
		}
		else if (yan < eps_y / 50)
		{
			h = h * 2;
			hh = h / 2;
			goto s1;
		}
		else
		{
			t0 = t0 + h;
			x = x_1;
		}

	}
	mat xx = join_rows(s_togo, x1);

	return xx.rows(0, i - 1);
}

double norm_2(vec x)
{
	double s = 0;
	for (uword i = 0; i < x.n_cols; i++)
	{
		s += x(i)*x(i);
	}
	return sqrt(s);
}

double rk3(double s_togo, double v)
{	
	double x = v;
	double k1, k2, k3, k4;
	uword i;
	double s = 0;
	vec sspan = linspace(s_togo, 0, 10000);
	double h = -s_togo / 9999;
	for (i = 0; i < 9999; i++)
	{
		s = sspan(i);
		k1 = dvds1(s, x);
		k2 = dvds1(s + h / 2, x + h * k1 / 2);
		k3 = dvds1(s + h / 2, x + h * k2 / 2);
		k4 = dvds1(s + h, x + h * k3);
 		x = x + (h / 6)*(k1 + 2 * k2 + 2 * k3 + k4);

	}
	return x;
}

mat rk4(double s_togo, vec x0)
{
	vec x = x0;
	vec k1(3), k2(3), k3(3), k4(3);
	uword i = 0;
	uword n = 500;
	vec sspan = linspace(s_togo, 0, n);
	double h = -s_togo / (n - 1);
	mat ex2(n, 5);
	double s = 0;
	double sigma0 = sig0_a / 180 * pi;
	double sigma;
	double vm = 5600 / Vc;
	double vf = Vf / Vc;
	double v = 0;
	for (i = 0; i < n-1; i++)
	{
		s = sspan(i);
		v = x(1);
		ex2(i, 0) = s;
		ex2.submat(i, 1, i, 3) = x.st();
		
		if (v > vm)
		{
			sigma = limit(sigma0 + (sigma_mid - sigma0) / (vm - vv0)*(v - vv0), v);
		}
		else
		{
			sigma = limit(sigmaf + (sigma_mid - sigmaf) / (vm - vf)*(v - vf), v);
		}
		ex2(i, 4) = sigma;
		k1 = dvds2(s, x,sigma);
		k2 = dvds2(s + h / 2, x + h * k1 / 2, sigma);
		k3 = dvds2(s + h / 2, x + h * k2 / 2, sigma);
		k4 = dvds2(s + h, x + h * k3, sigma);
		x = x + (h / 6)*(k1 + 2 * k2 + 2 * k3 + k4);

	}
	ex2(i, 0) = sspan(i);
	ex2.submat(i, 1, i, 3) = x.st();
	return ex2;
}

mat rk5(double s_togo, vec x1)
{
	uword num = 250;
	double v = x1(4);
	double x = v;
	double k1, k2, k3, k4;
	uword i;
	double s = 0;
	vec sspan = linspace(s_togo, 0, num);
	mat ex3(num, 4);
	double h = -s_togo / (num-1);
	double sigma0 = sig0_a / 180 * pi;
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

	for (i = 0; i < (num-1); i++)
	{
		ex3(i, 0) = sspan(i);
		ex3(i, 1) = x;
		if (x > vm)
		{
			sigma = limit(sigma0 + (sigma_mid - sigma0) / (vm - vv0)*(x - vv0), x);
		}
		else
		{
			sigma = limit(sigmaf + (sigma_mid - sigmaf) / (vm - vf)*(x - vf), x);
		}
		ex3(i, 2) = newton_r_s(x, sigma);
		s = sspan(i);
		k1 = dvds1(s, x);
		k2 = dvds1(s + h / 2, x + h * k1 / 2);
		k3 = dvds1(s + h / 2, x + h * k2 / 2);
		k4 = dvds1(s + h, x + h * k3);
		x = x + (h / 6)*(k1 + 2 * k2 + 2 * k3 + k4);

	}
	if (x > vm)
	{
		sigma = limit(sigma0 + (sigma_mid - sigma0) / (vm - vv0)*(x - vv0), x);
	}
	else
	{
		sigma = limit(sigmaf + (sigma_mid - sigmaf) / (vm - vf)*(x - vf), x);
	}
	ex3(i, 0) = sspan(i);
	ex3(i, 1) = x;
	ex3(i, 2) = newton_r_s(x, sigma);
	ex3(0, 3) = x1(5);
	for (i = 0; i < num - 1; i++)
	{
		ex3(i + 1, 3) = atan((ex3(i + 1, 2) - ex3(i, 2)) / (-h));
	}

	return ex3;
}

double newton_r_s( double v, double sigma)
{
	double f = 0;
	double x0 = 1.02;
	double fdot = 0;
	double x1 = 1.01;
	const double eps = 1e-10;
	double Cl = 0;
	vec aero = calcu_aero(x0, v);
	Cl = aero(2);
	while (fabs(x1 - x0) / fabs(x1) >= eps)
	{
		x0 = x1;
		f = 1 / x0 / x0 - v * v / x0 - 1.0 / 2 * rho0*exp(-(x0*Re - Re) / hs)*v*v*Cl*S / m * Re *cos(sigma);
		fdot = -2 / pow(x0, 3) + v * v / pow(x0, 2) + 1.0 / 2 * rho0*exp(-(x0*Re - Re) / hs)*v*v*Cl*S / m * Re * Re / hs * cos(sigma);
		x1 = x0 - f / fdot;
	}
	return x1;
}


mat rk(vec x0)
{
	uword n = 6001;
	vec tspan = linspace(0, 3, n);
	double h = 3.0 / (n - 1);
	vec k1(6), k2(6), k3(6), k4(6);
	uword i;
	mat xx(n,7);

	vec x = x0;
	double t = 0;

	mat Ref;
	Ref.load("E:/zaixianguihua/pre_corre/Ref.data");
	double ma;
	double alpha;
	double sigma;
	double sigma0 = sig0_a / 180 * pi;
	double vm = 5600 / Vc;
	double vf = Vf / Vc;
	
	vec sigma_x(n);
	sigma_x(0) = sig0;

	double s_togo;
	double dsigma;
	double dalpha;
	vec re(9);

	for (i = 0; i < n; i++)
	{
		t = tspan(i);
		xx.submat(i, 1, i, 6) = x.t();
		xx(i, 0) = t;
		if (x(3)*Vc < Vf)
		{
			break;
		}

		ma = x(3) * Vc / 340;
		if (ma >= 15)
		{
			alpha = 45 * pi / 180;
		}
		else
		{
			alpha = (45 - 0.21*(ma - 15)*(ma - 15)) * pi / 180;
		}
		if (x(3) > vv0)
		{
			sigma = sig0;
		}
		else
		{
			if (x(3) > vm)
			{
				sigma = limit(sigma0 + (sigma_mid - sigma0) / (vm - vv0)*(x(3) - vv0), x(3));
			}
			else
			{
				sigma = limit(sigmaf + (sigma_mid - sigmaf) / (vm - vf)*(x(3) - vf), x(3));
			}
		}

		
		s_togo = acos(cos(x(2))*cos(phif)*cos(x(1) - thetaf) + sin(x(2))*sin(phif));
		re = interp1(s_togo, Ref);
		dsigma = d_mag(re(3)*(x(0) - re(0)) + re(4)*(x(3) - re(1)) + re(5)*(x(4) - re(2)));
		dalpha = d_mag(re(6)*(x(0) - re(0)) + re(7)*(x(3) - re(1)) + re(8)*(x(4) - re(2)));
		//sigma = limit(sigma + dsigma, x(3));
		sigma = sigma + dsigma;
		alpha = alpha + dalpha;
		

		
		if (i > 0)
		{
			sigma_x(i) = sigma * sign_decide(sign(sigma_x(i - 1)), x);
		}


		k1 = dxdt2(t, x, sigma_x(i), alpha);
		k2 = dxdt2(t + h / 2, x + h * k1 / 2, sigma_x(i), alpha);
		k3 = dxdt2(t + h / 2, x + h * k2 / 2, sigma_x(i), alpha);
		k4 = dxdt2(t + h, x + h * k3, sigma_x(i), alpha);
		x= x + (h / 6)*(k1 + 2 * k2 + 2 * k3 + k4);


	}

	vec sigma1 = sigma_x.rows(0, i);
	sigma1.save("E:/zaixianguihua/pre_corre/sigma.data", raw_ascii);
	
	return xx.rows(0, i);
}



double sign_decide(double sign0, vec x)
{

	double delt_sigma1 = 10.0 / 180 * pi;
	double delt_sigma2 = 20.0 / 180 * pi;
	double delt_sigma3 = 5.0 / 180 * pi;
	double r = x(0);
	double theta = x(1);
	double phi = x(2);
	double v = x(3);
	double psi = x(5);
	double psi_t;
	if (phif - phi > 0)
		psi_t = asin(cos(phif)*sin(thetaf - theta) / sin(acos(cos(phi)*cos(phif)*cos(theta - thetaf) + sin(phi)*sin(phif))));
	else
		psi_t = pi - asin(cos(phif)*sin(thetaf - theta) / sin(acos(cos(phi)*cos(phif)*cos(theta - thetaf) + sin(phi)*sin(phif))));

	double delt_sigma = psi - psi_t;

	double sign;
	double V = v * Vc;
	double delt;
	if (V >= 5000)
	{
		if (delt_sigma > delt_sigma1)
			sign = -1;
		else if (delt_sigma <= delt_sigma1 && delt_sigma >= -delt_sigma1)
			sign = sign0;
		else
			sign = 1;
	}
	else if (4000 <= V && V < 5000)
	{
		if (delt_sigma > delt_sigma2)
			sign = -1;
		else if (delt_sigma <= delt_sigma2 && delt_sigma >= -delt_sigma2)
			sign = sign0;
		else
			sign = 1;
	}
	else if (3000 <= V && V < 4000)
	{
		delt = delt_sigma3 + (delt_sigma2 - delt_sigma3) / (4000 - 3000)*(V - 3000);
		if (delt_sigma > delt)
			sign = -1;
		else if (delt_sigma <= delt && delt_sigma >= -delt)
			sign = sign0;
		else
			sign = 1;
	}
	else
	{
		if (delt_sigma > delt_sigma3)
			sign = -1;
		else if (delt_sigma <= delt_sigma3 && delt_sigma >= -delt_sigma3)
			sign = sign0;
		else
			sign = 1;
	}

	return sign;
 }