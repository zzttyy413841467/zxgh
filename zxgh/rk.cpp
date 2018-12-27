#include "pch.h"

mat rk1(vec espan, vec x0)//≥ı º∂ŒπÊªÆ
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
	vec tspan(n);
	vec k1(6), k2(6), k3(6), k4(6);
	uword i;
	vec xxx(2);
	vec x_1(6), x_2(6);
	const double eps_y = 1e-5;
	double yan=0;
	double e = 0;
	for (i = 0; i < n; i++)
	{
		x1.row(i) = x.st();
		tspan(i) = t0;
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
	mat xx(n, 1 + x0.n_elem);
	xx = join_rows(tspan, x1);
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