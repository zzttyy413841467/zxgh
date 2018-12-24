#include "pch.h"

mat rk1(vec espan, vec x0)
{
	uword n = espan.n_elem;
	double ef = espan(n-1);
	double e0 = espan(0);
	double h = (ef - e0) / (n - 1);
	double e = 0;
	int result;
	vec x = x0;
	mat x1(n, x0.n_elem);
	vec k1(2), k2(2), k3(2), k4(5);
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
		k1 = f1(e, x);
		k2 = f1(e + h / 2, x + h * k1 / 2);
		k3 = f1(e + h / 2, x + h * k2 / 2);
		k4 = f1(e + h, x + h * k3);
		x = x + (h / 6)*(k1 + 2 * k2 + 2 * k3 + k4);
	}
	mat xx(n, 1 + x0.n_elem);
	xx = join_rows(espan, x1);
	return xx.rows(0, i-1);
}