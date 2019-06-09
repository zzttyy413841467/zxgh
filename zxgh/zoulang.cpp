#include "pch.h"

//计算hv走廊
void onLineTrack::zoulang()
{
	getTime();
	cout << "计算再入走廊" << endl;
	vec V = linspace(850, 7800, 140);
	vec aero(4);
	mat h_bianjie(140, 4);
	for (uword i = 0; i < 140; i++)
	{
		aero = calcu_aero(1, V(i) / Vc);
		h_bianjie(i, 0) = 2 * hs*log(sqrt(1.225)*pow(V(i), 3.15)*9.4369e-5 / (Qdot_max+50000));
		h_bianjie(i, 1) = hs * log(1.225*pow(V(i), 2) / 2 / q_max);
	
		h_bianjie(i, 2) = hs * log(1.225*S*V(i)*V(i)*sqrt(aero(2)*aero(2) + aero(3)*aero(3)) / 2 / n_max / g0 / m);	//法向过载
		h_bianjie(i, 3) = newton_h(1, V(i) / Vc, aero(2));
	}
	h_bianjie.save("E:/zaixianguihua/pre_corre/h_bianjie.data", raw_ascii);
	
}

//计算倾侧角边界
void onLineTrack::sigma_limit()
{
	getTime();
	cout << "计算倾侧角走廊" << endl;
	double K1 = 0;
	vec v = linspace(850/Vc, 7800/Vc, 140);
	vec aero(4);
	mat sigma_bianjie(140, 4);
	for (uword i = 0; i < 140; i++)
	{
		aero = calcu_aero(1.01, v(i));
		K1 = Re * aero(2) * S / 2 / m;
		sigma_bianjie(i, 0) = acos(k_q*k_q*(1 - v(i) * v(i))*pow(v(i), 4.3) / (K1*(Qdot_max + 0)*(Qdot_max + 0))) * 180 / pi;
		sigma_bianjie(i, 1) = acos(Vc*Vc*(1 - v(i) * v(i)) / (2 * K1*q_max)) * 180 / pi;
		
		sigma_bianjie(i, 2) = acos((1 - v(i) * v(i)) / n_max*(sqrt(aero(2) * aero(2) + aero(3) * aero(3)) / aero(2))) * 180 / pi;//法向过载
		sigma_bianjie(i, 3) = 15.0;
	}
	sigma_bianjie.save("E:/zaixianguihua/pre_corre/sigma_bianjie.data", raw_ascii);
}

//牛顿法(走廊())
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
	return x1 * Re - Re;
}