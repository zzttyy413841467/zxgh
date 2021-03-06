﻿#include "pch.h"

//printf("File    Fame: %s\n", __FILE__);      //文件名
//printf("Present Line: %d\n", __LINE__);      //所在行
//printf("Present Function: %s\n", __func__);  //函数名

//全局变量
double sig0_a = 0;
vector<Noflyzone> nfz;
double sig0 = 0;
double delta_psi;
double sigma_mid;
double vv0;
int item = 0;
vec sigma_qeg(501);
double sigmaf;
double vvf;
double Qmax;
double qmax;
double nmax;
mat waypoint=zeros<mat>(10,2);

//主函数
int main(int argc, char * argv[])
{
	Noflyzone nfz1(29.0 / 180.0 * pi, 9 / 180.0 * pi, 3.4 / 180.0 * pi, 2);
	Noflyzone nfz2(38.0 / 180.0 * pi, 18 / 180.0 * pi, 2.5 / 180.0 * pi, 1);
	Noflyzone nfz3(43.0 / 180.0 * pi, 13.5 / 180.0 * pi, 3 / 180.0 * pi, 2);
	//nfz.push_back(nfz1);
	nfz.push_back(nfz2);
	//nfz.push_back(nfz3);

	//Noflyzone nfz1(36.0 / 180.0 * pi, 2.5 / 180.0 * pi, 4 / 180.0 * pi, 1);
	//Noflyzone nfz2(48.0 / 180.0 * pi, -4.5 / 180.0 * pi, 5 / 180.0 * pi, 2);
	//Noflyzone nfz3(47.0 / 180.0 * pi, 20.5 / 180.0 * pi, 3 / 180.0 * pi, 1);
	//nfz.push_back(nfz1);
	//nfz.push_back(nfz2);
	////nfz.push_back(nfz3);

	for (int i = 0; i < nfz.size() - 1; i++)
	{
		waypoint.row(i) = newton_waypoint(nfz[i], nfz[i + 1]).t();
	}
	//把初始值封装成类，计算弹道作为类的方法
	double H0 = 121900;
	double V0 = 7630;
	double Hf = 26500;
	double Vf = 900;
	double r0 = (Re + H0) / Re;
	double v0 = V0 / Vc;
	double theta0 = 0;
	double phi0 = 0;
	double gamma0 = -0.5 * pi / 180;
	//double psi0 = 65 * pi / 180;
	double psi0 = 65 * pi / 180;
	double tau0 = 0;
	double rf = (Re + Hf) / Re;
	double vf = Vf / Vc;

	vvf = vf;

	double Qdot_max = 790000;
	double q_max = 1.5e4;
	double n_max = 2.5;
	double theta_f = 58 * pi / 180;
	double phi_f = 20 * pi / 180;
	vector<double> state0 = { r0,theta0,phi0,v0,gamma0,psi0 };
	vector<double> statef = { rf,theta_f,phi_f,vf };

	
	onLineTrack online_example(state0, statef, Qdot_max, q_max, n_max);

	online_example.zoulang();//再入走廊
	online_example.sigma_limit();

	online_example.calcu_online_track();

	mat xx = rk(state0, statef);
	getTime();
	cout << "滑翔段结束" << endl;
	xx.save("E:/zaixianguihua/pre_corre/xx.data", raw_ascii);
	

	//offLineTrack offline_example(state0, statef, Qdot_max, q_max, n_max);
	//offline_example.calcu_offline_track();


	return 0;
}

void onLineTrack::calcu_online_track()
{
	getTime();
	cout << "开始计算弹道" << endl;
	double thetaf = stateVaraible_f[1];
	double phif = stateVaraible_f[2];
	double rf = stateVaraible_f[0];
	double vf = stateVaraible_f[3];
	double theta0 = stateVaraible_0[1];
	double phi0 = stateVaraible_0[2];

	Qmax = Qdot_max;
	qmax = q_max;
	nmax = n_max;

	vec rv = calcu_aero(rf, vf);
	sigmaf = acos((1 / rf / rf - vf * vf / rf) / rv(0));
	double psi_t = 0;
	if (phif - phi0 > 0)
	{
		psi_t = asin(cos(phif)*sin(thetaf - theta0) / sin(acos(cos(phi0)*cos(phif)*cos(theta0 - thetaf) + sin(phi0)*sin(phif))));
	}
	else
	{
		psi_t = pi - asin(cos(phif)*sin(thetaf - theta0) / sin(acos(cos(phi0)*cos(phif)*cos(theta0 - thetaf) + sin(phi0)*sin(phif))));
	}

	calcu_sigma0();

	sig0 = -sigma0 * sign(stateVaraible_0[5] - psi_t);

	//vec xxx0 = stateVaraible_0;
	mat ex1 = rk2(stateVaraible_0, stateVaraible_f);

	uword ex1_n_r = ex1.n_rows;

	ex1.save("E:/zaixianguihua/pre_corre/ex1.data", raw_ascii);

	vec  x1 = ex1.row(ex1_n_r - 1).st();
	vv0 = x1(4);

	double sigma0_qeg;
	rv = calcu_aero(x1(1), x1(4));
	sigma0_qeg = acos((1 / x1(1) / x1(1) - x1(4) * x1(4) / x1(1)) / rv(0));

	if (sig0_a < 180 / pi * sigma0_qeg)
	{
		sig0_a = 180 / pi * sigma0_qeg;
	}

	calcu_qeg(x1);

	vec x10 = { x1(1),x1(4),x1(5) };

	double theta1 = x1(2);
	double phi1 = x1(3);
	double s_togo = acos(cos(phi1)*cos(phif)*cos(theta1 - thetaf) + sin(phi1)*sin(phif));


	mat ex2 = rk5(s_togo, x1);                                  //弹道倾角不连续

	//double e00 = 1 / x1(1) - x1(4)*x1(4) / 2;
	//vec x100 = { x1(1),x1(5),s_togo };
	//mat ex2 = rk6(e00, x100);

	ex2.save("E:/zaixianguihua/pre_corre/ex2.data", raw_ascii);

	//calcu_waypoint(xxx0, nfz1, nfz2);
	calcu_K(ex1, ex2);

}

void onLineTrack::calcu_sigma0()
{
	getTime();
	cout << "计算初始段倾侧角" << endl;
	double r0 = stateVaraible_0[0];
	double v0 = stateVaraible_0[3];
	double rf = stateVaraible_f[0];
	double vf = stateVaraible_f[3];
	double e0 = 1 / r0 - v0 * v0 / 2;
	double ef = 1 / rf - vf * vf / 2;
	vec espan = linspace(e0, ef, 2000);

	double gamma0 = stateVaraible_0[4];//初始弹道倾角
	vec x0 = { r0,gamma0 };
	vec sigma = zeros(2000, 1);
	vec qdot = zeros(2000, 1);
	sigma(0) = 30 * pi / 180;
	sigma(1) = 35 * pi / 180;
	uword i;
	mat ex;
	uword num_r = 0;
	double rho = 0;
	double v = 0;
	double qdot_delta = 0;
	for (i = 0; i < 1999; i++)
	{
		sig0_a = sigma(i);
		ex = rk1(espan, x0);
		num_r = ex.n_rows;
		ef = ex(num_r - 1, 0);
		rf = ex(num_r - 1, 1);
		rho = rho0 * exp(-(rf*Re - Re) / hs);
		v = sqrt(2 * (1 / rf - ef));
		qdot(i) = k_q * sqrt(rho)*pow(v, 3.15);
		qdot_delta = qdot(i) - Qdot_max + 100000;            //
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

////初始段仿真结束条件
//int trans(double e, vec x)
//{
//	int result = 0;
//	double r = x(0);
//	double gamma = x(1);
//	double v = sqrt(2 * (1 / r - e));
//	vec aero=calcu_aero(r, v);
//	double L = aero(0);
//	double D = aero(1);
//	double Cl = aero(2);
//	double drdv = v * sin(gamma) / (-D - sin(gamma) / r / r);
//	double drdt = v * sin(gamma);
//	double r0 = 1;
//	double r1 = newton_r(r0, v, Cl);
//	double rho1 = rho0 * exp(-(r1 * Re - Re) / hs);
//
//	double drhodr = -Re / hs * rho1;
//
//	double K = Re * S * Cl / 2 / m;
//
//	double drdvqeg = (-K * rho1 * 2 * v - 2 * v / r1) / (K * v * v * drhodr - v * v / r1 / r1 + 2 / pow(r1, 3));
//	//三个公式计算的结果一样
//	/*double drdvqeg = -2 * v * (K* rho1 * pow(r1, 3) + pow(r1, 2)) / (K*pow(r1, 3)*drhodr*pow(v, 2) - pow(v, 2)*r1 + 2);
//	double drdvqeg = -2 * v * pow((K* rho1 * pow(r1, 2) + r1), 2) / (K * pow(r1, 2)*drhodr + 2 * K*rho1*r1 + 1);*/
//
//	double drdvq = 6.3*hs / Re / v;
//
//	const double delta = 1e-2;
//	//if (fabs(drdv - drdvqeg)/fabs(drdv) < delta)
//	if (fabs(drdt) < 5e-3)
//	{
//		result = 1;
//	}
//	else
//	{
//		result = 0;
//	}
//	return result;
//}

//备份trans()
int trans(double e, vec x)
{
	int result = 0;
	double r = x(0);
	double gamma = x(1);
	double v = sqrt(2 * (1 / r - e));
	vec aero = calcu_aero(r, v);
	double L = aero(0);
	double D = aero(1);
	double Cl = aero(2);
	double drdv = v * sin(gamma) / (-D - sin(gamma) / r / r);
	double r0 = 1;
	double r1 = newton_r(r0, v, Cl);
	double rho1 = rho0 * exp(-(r1 * Re - Re) / hs);

	double drhodr = -Re / hs * rho1;

	double K = Re * S * Cl / 2 / m;

	double drdvq = -6.3*rho1 / (drhodr*v);
	//double drdvq = -2 * v * (K* rho1 * pow(r1, 3) + pow(r1, 2)) / (K*pow(r1, 3)*drhodr*pow(v, 2) - pow(v, 2)*r1 + 2);
	//double drdvq = -2 * v * pow((K* rho1 * pow(r1, 2) + r1), 2) / (K * pow(r1, 2)*drhodr + 2 * K*rho1*r1 + 1);

	/*double sigma0_qeg = acos((1 / r / r - v * v / r) / L);
	double drdvqeg=(-K * rho1 * 2 * v *cos(sigma0_qeg) - 2 * v / r) / (K * cos(sigma0_qeg) * v * v * drhodr - v * v / r / r + 2 / pow(r, 3));*/

	double drdvqeg = (-K * rho1 * 2 * v - 2 * v / r1) / (K * v * v * drhodr - v * v / r1 / r1 + 2 / pow(r1, 3));

	const double delta = 0.001;
	if (fabs(drdv - drdvqeg) < delta)
	{
		result = 1;
	}
	else
	{
		result = 0;
	}
	return result;
}



//牛顿法(trans()函数中)
double newton_r(double x0, double v, double Cl)
{
	double f = 0;
	double fdot = 0;
	double x1 = 1.01;
	const double eps = 1e-12;
	double sigma0 = sig0_a * pi / 180;
	while (fabs(x1 - x0) / fabs(x1) >= eps)
	{
		x0 = x1;
		f = 1 / x0 / x0 - v * v / x0 - 1.0 / 2 * rho0*exp(-(x0*Re - Re) / hs)*v*v*Cl*S / m * Re ;
		fdot = -2 / pow(x0, 3) + v * v / pow(x0, 2) + 1.0 / 2 * rho0*exp(-(x0*Re - Re) / hs)*v*v*Cl*S / m * Re*Re / hs ;
		x1 = x0 - f / fdot;
	}
	return x1;
}

//限幅
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

//用割线法在终端速度满足约束的情况下求取滑翔段的倾侧角
void onLineTrack::calcu_qeg(vec x1)     
{
	getTime();
	cout << "计算滑翔段倾侧角" << endl;
	double thetaf = stateVaraible_f[1];
	double phif = stateVaraible_f[2];
	double theta = x1(2);
	double phi = x1(3);
	double psi = x1(6);
	double s_togo = acos(cos(phi)*cos(phif)*cos(theta - thetaf) + sin(phi)*sin(phif));
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

	double vf = stateVaraible_f[3];
	uword n = 1000;
	vec sigma_m(n);
	vec v_f(n);
	sigma_m(0) = 45 * pi / 180;
	sigma_m(1) = 60 * pi / 180;
	uword i;
	double delta_v;
	double mid;
	for (i = 0; i < 999; i++)
	{
		sigma_mid = sigma_m(i);
		v_f(i) = rk3(s_togo, vv0);

		delta_v = v_f(i) - vf;
		if (fabs(delta_v) < 1e-12)
		{
			break;
		}

		if (i > 0)
		{
			mid = sigma_m(i) - (sigma_m(i) - sigma_m(i - 1)) / (v_f(i) - v_f(i - 1))*delta_v;
			if (mid > 85.0 / 180 * pi)
				mid = 85.0 / 180 * pi;
			sigma_m(i + 1) = mid;
		}
	}
}

////用割线法在终端速度满足约束的情况下求取滑翔段的倾侧角（自变量为e）
//void calcu_qeg1(vec x1)     
//{
//
//	double r0 = x1(1);
//	double v0 = x1(4);
//	double gamma0 = x1(5);
//	double theta = x1(2);
//	double phi = x1(3);
//	double psi = x1(6);
//	double s_togo = acos(cos(phi)*cos(phif)*cos(theta - thetaf) + sin(phi)*sin(phif));
//	double e0 = 1 / r0 - v0 * v0 / 2;
//	vec x10 = { x1(1),x1(5),s_togo };
//	uword n = 1000;
//	vec sigma_m(n);
//	vec s_f(n);
//	sigma_m(0) = 45 * pi / 180;
//	sigma_m(1) = 60 * pi / 180;
//	uword i;
//
//	for (i = 0; i < n - 1; ++i) 
//	{
//		sigma_mid = sigma_m(i);
//
//		s_f(i) = rk4(e0, x10);
//
//		if (fabs(s_f(i)) < 1e-12)
//		{
//		break;
//		}
//		if (i > 0)
//		{
//			sigma_m(i + 1) = sigma_m(i) - (sigma_m(i) - sigma_m(i - 1)) / (s_f(i) - s_f(i - 1))*s_f(i);
//		}
//
//	}
//
//
//}

//倾侧角走廊约束
double limit(double sigma, double v)
{

	double sigma0 = 15 * pi / 180;
	double r = 1.01;
	vec aero = calcu_aero(r, v);
	double Cl = aero(2);
	double Cd = aero(3);
	double K1 = Re * Cl * S / 2 / m;

	double sigmaqdot = acos(k_q * k_q * (1 - v * v) * pow(v, 4.3) / (K1 * Qmax * Qmax));
	double sigmaq = acos(Vc * Vc * (1 - v * v) / (2 * K1 * qmax));
	double sigman = acos((1 - v * v) / nmax*(sqrt(Cl * Cl + Cd * Cd) / Cl));

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

//LQR计算反馈增益
void calcu_K(mat ex1, mat ex2)
{
	getTime();
	cout << "计算反馈增益矩阵" << endl;
	mat ex2_1 = ex2.rows(1, ex2.n_rows - 1);
	vec s_togo = join_cols(ex1.col(0), ex2_1.col(0));
	vec r = join_cols(ex1.col(1), ex2_1.col(2));
	vec v = join_cols(ex1.col(4), ex2_1.col(1));
	vec gamma = join_cols(ex1.col(5), ex2_1.col(3));
	uword n = r.n_elem;
	uword i = 0;
	vec sigma(n);
	//
	double sigma0 = sig0_a / 180 * pi;
	double vm = 4200 / Vc;
	double vf = vvf;
	for (i = 0; i < n; i++)
	{
		if (v(i) > vv0)
		{
			sigma(i) = fabs(sig0);
		}
		else
		{
			if (v(i) > vm)
			{
				sigma(i) = limit(sigma0 + (sigma_mid - sigma0) / (vm - vv0)*(v(i) - vv0), v(i));
			}
			else
			{
				sigma(i) = limit(sigmaf + (sigma_mid - sigmaf) / (vm - vf)*(v(i) - vf), v(i));
			}
		}
	}
	
	sigma.save("E:/zaixianguihua/pre_corre/sigma.data", raw_ascii);

	vec ma = v * Vc / 340;
	vec alpha(n);
	for (i = 0; i < ma.n_elem; i++)
	{
		if (ma(i) >= 15)
		{
			alpha(i) = 45 * pi / 180;
		}
		else
		{
			alpha(i) = (45 - 0.21*(ma(i) - 15)*(ma(i) - 15)) * pi / 180;
		}
	}
	vec Cl = cl0 + cl1 * alpha + cl2 * alpha % alpha;
	vec Cd = cd0 + cd1 * Cl + cd2 * Cl % Cl;
	vec h = r * Re - Re;
	vec rho = rho0 * exp(-h / hs);
	vec q = 1.0 / 2 * rho % pow(v, 2);
	vec L = q % Cl * S / m * Re;
	vec D = q % Cd * S / m * Re;
	
	//vec Cl_v;
	//vec Cd_v;
	vec Cl_a = cl1 + 2 * cl2*alpha;
	vec Cd_a = (cd1 + 2 * cd2*Cl) % Cl_a;

	vec L_r = -Re / hs * L;
	//vec L_v = 2 * L / v + L % Cl_v / Cl;
	vec L_v = 2 * L / v;
	vec L_a = L % Cl_a / Cl;
	vec D_r = -Re / hs * D;
	//vec D_v = 2 * D / v + D % Cd_v / Cd;
	vec D_v = 2 * D / v;
	vec D_a = D % Cd_a / Cd;

	
	cube A(3, 3, n);
	A.subcube(0, 0, 0, 0, 0, n - 1) = -tan(gamma);
	A.subcube(0, 1, 0, 0, 1, n - 1) = zeros(n,1);
	A.subcube(0, 2, 0, 0, 2, n - 1) = -r / pow(cos(gamma), 2);
	A.subcube(1, 0, 0, 1, 0, n - 1) = (D_r%r + D - sin(gamma) / pow(r, 2)) / (v%cos(gamma));
	A.subcube(1, 1, 0, 1, 1, n - 1) = (D_v%r%v - D % r - sin(gamma) / r) / (v%v%cos(gamma));
	A.subcube(1, 2, 0, 1, 2, n - 1) = (1 / r + D % r%sin(gamma)) / (v%pow(cos(gamma), 2));
	A.subcube(2, 0, 0, 2, 0, n - 1) = -((L_r%r + L) % cos(sigma) / cos(gamma) + 1 / pow(r, 2)) / pow(v, 2);
	A.subcube(2, 1, 0, 2, 1, n - 1) = ((2 * L - L_v % v) % r%cos(sigma) - 2 * cos(gamma) / r) / (pow(v, 3) % cos(gamma));
	A.subcube(2, 2, 0, 2, 2, n - 1) = -L % r%cos(sigma) % sin(gamma) / (pow(v, 2) % pow(cos(gamma), 2));

	cube B(3, 2, n);
	B.subcube(0, 0, 0, 0, 0, n - 1) = zeros(n, 1);
	B.subcube(0, 1, 0, 0, 1, n - 1) = zeros(n, 1);
	B.subcube(1, 0, 0, 1, 0, n - 1) = zeros(n, 1);
	B.subcube(1, 1, 0, 1, 1, n - 1) = D_a % r / (v%cos(gamma));
	B.subcube(2, 0, 0, 2, 0, n - 1) = L % r%sin(sigma) / (pow(v, 2) % cos(gamma));
	B.subcube(2, 1, 0, 2, 1, n - 1) = -L_a % r%cos(sigma) / (pow(v, 2) % cos(gamma));
	cube K(2, 3, n);

#pragma omp parallel for
	for (int k = 0; k < n; k++)
	{
		K.slice(k)= lqr_iter(A.slice(k), B.slice(k));
	}
	
	//K.save("E:/zaixianguihua/pre_corre/K.data", raw_ascii);
	mat Ref(n, 10);
	Ref.col(0) = s_togo;
	Ref.col(1) = r;
	Ref.col(2) = v;
	Ref.col(3) = gamma;
	vec x = K.subcube(0, 0, 0, 0, 0, n - 1);
	Ref.col(4) = x;
	x = K.subcube(0, 1, 0, 0, 1, n - 1);
	Ref.col(5) = x;
	x = K.subcube(0, 2, 0, 0, 2, n - 1);
	Ref.col(6) = x;
	x = K.subcube(1, 0, 0, 1, 0, n - 1);
	Ref.col(7) = x;
	x = K.subcube(1, 1, 0, 1, 1, n - 1);
	Ref.col(8) = x;
	x = K.subcube(1, 2, 0, 1, 2, n - 1);
	Ref.col(9) = x;
	Ref.save("E:/zaixianguihua/pre_corre/Ref.data",raw_ascii);
}

//lqr迭代子程序
mat lqr_iter(mat A, mat B)
{
	double dr_max = 9000 / Re;
	double dv_max = 1000 / Vc;
	mat Q = { {1,0,0},{0,pow(dr_max,2) / pow(dv_max,2),0},{0,0,pow(dr_max,2) / 0.042} };
	mat R = { {pow(dr_max,2) / 0.01,0},{0,pow(dr_max,2) / 0.01} };
	mat P(3,3);
	mat I = eye(3, 3);
	mat C = inv(I - A);
	mat D = C.t();
	mat E = C * (I + A);
	mat G = 2 * C *C* B;
	mat H = R + B.t() * D * Q * C * B;
	mat W = Q * C * B;
	mat P0 = I;
	mat P1 = zeros(3, 3);
	mat M(3,2);
	while (norm(P1 - P0, 1) > 1e-4)
	{
		P0 = P1;
		M = E.t()*P1*G + W;
		P1 = E.t()*P1*E - M * inv(G.t()*P1*G + H)*M.t() + Q;
	}
	P = 2 * D * P1 * C;
	mat K = inv(R)*B.t()*P;
	return K;
}

void getTime()
{
	time_t t = time(0);
	tm a;
	localtime_s(&a, &t);
	printf("%02d:%02d:%02d\t", a.tm_hour, a.tm_min, a.tm_sec);

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//离线部分


void offLineTrack::calcu_offline_track ()
{
	getTime();
	cout << "计算离线轨迹" << endl;
	
	double r0 = stateVaraible_0[0];
	double v0 = stateVaraible_0[3];
	double theta0 = stateVaraible_0[1];
	double phi0 = stateVaraible_0[2];
	double rf = stateVaraible_f[0];
	double vf = stateVaraible_f[3];
	double thetaf = stateVaraible_f[1];
	double phif = stateVaraible_f[2];

	double s_togo = acos(cos(phi0)*cos(phif)*cos(theta0 - thetaf) + sin(phi0)*sin(phif));

	double e0 = 1 / r0 - v0 * v0 / 2;
	double ef = 1 / rf - vf * vf / 2;
	vec espan = linspace(e0, ef, 5000);

	double gamma0 = stateVaraible_0[4];//初始弹道倾角
	vec x0 = { r0,gamma0,s_togo};

	vec sigma = zeros(500, 1);
	sigma(0) = -20 * pi / 180;
	sigma(1) = 10 * pi / 180;

	uword i = 0;

	double psi_t = 0;
	if (phif - phi0 > 0)
	{
		psi_t = asin(cos(phif)*sin(thetaf - theta0) / sin(s_togo));
	}
	else
	{
		psi_t = pi - asin(cos(phif)*sin(thetaf - theta0) / sin(s_togo));
	}

	vec s_togo_rest(500);//当前弹道的航程与需要航程的差
	double sigma_init;
	
	for (i = 0; i < 500; i++)
	{
		sig0_a = sigma(i);

		//计算初始段弹道及航程

		//此处rk1()可以重载，输出为最后一组纵向参数
		mat ex = rk_off_init(espan, x0);

		ex.save("E:/zaixianguihua/pre_corre/off_ex.data", raw_ascii);

		uword num_r = ex.n_rows;
		double et = ex(num_r - 1, 0);
		double rt = ex(num_r - 1, 1);
		double vt = sqrt(2 * (1 / rt - et));
		double gammat = ex(num_r - 1, 2);
		double s_togot = ex(num_r - 1, 3);

		//解析规划滑翔段弹道
		vec aerot = calcu_aero(rt, vt);
		double Dt = aerot(1);
		vec x_t = { rt,vt, gammat };
		/*x_t.print();*/
		//Hermite插值
		double vm = (vt + vf) / 2;
		vec v_h = { vt,vm,vf };       //插值节点
		double rm = calcu_mid_r(vm);
		vec r_h = { rt,rm,rf };
		double drdvt = vt * sin(gammat) / (-Dt - sin(gammat) / rt / rt);

		vec aerof = calcu_aero(rf, vf);
		double Cl = aerof(2);
		double r1 = newton_r(1.012, vf, Cl);
		double rho1 = rho0 * exp(-(r1 * Re - Re) / hs);
		double drhodr = -Re / hs * rho1;
		double K = Re * S * Cl / 2 / m;
		double drdvqegf = (-K * rho1 * 2 * vf - 2 * vf / r1) / (K * vf * vf * drhodr - vf * vf / r1 / r1 + 2 / pow(r1, 3));
		double drdvqf = 2 * hs / Re / vf;
		double drdvf = 0.026;

		vec rdot_h = { drdvt,drdvf };
		mat A = { {pow(vt,4),pow(vt,3),pow(vt,2),vt,1},{pow(vm,4),pow(vm,3),pow(vm,2),vm,1},{pow(vf,4),pow(vf,3),pow(vf,2),vf,1},{4 * pow(vt,3),3 * pow(vt, 2),2 * vt,1,0},{4 * pow(vf,3),3 * pow(vf, 2),2 * vf,1,0} };
		vec b = { rt,rm,rf,drdvt,drdvf };
		vec xishu = solve(A, b);
		xishu.save("E:/zaixianguihua/pre_corre/xishu.data", raw_ascii);

		//计算滑翔段弹道参数及预测航程
		uword n_qeg = 501;
		vec v_qeg = linspace(vt, vf, n_qeg);
		vec r_qeg = xishu(0)*pow(v_qeg, 4) + xishu(1)*pow(v_qeg, 3) + xishu(2)*pow(v_qeg, 2) + xishu(3)*v_qeg + xishu(4);

		vec drdv_qeg = 4 * xishu(0)*pow(v_qeg, 3) + 3 * xishu(1)*pow(v_qeg, 2) + 2 * xishu(2)*v_qeg + xishu(3);
		vec D_qeg(n_qeg);
		for (uword j = 0; j < n_qeg; j++)
		{
			vec aero_q = calcu_aero(r_qeg(j), v_qeg(j));
			D_qeg(j) = aero_q(1);
		}
		vec gamma_qeg = asin(-drdv_qeg % D_qeg / (v_qeg + drdv_qeg / pow(r_qeg, 2)));
		mat ex_qeg(n_qeg, 3);
		ex_qeg.col(0) = r_qeg;
		ex_qeg.col(1) = v_qeg;
		ex_qeg.col(2) = gamma_qeg;

		ex_qeg.save("E:/zaixianguihua/pre_corre/off_ex_qeg.data", raw_ascii);

		vec s_togodot = drdv_qeg / (-r_qeg % tan(gamma_qeg));
		s_togodot.save("E:/zaixianguihua/pre_corre/sdot.data", raw_ascii);

		double s_togo_qeg = 0;
		double h = (vf - vt) / (n_qeg - 1);
		for (uword j = 0; j < n_qeg / 2; j++)
		{
			s_togo_qeg += (s_togodot[2 * j] + 4 * s_togodot[2 * j + 1] + s_togodot[2 * j + 2]) * h / 3;
		}
		
		//计算航程误差
		s_togo_rest[i] = s_togot + s_togo_qeg;
		//
		//记录2019.6.6

		//
		if (fabs(s_togo_rest[i]) < 1e-3)
		{
			sigma_init = sigma[i];
			ex.save("E:/zaixianguihua/pre_corre/off_ex.data", raw_ascii);
			xishu.save("E:/zaixianguihua/pre_corre/xishu.data", raw_ascii);
			ex_qeg.save("E:/zaixianguihua/pre_corre/off_ex_qeg.data", raw_ascii);
			break;
		}
		else
		{
			if (i > 0)
			{
				sigma[i + 1] = sigma[i] - (sigma[i] - sigma[i - 1]) / (s_togo_rest[i] - s_togo_rest[i - 1])*s_togo_rest[i];
			}
			if (sigma[i + 1] < 0)
			{
				sigma[i + 1] = 0;
			}
			if (sigma[i + 1] > 40 * pi /180)
			{
				sigma[i + 1] = 40 * pi / 180;
			}
		}


	}



}



double offLineTrack::calcu_mid_r(double v)
{
	vec aero(4);

	aero = calcu_aero(1, v);
	double hqdot = 2 * hs*log(sqrt(1.225)*pow(v*Vc, 3.15)*9.4369e-5 / (Qdot_max ));
	double hq = hs * log(1.225*pow(v*Vc, 2) / 2 / q_max);
	double V = v * Vc;
	double hn = hs * log(1.225*S*V*V*sqrt(aero(2)*aero(2) + aero(3)*aero(3)) / 2 / n_max / g0 / m);	//法向过载
	double hqeg = newton_h(1, v, aero(2));
	
	return (Re + (max(max(hq, hn),hqdot) + 4*hqeg) / 5) / Re;

}