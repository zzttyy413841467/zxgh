#include "pch.h"

//初始段规划
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

//计算初始段参考轨迹
mat rk2(vec x0, vec statef)
{
	double thetaf = statef[1];
	double phif = statef[2];
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

//2范数
double norm_2(vec x)
{
	double s = 0;
	for (uword i = 0; i < x.n_cols; i++)
	{
		s += x(i)*x(i);
	}
	return sqrt(s);
}

//滑翔段规划(单)
double rk3(double s_togo, double v)
{	
	double x = v;
	double k1, k2, k3, k4;
	uword i;
	double s = 0;
	vec sspan = linspace(s_togo, 0, 3000);
	double h = -s_togo / 2999;
	for (i = 0; i < 2999; i++)
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

////滑翔段规划(e为自变量)
//double rk4(double e0, vec x0)
//{
//	vec x = x0;
//	vec k1(3), k2(3), k3(3), k4(3);
//	uword i = 0;
//	uword n = 500;
//	double ef = 1.0 / ((Re + Hf) / Re) - pow(Vf / Vc, 2) / 2;
//	vec espan = linspace(e0, ef, n);
//	double h = (ef - e0) / (n - 1);
//	double e = 0;
//	double sigma0 = sig0_a / 180 * pi;
//	double sigma;
//	double vm = 4200 / Vc;
//	double vf = Vf / Vc;
//	double v = 0;
//	for (i = 0; i < n - 1; i++)
//	{
//		e = espan(i);
//		v = sqrt(2 * (1 / x(0) - e));
//		
//		if (v > vm)
//		{
//			sigma = limit(sigma0 + (sigma_mid - sigma0) / (vm - vv0)*(v - vv0), v);
//		}
//		else
//		{
//			sigma = limit(sigmaf + (sigma_mid - sigmaf) / (vm - vf)*(v - vf), v);
//		}
//		k1 = dxde2(e, x, sigma);
//		k2 = dxde2(e + h / 2, x + h * k1 / 2, sigma);
//		k3 = dxde2(e + h / 2, x + h * k2 / 2, sigma);
//		k4 = dxde2(e + h, x + h * k3, sigma);
//		x = x + (h / 6)*(k1 + 2 * k2 + 2 * k3 + k4);
//
//	}
//
//	return x(2);
//}

//滑翔段仿真
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
	//**********************************
	for (i = 0; i < num - 1; i++)
	{
		ex3(i + 1, 3) = atan((ex3(i + 1, 2) - ex3(i, 2)) / (-h * ex3(i, 2)));
	}
	//**********************************
	return ex3;
}


//mat rk6(double e0, vec x0)
//{
//	vec x = x0;
//	vec k1(3), k2(3), k3(3), k4(3);
//	uword i = 0;
//	uword n = 500;
//	double ef = 1.0 / ((Re + Hf) / Re) - pow(Vf / Vc, 2) / 2;
//	vec espan = linspace(e0, ef, n);
//	double h = (ef - e0) / (n - 1);
//	double e = 0;
//	double sigma0 = sig0_a / 180 * pi;
//	double sigma;
//	double vm = 4200 / Vc;
//	double vf = Vf / Vc;
//	double v = 0;
//	mat ex2(n, 6);
//	for (i = 0; i < n - 1; i++)
//	{
//		e = espan(i);
//		v = sqrt(2 * (1 / x(0) - e));
//		ex2(i, 5) = e;
//		ex2(i, 0) = x(2);
//		ex2(i, 1) = x(0);
//		ex2(i, 3) = x(1);
//		ex2(i, 2) = v;
//
//
//		if (v > vm)
//		{
//			sigma = limit(sigma0 + (sigma_mid - sigma0) / (vm - vv0)*(v - vv0), v);
//		}
//		else
//		{
//			sigma = limit(sigmaf + (sigma_mid - sigmaf) / (vm - vf)*(v - vf), v);
//		}
//		ex2(i, 4) = sigma;
//		k1 = dxde2(e, x, sigma);
//		k2 = dxde2(e + h / 2, x + h * k1 / 2, sigma);
//		k3 = dxde2(e + h / 2, x + h * k2 / 2, sigma);
//		k4 = dxde2(e + h, x + h * k3, sigma);
//		x = x + (h / 6)*(k1 + 2 * k2 + 2 * k3 + k4);
//
//	}
//	e = espan(i);
//	v = sqrt(2 * (1 / x(0) - e));
//	ex2(i, 5) = e;
//	ex2(i, 0) = x(2);
//	ex2(i, 1) = x(0);
//	ex2(i, 3) = x(1);
//	ex2(i, 2) = v;
//
//	return ex2;
//}


//牛顿法(rk5())
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

//全段仿真
mat rk(vec x0,vec statef)
{
	getTime();
	cout << "再入仿真开始" << endl;
	double thetaf = statef[1];
	double phif = statef[2];
	uword n = 6001;
	vec tspan = linspace(0, 3, n);
	double h = 3.0 / (n - 1);
	vec k1(6), k2(6), k3(6), k4(6);
	uword i;
	mat xx(n,7);
	//有效导航比
	const double N_P = 4;

	vec x = x0;
	double t = 0;

	mat Ref;
	Ref.load("E:/zaixianguihua/pre_corre/Ref.data");
	double ma;
	double alpha;
	double sigma;
	double sigma0 = sig0_a / 180 * pi;
	double vm = 4200 / Vc;
	double vf = statef[3];
	
	vec sigma_x(n);
	vec alpha_x(n);
	sigma_x(0) = sig0;

	double s_togo;
	double dsigma;
	double dalpha;
	vec re(9);


	for (i = 0; i < n; i++)
	{
		//中间某时刻突然发现一个禁飞区(bool)
		//if (i == 1500)
		//{
		//	getTime();
		//	cout << "加入新的禁飞区" << endl;
		//	//Noflyzone nfz4(43.0 / 180.0 * pi, 13.5 / 180.0 * pi, 3 / 180.0 * pi, 2);
		//	Noflyzone nfz4(47.0 / 180.0 * pi, 20.5 / 180.0 * pi, 3 / 180.0 * pi, 1);
		//	nfz.push_back(nfz4);
		//	waypoint.row(nfz.size() - 2) = newton_waypoint(nfz[nfz.size() - 2], nfz.back()).t();
		//}
		t = tspan(i);
		xx.submat(i, 1, i, 6) = x.t();
		xx(i, 0) = t;
		if (x(3) < vf)
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
			sigma = fabs(sig0);
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

		
		s_togo = acos(cos(x(2)) * cos(phif) * cos(x(1) - thetaf) + sin(x(2)) * sin(phif));
		if (Re*s_togo < 50000)
		{
			//判断是否进入末制导阶段
			break;
		}
		re = interp2(s_togo, Ref);
		dsigma = d_mag(re(3)*(x(0) - re(0)) + re(4)*(x(3) - re(1)) + re(5)*(x(4) - re(2)));
		dalpha = d_mag(re(6)*(x(0) - re(0)) + re(7)*(x(3) - re(1)) + re(8)*(x(4) - re(2)));
		//sigma = limit(sigma + dsigma, x(3));

		sigma = sigma + dsigma;
		alpha = alpha + dalpha;
		
		alpha_x(i) = alpha;
		if (i > 0)
		{
			sigma_x(i) = sigma * sign_decide(sign(sigma_x(i - 1)), x, statef);
		}

		k1 = dxdt2(t, x, sigma_x(i), alpha);
		k2 = dxdt2(t + h / 2, x + h * k1 / 2, sigma_x(i), alpha);
		k3 = dxdt2(t + h / 2, x + h * k2 / 2, sigma_x(i), alpha);
		k4 = dxdt2(t + h, x + h * k3, sigma_x(i), alpha);
		x= x + (h / 6)*(k1 + 2 * k2 + 2 * k3 + k4);

	}

	double R_TM;
	double R_TM1;
	double R_TM2;
	double R_TM3;
	double V_T;
	double V_TM1;
	double V_TM2;
	double V_TM3;
	double V_close;
	double r;
	double phi;
	double theta;
	double gamma;
	double psi_t;
	double psi;
	double lambdadot;

	double betadot;

	double n1;
	double u0;
	double v;
	double m_g = 2;
	double n_g = 1;
	double n_vertical;
	double n_lateral;
	double gammaf = -60 / 180 * pi;
	double alpha1 = 0;
	double R_TM_old = Re;
	//末制导
	for (; i < n; i++)
	{
		t = tspan(i);
		xx.submat(i, 1, i, 6) = x.t();
		xx(i, 0) = t;

		v = x(3);
		ma = v * Vc / 340;
		r = x(0);
		theta = x(1);
		phi = x(2);
		psi = x(5);
		gamma = x(4);
		s_togo = acos(cos(phi)*cos(phif)*cos(theta - thetaf) + sin(phi)*sin(phif));

		R_TM = Re * s_togo;
		double R_TMa = 2 * Re * sin(s_togo / 2);

		if (R_TM > R_TM_old || r < 1)
		{
			break;
		}
		/*if (theta > thetaf)
		{
			break;
		}*/

		if (phif - phi > 0)
			psi_t = asin(cos(phif)*sin(thetaf - theta) / sin(acos(cos(phi)*cos(phif)*cos(theta - thetaf) + sin(phi)*sin(phif))));
		else
			psi_t = pi - asin(cos(phif)*sin(thetaf - theta) / sin(acos(cos(phi)*cos(phif)*cos(theta - thetaf) + sin(phi)*sin(phif))));

		R_TM_old = R_TM;
		R_TM1 = Re * (phif - phi);

		double R_TM1a = R_TMa * cos(psi_t);
		
		R_TM2 = Re * cos(phif) * (thetaf - theta);
		
		double R_TM2a = R_TMa * sin(psi_t);
		
		R_TM3 = -Re * (r - 1);
		V_T = v * Vc*cos(gamma);
		V_TM1 = -V_T * cos(psi);
		V_TM2 = -V_T * sin(psi);
		V_TM3 = v * Vc*sin(gamma);
		lambdadot = (R_TM1*V_TM2 - R_TM2 * V_TM1) / pow(R_TM, 2);
		double lambdadota = (R_TM1a*V_TM2 - R_TM2a * V_TM1) / pow(R_TMa, 2);
		betadot = (V_TM3* pow(R_TMa, 2) - R_TM3 * (R_TM1a*V_TM1 + R_TM2a * V_TM2)) / R_TMa / (pow(R_TMa, 2) + pow(R_TM3, 2));
		double beta = fabs(atan(R_TM3 / R_TMa));
		
		double aa = psi_t - psi;
		V_close = V_T * cos(psi_t - psi);
		double v_close = -(R_TM1a*V_TM1 + R_TM2a * V_TM2) / R_TMa;
	/*	if (v_close<0)
		{

			break;
		}*/
		
		double V_close = (R_TMa * v_close + R_TM3 * V_TM3) / sqrt(pow(R_TMa, 2) + pow(R_TM3, 2));

		n1 = N_P * v_close * lambdadot;
		double n2 = 22 * V_close * betadot + 5 * (gamma - gammaf) * V_close / (R_TM / v_close);

		u0 = gamma * (m_g + n_g + 3) / R_TMa - gammaf * (m_g + 1)*(n_g + 1) / R_TMa - (r*Re - Re) * (m_g + 2)*(n_g + 2) / pow(R_TMa, 2);
		//n_vertical = (u0*pow(v*Vc, 2) - (pow(v*Vc, 2) / (r*Re) - g0))*cos(gamma);
		n_vertical = (n2 - (pow(v*Vc, 2) / (r*Re) - g0))*cos(gamma);
		n_lateral = (-n1 - pow(v*Vc, 2) / (r*Re)*cos(gamma)*sin(psi)*tan(phi))*cos(gamma);

		sigma = atan(n_lateral / n_vertical);

		double n_control = sqrt(pow(n_vertical, 2) + pow(n_lateral, 2));
		alpha1=calcu_terminal_alpha(n_control, x);

		//需要修改
		if (i > 0)
		{
			//sigma_x(i) = sigma * sign_decide(sign(sigma_x(i - 1)), x, statef);
			sigma_x(i) = sigma;
		}

		k1 = dxdt2(t, x, sigma_x(i), alpha1);
		k2 = dxdt2(t + h / 2, x + h * k1 / 2, sigma_x(i), alpha1);
		k3 = dxdt2(t + h / 2, x + h * k2 / 2, sigma_x(i), alpha1);
		k4 = dxdt2(t + h, x + h * k3, sigma_x(i), alpha1);
		x = x + (h / 6)*(k1 + 2 * k2 + 2 * k3 + k4);

	}

	vec sigma1 = sigma_x.rows(0, i - 1);
	sigma1.save("E:/zaixianguihua/pre_corre/sigma.data", raw_ascii);
	alpha_x.save("E:/zaixianguihua/pre_corre/alpha.data", raw_ascii);
	return xx.rows(0, i - 1);
}

double calcu_terminal_alpha(double n, vec x)
{
	double v = x(3);
	double r = x(0);
	double theta = x(1);
	double phi = x(2);
	double psi = x(5);
	double gamma = x(4);

	double V = v * Vc;
	double rho = rho0 * exp(-Re * (r - 1) / hs);
	double q = 1.0 / 2 * rho*pow(V, 2);
	double Cl_p = n * m / q / S;
	double a = cl2;
	double b = cl1;
	double c = cl0 - Cl_p;

	double alpha = (-b + sqrt(b*b - 4 * a*c)) / 2 / a;
	if (alpha > 60.0 / 180 * pi)
	{
		alpha = 60.0 / 180 * pi;
	}
	return alpha;
}

double sign_decide(double sign0, vec x, vec statef)
{
	int n = (int)nfz.size();

	double thetaf;
	double phif;

	thetaf = statef[1];
	phif = statef[2];

	double delt_sigma = 7.0 / 180 * pi;
	double delt_nz = 4.0 / 180 * pi;

	double r = x(0);
	double theta = x(1);
	double phi = x(2);
	double v = x(3);
	double psi = x(5);
	double psi_t;

	/*for (int i = 0; i < n - 2; i++)
	{
		if (waypoint(i, 0) > theta)
		{
			thetaf = waypoint(i, 0);
			phif = waypoint(i, 1);
			break;
		}
	}路径点*/

	if (phif - phi > 0)
		psi_t = asin(cos(phif)*sin(thetaf - theta) / sin(acos(cos(phi)*cos(phif)*cos(theta - thetaf) + sin(phi)*sin(phif))));
	else
		psi_t = pi - asin(cos(phif)*sin(thetaf - theta) / sin(acos(cos(phi)*cos(phif)*cos(theta - thetaf) + sin(phi)*sin(phif))));

	

	double sign;

	double psi_min;
	double psi_max;

	psi_min = psi_t - delt_sigma;
	psi_max = psi_t + delt_sigma;

	double psi_nz_center;
	double psi_nz_tan;
	double distance;
	for (int i = n - 1; i >= 0; i--)
	{
		if (nfz[i].theta + 1 / 5.0 * nfz[i].radius < theta)
		{
			continue;
		}

		if (nfz[i].location == 1)
		{
			distance = acos(cos(phi)*cos(nfz[i].phi)*cos(theta - nfz[i].theta) + sin(phi)*sin(nfz[i].phi));
			psi_nz_center = asin(cos(nfz[i].phi)*sin(nfz[i].theta - theta) / sin(distance));
			psi_nz_tan = psi_nz_center + asin(nfz[i].radius / distance);

			psi_min = max(psi_min, psi_nz_tan);
			psi_max = max(psi_max, psi_min + delt_nz);
		}
		else if (nfz[i].location == 2)
		{
			distance = acos(cos(phi)*cos(nfz[i].phi)*cos(theta - nfz[i].theta) + sin(phi)*sin(nfz[i].phi));
			psi_nz_center = pi - asin(cos(nfz[i].phi)*sin(nfz[i].theta - theta) / sin(distance));
			psi_nz_tan = psi_nz_center - asin(nfz[i].radius / distance);

			psi_max = min(psi_max, psi_nz_tan);
			psi_min = min(psi_min, psi_max - delt_nz);
		}
	}

	if (psi < psi_min)
		sign = 1;
	else if (psi <= psi_max && psi >= psi_min)
		sign = sign0;
	else
		sign = -1;


	return sign;
}

double max(double a, double b)
{
	if (a < b)
		a = b;
	return a;
}

double min(double a, double b)
{
	if (a > b)
		a = b;
	return a;
}

mat rk_off_init(vec espan, vec x0)
{
	uword n = espan.n_elem;
	double ef = espan(n - 1);
	double e0 = espan(0);
	double h = (ef - e0) / (n - 1);
	double e = 0;
	int result;
	vec x = x0;
	mat x1(n, x0.n_elem);
	vec k1(3), k2(3), k3(3), k4(3);
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
		k1 = dxde_off_init(e, x);
		k2 = dxde_off_init(e + h / 2, x + h * k1 / 2);
		k3 = dxde_off_init(e + h / 2, x + h * k2 / 2);
		k4 = dxde_off_init(e + h, x + h * k3);
		x = x + (h / 6)*(k1 + 2 * k2 + 2 * k3 + k4);
	}
	mat xx(n, 1 + x0.n_elem);
	xx = join_rows(espan, x1);
	return xx.rows(0, i - 1);
}




/***********************************************************************************************************************************************************/
//路径点策略，没完。。
/***********************************************************************************************************************************************************/


vec newton_waypoint(Noflyzone nfz1, Noflyzone nfz2)
{
	vec pos(2);
	pos(0) = (nfz1.theta + nfz2.theta) / 2;
	pos(1) = (nfz1.phi + nfz2.phi) / 2;
	
	vec delta_pos = { 1e-3, 1e-3 };
	vec Fx(2);
	mat Fdotx(2, 2);

	const double eps = 1e-8;
	while (norm(delta_pos,2) / norm(pos,2) >= eps)
	{
		pos += delta_pos;
		Fx = F_x(pos, nfz1, nfz2);
		Fdotx = dF_dx(pos, nfz1, nfz2);
		delta_pos = solve(Fdotx, - Fx);
	}
	return pos;
}

vec F_x(vec x, Noflyzone nfz1, Noflyzone nfz2)
{
	double distance;
	double distance1;
	double distance2;
	if (nfz1.radius == nfz2.radius)
	{
		distance = acos(cos(nfz1.phi)*cos(nfz2.phi)*cos(nfz1.theta - nfz2.theta) + sin(nfz1.phi)*sin(nfz2.phi));
		distance1 = distance / 2;
		distance2 = distance1;
	}
	else
	{
		distance = acos(cos(nfz1.phi)*cos(nfz2.phi)*cos(nfz1.theta - nfz2.theta) + sin(nfz1.phi)*sin(nfz2.phi));
		distance1 = newton_distance(nfz1.radius, nfz2.radius, distance);
		distance2 = distance - distance1;
	}
	vec result(2);
	result(0) = cos(nfz1.phi)*cos(x(1))*cos(nfz1.theta - x(0)) + sin(nfz1.phi)*sin(x(1)) - cos(distance1);
	result(1) = cos(nfz2.phi)*cos(x(1))*cos(nfz2.theta - x(0)) + sin(nfz2.phi)*sin(x(1)) - cos(distance2);


	return result;
}

mat dF_dx(vec x, Noflyzone nfz1, Noflyzone nfz2)
{
	mat result(2, 2);
	result(0, 0) = cos(nfz1.phi)*cos(x(1))*sin(nfz1.theta - x(0));
	result(0, 1) = -cos(nfz1.phi)*sin(x(1))*cos(nfz1.theta - x(0)) + sin(nfz1.phi)*cos(x(1));
	result(1, 0) = cos(nfz2.phi)*cos(x(1))*sin(nfz2.theta - x(0));
	result(1, 1) = -cos(nfz2.phi)*sin(x(1))*cos(nfz2.theta - x(0)) + sin(nfz2.phi)*cos(x(1));
	return result;
}

double newton_distance(double r1,double r2,double distance)
{
	double x = distance / 2;
	double x0 = distance / 16;
	double eps = 1e-10;
	double f;
	double fdot;
	while (fabs(x - x0) / fabs(x) >= eps)
	{
		x0 = x;
		f = sin(r2)*sin(x) - sin(r1)*sin(distance - x);
		fdot = sin(r2)*cos(x) + sin(r1)*cos(distance - x);
		x = x - f / fdot;
	}
	return x;
}


vec rk_waypoint(vec x0, double theta, double t_rev, int flag)
{
	uword n = 10001;
	vec tspan = linspace(0, 3, n);
	double h = 3.0 / (n - 1);
	vec k1(6), k2(6), k3(6), k4(6);

	vec x = x0;
	double t = 0;

	double ma;
	double alpha;
	double sigma;
	double sigma0 = sig0_a / 180 * pi;
	double vm = 4200 / Vc;
	double vf = vvf;

	vec re(7);

	for (auto t:tspan)
	{

		//////
		if (x(1) > theta)
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

		if (flag == 1)
		{
			if (t > t_rev)
			{
				sigma = -sigma;
			}
		}
		else if (flag == 2)
		{
			if (t < t_rev)
			{
				sigma = -sigma;
			}
		}
		else
		{
			cout << "error" << endl;
		}
		

		k1 = dxdt2(t, x, sigma, alpha);
		k2 = dxdt2(t + h / 2, x + h * k1 / 2, sigma, alpha);
		k3 = dxdt2(t + h / 2, x + h * k2 / 2, sigma, alpha);
		k4 = dxdt2(t + h, x + h * k3, sigma, alpha);
		x = x + (h / 6)*(k1 + 2 * k2 + 2 * k3 + k4);

	}
	re(0) = t;
	re.rows(1, 6) = x;
	return re;
}

void calcu_waypoint(vec x0, Noflyzone nfz1, Noflyzone nfz2)
{
	vec waypoint_pos = newton_waypoint(nfz1, nfz2);

	double theta_waypoint = waypoint_pos(0);
	double phi_waypoint = waypoint_pos(1);
	vec re0 = rk_waypoint(x0, theta_waypoint, 3, 1);
	vec t_rev(1000);
	t_rev(0) = 1.0 / 2 * re0(0);
	t_rev(1) = t_rev(0) + 0.01;
	uword i;
	vec re(7);
	vec phi_theta(1000);

	double t_a;
	vec psi_wp_m(2);
	promise<double> promiseObj;
	future<double> futureObj = promiseObj.get_future();
	thread th(th_calcu_reasonable_angle, &promiseObj, x0, waypoint_pos);
	psi_wp_m(1) = futureObj.get();
	th.detach();
	for (i = 0; i < 999; i++)
	{
		re = rk_waypoint(x0, theta_waypoint, t_rev(i), 1);
		phi_theta(i) = re(3);

		if (fabs(phi_waypoint - phi_theta(i))/fabs(phi_waypoint) < 1e-5)
		{
			break;
		}
		if (i > 0)
		{
			t_a = t_rev(i) - (t_rev(i) - t_rev(i - 1)) / (phi_theta(i) - phi_theta(i - 1))*(phi_theta(i) - phi_waypoint);
			if (t_a < 0)
				t_a = 0.01;
			else if (t_a > re0(0))
				t_a = re0(0) - 0.01;
			t_rev(i + 1) = t_a;
		}
	}
	psi_wp_m(0) = re(6);
	/*for (i = 0; i < 999; i++)
	{
		re = rk_waypoint(x0, theta_waypoint, t_rev(i), 2);
		phi_theta(i) = re(3);

		if (fabs(phi_waypoint - phi_theta(i)) / fabs(phi_waypoint) < 1e-5)
		{
			break;
		}
		if (i > 0)
		{
			t_a = t_rev(i) - (t_rev(i) - t_rev(i - 1)) / (phi_theta(i) - phi_theta(i - 1))*(phi_theta(i) - phi_waypoint);
			if (t_a < 0)
				t_a = 0.01;
			else if (t_a > re0(0))
				t_a = re0(0) - 0.01;
			t_rev(i + 1) = t_a;
		}
	}
	psi_wp_m(1) = re(6);*/
	
	double psi_nz_center;
	vec psi_nz_tan(2);
	double distance;
	if (nfz1.location == 1)
	{
		distance = acos(cos(phi_waypoint)*cos(nfz1.phi)*cos(theta_waypoint - nfz1.theta) + sin(phi_waypoint)*sin(nfz1.phi));
		psi_nz_center = asin(cos(nfz1.phi)*sin(nfz1.theta - theta_waypoint) / sin(distance));
		psi_nz_tan(0) = psi_nz_center + asin(nfz1.radius / distance);
	}
	else if (nfz1.location == 2)
	{
		distance = acos(cos(phi_waypoint)*cos(nfz1.phi)*cos(theta_waypoint - nfz1.theta) + sin(phi_waypoint)*sin(nfz1.phi));
		psi_nz_center = pi - asin(cos(nfz1.phi)*sin(nfz1.theta - theta_waypoint) / sin(distance));
		psi_nz_tan(0) = psi_nz_center - asin(nfz1.radius / distance);
	}
	if (nfz2.location == 1)
	{
		distance = acos(cos(phi_waypoint)*cos(nfz2.phi)*cos(theta_waypoint - nfz2.theta) + sin(phi_waypoint)*sin(nfz2.phi));
		psi_nz_center = asin(cos(nfz2.phi)*sin(nfz2.theta - theta_waypoint) / sin(distance));
		psi_nz_tan(1) = psi_nz_center + asin(nfz2.radius / distance);
	}
	else if (nfz2.location == 2)
	{
		distance = acos(cos(phi_waypoint)*cos(nfz2.phi)*cos(theta_waypoint - nfz2.theta) + sin(phi_waypoint)*sin(nfz2.phi));
		psi_nz_center = pi - asin(cos(nfz2.phi)*sin(nfz2.theta - theta_waypoint) / sin(distance));
		psi_nz_tan(1) = psi_nz_center - asin(nfz2.radius / distance);
	}
	vec qujian(2);
	if (min(psi_wp_m) > max(psi_nz_tan) || min(psi_nz_tan) > max(psi_wp_m))
	{
		cout << "can't arrive the waypoint " << endl;
	}
	else
	{
		qujian(0) = max(min(psi_wp_m), min(psi_nz_tan));
		qujian(1) = min(max(psi_wp_m), max(psi_nz_tan));
	}
	double psi_waypoint = (qujian(0) + qujian(1)) / 2;
	cout << psi_waypoint << endl;



}

void th_calcu_reasonable_angle(promise<double> *promObj, vec x0, vec waypoint_pos)
{
	double theta_waypoint = waypoint_pos(0);
	double phi_waypoint = waypoint_pos(1);
	vec re0 = rk_waypoint(x0, theta_waypoint, 3, 2);
	double t_rev0 = re0(0);
	vec t_rev(1000);
	t_rev(0) = 2.0 / 3 * t_rev0;
	t_rev(1) = t_rev(0) - 0.01;
	uword i;
	vec re(7);
	vec phi_theta(1000);

	double t_a;
	for (i = 0; i < 999; i++)
	{
		re = rk_waypoint(x0, theta_waypoint, t_rev(i), 2);
		phi_theta(i) = re(3);

		if (fabs(phi_waypoint - phi_theta(i)) / fabs(phi_waypoint) < 1e-5)
		{
			break;
		}
		if (i > 0)
		{
			t_a = t_rev(i) - (t_rev(i) - t_rev(i - 1)) / (phi_theta(i) - phi_theta(i - 1))*(phi_theta(i) - phi_waypoint);
			if (t_a < 0)
				t_a = 0.01;
			else if (t_a > t_rev0)
				t_a = t_rev0 - 0.01;
			t_rev(i + 1) = t_a;
		}
	}
	promObj->set_value(re(6));

}

/***********************************************************************************************************************************************************/
/***********************************************************************************************************************************************************/