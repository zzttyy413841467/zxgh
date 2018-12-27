
#ifndef PCH_H
#define PCH_H

// TODO: 添加要在此处预编译的标头
#include <iostream>
#include <armadillo>
#include <cmath>
#include <iomanip>
#include <fstream>

using namespace arma;
using namespace std;

//常量
const double pi = acos(-1);
const double m = 35828;
const double S = 149.4;
const double q_max = 790000;
const double cl1 = 0.016292*(180 / pi);
const double cl2 = 0.00026024*(180 / pi)*(180 / pi);
const double cl0 = -0.041065;
const double cd1 = -0.03026;
const double cd2 = 0.86495;
const double cd0 = 0.080505;
const double rho0 = 1.225;
const double g0 = 9.81;
const double Re = 6371393;
const double mu = g0 * Re *Re;
const double hs = 7320;
const double Vc = sqrt(g0*Re);
const double k_q = 9.4369e-5*pow(Vc, 3.15);

const double H0 = 121900;
const double V0 = 7630;
const double Hf = 26500;
const double Vf = 900;

//函数声明
double calcu_sigma0();
vec calcu_aero(double, double);//输入无量纲地心距和速度，输出无量纲气动力
int trans(double e, vec x);
double newton_r(double x0, double v, double Cl);
double magnitude(double x);


vec dxde1(double e, vec x);
vec dxdt1(double t, vec x);

mat rk1(vec espan, vec x0);
mat rk2(vec x0);
double norm_2(vec x);


void zoulang();
double newton_h(double x0, double v, double Cl);

//全局变量声明
extern double sig0_a;//初始段规划
extern double sig0;//初始段仿真


//禁飞区
struct Noflyzone
{
	double theta;
	double phi;
	double radius;
};
extern Noflyzone nfz1, nfz2;

#endif //PCH_H
