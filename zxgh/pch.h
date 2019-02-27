
#ifndef PCH_H
#define PCH_H

// TODO: 添加要在此处预编译的标头
#include <iostream>
#include <armadillo>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <vector>

#include <thread>
#include <future>

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
const double tc = sqrt(Re / g0);

const double H0 = 121900;
const double V0 = 7630;
const double Hf = 26500;
const double Vf = 900;


const double thetaf = 55 * pi / 180;
const double phif = 20 * pi / 180;

//禁飞区
typedef struct tagNoflyzone
{
	double theta;
	double phi;
	double radius;
	int location;      //1为禁飞区在视线上，2为禁飞区在视线下
}Noflyzone;
extern vector<Noflyzone> nfz;


//函数声明
double calcu_sigma0();//初始段倾侧角
vec calcu_aero(double, double);//输入无量纲地心距和速度，输出无量纲气动力
int trans(double e, vec x);
double newton_r(double x0, double v, double Cl);
double magnitude(double x);
void calcu_qeg(vec x1);
double limit(double sigma, double v);
void calcu_K(mat ex1, mat ex2);
mat lqr_iter(mat A, mat B);

vec dxde1(double e, vec x);
vec dxdt1(double t, vec x);
double dvds1(double s, double v);
vec dvds2(double s, vec x,double sigma);
vec dxdt2(double t, vec x, double sigma, double alpha);

mat rk1(vec espan, vec x0);
mat rk2(vec x0);
double norm_2(vec x);
double rk3(double s_togo, double v);
mat rk4(double s_togo, vec x0);
mat rk5(double s_togo, vec x1);
double newton_r_s(double v, double sigma);
mat rk(vec x0);
vec interp2(double s_togo, mat Ref);
double d_mag(double theta);
double sign_decide(double sign0, vec x);
double max(double a, double b);
double min(double a, double b);
vec newton_waypoint(Noflyzone nfz1, Noflyzone nfz2);
vec F_x(vec x, Noflyzone nfz1, Noflyzone nfz2);
mat dF_dx(vec x, Noflyzone nfz1, Noflyzone nfz2);
double newton_distance(double r1, double r2, double distance);
vec rk_waypoint(vec x0, double theta, double t_rev, int flag);
void calcu_waypoint(vec x0, Noflyzone nfz1, Noflyzone nfz2);

void th_calcu_reasonable_angle(promise<double> *promObj, vec x0, vec waypoint_pos);

void zoulang();
double newton_h(double x0, double v, double Cl);
void sigma_limit();

//全局变量声明
extern double sig0_a;//初始段规划
extern double sig0;//初始段仿真
extern double delta_psi;
extern double sigma_mid;
extern double vv0;
extern double sigmaf;

extern int item;
extern vec sigma_qeg;


#endif //PCH_H
