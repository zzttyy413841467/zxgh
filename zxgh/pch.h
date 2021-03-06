﻿
#ifndef PCH_H
#define PCH_H

// TODO: 添加要在此处预编译的标头
#include <iostream>
#include <armadillo>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <vector>
#include <time.h>

#include <thread>
#include <future>

using namespace arma;
using namespace std;

//常量
const double pi = acos(-1);
const double m = 35828;
const double S = 149.4;
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

//全局变量声明
extern double sig0_a;//初始段规划
extern double sig0;//初始段仿真
extern double delta_psi;
extern double sigma_mid;
extern double vv0;
extern double sigmaf;
extern double vvf;
extern double Qmax;
extern double qmax;
extern double nmax;

extern int item;
extern vec sigma_qeg;

extern mat waypoint;

void getTime();


//禁飞区
struct Noflyzone
{
public:
	double theta;
	double phi;
	double radius;
	int location;      //1为禁飞区在视线上，2为禁飞区在视线下
	Noflyzone() = default;
	Noflyzone(double theta1, double phi1, double radius1, int location1) : 
		theta(theta1),phi(phi1),radius(radius1),location(location1){}
	~Noflyzone() {};
};
extern vector<Noflyzone> nfz;

class onLineTrack
{
private:
	vector<double> stateVaraible_0;//顺序为r0,theta0,phi0,v0,gamma0,psi0
	vector<double> stateVaraible_f;//与上面一样，但只有rf,thetaf,phif,vf
	double Qdot_max;//热流密度
	double q_max;//动压
	double n_max;//过载
public:
	double sigma0;//初始段倾侧角

public:
	onLineTrack() = default;
	onLineTrack(vector<double> &state_0, vector<double> &state_f, double Qdot_max1, double q_max1, double n_max1) :
		stateVaraible_0(state_0), stateVaraible_f(state_f), Qdot_max(Qdot_max1), q_max(q_max1), n_max(n_max1) {
		getTime();
		cout << "创建在线轨迹对象" << endl;
	}
	~onLineTrack() { 
		getTime();
		cout << "在线轨迹对象销毁" << endl; }
	void calcu_sigma0();//计算初始段倾侧角
	void zoulang();//hv走廊
	void sigma_limit();//倾侧角走廊
	void calcu_qeg(vec x1);//滑翔段
	void calcu_online_track();//计算弹道
};

class offLineTrack
{
private:
	vector<double> stateVaraible_0;//顺序为r0,theta0,phi0,v0,gamma0,psi0
	vector<double> stateVaraible_f;//与上面一样，但只有rf,thetaf,phif,vf
	double Qdot_max;//热流密度
	double q_max;//动压
	double n_max;//过载
public:
	double sigma0;//初始段倾侧角

public:
	offLineTrack() = default;
	offLineTrack(vector<double> &state_0, vector<double> &state_f, double Qdot_max1, double q_max1, double n_max1) :
		stateVaraible_0(state_0), stateVaraible_f(state_f), Qdot_max(Qdot_max1), q_max(q_max1), n_max(n_max1) {
		getTime();
		cout << "创建离线轨迹对象" << endl;
	}
	~offLineTrack() {
		getTime();
		cout << "离线轨迹对象销毁" << endl;
	}
	void calcu_offline_track();
	double calcu_mid_r(double v);//计算再入走廊内部给定速度的r

};


//函数声明

vec calcu_aero(double, double);//输入无量纲地心距和速度，输出无量纲气动力
int trans(double e, vec x);
double newton_r(double x0, double v, double Cl);
double magnitude(double x);
double limit(double sigma, double v);
void calcu_K(mat ex1, mat ex2);
mat lqr_iter(mat A, mat B);
//void calcu_qeg1(vec x1);//


vec dxde1(double e, vec x);
vec dxdt1(double t, vec x);
double dvds1(double s, double v);
vec dvds2(double s, vec x,double sigma);
vec dxdt2(double t, vec x, double sigma, double alpha);
vec dxde2(double e, vec x, double sigma);
vec dxde_off_init(double e, vec x);

mat rk1(vec espan, vec x0);
mat rk2(vec , vec);
double norm_2(vec x);
double rk3(double s_togo, double v);
//double rk4(double e0, vec x0);
mat rk5(double s_togo, vec x1);
double newton_r_s(double v, double sigma);
mat rk(vec x0, vec statef);
vec interp2(double s_togo, mat Ref);
double d_mag(double theta);
double sign_decide(double sign0, vec x, vec statef);
double max(double a, double b);
double min(double a, double b);
vec newton_waypoint(Noflyzone nfz1, Noflyzone nfz2);
vec F_x(vec x, Noflyzone nfz1, Noflyzone nfz2);
mat dF_dx(vec x, Noflyzone nfz1, Noflyzone nfz2);
double newton_distance(double r1, double r2, double distance);
vec rk_waypoint(vec x0, double theta, double t_rev, int flag);
void calcu_waypoint(vec x0, Noflyzone nfz1, Noflyzone nfz2);
//mat rk6(double e0, vec x0);
mat rk_off_init(vec espan, vec x0);
double calcu_terminal_alpha(double n, vec x);

void th_calcu_reasonable_angle(promise<double> *promObj, vec x0, vec waypoint_pos);

double newton_h(double x0, double v, double Cl);


uword find_index(vec x);

#endif //PCH_H
