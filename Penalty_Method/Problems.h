#include <bits/stdc++.h>
using namespace std;
#define pi 3.14159265359

// bracket operator
double bracket(double gval){
	if(gval<0) return pow(gval,2);
	else return 0;
}
double g1(vector<double> &x){
	int n = x.size();
	double h = -1, res = 0;
	for(int i = 0; i < n; i++)h+=x[i]*x[i];
	// y[i] = x[i] x[i]>=0
	for(int i = 0; i < n; i++)res += bracket(x[i]) + bracket(1-x[i]);
	res+=pow(h,2);
	return res;
}
double f1(vector<double> &x){
	int n = x.size();
	double res = 1;
	for(int i = 0; i < n; i++)res*=x[i];
	return -res*pow(sqrt((double)n), ((double)n));
}

double g2(vector<double> &x){
	double g[2], h[4], res = 0;
	g[0] = 1-(pow(x[0]-5, 2) + pow(x[1]-5, 2))/100;
	g[1] = (pow(x[0]-5, 2) + pow(x[1]-5, 2))/82.81-1;
	h[0] = (x[0]-13)/100;
	h[1] = (100-x[0])/100;
	h[2] = (x[1])/100;
	h[3] = (100-x[1])/100;
	for(int i = 0; i < 2; i++)res+=bracket(g[i]);
	for(int i = 0; i < 4; i++)res+=bracket(h[i]);
	return res;
}
double f2(vector<double> &x){
	return pow(x[0]-10, 3) + pow(x[1]-20, 3);
}
double g3(vector<double> &x){
	double g[2], h[4], res = 0;
	g[0] = x[1]-1-x[0]*x[0];
	g[1] = x[0] - 1 - pow(x[1]-4, 2);
	h[0] = x[0]/10;
	h[1] = 1-x[0]/10;
	h[2] = x[1]/10;
	h[3] = 1-x[1]/10;
	for(int i = 0; i < 2; i++)res+=bracket(g[i]);
	for(int i = 0; i < 4; i++)res+=bracket(h[i]);
	return res;
}
double f3(vector<double> &x){
	return -(pow(sin(2*pi*x[0]), 3)/(pow(x[0], 3))*sin(2*pi*x[1])/(x[0]+x[1]));
}
double g4(vector<double> &x){
	double g[2], h[4], res = 0;
	g[0] = 1-0.0025*(x[3]+x[5]);
	g[1] = 1-0.0025*(-x[3]+x[4]+x[6]);
	g[2] = 1-0.01*(-x[4]+x[7]);
	g[3] = -(100*x[0] - x[0]*x[5] + 833.33252*x[3] - 83333.333);
	g[4] = -(x[1]*x[3]-x[1]*x[6]-1250*x[3]+1250*x[4]);
	g[5] = -(x[2]*x[4]-x[2]*x[7]-2500*x[4]+1250000);
	h[0] = x[0]/10;
	h[1] = 1-x[0]/10;
	h[2] = x[1]/10;
	h[3] = 1-x[1]/10;
	for(int i = 0; i < 2; i++)res+=bracket(g[i]);
	for(int i = 0; i < 4; i++)res+=bracket(h[i]);
	return res;
}
double f4(vector<double> &x){
	return x[0]+x[1]+x[3];
}

double g5(vector<double> &x){
	double h[3], y[10], res = 0;
	h[0] = -10;
	for(int i = 0; i < 5; i++)h[0]+=x[i]*x[i];
	h[1] = x[1]*x[2]-5*x[3]*x[4];
	h[2] = pow(x[0], 3) + pow(x[1], 3) + 1;
	y[0] = 1-x[0]/2.3;
	y[1] = 1-x[1]/2.3;
	y[2] = 1-x[2]/3.2;
	y[3] = 1-x[3]/3.2;
	y[4] = 1-x[4]/3.2;
	y[5] = x[0]/2.3+1;
	y[6] = x[1]/2.3+1;
	y[7] = x[2]/3.2+1;
	y[8] = x[3]/3.2+1;
	y[9] = x[4]/3.2+1;
	for(int i = 0; i < 2; i++)res+=bracket(y[i]);
	res += pow(h[0],2)+pow(h[1],2)+pow(h[2],2);
	return res;
}
double f5(vector<double> &x){
	return exp(x[0]*x[1]*x[2]*x[3]*x[4]);
}

double g6(vector<double> &x){
	double g, h[2], res = 0;
	g = (pow(x[0]-5,2)+x[1]*x[1])/26-1;
	h[0] = x[0];
	h[1] = x[1];
	res+=bracket(g);
	for(int i = 0; i < 2; i++)res+=bracket(h[i]);
	return res;
}
double f6(vector<double> &x){
	return pow(x[0]*x[0]+x[1]-11, 2) + pow(x[0]+x[1]*x[1]-7, 2);
}