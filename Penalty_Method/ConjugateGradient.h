#include <bits/stdc++.h>
#include "Univariate.h"						//Include the Univariate file
using namespace std;
// long avgCGItr = 0, avgfneval = 0;	//GLobal variable for no. of variables and no. of fn eval
// double minavg = 0;

double dot(vector<double> &v1, vector<double> &v2){			//Dot Product of v1 & v2
	double dotprod = 0;
	int n = v1.size();
	for(int i = 0; i < n; i++){
		dotprod += v1[i]*v2[i];
	}
	return dotprod;
}
double norm(vector<double> &v){								//Norm of v
	double res = 0;
	int n = v.size();
	for(int i = 0; i < n; i++){
		res += v[i]*v[i];
	}
	return sqrt(res);
}
double norm2(vector<double> &v2, vector<double> &v1){		//Norm of v1-v2
	double res = 0;
	int n = v1.size();
	for(int i = 0; i < n; i++){
		res += pow(v2[i]-v1[i], 2);
	}
	return sqrt(res);
}
void unit(vector<double> &v){								//To obtain a unit vector in direction v
	double mod = norm(v);
	int n = v.size();
	for(int i = 0; i < n; i++)v[i]/=mod;
}

class ConjugateGradient
{
	double (*func)(vector<double>&), eps1, eps2, eps3;		//Class member, func stores the address of function passed
	int n, restart, k;
	vector<double> &A, &B;
private:
	vector<double> grad(vector<double> &x){					//Computes the grad of x using central diff. method
		double delta = 0.0001, f1, f2;
		vector<double> grad(n);
		for(int i = 0; i < n; i++){
			x[i] = x[i] - delta;   f1 = func(x);
			x[i] = x[i] + 2*delta; f2 = func(x);
			grad[i] = (f2 - f1)/(2*delta);
			x[i] = x[i] - delta;
		}
		return grad;
	}
public:
	ConjugateGradient(double (*func)(vector<double>&), vector<double> &A, vector<double> &B):
	A(A), B(B){
		this->func = func;
		this->eps1 = 1e-6, this->eps2 = 1e-8, this->eps3 = 1e-10;
		this->n = A.size();
	}
	vector<double> optima(vector<double> &x){				//Find the optima from initial pt. x
		k = 0, restart = 0;									//Initialise k, restart, cgfneval
		return cGmin(x);									//Calls the cGmin function
	}
private:
	vector<double> cGmin(vector<double> &x){				//CG method invoked by optima function to find the optima
		vector<double> g(n), gn(n);
		vector<double> s(n), sn(n);
		vector<double> xn(n);
		g = grad(x);
		for(int i = 0; i < n; i++) s[i] = -g[i];
		unit(s);
		{
			Univariate uni_gx(func, A, B, x, s, eps1);					//Performing a univariate along grad(x)
			xn = uni_gx.optima();
			gn = grad(xn);
		}
		double d;
		while(norm(gn) > eps3 && norm2(xn,x)/norm(x) > eps2){			//Termination Condition
			double r = pow((norm(gn)/norm(g)),2);						//square ratio of norm(gn)/norm(g)
			//sn = -gn + r*s;
			for(int i = 0; i < n; i++) sn[i] = -gn[i] + r*s[i];			//Setting the new search direction
			x = xn; g = gn;												//Storing the previously calculated x, grad(x)
			unit(sn);													//Getting a unit vector in direction sn
			d = dot(s,sn);												//Dot product of s & sn

			if(d>0.999 && restart <= 100){restart++;return cGmin(xn);}	//Check if too many restarts are occuring
			else if(d>0.99 && restart>100 && norm(gn) < eps3 && norm2(xn,x)/norm(x) < eps2){return xn;}
			else if(d>0.99 && restart>100){break;}

			Univariate uni_gx(func, A, B, x, sn, eps1);								//Performing a univariate along sn
			xn = uni_gx.optima();
			gn = grad(xn);
			k++;
		}
		return xn;
	}
	
};