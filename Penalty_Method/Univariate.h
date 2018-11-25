#include <bits/stdc++.h>
using namespace std;
typedef pair<double, double> ip;

double del(double x){
	if(abs(x) > 0.0001) return 0.0001*abs(x);
	else return 0.0001;
}

class Univariate
{
	double (*fx)(vector<double>&);
	int univariatefneval, n;
	vector<double> &xo, &s;
	double epsilon, epsilon1, low_bound, upper_bound;
public:
	Univariate(double (*fx)(vector<double>&), vector<double> &A, vector<double> &B, vector<double> &xo, vector<double> &s, double epsilon):
	xo(xo), s(s){
		this->fx = fx;
		this->epsilon = epsilon;
		epsilon1 = epsilon;
		n = xo.size();
		low_bound = DBL_MIN, upper_bound = DBL_MAX;
		double l, h;
		for(int i = 0; i < n; i++){					//Finding the bounds for alpha as it is the new variable
			l = (A[i]-xo[i])/s[i];					//Lower bound
			h = (B[i]-xo[i])/s[i];					//Upper bound
			if(l>h)swap(l, h);
			low_bound = max(low_bound, l);
			upper_bound = min(upper_bound, h);
		}
		univariatefneval = 0;
	}

	vector<double> optima(){
		ip bp = bpphase(low_bound, upper_bound);		//First BoundedPhase method is called to shorten the range
		ip gs = goldenSection(bp.first, bp.second);		//Then golden dection is used to optima
		double optima_alpha = (gs.first + gs.second)/2;

		vector<double> uni_optima(n);
		for(int i = 0; i < n; i++){
			uni_optima[i] = xo[i] + optima_alpha*s[i];
		}
		return uni_optima;
	}
	
private:
	double func(double alpha){
		vector<double> x(n);
		for(int i = 0; i < n; i++)x[i] = xo[i] + alpha*s[i];
		univariatefneval++;
		return fx(x);
	}

	ip bpphase(double a, double b){
		//cout<<"\tEntering BP phase...\n";
		double alpha = (a - b) * ((double) rand() / (double) RAND_MAX) + b;
		double f1, f, f3, del = 0.001;
		int i = 0;
		while(i++ < 10){
		    f1 = func(alpha-del), f = func(alpha), f3 = func(alpha+del);
		    if(f1 > f && f > f3) break;
		    if(f1 < f && f < f3) {del = -del; break;}
		    alpha = (a - b) * ((double) rand() / (double) RAND_MAX) + b;
		}
		if(i == 10 && f1<f && f>f3)return {a,b};
		int k = 0;
		double alphan = alpha + del, fn = func(alphan);
		while(fn < f && a<=alphan && alphan<=b){
			alpha = alphan; f = fn;
			alphan = alpha + pow(2,k)*del;
			fn = func(alphan);
			k++;
		}
		alpha = alpha - pow(2,k-1)*del;
		ip p = {min(alpha, alphan), max(alpha, alphan)};
		ip bp = {max(p.first, a), min(p.second, b)};
		return bp;
	}

	double alphaX(double &w, double &a, double &b){return  w*(b-a) + a;}

	ip goldenSection(double a, double b){
		double aw = 0, bw = 1, Lw = 1;
		double w1, w2, f1, f2;
		w1 = aw + 0.618*Lw; f1 = func(alphaX(w1, a, b));
		w2 = bw - 0.618*Lw; f2 = func(alphaX(w2, a, b));
		int k=1;
		while(abs(Lw) > epsilon1){
			if(f1<f2)aw=w2;
			else if(f1>f2)bw = w1;
			else aw = w2, bw = w1;
			Lw = bw-aw;
			w1 = aw + 0.618*Lw; f1 = func(alphaX(w1, a, b));
			w2 = bw - 0.618*Lw; f2 = func(alphaX(w2, a, b));
			k = k+1;
		}
		double x1 = alphaX(w2, a, b), x2 = alphaX(w1, a, b);
		return {min(x1, x2), max(x1, x2)};
	}
};



	// double NewtonRaphson(double a, double b){
	// 	//cout<<"Entering Newton Raphson...\n";
	// 	double alpha, fialpha, fiialpha, delta, alphan;

	// 	delta = del(alpha);
		
	// 	alpha = (a+b)/2; 			 					//Initial Guess
	// 	double f1 = func(alpha-delta), f = func(alpha), f3 = func(alpha+delta);
	// 	fialpha = (f3 - f1)/(2 * delta);				 					//Compute F'(alpha)
	// 	//printf("Initial Point = %2.5f\nF(alpha) = %2.5f\t\tF'(alpha) = %2.5f\n", alpha, func(alpha), fialpha);

	// 	int itr = 0, max_itr = 100;
	// 	while(itr++ <= max_itr && abs(fialpha)>epsilon){					//If fi(alphan) < epsilon, Then Terminate
	// 		alphan = alpha - f/fialpha;										//Compute New Point
	// 		if(alphan<a || b<alphan)return alpha;
	// 		alpha = alphan;
	// 		delta = del(alpha);
	//         f1 = func(alpha-delta), f = func(alpha), f3 = func(alpha+delta);
	// 		fialpha = (f3 - f1)/(2 * delta);								//Compute F'(alpha)
	// 	}
	// 	//cout<<"Newton Raphson iterations required = "<<itr<<endl;
	// 	//printf("Final Optima Point = %2.10f\n", alpha);
	// 	return alpha;
	// }