#include <bits/stdc++.h>
#include "ConjugateGradient.h"
#include "Problems.h"
using namespace std;

// Global variables used by ConstraintPenaltyMethod
double (*f)(vector<double>&), (*g)(vector<double>&), r; 
int fneval, n;
long avgItr = 0, avgfneval = 0;
class ConstraintPenaltyMethod
{
	double eps;		//Class member, f stores the address of objective function passed and g stores the constraint functions
	vector<double> &A, &B; // Variable bounds used by ConjugateGradient(univariate search)
public:
	ConstraintPenaltyMethod(double (*fn)(vector<double>&), double (*gn)(vector<double>&), vector<double> &A, 
		vector<double> &B):A(A), B(B){
		f = fn, g = gn;
		n = A.size();
		eps = 0.001;
	}
	static double func(vector<double> &x){      //Modified func
		fneval++;
		// return modified func value with penalty
		return f(x) + r*g(x);
	}
	vector<double> optima(vector<double> x){
		int itr = 0;
		fneval = 0, r = 0.1;
		vector<double> xn(x.size());
		double fvalue, fvaluenew;
		fvalue = func(x);
		// create an object of class ConjugateGradient and call optima
		ConjugateGradient cg(&func, A, B);
		xn = cg.optima(x);

		fvaluenew = func(xn);
		while(itr<100 && abs(fvaluenew - fvalue)>eps){ //Check Termination condition
			r = 5*r;								   //Update the penalty constant
			x = xn;									   //Update x
			xn = cg.optima(x);						   //Call conjugate gradient
			itr++;
		}
		x = xn;
		cout<<"\tPenalty Method Iterations      = "<<itr<<endl;
		cout<<"\tFunction evaluations           = "<<fneval<<endl;
		cout<<"\tFunction evaluations/Iteration = "<<fneval/itr<<endl;
		// cout<<itr<<"\t";
		// cout<<fneval<<"\t";
		// cout<<fneval/itr<<"\t";
		avgItr += itr;
		avgfneval+=fneval;
		return x;
	}

};

int main()
{
	// cout<<std::fixed;
	// cout<<std::setprecision(6);
	#ifndef ONLINE_JUDGE
	    // for getting input from input.talphat
	    freopen("./Problems/input5.txt", "r", stdin);
	    // for writing output to output.talphat
	    freopen("./Problems/outputC5.csv", "w", stdout);
	#endif

	int q,P;
	double (*f)(vector<double>&), (*g)(vector<double>&);
	cin>>q;							//Input Function No.
	cin>>n;							//Input no. of variables
	vector<double> x(n), A(n), B(n);
	for(int i = 0; i < n; i++){     //Input Domain Bounds
		cin>>A[i]>>B[i];
	}
	switch(q){						//Setting the refrence of the function acc to Func No. to f
		case 1:
			f = &f1; g = &g1;break;
		case 2:
			f = &f2; g = &g2;break;
		case 3:
			f = &f3; g = &g3;break;
		case 4:
			f = &f4; g = &g4;break;
		case 5:
			f = &f5; g = &g5;break;
	}
	// create an object of class ConstraintPenaltyMethod
	ConstraintPenaltyMethod fp(f, g, A, B);
	cin>>P;
	double minavg = 0;
	for(int p = 0; p < P; p++){
		// get the starting point
		cout<<"Initial pt.\t\t";
		for (int i = 0; i < n; i++){
			cin>>x[i];
			cout<<x[i]<<"\t";
		}
		cout<<endl;
		// call optima Bracket Penalty Method
		vector<double> res = fp.optima(x);
		double fmin = f(res);
		for(int i = 0; i < n; i++) cout<<res[i]<<" ";
		cout<<"\tMin. value of fn = "<<fmin<<"\n\n";// << std::scientific;
		// cout<<fmin<<"\t";
		cout<<"\n";
		minavg += fmin;
	}
	cout<<"Average Bracket Penalty Itr = "<<avgItr/P<<endl;
	cout<<"Average function eval       = "<<avgfneval/P<<endl;
	cout<<"Average final minima obtained       = "<<minavg/P<<endl;
	// cout<<avgItr/P<<"\n";
	// cout<<avgfneval/P<<"\n";
	// cout<<minavg/P<<"\n";

	return 0;
}