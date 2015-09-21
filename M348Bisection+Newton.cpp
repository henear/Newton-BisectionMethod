//============================================================================
// Name        : M348Bisection+Newton.cpp
// Author      : Haocheng An
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

///////////////////////////////////////////////////////////////////////////////
// NEWTON'S METHOD
//
// state newton(double& x, double tolerance, int maxIteration)
//
// Inputs:
//   x              The initial guess at the solution.
//   tolerance      The convergence tolerance (must be > 0).
//   maxIteration   The maximum number of iterations that can be taken.
//   debug          Boolean to set debugging output.
// Outputs:
//   x              The solution.
// Return:
//   state          An error status code.
//     SUCCESS      Sucessful termination.
//     WONT_STOP    Error: Exceeded maximum number of iterations.
//     BAD_ITERATE  Error: The function had a vanishing derivative.
//
//  Remark: We assume we are given the two functions
//   f              The name of the function for which a root is sought.
//   df             The name of the derivative of the function.  The derivative
//                  of f must be computed by hand and coded correctly as df,
//                  or newton will not work!
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>
#include <math.h>
using namespace std;

enum state {SUCCESS=0, WONT_STOP, BAD_DATA, BAD_ITERATE};

double f(double x) { return x*(pow(x*x+1,-0.5)); }
double df(double x) { return pow(x*x+1,-1/2)-x*x*pow(x*x+1,-1.5); }
double sgn(double x){
	if(x>0){
		return 1;
	}
	else if(x<0){
		return -1;
	}
	else if(x==0){
		return 0;
	}else return x;
}
state newton(double& x, double tolerance, int maxIteration, int debug, double a, double b) {
  if(sgn(f(a))*sgn(f(b))>0.0){
	  return BAD_DATA;
  }
  double fa=f(a);
  double fb=f(b);

  if(a > b) {
      double c = a, fc = fa;
      a = b; fa = fb;
      b = c; fb = fc;
    }
  double dx=(b-a)/2;
  x=a+dx;

  if(debug) {
    cout << "Guess: x = " << x << endl;
  }

  for(int iteration = 1; iteration <= maxIteration; iteration++) {
    double dfx = df(x);
    while(dfx == 0.0) {
    	if(f(x)<0){
    		a=x;
    	}else{
    	b=x;
    	}
    	dx=(b-a)/2;
    	x=a+dx;
    }

    dx = -f(x)/dfx;
    x += dx;
    if(x<=a||x>=b){
    	dx = -f(x)/dfx;
    	    x += dx;
    }
    if(debug) {
      cout << "Iter " << iteration << ": x = " << x << ", dx = " << dx << endl;
    }

    if(fabs(dx) < tolerance||f(x)==0) return SUCCESS;
  }
  return WONT_STOP;
}

int main() {
  double root,tol,a,b;
  int maxIter,debug;
  state s;

  // Input

  cout << "Enter guesses of 2 roots, represented by a and b respectively: ";
  cin >> a;
  cin >> b;
  cout << "Enter tolerance,  max iteration, and debug flag: ";
  cin >> tol >> maxIter >> debug;

  // Solve

  s = newton(root, tol, maxIter, debug, a, b);

  // Report results

  switch(s) {
  case SUCCESS: {
    int prec = (int) (log10(1.0/tol) + log10(fabs(root))) + 1;
    cout << "The root is " << setprecision(prec) << root << endl;
    return 0;
  }
  case WONT_STOP:
    cout << "ERROR: Failed to converge in " << maxIter << " iterations!"
         << endl;
    break;
  case BAD_ITERATE:
    cout << "ERROR: Obtained a vanishing derivative!" << endl;
    break;
  default:
    cout << "ERROR: Coding error!" << endl;
  }
  return 1;
}
