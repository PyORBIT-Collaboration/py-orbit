//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    MathPolynomial.cc
//
// AUTHOR
//    T. Gorlov
//
// CREATED
//    06/28/2008
//
// DESCRIPTION
//    This class is a collection of the static methods to calculate 
//    different polynomials.
//
///////////////////////////////////////////////////////////////////////////
#include "orbit_mpi.hh"
#include "MathPolynomial.hh"

using namespace OrbitUtils;

int MathPolynomial::fact_n_max = 30;
int* MathPolynomial::fact_arr = MathPolynomial::getFactorialArr();

/** The method calculates the Hermite polynomial of order n and value x. */
double MathPolynomial::ReHermite(int n, double x){
	double a,b,ff;
	
	a=1.0;
	b=2.0*x;
	
	if(n==0) {
		ff=a; 
	}else if(n==1) {
		ff=b;
	} else if(n>=2)	{
		for(int i=2; i<=n;i++)	{
			ff=2.0*(x*b-(i-1)*a);
			a=b;
			b=ff;
		}
	}
	return ff;
}

/** The method calculates the complex Hermite polynomial of order n and value x. */
tcomplex MathPolynomial::ComplexHermite(int n, tcomplex x){
	tcomplex a,b,ff;
	
	a=1.0;
	b=2.0*x;
	
	if(n==0){
		ff=a;
	} else if(n==1) {
		ff=b;
	} else if(n>=2)	{
		for(int i=2; i<=n;i++)	{
			ff=2.0*(x*b-double(i-1)*a);
			a=b;
			b=ff;
		}
	}
	return ff;
}

/** The method calculates the factorial of n = n! */
int MathPolynomial::Factorial(int n){
	if(n >= 0 && n <= fact_n_max) return fact_arr[n];
  ORBIT_MPI_Finalize("pyORBIT Utils: MathPolynomial::Factorial  0<= n <= n_max. Stop.");
	return 0;
}

/** The method calculates the array with n! */
int* MathPolynomial::getFactorialArr(){
	fact_arr = new int[fact_n_max+1];
	fact_arr[0] = 1;
	for(int i=1; i<=fact_n_max; i++){
		fact_arr[i]*=i*fact_arr[i-1];
	}
	return fact_arr;
}

