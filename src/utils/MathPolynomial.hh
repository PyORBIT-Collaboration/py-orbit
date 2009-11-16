#ifndef MATHPOLYNOMIAL_HH_
#define MATHPOLYNOMIAL_HH_

//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    MathPolynomial.hh
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
#include <complex>
#include <cmath>
#include "tcomplex.hh"
#include <iostream>

namespace OrbitUtils{

/**    
  This class is a collection of the static methods to calculate 
  different polynomials.
*/

  class  MathPolynomial
	{
	public:

		/** The method calculates the Hermite polynomial of order n and value x. */
		static  double 	ReHermite(int n, double x);
		
		/** The method calculates the complex Hermite polynomial of order n and value x. */
		static  tcomplex 	ComplexHermite(int n, tcomplex x);
		
		/** The method calculates the factorial of n = n! Returns integer. */
		static int Factorial(int n);
		
		/** The method calculates the factorial of n = n! Returns double. */
		static double FactorialDouble(int n);
		
	private:
		
		/** The method calculates the array with n! */
		static int* getFactorialArr();
		static double* getFactorialDoubleArr();
		
		static int fact_n_max;
		static int* fact_arr;
		
		static int fact_double_n_max;
		static double* fact_double_arr;		
		
	};

};

#endif /*MATHPOLYNOMIAL_HH_*/
