//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    bessj0.hh
//
//    This subroutine calculates the First and Second Kind Bessel Function of
//    order 0,1, and n, for any real number x. The polynomial approximation by
//    series of Chebyshev polynomials is used for 0<x<8 and 0<8/x<1.
//    REFERENCES:
//    M.ABRAMOWITZ,I.A.STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS, 1965.
//    C.W.CLENSHAW, NATIONAL PHYSICAL LABORATORY MATHEMATICAL TABLES,
//    VOL.5, 1962.
//
///////////////////////////////////////////////////////////////////////////
#ifndef ORBIT_UTILS_BESSJ0_H
#define ORBIT_UTILS_BESSJ0_H

#include <cmath>

namespace OrbitUtils{
	
	double bessj0(double x);
	double bessj1(double x);
	double bessj(int n, double x);

	double bessi0(double x);
	double bessi1(double x);
	double bessi(int n, double x);
	
	double BSign(double X, double Y);

}

#endif

	 
