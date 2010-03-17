//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    gauss-legendre_points.hh
//
// CREATED
//    03/16/2010
//
// DESCRIPTION
//    Generates points and weights for gauss-legendre integration. 
///////////////////////////////////////////////////////////////////////////
#ifndef GAUSS_LEGENDRE_POINTS_H
#define GAUSS_LEGENDRE_POINTS_H

#include "OU_Function.hh"

namespace OrbitUtils{
	
#ifdef __cplusplus
extern "C" {
#endif	
		
double gauss_legendre_generator(int n, double a, double b, OrbitUtils::Function* fn);

#ifdef __cplusplus
}
#endif	

};


#endif
