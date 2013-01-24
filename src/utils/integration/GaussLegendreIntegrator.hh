//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    GaussLegendreIntegrator.hh
//
// AUTHOR
//    A. Shishlo
//
// CREATED
//    01/23/2013
//
// DESCRIPTION
//    The integrator for the Gauss-Legendre schema. 
//
///////////////////////////////////////////////////////////////////////////
#ifndef ORBIT_UTILS_GAUSS_LEGENDRE_INTEGRATOR_H
#define ORBIT_UTILS_GAUSS_LEGENDRE_INTEGRATOR_H

#include <iostream> 
#include <fstream>
#include <cstdlib>
#include <cmath>

#include "CppPyWrapper.hh"

#include "OU_Function.hh"
#include "OU_SplineCH.hh"

using namespace std;

namespace OrbitUtils{
		
	class  GaussLegendreIntegrator : public CppPyWrapper
	{
	public:
		//-----------------------------------------
		//the public methods of the GaussLegendreIntegrator class
		//-----------------------------------------
		
		/** This constructor creates an integrator with 1024 points and (0,1) limits */
		GaussLegendreIntegrator();
		
		/** This constructor creates an integrator with nPoints points and (0,1) limits */
		GaussLegendreIntegrator(int nPoints);
		
		/** This constructor creates an integrator with nPoints points and (x_from,x_to) limits */
		GaussLegendreIntegrator(int nPoints, double x_from, double x_to);
		
		virtual ~GaussLegendreIntegrator();
		
		/** It sets the number of integration points */
		void setnPoints(int nPoints);
		
		/** Returns the number of integration points */
		int getnPoints();
		
		/** Sets the integration limits */
		void setLimits(double x_from, double x_to);
		
		/** Returns the Function class instance with integration points and weights */
		Function* getPointsAndWeightFunc();
		
		/** Calculates integral for the Function instance */
		double integral(Function* func);
		
		/** Calculates integral for the cubic Hermite spline instance */
		double integral(SplineCH* spline);
		
	private:
		//------------------------------------------
		//the private members of the GaussLegendreIntegrator class
		//------------------------------------------
		
		//the number of integration points
		int n_int_points;
		
		//
		double x0,x1;
		
		//number of pair of (x,y)
		Function* pw_finc;
		
	};
	
}

#endif
