//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    SplineCH.hh
//
// AUTHOR
//    A. Shishlo
//
// CREATED
//    03/25/2010
//
// DESCRIPTION
//    A cubic Hermite spline. 
//    http://en.wikipedia.org/wiki/Cubic_Hermite_spline
//    y(t) = h00(t)*y0 + h10(t)*m0 + h01(t)*y1 + h11(t)*m1
//    t = (x-x0)/(x1-x0)
//    h00(t)=2*t^3 - 3*t^2 +1
//    h10(t)=t^3 - 2*t^2 +t
//    h01(t)=-2*t^3 + 3*t^2
//    h11(t)=t^3 - t^2
//    m_k = 0.5*((y[k+1]-y[k])/(x[k+1]-x[k])+(y[k]-y[k-1])/(x[k]-x[k-1]))
//    m0 = (y[1] - y[0])/(x[1]-x[0])
//    m[n-2] = (y[n-1] - y[n-2])/(x[n-1]-x[n-2])
//
///////////////////////////////////////////////////////////////////////////
#ifndef ORBIT_UTILS_SPLINE_CMP_H
#define ORBIT_UTILS_SPLINE_CMP_H

#include <iostream> 
#include <fstream>
#include <cstdlib>
#include <cmath>

#include "CppPyWrapper.hh"
#include "OU_Function.hh"

namespace OrbitUtils{
	
	class  SplineCH : public CppPyWrapper
	{
	public:
		//-----------------------------------------
		//the public methods of the SplineCH class
		//-----------------------------------------
		SplineCH();
		
		virtual ~SplineCH();
		
		int compile(OrbitUtils::Function* f);
		
		int getSize();
		
		double x(int ind);
		double y(int ind);
				
		double getY(double x);
		double getYP(double x);
		
		void print(std::ostream& Out);
		void print(const char* fileName);
		
	private:
		//------------------------------------------
		//the private methods of the SplineCH class
		//------------------------------------------
		void finalize(const char* message);
		
	private:
		//------------------------------------------
		//the private members of the SplineCH class
		//------------------------------------------

		//number of pair of (x,y)
		int size;
		
		//x and y array 
		double* x_arr;
		double* y_arr;
		
		//ranning derivetive m array
		double* m_arr;
		
		//MPI members
		int iMPIini; 
		int rank_MPI; 
		int size_MPI;
		
	};
	
}

#endif
