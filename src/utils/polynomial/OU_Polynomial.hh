//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    Polynomial.hh
//
// AUTHOR
//    A. Shishlo
//
// CREATED
//    01/31/2013
//
// DESCRIPTION
//    order = 0   y = coef[0]
//    order = 1   y = coef[0] + coef[1]*x
//    order = 2   y = coef[0] + coef[1]*x + coef[2]*x^2
//
///////////////////////////////////////////////////////////////////////////
#ifndef ORBIT_UTILS_POLYNOMIAL_H
#define ORBIT_UTILS_POLYNOMIAL_H

#include "CppPyWrapper.hh"

using namespace std;

namespace OrbitUtils{
		
	class  Polynomial : public CppPyWrapper
	{
	public:
		//-----------------------------------------
		//the public methods of the Polynomial class
		//-----------------------------------------
		Polynomial(int order_in);
		
		virtual ~Polynomial();
		
		void setOrder(int order_in);
		
		int getOrder();
		
		void setCoef(int index, double val);
		
		double getCoef(int index);
		
		double value(double x);
		
		double derivative(double x);
		
		void derivativeTo(Polynomial* derivP);
		
		void copyTo(Polynomial* p);
		
		double getMinX();
		double getMaxX();
		
		void setMinX(double x_min);
		void setMaxX(double x_max);
		
	private:
		//------------------------------------------
		//the private members of the Polynomial class
		//------------------------------------------
		
		int order;
		
		double* coef_arr;
		
		// minimal and maximal x that give a good accuracy 
		// for a polynomial representation of external data
		// By default are just negative and positive infinity.
		double x_min;
		double x_max;
		
	};
	
}

#endif
