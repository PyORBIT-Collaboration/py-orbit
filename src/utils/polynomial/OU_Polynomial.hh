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
//    Calculates the polynomial:
//    order = 0   y = coef[0]
//    order = 1   y = coef[0] + coef[1]*x^2
//    etc.
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
		
	private:
		//------------------------------------------
		//the private members of the Polynomial class
		//------------------------------------------
		
		int order;
		
		double* coef_arr;
		
	};
	
}

#endif
