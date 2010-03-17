//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    Function.hh
//
// AUTHOR
//    Y. Sato, A. Shishlo
//
// CREATED
//    12/31/2003
//
// DESCRIPTION
//    Specification and inline functions for a class that keeps
//    table y(x) and does some operation with the tables.
//    It is using linear interpolation.
//
///////////////////////////////////////////////////////////////////////////
#ifndef ORBIT_UTILS_FUNCTION_H
#define ORBIT_UTILS_FUNCTION_H

#include <iostream> 
#include <fstream>
#include <cstdlib>
#include <cmath>

#include "CppPyWrapper.hh"

using namespace std;

namespace OrbitUtils{
		
	class  Function : public CppPyWrapper
	{
	public:
		//-----------------------------------------
		//the public methods of the Function class
		//-----------------------------------------
		Function();
		
		virtual ~Function();
		
		void add(double x, double y);
		
		int getSize();
		
		double x(int ind);
		double y(int ind);
		
		double* xArr();
		double* yArr();
		
		double getMinX();
		double getMinY();
		double getMaxX();
		double getMaxY();
		
		void clean(); 
		
		double getY(double x);
		
		//this method should be used only for monotonic function
		// f(x1) < f(x2) if x1 < x2
		double getX(double y);
		
		//set the info variable info=1 const step info=0 non-const
		void setConstStep(int info);
		
		//return 1 if step on x is constant and 0 - otherwise
		int isStepConst();
		
		//sets the inverse function. The x-coordinates 
		//     of f_inv could be defined already
		void setInverse(Function* f_inv);
		
		void print(std::ostream& Out);
		void print(char* fileName);
		
		//auxilary method to create normalize cumulative function 
		//for probability distribution
		//with y_min = 0 and y_max = 1.0
		void normalize();
		
		
	private:
		//------------------------------------------
		//the private methods of the Function class
		//------------------------------------------
		void resize();
		void finalize(const char* message);
		
	private:
		//------------------------------------------
		//the private members of the Function class
		//------------------------------------------
		
		//-----------------------
		//info variables
		//-----------------------
		
		//inf_const_step = 1 if the step is const 0 - otherwise
		int inf_const_step;
		
		double x_step;
		
		//number of pair of (x,y)
		int size;
		
		int maxSize;
		int sizeChunk;
		
		//min and max values
		double xMin, xMax;
		double yMin, yMax;
		
		//x and y array 
		double* x_arr;
		double* y_arr;
		
		//MPI members
		int iMPIini; 
		int rank_MPI; 
		int size_MPI;
		
	};
	
}

#endif
