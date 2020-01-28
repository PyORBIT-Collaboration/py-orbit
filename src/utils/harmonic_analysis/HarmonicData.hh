//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    HarmonicData.cc
//
// AUTHOR
//    A. Shishlo
//
// CREATED
//    04/17/2018
//
// DESCRIPTION
//    Keeps the data y(x) for set of x(i), i=1,Ndata values where each x(i) 
//    is a phase between -180 to +180 deg.
//    It also has parameters of the fitting representation of this y(x) data:
//    y_fit(x) = A0 + A1*cos((PI/180.)*(x + phase_shift_1)) + 
//                    A2*cos((PI/180.)*(2*x + parse_shift_2)) +
//                    ...
//    with n_harm the maximal harmonic number in the y_fit function.
//    This class has a method getSumDiff2 that returns
//     sumDiff2 = sum((y(x(i))-y_fit(x(i)))**2, i= 1,Ndata)
//
//    The fast calculation of the sumDiff2 is the main purpose of this 
//    class. This value will be used for a fast preliminary analysis
//    of the phase scan data of RF cavities with one or few RF gaps.
//   
//    This class is used for the harmonics analysis of the scan data for 
//    the RF cavities when we measure the phase response from BPMs.
//
///////////////////////////////////////////////////////////////////////////
#ifndef ORBIT_UTILS_HARMONICDATA_H
#define ORBIT_UTILS_HARMONICDATA_H

#include "CppPyWrapper.hh"
#include "OU_Function.hh"

using namespace std;

namespace OrbitUtils{
		
	class  HarmonicData : public CppPyWrapper
	{
	public:
		//-----------------------------------------
		//the public methods of the HarmonicData class
		//-----------------------------------------
		
		/**
		Constructor of HarmonicData with order and Function with data.
		*/
		HarmonicData(int order_in, Function* inFunc);
		
		/**
		Copy Constructor of HarmonicData.
		*/		
    HarmonicData(HarmonicData* harmonicData);
		
		/**
		Destructor of HarmonicData.
		*/	
		virtual ~HarmonicData();
		
		/**
		  Sets order of the harmonic fit A0,A1,A2,..,An where order = 2*n+1.
		*/			
		void setOrder(int order_in);
		
		/**
		  Returns order of the harmonic fit A0,A1,A2,..,An where order = 2*n+1.
		*/				
		int getOrder();
		
		/**
		  Sets the new initial data for fiiting as new Function.
		*/			
		void setDataFunction(Function* inFunc);		
		
		/**
		  Sets A or phase parameter in param_arr = [A0,A1,phase1,A2,phase2, ...].
		*/		
		void setParam(int index, double val);
		
		/**
		  Returns A or phase parameter from param_arr = [A0,A1,phase1,A2,phase2, ...].
		*/		
		double getParam(int index);
		
		/**
		  Returns the number of points in the Y(x) input scan data.
		*/
		int dataSize();
		
		/**
		  Returns the Y value for data point with the indexX in the input scan data.
		*/		
		double valueY(int indexX);

		/**
		  Returns the error of Y value for data point with the indexX in the input scan data.
		*/		
		double valueErr(int indexX);	
		
		/**
		  Returns the X value for data point with the indexX in the input scan data.
		*/			
		double valueX(int indexX);
		
		/**
		  Returns the Y value of the harmonic fit.
		*/		
		double fitValueY(double x);
		
		/**
		  Sets all parameter in param_arr = [A0,A1,phase1,A2,phase2, ...] to zero.
		*/			
		void clean();
		
		/**
		  Returns the sum of square of difference between data and fit.
		*/			
		double sumDiff2();
		
		/**
		  Returns the Function with initial data points.
		*/			
		Function* getDataFunction();

	private:
		
		void init(int order_in, Function* inFunc);
		
		//------------------------------------------
		//the private members of the HarmonicData class
		//------------------------------------------
		
		int order;
		
		// param_arr = [A0,A1,phase1,A2,phase2, ...]
		double* param_arr;
		
		Function*  dataFunc;
		
	};
	
}

#endif
