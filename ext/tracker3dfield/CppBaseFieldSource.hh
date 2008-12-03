//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    PyBaseFieldSource.hh
//
// CREATED
//    04/21/2003
//
// DESCRIPTION
//    The base class for Python implementation of a field source. 
//    It should be sub-classed on Python level and implement 
//    getElectricField(x,y,z,t) and getMagneticField (x,y,z,t) methods.
//    The results of these methods will be available from the c++ level.
//    This is an example of embedding Python in C++ Orbit level.
//
//    ATTENTION: Using this class in real calculations is not wise! 
//               It is slow, because it delegates the field calculations
//               to the Python level. It is for prototyping and 
//               debugging only.
//
///////////////////////////////////////////////////////////////////////////
#ifndef CPPBASEFIELDSOURCE_HH_
#define CPPBASEFIELDSOURCE_HH_

#include "Python.h"

#include "BaseFieldSource.hh"

namespace OrbitUtils{
	
	class  CppBaseFieldSource: public BaseFieldSource
	{
		public:
		
			CppBaseFieldSource();
			~CppBaseFieldSource();
		
			void getElectricField(double x, double y, double z, double t, double& f_x, double& f_y, double& f_z);
			void getMagneticField(double x, double y, double z, double t, double& f_x, double& f_y, double& f_z);

	};
};




#endif /*CPPBASEFIELDSOURCE_HH_*/



