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
//    It should be sub-classed on Python level and implements 
//    getElectricMagneticField(x,y,z,t) method (returns (Ex,Ey,Ez,Bx,By,Bz)).
//    The results of these methods will be available from the c++ level.
//    This is an example of embedding Python in C++ Orbit level.
//
//    ATTENTION: Using this class in real calculations is not wise! 
//               It is slow, because it delegates the field calculations
//               to the Python level. It is only for prototyping and 
//               debugging only.
//
///////////////////////////////////////////////////////////////////////////
#ifndef PY_BASE_FIELD_SOURCE_H
#define PY_BASE_FIELD_SOURCE_H

#include "Python.h"

#include "BaseFieldSource.hh"

namespace OrbitUtils{
	
	
  /**
	   The base class for Python implementation of a field source.
		 It should be sub-classed on Python level and implements 
		 getElectricMagneticField(x,y,z,t) method (returns (Ex,Ey,Ez,Bx,By,Bz)).
		 The results of these methods will be available from the c++ level.
		 This is an example of embedding Python in C++ Orbit level.
		 ATTENTION: It is only for prototyping and debugging only. It is slow.
	*/
	
	
	class  PyBaseFieldSource: public BaseFieldSource
	{
		public:
		
			/** Constructor. Zero fields inside. */
			PyBaseFieldSource(PyObject* py_wrapperIn);
			
			/** Destructor. */
			~PyBaseFieldSource();
		
			/** Returns E and B as functions of the x,y,z and t-time. */
			void getElectricMagneticField(
				double x, double y, double z, double t, 
				double& fe_x, double& fe_y, double& fe_z,
				double& fm_x, double& fm_y, double& fm_z);

	};
};

#endif
