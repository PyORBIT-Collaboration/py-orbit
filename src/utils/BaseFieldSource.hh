//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    BaseFieldSource.hh
//
// CREATED
//    04/21/2003
//
// DESCRIPTION
//    The base class for electro-magnetic field source. 
///////////////////////////////////////////////////////////////////////////
#ifndef BASE_FIELD_SOURCE_H
#define BASE_FIELD_SOURCE_H

#include "CppPyWrapper.hh"

namespace OrbitUtils{
	
	
/**
  The base class for electro-magnetic field source. It should be sub-classed. The units
  for E and B are unknown at this point. They should be defined in 
  subclasses.	
*/
	
	class  BaseFieldSource: public CppPyWrapper
	{
	public:
		
		/** Constructor. It does nothing. */
		BaseFieldSource();
		
		/** Destructor. It does nothing. */
		virtual ~BaseFieldSource();
		
		/** Returns components of the electric and magnetic filds. */
		virtual void getElectricMagneticField(
			  double x, double y, double z, double t, 
				double& E_x, double& E_y, double& E_z,
				double& H_x, double& H_y, double& H_z);

	};
};

#endif
