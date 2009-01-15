//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    BaseFieldSource.hh
//
// CREATED
//    04/21/2003
//
// DESCRIPTION
//    The base class for field source. It should be sub-classed. The units
//    for E and B are unknown at this point. They should be defined in 
//    subclasses.
//
///////////////////////////////////////////////////////////////////////////
#ifndef BASE_FIELD_SOURCE_H
#define BASE_FIELD_SOURCE_H

#include "CppPyWrapper.hh"

namespace OrbitUtils{
	class  BaseFieldSource: public CppPyWrapper
	{
	public:
		
		BaseFieldSource();
		virtual ~BaseFieldSource();
		
		virtual void getElectricMagneticField(double x, double y, double z, double t, 
				double& E_x, double& E_y, double& E_z,
				double& H_x, double& H_y, double& H_z);

	};
};

#endif
