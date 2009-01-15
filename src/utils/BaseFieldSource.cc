//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    BaseFieldSource.cc
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


#include "BaseFieldSource.hh"

using namespace OrbitUtils;

BaseFieldSource::BaseFieldSource(): CppPyWrapper(NULL)
{
}

BaseFieldSource::~BaseFieldSource()
{
}

void BaseFieldSource::getElectricMagneticField(double x, double y, double z, double t, 
		double& E_x, double& E_y, double& E_z,
		double& H_x, double& H_y, double& H_z){
	
	E_x=0.;	E_y=0.;	E_z=0.;
	H_x=0.;	H_y=0.;	H_z=0.;
	
}


