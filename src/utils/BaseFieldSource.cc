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

void BaseFieldSource::getElectricField(double x, double y, double z, double t, double& f_x, double& f_y, double& f_z)
{
  f_x = 0.0; f_y = 0.0; f_z = 0.0; 
}

void BaseFieldSource::getMagneticField(double x, double y, double z, double t, double& f_x, double& f_y, double& f_z)
{
  f_x = 0.0; f_y = 0.0; f_z = 0.0; 
}

