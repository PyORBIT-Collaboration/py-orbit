//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    BaseFieldSource.cc
//
// AUTHOR
//    Y. Sato, A. Shishlo
//
// CREATED
//    03/31/2003
//
// MODIFIED
//    09/19/2006  -> by Perkett - for use with ORBIT_SH instead of ORBIT_MPI
//                ->    added time as one of the parameters in getElectricField()
//                ->    & getMagneticField() functions
// DESCRIPTION
//    The base class for field source. Can be sub-classed.
//
///////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////
//
// include files
//
/////////////////////////////////////////////////////////////////////////////

#include "BaseFieldSource.hh"

BaseFieldSource::BaseFieldSource()
{
}

BaseFieldSource::~BaseFieldSource()
{
}

void BaseFieldSource::getElectricField(double t, double x, double y, double z, double& f_x, double& f_y, double& f_z)
{
  f_x = 0.0; f_y = 0.0; f_z = 0.0; 
}

void BaseFieldSource::getMagneticField(double t, double x, double y, double z, double& f_x, double& f_y, double& f_z)
{
  f_x = 0.0; f_y = 0.0; f_z = 0.0; 
}

///////////////////////////////////////////////////////////////////////////
//
// END OF FILE
//
///////////////////////////////////////////////////////////////////////////
