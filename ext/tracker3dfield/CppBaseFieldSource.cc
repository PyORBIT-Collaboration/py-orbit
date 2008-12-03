//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    CppBaseFieldSource.cc
//
// CREATED
//    04/21/2003
//
// DESCRIPTION
//    The base class for Python implementation of a field source. 
//    It should be sub-classed on Python level and implements 
//    getElectricField(x,y,z,t) and getMagneticField (x,y,z,t) methods.
//    The results of these methods will be available from the c++ level.
//    This is an example of embedding Python in C++ Orbit level.
//
//
///////////////////////////////////////////////////////////////////////////
#include "CppBaseFieldSource.hh"

#include "orbit_mpi.hh"
#include <iostream>

using namespace OrbitUtils;

CppBaseFieldSource::CppBaseFieldSource()
{
}

CppBaseFieldSource::~CppBaseFieldSource()
{
}

void CppBaseFieldSource::getElectricField(double x, double y, double z, double t, double& f_x, double& f_y, double& f_z)
{	  
	 f_x = 0.0; f_y = 0.0; f_z = 0.0; 
}

void CppBaseFieldSource::getMagneticField(double x, double y, double z, double t, double& f_x, double& f_y, double& f_z)
{
	 f_x = 2.0; f_y = 0.0; f_z = 0.0; 
	 

}

