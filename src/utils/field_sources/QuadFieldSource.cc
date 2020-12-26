//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//   QuadFieldSource.cc
//
// AUTHOR
//    A. Shishlo
//
// CREATED
//    12/25/2020
//
// DESCRIPTION
//    A class is an implementation of ShiftedFieldSource class.
//    This class describes hard edge quad field. It needs the quad length and
//    gradient in [T/m]. The sign of the gradient describes what quad it will
//    be - horizontally or vertically focusing. The center of the quad is 
//    at (0,0,0).
//
///////////////////////////////////////////////////////////////////////////

#include "orbit_mpi.hh"
#include "BufferStore.hh"

#include <cstdlib>
#include <iostream>
#include <cmath>
#include <cfloat>

#include "ShiftedFieldSource.hh"
#include "QuadFieldSource.hh"

using namespace OrbitUtils;

QuadFieldSource::QuadFieldSource(): ShiftedFieldSource()
{
	length = 0.;
	gradient = 0.;
}

QuadFieldSource::~QuadFieldSource()
{
}

/** Sets the length of the quad */
void QuadFieldSource::setLength(double length)
{
	this->length = length;
}

/** Returns the length of the quad */
double QuadFieldSource::getLength()
{
	return length;
}

/** Sets the gradient of the quad in [T/m] */
void QuadFieldSource::setGradient(double gradient)
{
	this->gradient = gradient;
}

/** Returns the gradient of the quad in [T/m] */
double QuadFieldSource::getGradient()
{
	return gradient;
}

/** Returns inner components of the electric and magnetic filds. */
void QuadFieldSource::getInnerElectricMagneticField(
	  double x, double y, double z, double t, 
		double& E_x, double& E_y, double& E_z,
		double& H_x, double& H_y, double& H_z){

		E_x=0.;	E_y=0.;	E_z=0.;
		H_x=0.;	H_y=0.;	H_z=0.;

		if(z < -length/2 || z > length/2) { return; }
		
		H_x = gradient*y;
		H_y = -gradient*x;
}