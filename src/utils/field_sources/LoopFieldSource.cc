//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//   LoopFieldSource.cc
//
// AUTHOR
//    A. Shishlo
//
// CREATED
//    04/05/2022
//
// DESCRIPTION
//    A class is an implementation of ShiftedFieldSource class.
//    This class describes the 3D field from the loop with current I. 
//    It needs the value of the current and the radius of the loop.
//    The center of the loop is at (0,0,0) and the loop plane is x-y. 
//    For positive current value the field at the center is directed 
//    along the positive direction of the z-axis.
// 
//
///////////////////////////////////////////////////////////////////////////

#include "orbit_mpi.hh"
#include "BufferStore.hh"

#include <cstdlib>
#include <iostream>
#include <math.h>
#include <float.h>

#include "ShiftedFieldSource.hh"
#include "LoopFieldSource.hh"
#include "elliptint.hh"

#include "OrbitConst.hh"

using namespace OrbitUtils;

LoopFieldSource::LoopFieldSource(): ShiftedFieldSource()
{
	radius = 1.0;
	current = 0.;
	
	//maximal distance along z-axis
	z_max = DBL_MAX;
	
	//maximal distance from the center in x-y plane
	r_max = DBL_MAX;;	
	
	
}

LoopFieldSource::~LoopFieldSource()
{
}

/** Sets the radius of the loop */
void LoopFieldSource::setRadius(double radius)
{
	this->radius = radius;
}

/** Returns the radius of the loop */
double LoopFieldSource::getRadius()
{
	return radius;
}

/** Sets the current in [A] */
void LoopFieldSource::setCurrent(double current)
{
	this->current = current;
}

/** Returns the  current in [A] */
double LoopFieldSource::getCurrent()
{
	return current;
}

/** Sets the max distance in z direction */
void LoopFieldSource::setMaxZ(double z_max)
{
	this->z_max = z_max;
}

/** Returns the max distance in z direction */
double LoopFieldSource::getMaxZ()
{
	return z_max;
}

/** Sets the max distance in x-y plane */
void LoopFieldSource::setMaxR(double r_max)
{
	this->r_max = r_max;
}

/** Returns the max distance in x-y plane */
double LoopFieldSource::getMaxR()
{
	return r_max;
}

/** Returns inner components of the electric and magnetic fields. */
void LoopFieldSource::getInnerElectricMagneticField(
	double x, double y, double z, double t, 
	double& E_x, double& E_y, double& E_z,
	double& H_x, double& H_y, double& H_z)
{
	
	E_x=0.;	E_y=0.;	E_z=0.;
	H_x=0.;	H_y=0.;	H_z=0.;

	if(current == 0.) { return; }	
	double r2 = x*x + y*y;
	double r = sqrt(r2);

	if(r > r_max || fabs(z) > z_max) { return; }
	
	double radius2 = radius*radius;
	double z2 = z*z;
	double z_a_r_2 = z2 + pow(radius+r,2);
	double z_r_m_a_2 = z2 + pow(radius-r,2);
	double k_coeff2 = 4*r*radius/z_a_r_2;
	double k_coeff = sqrt(k_coeff2);
	
	//  Legendre elliptic integral of first and second kind 
	double E1 = EllipticalIntegrals::ellf(OrbitConst::PI/2,k_coeff);
	double E2 = EllipticalIntegrals::elle(OrbitConst::PI/2,k_coeff);

	double Bz = (radius2 - z2 - r2)*E2/z_r_m_a_2 +E1;
	
	double coeff =  OrbitConst::permeability*current/(2*OrbitConst::PI*sqrt(z_a_r_2));

	H_z = Bz*coeff;
	
	if(r != 0.){
		double Br = coeff*z*((z2 + r2 + radius2)*E2/z_r_m_a_2 - E1)/r;
		H_x = Br*x/r;
		H_y = Br*y/r;
	}
}