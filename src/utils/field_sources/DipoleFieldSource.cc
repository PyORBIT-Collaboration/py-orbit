//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//   DipoleFieldSource.cc
//
// AUTHOR
//    A. Shishlo
//
// CREATED
//    12/25/2020
//
// DESCRIPTION
//    A class is an implementation of ShiftedFieldSource class.
//    This class describes dipole magnetic field. It needs dipole sizes
//    [m] in three directions and the field strength in [T].
//    The center of the dipole is at [0,0,0]. 
//
///////////////////////////////////////////////////////////////////////////

#include "orbit_mpi.hh"
#include "BufferStore.hh"

#include <cstdlib>
#include <iostream>
#include <cmath>
#include <cfloat>

#include "ShiftedFieldSource.hh"
#include "DipoleFieldSource.hh"

using namespace OrbitUtils;

DipoleFieldSource::DipoleFieldSource(): ShiftedFieldSource()
{	
	sizeX = 0.;	
	sizeY = 0.;	
	sizeZ = 0.;	
	
	fieldX = 0.;	
	fieldY = 0.;	
	fieldZ = 0.;		
}

DipoleFieldSource::~DipoleFieldSource()
{
}

/** Sets the sizes of the dipole field */
void DipoleFieldSource::setSizes(double sizeX, double sizeY, double sizeZ)
{
	this->sizeX = sizeX;
	this->sizeY = sizeY;
	this->sizeZ = sizeZ;
}

/** Returns the sizes of the dipole field */
void DipoleFieldSource::getSizes(double& sizeX, double& sizeY, double& sizeZ)
{
	sizeX = this->sizeX;
	sizeY = this->sizeY;
	sizeZ = this->sizeZ;
} 

/** Sets the fields of the dipole in [T] */
void DipoleFieldSource::setFields(double fieldX, double fieldY, double fieldZ)
{
	this->fieldX = fieldX;
	this->fieldY = fieldY;
	this->fieldZ = fieldZ;
}

/** Returns the fields of the dipole in [T] */
void DipoleFieldSource::getFields(double& fieldX, double& fieldY, double& fieldZ)
{
	fieldX = this->fieldX;
	fieldY = this->fieldY;
	fieldZ = this->fieldZ;
} 

/** Returns inner components of the electric and magnetic filds. */
void DipoleFieldSource::getInnerElectricMagneticField(
	  double x, double y, double z, double t, 
		double& E_x, double& E_y, double& E_z,
		double& H_x, double& H_y, double& H_z){

		E_x=0.;	E_y=0.;	E_z=0.;
		H_x=0.;	H_y=0.;	H_z=0.;

		if(x < -sizeX/2 || x > sizeX/2) { return; }
		if(y < -sizeY/2 || y > sizeY/2) { return; }
		if(z < -sizeZ/2 || z > sizeZ/2) { return; }
		
		H_x = fieldX;
		H_y = fieldY;
		H_z = fieldZ;
}