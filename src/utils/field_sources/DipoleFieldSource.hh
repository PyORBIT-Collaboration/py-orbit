//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//   DipoleFieldSource.hh
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
#ifndef DIPOLE_FIELD_SOURCE_H
#define DIPOLE_FIELD_SOURCE_H

#include "Grid3D.hh"
#include "ShiftedFieldSource.hh"

namespace OrbitUtils{
	
	/** A class implements of BaseFiledSource class with magnetic fields of dipiole */
	
	class DipoleFieldSource : public ShiftedFieldSource
	{
		public:
		
			/** Constructor. */
			DipoleFieldSource();
			
			/** Destructor */
			~DipoleFieldSource();
			
			/** Sets the sizes of the dipole field */
			void setSizes(double sizeX, double sizeY, double sizeZ);
			
			
			/** Returns the sizes of the dipole field */
			void getSizes(double& sizeX, double& sizeY, double& sizeZ);
			
			
			/** Sets the fields of the dipole in [T] */
			void setFields(double fieldX, double fieldY, double fieldZ);
			
			
			/** Returns the fields of the dipole in [T] */
			void getFields(double& fieldX, double& fieldY, double& fieldZ);
			
			/** Returns inner components of the electric and magnetic filds. */
			virtual void getInnerElectricMagneticField(
			     double x, double y, double z, double t, 
				   double& E_x, double& E_y, double& E_z,
				   double& H_x, double& H_y, double& H_z);

	private:
		
		//sizes in [m]
		double sizeX;
		double sizeY;
		double sizeZ;
		
		// fields in [T]
		double fieldX;
		double fieldY;
		double fieldZ;
		
	};	
};

///////////////////////////////////////////////////////////////////////////
//
// END OF FILE
//
///////////////////////////////////////////////////////////////////////////

#endif
