//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//   QuadFieldSource.hh
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

#ifndef QUAD_FIELD_SOURCE_H
#define QUAD_FIELD_SOURCE_H

#include "Grid3D.hh"
#include "ShiftedFieldSource.hh"

namespace OrbitUtils{
	
	/** A class implements of BaseFiledSource class with magnetic fields of hard-edge quad */
	
	class QuadFieldSource : public ShiftedFieldSource
	{
		public:
		
			/** Constructor. */
			QuadFieldSource();
			
			/** Destructor */
			~QuadFieldSource();
			
			/** Sets the length of the quad */
			void setLength(double length);
			
			/** Returns the length of the quad */
			double getLength();
			
			/** Sets the gradient of the quad in [T/m] */
			void setGradient(double gradient);
			
			/** Returns the gradient of the quad in [T/m] */
			double getGradient();
			
			/** Returns inner components of the electric and magnetic filds. */
			virtual void getInnerElectricMagneticField(
			     double x, double y, double z, double t, 
				   double& E_x, double& E_y, double& E_z,
				   double& H_x, double& H_y, double& H_z);

	private:
		
		//length is [m]
		double length;
		
		// field gardient in [T/m]
		double gradient;
		
	};	
};

///////////////////////////////////////////////////////////////////////////
//
// END OF FILE
//
///////////////////////////////////////////////////////////////////////////

#endif
