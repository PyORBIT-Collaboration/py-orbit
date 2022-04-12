//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//   LoopFieldSource.hh
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

#ifndef CURRENT_LOOP_FIELD_SOURCE_H
#define CURRENT_LOOP_FIELD_SOURCE_H

#include "ShiftedFieldSource.hh"

namespace OrbitUtils{
	
	/** A class implements of BaseFiledSource class with magnetic fields of hard-edge quad */
	class LoopFieldSource : public ShiftedFieldSource
	{
	public:
		
		/** Constructor. */
		LoopFieldSource();
		
		/** Destructor */
		~LoopFieldSource();
		
		/** Sets the radius of the loop */
		void setRadius(double radius);
		
		/** Returns the radius of the loop */
		double getRadius();
		
		/** Sets the current in [A] */
		void setCurrent(double current);
		
		/** Returns the  current in [A]*/
		double getCurrent();
		
		/** Sets the max distance in z direction */
		void setMaxZ(double z_max);
		
		/** Returns the max distance in z direction */
		double getMaxZ();
		
		/** Sets the max distance in x-y plane */
		void setMaxR(double r_max);
		
		/** Returns the max distance in x-y plane */
		double getMaxR();
		
		
		/** Returns inner components of the electric and magnetic fields. */
		virtual void getInnerElectricMagneticField(
			double x, double y, double z, double t, 
			double& E_x, double& E_y, double& E_z,
			double& H_x, double& H_y, double& H_z);

	private:
		
		//radius of the loop [m]
		double radius;
		
		//current [A]
		double current;
		
		//maximal distance along z-axis
		double z_max;
		
		//maximal distance from the center in x-y plane
		double r_max;
		
	};	
};

///////////////////////////////////////////////////////////////////////////
//
// END OF FILE
//
///////////////////////////////////////////////////////////////////////////

#endif