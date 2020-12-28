//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//   MagnetFieldSourceGrid3D.hh
//
// AUTHOR
//    A. Shishlo
//
// CREATED
//    12/23/2020
//
// DESCRIPTION
//    A class is an implementation of ShiftedFieldSource class.
//    This class contains 3 Grid3D instances with Bx,By,Bz magnetic fields in [T]
//    The coordTransform Matrix is (4x4) matrix object that translates the
//    spacial coordinates from external coordinates to the inner coordinates
//    of each Grid3D. The last raw in coordTransform Matrix is [0,0,0,1], and
//    The last column is a [a,b,c,1] where [-a,-b,-c] coordinates of the Grid3D
//    origin in the external coordinates system.
//    The vector (Bx,By,Bz) transforms back to the external coordinate
//    system by coordTransformBack (3x3) matrix.
//
///////////////////////////////////////////////////////////////////////////

#ifndef MAGNET_FIELD_SOURCE_GRID3D_H
#define MAGNET_FIELD_SOURCE_GRID3D_H

#include "Grid3D.hh"
#include "ShiftedFieldSource.hh"

namespace OrbitUtils{

	/** A class implements of BaseFiledSource class with magnetic fields in Grid3D.*/

	class MagnetFieldSourceGrid3D : public ShiftedFieldSource
	{
		public:

			/** Constructor. */
			MagnetFieldSourceGrid3D(Grid3D* BxGrid, Grid3D* ByGrid, Grid3D* BzGrid);

			/** Destructor */
			~MagnetFieldSourceGrid3D();

			/** Sets symmetry properties in Grid3D fields along x,y,z axises  */
			void setSymmetry(int symmetry_x, int symmetry_y, int symmetry_z);

			/** Returns symmetry properties in Grid3D fields along x,y,z axises  */
			void getSymmetry(int& symmetry_x, int& symmetry_y, int& symmetry_z);

			/** Sets signs for fields in different quadrants that defined by signs of  signX, signY, signZ */
			void setFieldSignsForQuadrants(int signX, int signY, int signZ, int signBx, int signBy, int signBz);

			/** Returnss signs for fields in different quadrants that defined by signs of  signX, signY, signZ */
			void getFieldSignsForQuadrants(int signX, int signY, int signZ, int& signBx, int& signBy, int& signBz);

			/** Sets the scaling coefficient for inner fields in Grid3D instances */
			void setFieldCoeff(double field_coeff);

			/** Returns the scaling coefficient for inner fields in Grid3D instances */
			double getFieldCoeff();

			/** Returns inner components of the electric and magnetic filds. */
			virtual void getInnerElectricMagneticField(
			     double x, double y, double z, double t,
				   double& E_x, double& E_y, double& E_z,
				   double& H_x, double& H_y, double& H_z);

	private:
		
		Grid3D* BxGrid;
		Grid3D* ByGrid;
		Grid3D* BzGrid;
		
		int symmetry_x;
		int symmetry_y;
		int symmetry_z;
		
		int** field_sign_arr;
		
		//used to correct Bx,By,Bz fields if we have the same 
		//distributions for Grid3D (identical magnets) with different fields
		//It will save the memory for us. 
		double field_coeff;
		
	};	
};

///////////////////////////////////////////////////////////////////////////
//
// END OF FILE
//
///////////////////////////////////////////////////////////////////////////

#endif
