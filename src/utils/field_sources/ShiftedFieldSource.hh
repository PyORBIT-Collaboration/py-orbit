//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//   ShiftedFieldSource.hh
//
// AUTHOR
//    A. Shishlo
//
// CREATED
//    12/23/2020
//
// DESCRIPTION
//    A class is an implementation of BaseFiledSource class.
//    This class describe the instance of the Field Source (FS) class which is
//    shifted and rotated with respect of the center (0,0,0) original 
//    coordinates.
//    The coordTransform Matrix is (4x4) matrix object that translates the 
//    spacial coordinates from external coordinates to the inner coordinates
//    of this FS. The last raw in coordTransform Matrix is [0,0,0,1], and
//    The last column is a [a,b,c,1] where [-a,-b,-c] coordinates of this FS
//    origin in the external coordinates system.
//    The vectors (Ex,Ey,Ez) and (Bx,By,Bz) are transformed back to 
//    the external coordinate system by coordTransformBack (3x3) matrix.
//
///////////////////////////////////////////////////////////////////////////

#ifndef SHIFTED_FIELD_SOURCE_H
#define SHIFTED_FIELD_SOURCE_H

#include "orbit_mpi.hh"

#include <cstdlib>
#include <iostream>
#include <cmath>
#include <cfloat>

#include "BaseFieldSource.hh"
#include "Matrix.hh"
#include "PhaseVector.hh"

namespace OrbitUtils{
	
	/** A class implements of BaseFiledSource class for shifted coordinate system.*/
	
	class ShiftedFieldSource : public BaseFieldSource
	{
		public:
		
			/** Constructor. */
			ShiftedFieldSource();
			
			/** Destructor */
			virtual ~ShiftedFieldSource();
			
			/** Returns components of the electric and magnetic fields. */
			void getElectricMagneticField(
			     double x, double y, double z, double t, 
				   double& E_x, double& E_y, double& E_z,
				   double& H_x, double& H_y, double& H_z);
			
			/** 
			    Returns components of the electric and magnetic fields in 
			    the inner coordinate system. These components will be 
			    transformed to the external system in getElectricMagneticField
			    method.
			    This is a virtual method, and subclasses should implement this 
			    abstarct method.
			 */
			virtual void getInnerElectricMagneticField(
			     double x, double y, double z, double t, 
				   double& E_x, double& E_y, double& E_z,
				   double& H_x, double& H_y, double& H_z);			
			
			/** 
			     Returns coordinates transformation matrix 4x4.
			     It includes rotation 3x3 matrix and origin shift.
			 */
			Matrix* getCoordsTransformMatrix();
			
			/** 
			      Sets coordinates transformation matrix 4x4.
			      It includes rotation 3x3 matrix and origin shift.
			 */
			void setCoordsTransformMatrix(Matrix* coordTransformM4x4);
			
			
		protected:
			
			//Matrix for transformation from external to shifted system
			Matrix* coordTransformM4x4;
			
			//-----------------------------------
			//coordinate vectors 
			//-----------------------------------
			//coordVectExt - from external system
			PhaseVector* coordVectExt;
			//coordVectInn - inside shifted system
			PhaseVector* coordVectInn;
			
			//Fileds (E,H) vectors
			PhaseVector* fieldVectExt;
			PhaseVector* fieldVectInn;
			
			
			//Matrix for transformation E and B from shifted to external system
			Matrix* coordTransformM3x3;	
			
	};
};

///////////////////////////////////////////////////////////////////////////
//
// END OF FILE
//
///////////////////////////////////////////////////////////////////////////

#endif