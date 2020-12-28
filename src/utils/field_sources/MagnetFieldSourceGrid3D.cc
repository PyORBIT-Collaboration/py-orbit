//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//   MagnetFieldSourceGrid3D.cc
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

#include "orbit_mpi.hh"
#include "BufferStore.hh"

#include <cstdlib>
#include <iostream>
#include <cmath>
#include <cfloat>

#include "ShiftedFieldSource.hh"
#include "MagnetFieldSourceGrid3D.hh"

using namespace OrbitUtils;

MagnetFieldSourceGrid3D::MagnetFieldSourceGrid3D(Grid3D* BxGrid_In, Grid3D* ByGrid_In, Grid3D* BzGrid_In): ShiftedFieldSource()
{
	BxGrid = BxGrid_In;
	ByGrid = ByGrid_In;
	BzGrid = BzGrid_In;
	
	//if symmetry for axis != 0 (e.g. 1) only abs(coord_value) will be used
	//to get information from Grid3D objects
	symmetry_x = 0;
	symmetry_y = 0;
	symmetry_z = 0;
	
	// there are 8 possible quadrant (positive also include "equal to 0" case)
	// #1  |  #0                      I   #5  |  #4
	//-----------  x-y plane for z>0  I  ----------  x-y plane for z<0 
	// #3  |  #2                      I   #7  |  #6
	//------------------------------------------------------------------------
	// #0  -  x > 0   y > 0   z >0       
	// #1  -  x < 0   y > 0   z >0
	// #2  -  x > 0   y < 0   z >0
	// #3  -  x < 0   y < 0   z >0
	// #4  -  x > 0   y > 0   z <0
	// #5  -  x < 0   y > 0   z <0
	// #6  -  x > 0   y < 0   z <0
	// #7  -  x < 0   y < 0   z <0
	//-------------------------------------------------------------------------
	// field_sign_arr[#ind][axis_ind] - multiplier for Bx,By, or Bz (defined by axis_ind)
	// ind_x = 0; ind_y = 0; ind_z = 0;
	// if x < 0 -> ind_x = 1
	// if y < 0 -> ind_y = 1
	// if z < 0 -> ind_z = 1
	// #ind = ind_x + 2*ind_y + 4*ind_z

	field_sign_arr = new int*[8];
	for(int i = 0; i < 8; i++){
		field_sign_arr[i] = new int[3];
		field_sign_arr[i][0] = 1;
		field_sign_arr[i][1] = 1;
		field_sign_arr[i][2] = 1;
	}
	
	//used to correct Bx,By,Bz fields if we have the same 
	//distributions for Grid3D (identical magnets) with different fields
	//It will save the memory for us. 
	field_coeff = 1.;	
	
}

MagnetFieldSourceGrid3D::~MagnetFieldSourceGrid3D()
{

	if(BxGrid->getPyWrapper() == NULL){
		delete BxGrid;
	}
	else {
		Py_XDECREF(BxGrid->getPyWrapper());
	}	

	if(ByGrid->getPyWrapper() == NULL){
		delete ByGrid;
	}
	else {
		Py_XDECREF(ByGrid->getPyWrapper());
	}

	if(BzGrid->getPyWrapper() == NULL){
		delete BzGrid;
	}
	else {
		Py_XDECREF(BzGrid->getPyWrapper());
	}
	
	for(int i = 0; i < 8; i++){
		delete [] field_sign_arr[i];
	}
	delete [] field_sign_arr;
	
}


/** Sets symmetry properties in Grid3D fields along x,y,z axises  */
void MagnetFieldSourceGrid3D::setSymmetry(int symmetry_x, int symmetry_y, int symmetry_z)
{
	this->symmetry_x = symmetry_x;
	this->symmetry_y = symmetry_y;
	this->symmetry_z = symmetry_z;
}

/** Returns symmetry properties in Grid3D fields along x,y,z axises  */
void MagnetFieldSourceGrid3D::getSymmetry(int& symmetry_x, int& symmetry_y, int& symmetry_z){
	symmetry_x = this->symmetry_x;
	symmetry_y = this->symmetry_y;
	symmetry_z = this->symmetry_z;
}              

/** Sets signs for fields in different quadrants that defined by signs of  signX, signY, signZ */
void MagnetFieldSourceGrid3D::setFieldSignsForQuadrants(int signX, int signY, int signZ, int signBx, int signBy, int signBz){
	
    int ind_x = 0;
    int ind_y = 0;
    int ind_z = 0;
    
    if(signX < 0){ ind_x = 1; }
    if(signY < 0){ ind_y = 1; }
    if(signZ < 0){ ind_z = 1; }
    
    int quadrant_ind = ind_x + 2*ind_y + 4*ind_z;
    int* axis_sign_for_fields_arr = field_sign_arr[quadrant_ind];
    
    axis_sign_for_fields_arr[0] = 1;
    axis_sign_for_fields_arr[1] = 1;
    axis_sign_for_fields_arr[2] = 1;
    
    if(signBx < 0) { axis_sign_for_fields_arr[0] = -1; }
    if(signBy < 0) { axis_sign_for_fields_arr[1] = -1; }
    if(signBz < 0) { axis_sign_for_fields_arr[2] = -1; }
}

/** Sets the scaling coefficient for inner fields in Grid3D instances */
void MagnetFieldSourceGrid3D::setFieldCoeff(double field_coeff)
{
	this->field_coeff = field_coeff;
}

/** Returns the scaling coefficient for inner fields in Grid3D instances */
double MagnetFieldSourceGrid3D::getFieldCoeff()
{
	return field_coeff;
}

/** Returns signs for fields in different quadrants that defined by signs of  signX, signY, signZ */
void MagnetFieldSourceGrid3D::getFieldSignsForQuadrants(int signX, int signY, int signZ, int& signBx, int& signBy, int& signBz){
	
    int ind_x = 0;
    int ind_y = 0;
    int ind_z = 0;
    
    if(signX < 0){ ind_x = 1; }
    if(signY < 0){ ind_y = 1; }
    if(signZ < 0){ ind_z = 1; }
    
    int quadrant_ind = ind_x + 2*ind_y + 4*ind_z;
    int* axis_sign_for_fields_arr = field_sign_arr[quadrant_ind];
    
    signBx = axis_sign_for_fields_arr[0];
    signBy = axis_sign_for_fields_arr[1];
    signBz = axis_sign_for_fields_arr[2];
}


/** Returns inner components of the electric and magnetic filds. */
void MagnetFieldSourceGrid3D::getInnerElectricMagneticField(
	  double x, double y, double z, double t, 
		double& E_x, double& E_y, double& E_z,
		double& H_x, double& H_y, double& H_z){

		E_x=0.;	E_y=0.;	E_z=0.;
		H_x=0.;	H_y=0.;	H_z=0.;

    double x_inn = x;
    double y_inn = y;
    double z_inn = z;
    if(symmetry_x > 0){ x_inn = abs(x_inn); }
    if(symmetry_y > 0){ y_inn = abs(y_inn); }
    if(symmetry_z > 0){ z_inn = abs(z_inn); }
    
    if(x_inn > BxGrid->getMaxX()){ return; }
    if(x_inn < BxGrid->getMinX()){ return; }
    
    if(y_inn > ByGrid->getMaxY()){ return; }
    if(y_inn < ByGrid->getMinY()){ return; }
     
    if(z_inn > BzGrid->getMaxZ()){ return; }
    if(z_inn < BzGrid->getMinZ()){ return; }
     
    int ind_x = 0;
    int ind_y = 0;
    int ind_z = 0;
    if(x < 0.){ ind_x = 1; }
    if(y < 0.){ ind_y = 1; }
    if(z < 0.){ ind_z = 1; }
    
    int quadrant_ind = ind_x + 2*ind_y + 4*ind_z;
    int* axis_sign_for_fields_arr = field_sign_arr[quadrant_ind];
    
    H_x = field_coeff*axis_sign_for_fields_arr[0]*BxGrid->getValue(x_inn,y_inn,z_inn);
    H_y = field_coeff*axis_sign_for_fields_arr[1]*ByGrid->getValue(x_inn,y_inn,z_inn);
    H_z = field_coeff*axis_sign_for_fields_arr[2]*BzGrid->getValue(x_inn,y_inn,z_inn);
}