//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    FieldSourceContainer.cc
//
// CREATED
//    06/27/2008
//
// AUTHOR
//    T. Gorlov
//
// DESCRIPTION
//    The container for instances of the BaseFieldSource class.
//
///////////////////////////////////////////////////////////////////////////
#include "orbit_mpi.hh"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>

#include "FieldSourceContainer.hh"


using namespace OrbitUtils;

FieldSourceContainer::FieldSourceContainer():BaseFieldSource() 
{
}

FieldSourceContainer::~FieldSourceContainer(){
	
	for (int i=0;i<ref.size();i++){ 
		if(ref[i]->getPyWrapper() == NULL){
			delete ref[i];
		} else {
			Py_XDECREF(ref[i]->getPyWrapper());
		}
	}

}

/** Adds the instance of the  ExternalEffects class to the container. */
void FieldSourceContainer::AddFieldSource(BaseFieldSource* fs)	{

	if(fs->getPyWrapper() != NULL){
		Py_INCREF(fs->getPyWrapper());
	}

	ref.push_back(fs);
	
}

/** Adds the instance of the  ExternalEffects class to the container. */
void FieldSourceContainer::getElectricMagneticField(double x, double y, double z, double t, 
				double& E_x, double& E_y, double& E_z,
				double& H_x, double& H_y, double& H_z)	{
	
	double Ex, Ey, Ez, Hx, Hy, Hz;
	
	E_x=0;
	E_y=0;
	E_z=0;
	H_x=0;
	H_y=0;
	H_z=0;
	
	for (int i=0;i<ref.size();i++){
		ref[i]->getElectricMagneticField(x, y, z, t, Ex, Ey, Ez, Hx, Hy, Hz);
		
		E_x+=Ex;
		E_y+=Ey;
		E_z+=Ez;
		
		H_x+=Hx;
		H_y+=Hy;
		H_z+=Hz;
	
	}
	

	
}



























