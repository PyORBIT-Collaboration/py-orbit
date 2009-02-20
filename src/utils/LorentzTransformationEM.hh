#ifndef LORENTZ_TRANSFORMATION_EM_
#define LORENTZ_TRANSFORMATION_EM_

//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    LorentzTransformationEM.hh
//
// AUTHOR
//    T. Gorlov
//
// CREATED
//    06/28/2008
//
// DESCRIPTION
//    This class provides Lorentz transformations for the electromagnetic field
//    from the laboratory frame to the particle rest frame.
//    mass - mass of the particle in GeV
//    px,py,pz - momentum of the particle the lab frame in GeV/c
//    E_x,E_y,E_z - components of the electric field V/m (parameters are replaced in place) 
//    B_x,B_y,B_z - components of the magnetic field [T] (parameters are replaced in place)   
//
//    OrbitConst::c in [m/sec]
///////////////////////////////////////////////////////////////////////////
#include "tcomplex.hh"

class  LorentzTransformationEM
{
public:
  static  void 	transform(double mass,
		                      double px, double py, double pz,
													double& E_x, double& E_y, double& E_z,
													double& B_x, double& B_y, double& B_z);
  
  static  void 	complex_transform(double mass,
  		                      double px, double py, double pz,
  		                    tcomplex& E_x, tcomplex& E_y, tcomplex& E_z,
  		                    tcomplex& B_x, tcomplex& B_y, tcomplex& B_z);
  
};

#endif /*LORENTZ_TRANSFORMATION_EM_*/


