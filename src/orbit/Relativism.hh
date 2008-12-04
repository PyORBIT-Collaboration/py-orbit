#ifndef RELATIVISM_HH_
#define RELATIVISM_HH_

//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    LorentzTransformationEM.hh
//
// AUTHOR
//    T. Gorlov
//
// CREATED
//    06/28/2005
//
// DESCRIPTION
//    This class provides Lorentz transformations for the electromagnetic field
//    from the laboratory frame to the particle rest frame.
//    mass - mass of the particle in GeV
//    px,py,pz - momentum of the particle the lab frame in GeV/c
//    E_x,E_y,E_z - components of the electric field V/m (parameters are replaced in place) 
//    B_x,B_y,B_z - components of the magnetic field [T] (parameters are replaced in place)   
//
///////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////
//
// include files
//
/////////////////////////////////////////////////////////////////////////////


class  LorentzTransformationEM
{
public:
  static  void 	transform(double mass,
		                      double px, double py, double pz,
													double& E_x, double& E_y, double& E_z,
													double& B_x, double& B_y, double& B_z);	
};

///////////////////////////////////////////////////////////////////////////
//
// END OF FILE
//
///////////////////////////////////////////////////////////////////////////



#endif /*RELATIVISM_HH_*/


