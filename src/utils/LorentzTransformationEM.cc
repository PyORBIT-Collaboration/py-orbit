//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    LorentzTransformationEM.cc
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
//
///////////////////////////////////////////////////////////////////////////

#include "LorentzTransformationEM.hh"
#include "OrbitConst.hh"
#include <cmath>

using namespace OrbitUtils;

void 	LorentzTransformationEM::transform(double mass, 
	                                      double px, double py, double pz,
	                                      double& E_x, double& E_y, double& E_z,
																				double& B_x, double& B_y, double& B_z)
{
	
	double mp2= mass*mass;
	double p2=px*px+py*py+pz*pz;
	double E=sqrt(mp2+p2);
	double gamma=E/mass;
	double _b=(OrbitConst::c)/E;
	double vx=_b*px;
	double vy=_b*py;
	double vz=_b*pz;
	
	double c2=(OrbitConst::c)*(OrbitConst::c);
	double k=gamma*gamma/(c2*(1+gamma));
	

	double EK=(E_x*vx + E_y*vy + E_z*vz)*k;
	double BK=(B_x*vx + B_y*vy + B_z*vz)*k;
	
	double BV1=E_x+B_z*vy-B_y*vz;
	double BV2=E_y+B_x*vz-B_z*vx;
	double BV3=E_z+B_y*vx-B_x*vy;
	
	double EV1=B_x+(E_y*vz-E_z*vy)/c2;
	double EV2=B_y+(E_z*vx-E_x*vz)/c2;
	double EV3=B_z+(E_x*vy-E_y*vx)/c2;
	
	
	// Electric field in the frame of particle	
	E_x=gamma*BV1-vx*EK;
	E_y=gamma*BV2-vy*EK;
	E_z=gamma*BV3-vz*EK;
	
	// Magnetic field in the frame of particle	(ignored in case of useless for faster computations)
	B_x=gamma*EV1-vx*BK;
	B_y=gamma*EV2-vy*BK;
	B_z=gamma*EV3-vz*BK;	
}


void 	LorentzTransformationEM::complex_transform(double mass, double px, double py, double pz,
        tcomplex& E_x, tcomplex& E_y, tcomplex& E_z,
        tcomplex& B_x, tcomplex& B_y, tcomplex& B_z)
{
	
	double mp2= mass*mass;
	double p2=px*px+py*py+pz*pz;
	double E=sqrt(mp2+p2);
	double gamma=E/mass;
	double _b=(OrbitConst::c)/E;
	double vx=_b*px;
	double vy=_b*py;
	double vz=_b*pz;
	
	double c2=pow(OrbitConst::c,2);
	double k=gamma*gamma/(c2*(1+gamma));
	

tcomplex EK=(E_x*vx + E_y*vy + E_z*vz)*k;
tcomplex BK=(B_x*vx + B_y*vy + B_z*vz)*k;
	
tcomplex BV1=E_x+B_z*vy-B_y*vz;
tcomplex BV2=E_y+B_x*vz-B_z*vx;
tcomplex BV3=E_z+B_y*vx-B_x*vy;
	
tcomplex EV1=B_x+(E_y*vz-E_z*vy)/c2;
tcomplex EV2=B_y+(E_z*vx-E_x*vz)/c2;
tcomplex EV3=B_z+(E_x*vy-E_y*vx)/c2;
	
	
	// Electric field in the frame of particle	
	E_x=gamma*BV1-vx*EK;
	E_y=gamma*BV2-vy*EK;
	E_z=gamma*BV3-vz*EK;
	
	// Magnetic field in the frame of particle	
	B_x=gamma*EV1-vx*BK;
	B_y=gamma*EV2-vy*BK;
	B_z=gamma*EV3-vz*BK;	
}




void 	LorentzTransformationEM::complex_electric_transform(double mass, double px, double py, double pz,
        tcomplex& E_x, tcomplex& E_y, tcomplex& E_z,
        tcomplex B_x, tcomplex B_y, tcomplex B_z)
{
	double mp2= mass*mass;
	double p2=px*px+py*py+pz*pz;
	double E=sqrt(mp2+p2);
	double gamma=E/mass;
	double _b=(OrbitConst::c)/E;
	double vx=_b*px;
	double vy=_b*py;
	double vz=_b*pz;
	
	double c2=pow(OrbitConst::c,2);
	double k=gamma*gamma/(c2*(1+gamma));
	

tcomplex EK=(E_x*vx + E_y*vy + E_z*vz)*k;
	
tcomplex BV1=E_x+B_z*vy-B_y*vz;
tcomplex BV2=E_y+B_x*vz-B_z*vx;
tcomplex BV3=E_z+B_y*vx-B_x*vy;
	
		
	// Electric field in the frame of particle	
	E_x=gamma*BV1-vx*EK;
	E_y=gamma*BV2-vy*EK;
	E_z=gamma*BV3-vz*EK;
	
	

}



























