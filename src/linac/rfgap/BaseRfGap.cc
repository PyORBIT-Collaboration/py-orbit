//This class repersents a simple RF gap. For this RF gap we know the E0TL parameter only.

#include "BaseRfGap.hh"
#include "ParticleMacroSize.hh"

#include <iostream>
#include <cmath>

#include "Bunch.hh"
#include "bessel.hh"
#include "OrbitConst.hh"

using namespace OrbitUtils;

// Constructor
BaseRfGap::BaseRfGap(): CppPyWrapper(NULL)
{
}

// Destructor
BaseRfGap::~BaseRfGap()
{
}

/** Tracks the Bunch trough the RF gap. */	
void BaseRfGap::trackBunch(Bunch* bunch, double frequency, double E0TL, double phase){
	bunch->compress();
	SyncPart* syncPart = bunch->getSyncPart();
	double gamma = syncPart->getGamma();
	double beta = syncPart->getBeta();
	double mass = bunch->getMass();
	double charge = bunch->getCharge();
	double eKin_in = syncPart->getEnergy();
	double chargeE0TLsin = charge*E0TL*sin(phase);	
	double delta_eKin = charge*E0TL*cos(phase);
	//calculate params in the middle of the gap
	syncPart->setMomentum(syncPart->energyToMomentum(eKin_in + delta_eKin/2.0));
	double gamma_gap = syncPart->getGamma();
	double beta_gap = syncPart->getBeta();	
	//now move to the end of the gap
	double eKin_out = eKin_in + delta_eKin;
	syncPart->setMomentum(syncPart->energyToMomentum(eKin_out));	
	//the base RF gap is simple - no phase correction. The delta time in seconds
	double delta_time = 0.;	
	syncPart->setTime(syncPart->getTime() + delta_time);	
	double gamma_out =	syncPart->getGamma();
	double beta_out = syncPart->getBeta();	
	double prime_coeff = (beta*gamma)/(beta_out*gamma_out);
	//std::cout <<" debug delta_eKin="<< delta_eKin*1.0e+3 <<"    prime_coeff ="<<prime_coeff <<std::endl;
	//wave momentum
	double k = 2.0*OrbitConst::PI*frequency/OrbitConst::c;
	double phase_time_coeff = k/beta;
	double kr = k/(gamma*beta);
	//transverse focusing coeff
	double cappa = - charge*E0TL*k/(2.0*mass*beta_gap*beta_gap*beta_out*gamma_gap*gamma_gap*gamma_out);
	//std::cout <<" debug E0TL="<< E0TL << "  kr="<<kr<<"  cappa="<<cappa<< std::endl;	
	double x,y,r,rp, d_phi;
	double d_rp = cappa*sin(phase);
	//std::cout<<"debug BaseRfGap::trackBunch eKin_in=     "<< eKin_in*1.e+3 <<"     xp_coeff=      "<<prime_coeff<<"      x_coeff=     "<<d_rp<<"   phase="<< ((phase*180/OrbitConst::PI)+180.)<<std::endl;
  //std::cout<<"debug BaseRfGap::trackBunch eKin_in=     "<< eKin_in*1.e+3 <<"  chargeE0TLsin="<<chargeE0TLsin<<" phase_time_coeff="<<phase_time_coeff<<std::endl;
	double I0, I1;
	for(int i = 0, n = bunch->getSize(); i < n; i++){
		x = bunch->x(i);
		y = bunch->y(i);
		r = sqrt(x*x + y*y);
    I0 = bessi0(kr*r);
		I1 = bessi1(kr*r);
		//longitudinal-energy part
		d_phi = - bunch->z(i)*phase_time_coeff;
		bunch->z(i) = bunch->z(i)*beta_out/beta;
		//bunch->dE(i) = bunch->dE(i) + charge*E0TL*cos(phase + d_phi)*I0 - delta_eKin;
		//bunch->dE(i) = (bunch->dE(i)/beta - chargeE0TLsin*d_phi/beta_gap)*beta_out;
		bunch->dE(i) =bunch->dE(i)  - chargeE0TLsin*d_phi;
		//transverse focusing 
		//d_rp = cappa*sin(phase + d_phi)*2.0*I1/kr;		
		bunch->xp(i) = bunch->xp(i)*prime_coeff + d_rp*x;
		bunch->yp(i) = bunch->yp(i)*prime_coeff + d_rp*y;
		//bunch->xp(i) = bunch->xp(i)*prime_coeff;
		//bunch->yp(i) = bunch->yp(i)*prime_coeff;			
	}	

}


