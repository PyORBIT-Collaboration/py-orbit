/**
This class represents a simplified RF gap. It acts on the coordinates like a transport matrix. 
There are no nonlinear effects. It should be analog of RF Gap of XAL online model or Trace3D.
For this RF gap we know the E0TL, frequency, and phase only.
*/

#include "MatrixRfGap.hh"
#include "ParticleMacroSize.hh"

#include <iostream>
#include <cmath>

#include "Bunch.hh"
#include "bessel.hh"
#include "OrbitConst.hh"

using namespace OrbitUtils;

// Constructor
MatrixRfGap::MatrixRfGap(): CppPyWrapper(NULL)
{
}

// Destructor
MatrixRfGap::~MatrixRfGap()
{
}

/** Tracks the Bunch trough the RF gap. */	
void MatrixRfGap::trackBunch(Bunch* bunch, double frequency, double E0TL, double phase){
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
	//wave momentum
	double k = 2.0*OrbitConst::PI*frequency/OrbitConst::c;
	double phase_time_coeff = k/beta;
	//transverse focusing coeff
	double cappa = - charge*E0TL*k/(2.0*mass*beta_gap*beta_gap*beta_out*gamma_gap*gamma_gap*gamma_out);
	double d_rp = cappa*sin(phase);
	for(int i = 0, n = bunch->getSize(); i < n; i++){
		//longitudinal-energy part
		bunch->dE(i) =bunch->dE(i)  - chargeE0TLsin*phase_time_coeff*bunch->z(i);		
		bunch->z(i) = bunch->z(i)*beta_out/beta;
		//transverse focusing 
		bunch->xp(i) = bunch->xp(i)*prime_coeff + d_rp*bunch->x(i);
		bunch->yp(i) = bunch->yp(i)*prime_coeff + d_rp*bunch->y(i);		
	}	

}


