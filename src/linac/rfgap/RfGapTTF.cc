/**
   This class represents a Parmila type RF gap. It acts on the coordinates 
   of the particle by using the transit time factors. The model includes  
   non-linearity in transverse direction.
	 TTFs (T,S,Tp,Sp) are funcftions of the kappa variable = 2*pi*f/(c*beta).

   The description of the models can be found in
   A. Shishlo, J. Holmes, 
   "Physical Models for Particle Tracking Simulations in the RF Gap", 
   ORNL Tech. Note ORNL/TM-2015/247, June 2015	 	 
*/

#include <iostream>
#include <cmath>

#include "Bunch.hh"
#include "bessel.hh"
#include "OrbitConst.hh"
#include "RfGapTTF.hh"

using namespace OrbitUtils;


// Constructor
RfGapTTF::RfGapTTF(): CppPyWrapper(NULL)
{
}

// Destructor
RfGapTTF::~RfGapTTF()
{
}

/** Tracks the Bunch trough the RF gap as whole using the TTF T and S. 
    There are a lot of possibilities to modify this algorithm, but for now
		we believe that all changes will give a small difference in results if
		the energy spread in the bunch and the energy gain in the gap are small
		relative to the initial energy.
		A curious user could try to redefine how we calculate the energy gain,
		do we use the same TTF for all particles, recursive procedures etc.
		If you find something interesting, please let us know.
*/	
void RfGapTTF::trackBunch(Bunch* bunch, double frequency, double E0L, double phase,
	OrbitUtils::Polynomial* Tttf,
	OrbitUtils::Polynomial* Sttf,
	OrbitUtils::Polynomial* Tpttf,
	OrbitUtils::Polynomial* Spttf){
	// E0L, energy, and mass of particles in GeV
	// RF frequency is in Hz
	// RF phase in radians		
	bunch->compress();
	SyncPart* syncPart = bunch->getSyncPart();
	double gamma_in = syncPart->getGamma();
	double beta_in = syncPart->getBeta();
	double mass = bunch->getMass();
	double charge = bunch->getCharge();
	double eKin_in = syncPart->getEnergy();
	double kappa_in = 2.0*OrbitConst::PI*frequency/(OrbitConst::c*beta_in);
	double delta_eKin = charge*E0L*(Tttf->value(kappa_in)*cos(phase) - Sttf->value(kappa_in)*sin(phase));
	//calculate params in the middle of the gap
	syncPart->setMomentum(syncPart->energyToMomentum(eKin_in + delta_eKin/2.0));
	double gamma_gap = syncPart->getGamma();
	double beta_gap = syncPart->getBeta();	
	double kappa_gap = 2.0*OrbitConst::PI*frequency/(OrbitConst::c*beta_gap);
	// T,S,Tp,Sp for kappa = kappa_gap, we assume a small energy spread
  double ttf_t = Tttf->value(kappa_gap);
  double ttf_s = Sttf->value(kappa_gap);
  double ttf_tp = Tpttf->value(kappa_gap);
  double ttf_sp = Spttf->value(kappa_gap);		
	//the TTF RF gap has the phase correction to simplectic tracking. The delta time in seconds
	double delta_phase = charge*E0L*kappa_gap*(ttf_tp*sin(phase) + ttf_sp*cos(phase))
	                     /(mass*beta_gap*beta_gap*gamma_gap*gamma_gap*gamma_gap);							 
	double delta_time = delta_phase/(2.0*OrbitConst::PI*frequency);	
	syncPart->setTime(syncPart->getTime() + delta_time);
	//now move to the end of the gap
	double eKin_out = eKin_in + delta_eKin;
	syncPart->setMomentum(syncPart->energyToMomentum(eKin_out));	
	double gamma_out = syncPart->getGamma();
	double beta_out = syncPart->getBeta();	
	double kappa_out = 2.0*OrbitConst::PI*frequency/(OrbitConst::c*beta_out);
	//(wave momentum)/beta
	double Kr = kappa_gap/gamma_gap;
	double kappa_Kr = kappa_gap/Kr;
	//phase coeff
	double phase_coeff = charge*E0L*kappa_gap/(mass*beta_gap*beta_gap*gamma_gap*gamma_gap*gamma_gap);
  //transverse coeff
	double trans_coeff = charge*E0L/(mass*beta_gap*beta_gap*gamma_gap*gamma_gap);
	double prime_coeff = (beta_in * gamma_in)/(beta_out * gamma_out); 
	double x, y, r, I0,I1, phase_in , phase_out, phase_rf, d_rp;
	double sin_phRf, cos_phRf;
	for(int i = 0, n = bunch->getSize(); i < n; i++){
    x = bunch->x(i);
    y = bunch->y(i);
    r = sqrt(x * x + y * y);
    I0 = bessi0(Kr * r);
    I1 = bessi1(Kr * r);		
		phase_in = bunch->z(i)*kappa_in;
		phase_rf = phase - phase_in;	
		sin_phRf = sin(phase_rf);
		cos_phRf = cos(phase_rf);
		//longitudinal-energy part
		bunch->dE(i) = bunch->dE(i)  + charge*E0L*I0*(ttf_t*cos_phRf - ttf_s*sin_phRf) - delta_eKin;	
		phase_out = phase_in + phase_coeff*(I0*(ttf_tp*sin_phRf + ttf_sp*cos_phRf) +
			                     r*kappa_Kr*I1*(ttf_t*sin_phRf + ttf_s*cos_phRf));
		bunch->z(i) = phase_out/kappa_out;
		//transverse focusing 
		if(r == 0.){
			d_rp = 0.;
		}
		else{
			d_rp = - trans_coeff*I1*(ttf_t*sin_phRf + ttf_s*cos_phRf)/r;
		}
		bunch->xp(i) = bunch->xp(i)*prime_coeff + d_rp*bunch->x(i);
		bunch->yp(i) = bunch->yp(i)*prime_coeff + d_rp*bunch->y(i);		
	}
}
