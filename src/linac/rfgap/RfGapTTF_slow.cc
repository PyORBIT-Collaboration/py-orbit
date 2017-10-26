/**
   This class represents a Parmila type RF gap. It acts on the coordinates 
   of the particle by using the transit time factors. The model includes  
   non-linearity in transverse direction.
	 TTFs (T,S,Tp,Sp) are funcftions of the kappa variable = 2*pi*f/(c*beta).

   The description of the models can be found in
   A. Shishlo, J. Holmes, 
   "Physical Models for Particle Tracking Simulations in the RF Gap", 
   ORNL Tech. Note ORNL/TM-2015/247, June 2015	
   
   The arrival time is defined by dt = w*z/(beta*c)
   instead of dt = w*z/(beta_synch*c) in RfGapTTF.cc file   
 	 
   The transformation coefficients are defined for each particle separately,
   so the name of the class has the word slow in the na  
   
*/

#include <iostream>
#include <cmath>

#include "Bunch.hh"
#include "bessel.hh"
#include "OrbitConst.hh"
#include "RfGapTTF_slow.hh"

using namespace OrbitUtils;


// Constructor
RfGapTTF_slow::RfGapTTF_slow(): CppPyWrapper(NULL)
{
}

// Destructor
RfGapTTF_slow::~RfGapTTF_slow()
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
void RfGapTTF_slow::trackBunch(Bunch* bunch, double frequency, double E0L, double phase,
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
  double p_s = syncPart->getMomentum();
  double p2_s = p_s*p_s;	
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
	delta_eKin = charge*E0L*(ttf_t*cos(phase) - ttf_s*sin(phase));
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
	double prime_coeff = (beta_in * gamma_in)/(beta_out * gamma_out); 
	double x, y, r, I0,I1, phase_in , phase_out, phase_rf, d_rp;
	double sin_phRf, cos_phRf;
	double kappa, trans_coeff;
	double Ekin, dE, p2, p_z2, beta_z, gamma_z, beta_z_out, gamma_z_out;
	double xp,yp;
	for(int i = 0, n = bunch->getSize(); i < n; i++){
		
		dE = bunch->dE(i);
		xp = bunch->xp(i);
		yp = bunch->yp(i);		
	    
 		Ekin = eKin_in + dE;
		p2 = Ekin*(Ekin+2.0*mass);
		p_z2 = p2 - (xp*xp + yp*yp)*p2_s;
		beta_z = sqrt(p_z2)/(Ekin+mass);
		gamma_z = 1.0/sqrt(1.0-beta_z*beta_z);
    kappa = 2.0*OrbitConst::PI*frequency/(OrbitConst::c*beta_z);

		Kr = kappa/gamma_z;	
    kappa_Kr = gamma_z;		
		
    x = bunch->x(i);
    y = bunch->y(i);
    r = sqrt(x * x + y * y);
    I0 = bessi0(Kr * r);
    I1 = bessi1(Kr * r);	

    trans_coeff = charge*E0L/(mass*beta_z*beta_z*gamma_z*gamma_z);
    phase_coeff = charge*E0L*kappa/(mass*beta_z*beta_z*gamma_z*gamma_z*gamma_z);
    
		phase_in = bunch->z(i)*kappa;
		phase_rf = phase - phase_in;	
		sin_phRf = sin(phase_rf);
		cos_phRf = cos(phase_rf);
		ttf_t = Tttf->value(kappa);
		ttf_s = Sttf->value(kappa);
		ttf_tp = Tpttf->value(kappa);
		ttf_sp = Spttf->value(kappa);
				
		//longitudinal-energy part
		bunch->dE(i) = bunch->dE(i)  + charge*E0L*I0*(ttf_t*cos_phRf - ttf_s*sin_phRf) - delta_eKin;	
		phase_out = phase_in + phase_coeff*(I0*(ttf_tp*sin_phRf + ttf_sp*cos_phRf) +
			                     r*kappa_Kr*I1*(ttf_t*sin_phRf + ttf_s*cos_phRf));
		phase_out -= delta_phase;
		
    dE = bunch->dE(i);
		Ekin = eKin_out + dE;
		p2 = Ekin*(Ekin+2.0*mass);
		p_z2 = p2 - (xp*xp + yp*yp)*p2_s;
		beta_z_out = sqrt(p_z2)/(Ekin+mass);
		kappa = 2.0*OrbitConst::PI*frequency/(OrbitConst::c*beta_z_out);
		
		//transverse focusing 
		if(r == 0.){
			d_rp = 0.;
		}
		else{
			d_rp = - trans_coeff*I1*(ttf_t*sin_phRf + ttf_s*cos_phRf)/r;
		}		
		
		bunch->xp(i) = bunch->xp(i)*prime_coeff + d_rp*bunch->x(i);
		bunch->yp(i) = bunch->yp(i)*prime_coeff + d_rp*bunch->y(i);
		
		bunch->z(i) = phase_out/kappa;
	}
}
