/**
This class represents a simple RF gap.
For this RF gap model we need the E0TL parameter only.

The description of the models can be found in
A. Shishlo, J. Holmes, 
"Physical Models for Particle Tracking Simulations in the RF Gap", 
ORNL Tech. Note ORNL/TM-2015/247, June 2015

The arrival time is defined by dt = w*z/(beta*c)
instead of dt = w*z/(beta_synch*c) in BaseRfGap.cc file

The transformation coefficients are defined for each particle separately,
so the name of the class has the word slow in the name.
   
*/

#include "BaseRfGap_slow.hh"
#include "ParticleMacroSize.hh"

#include <iostream>
#include <cmath>

#include "Bunch.hh"
#include "bessel.hh"
#include "OrbitConst.hh"

using namespace OrbitUtils;

// Constructor
BaseRfGap_slow::BaseRfGap_slow(): CppPyWrapper(NULL)
{
}

// Destructor
BaseRfGap_slow::~BaseRfGap_slow()
{
}

/** 
 Tracks the Bunch through the RF gap. This type of model includes the transverse 
 non-linearity.
*/

void BaseRfGap_slow::trackBunch(Bunch* bunch, double frequency, double E0TL, double phase)
{
	// E0TL is a maximal energy gain in the gap. It is in GeV.
	// RF frequency is in Hz
	// RF phase in radians
  bunch->compress();
  SyncPart* syncPart   = bunch->getSyncPart();
  double gamma         = syncPart->getGamma();
  double beta          = syncPart->getBeta();
  double mass          = bunch->getMass();
  double charge        = bunch->getCharge();
  double eKin_in       = syncPart->getEnergy();
  double p_s           = syncPart->getMomentum();
  double p2_s          = p_s*p_s;
  //double chargeE0TLsin = charge * E0TL * sin(phase);	
  double delta_eKin    = charge * E0TL * cos(phase);

  //calculate params in the middle of the gap
  syncPart->setMomentum(syncPart->energyToMomentum(eKin_in + delta_eKin/2.0));
  double gamma_gap = syncPart->getGamma();
  double beta_gap  = syncPart->getBeta();

  //now move to the end of the gap
  double eKin_out  = eKin_in + delta_eKin;
  syncPart->setMomentum(syncPart->energyToMomentum(eKin_out));

  //the base RF gap is simple - no phase correction. The delta time in seconds
  double delta_time  = 0.;
  syncPart->setTime(syncPart->getTime() + delta_time);
  double gamma_out   = syncPart->getGamma();
  double beta_out    = syncPart->getBeta();
  double prime_coeff = (beta * gamma)/(beta_out * gamma_out);

  //wave momentum
  double k = 2.0 * OrbitConst::PI*frequency / OrbitConst::c;
  double x, y, r, rp, d_phi;
	//the linear part - implemented in MatrixRFGap
  double kappa, d_rp, kr;
  double I0, I1;
  double dE, xp, yp, Ekin, p2, p_z2, beta_z, gamma_z,beta_z_out,gamma_z_out;
  for(int i = 0, n = bunch->getSize(); i < n; i++)
  {
  	
		dE = bunch->dE(i);
		xp = bunch->xp(i);
		yp = bunch->yp(i);
		
		Ekin = eKin_in + dE;
		p2 = Ekin*(Ekin+2.0*mass);
		p_z2 = p2 - (xp*xp + yp*yp)*p2_s;
		beta_z = sqrt(p_z2)/(Ekin+mass);
		gamma_z = 1.0/sqrt(1.0-beta_z*beta_z);

		kr = k / (gamma_z * beta_z);
		
    x = bunch->x(i);
    y = bunch->y(i);
    r = sqrt(x * x + y * y);
    I0 = bessi0(kr * r);
    I1 = bessi1(kr * r);
    //longitudinal-energy part
    d_phi = - bunch->z(i)*k/beta_z;
    bunch->z(i) = bunch->z(i)*beta_out/beta_z;
    bunch->dE(i) = bunch->dE(i) + charge * E0TL * cos(phase + d_phi) * I0 - delta_eKin;
    //the linear part - implemented in MatrixRFGap
    //bunch->dE(i) = bunch->dE(i) - chargeE0TLsin * d_phi;
    dE = bunch->dE(i);
		Ekin = eKin_out + dE;
		p2 = Ekin*(Ekin+2.0*mass);
		p_z2 = p2 - (xp*xp + yp*yp)*p2_s;
		beta_z_out = sqrt(p_z2)/(Ekin+mass);
		gamma_z_out = 1.0/sqrt(1.0-beta_z_out*beta_z_out);
		kappa = -charge * E0TL * k /(2.0 * mass * beta_z * beta_z * beta_z_out * gamma_z * gamma_z * gamma_z_out);
    //transverse focusing
		if(r == 0.){
			d_rp = 0.;
		}
		else{    
			d_rp = kappa *sin(phase + d_phi) * 2.0 * I1 / (kr * r);
		}
		//the transverse non-linear part is removed
		//d_rp = kappa *sin(phase + d_phi);
    bunch->xp(i) = xp * prime_coeff + d_rp * x;
    bunch->yp(i) = yp * prime_coeff + d_rp * y;
  }
}
