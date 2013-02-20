// This class represents a simple RF gap.
// For this RF gap we know the E0TL parameter only.

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

/** 
 Tracks the Bunch through the RF gap. This type of model includes the transverse 
 non-linearity.
*/

void BaseRfGap::trackBunch(Bunch* bunch, double frequency, double ampl,
                           double E0TL , double phase)
{
	//we do not use the ampl amplitude of the field in the gap
	//it can be used in other models of the gap
  bunch->compress();
  SyncPart* syncPart   = bunch->getSyncPart();
  double gamma         = syncPart->getGamma();
  double beta          = syncPart->getBeta();
  double mass          = bunch->getMass();
  double charge        = bunch->getCharge();
  double eKin_in       = syncPart->getEnergy();
  double chargeE0TLsin = charge * E0TL * sin(phase);	
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
  double k                = 2.0 * OrbitConst::PI*frequency / OrbitConst::c;
  double phase_time_coeff = k / beta;
  double kr               = k / (gamma * beta);

  //transverse focusing coeff
  double kappa = -charge * E0TL * k /
                 (2.0 * mass * beta_gap * beta_gap * beta_out *
                  gamma_gap * gamma_gap * gamma_out);
  double x, y, r, rp, d_phi;
	//the linear part - implemented in MatrixRFGap
  double d_rp = kappa *sin(phase);
  double I0, I1;
  for(int i = 0, n = bunch->getSize(); i < n; i++)
  {
    x = bunch->x(i);
    y = bunch->y(i);
    r = sqrt(x * x + y * y);
    I0 = bessi0(kr * r);
    I1 = bessi1(kr * r);
    //longitudinal-energy part
    d_phi = - bunch->z(i)*phase_time_coeff;
    bunch->z(i) = bunch->z(i) * beta_out / beta;
    bunch->dE(i) = bunch->dE(i) + charge * E0TL * cos(phase + d_phi) * I0 - delta_eKin;
    //the linear part - implemented in MatrixRFGap
    //bunch->dE(i) = bunch->dE(i) - chargeE0TLsin * d_phi;
    //transverse focusing
    d_rp = kappa *sin(phase + d_phi) * 2.0 * I1 / (kr * r);
		//the transverse non-linear part is removed
		//d_rp = kappa *sin(phase + d_phi);
    bunch->xp(i) = bunch->xp(i) * prime_coeff + d_rp * x;
    bunch->yp(i) = bunch->yp(i) * prime_coeff + d_rp * y;
  }
}
