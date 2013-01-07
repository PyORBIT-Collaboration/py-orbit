#include "Frequency_Cav.hh"
#include "ParticleMacroSize.hh"

#include <iostream>
#include <cmath>

#include "Bunch.hh"
#include "bessel.hh"
#include "OrbitConst.hh"

using namespace OrbitUtils;

// Constructor
Frequency_Cav::Frequency_Cav(double RFFreq,
                             double RFE0TL,
                             double RFPhase): CppPyWrapper(NULL)
{
  _RFFreq  = RFFreq;
  _RFE0TL  = RFE0TL;
  _RFPhase = RFPhase;
}

// Destructor
Frequency_Cav::~Frequency_Cav()
{
}

void Frequency_Cav::setRFFreq(double RFFreq)
{
  _RFFreq = RFFreq;
}

double Frequency_Cav::getRFFreq()
{
  return _RFFreq;
}

void Frequency_Cav::setRFE0TL(double RFE0TL)
{
  _RFE0TL = RFE0TL;
}

double Frequency_Cav::getRFE0TL()
{
   return _RFE0TL;
}

void Frequency_Cav::setRFPhase(double RFPhase)
{
  _RFPhase = RFPhase;
}

double Frequency_Cav::getRFPhase()
{
  return _RFPhase;
}

void Frequency_Cav::trackBunch(Bunch* bunch)
{
  double RFFreq  = _RFFreq;
  double RFE0TL  = _RFE0TL;
  double RFPhase = OrbitConst::PI * _RFPhase / 180.0;

  bunch->compress();
  SyncPart* syncPart   = bunch->getSyncPart();

  // Bunch characteristics
  double mass          = bunch->getMass();
  double charge        = bunch->getCharge();

  // Synchronous particle acceleration
  double chargeE0TLsin = charge * RFE0TL * sin(RFPhase);
  double delta_eK      = charge * RFE0TL * cos(RFPhase);

  // Synchronous particle parameters at the start of the gap
  double eK_in         = syncPart->getEnergy();
  double gamma_in      = syncPart->getGamma();
  double beta_in       = syncPart->getBeta();

  // Synchronous particle parameters in the middle of the gap
  double eK_gap        = eK_in + delta_eK / 2.0;
  syncPart->setMomentum(syncPart->energyToMomentum(eK_gap));
  double gamma_gap     = syncPart->getGamma();
  double beta_gap      = syncPart->getBeta();

  // Synchronous particle parameters at the end of the gap
  double eK_out        = eK_in + delta_eK;
  syncPart->setMomentum(syncPart->energyToMomentum(eK_out));
  double gamma_out     = syncPart->getGamma();
  double beta_out      = syncPart->getBeta();

  double adbtc_dmp     = (beta_in * gamma_in)/(beta_out * gamma_out);

  // The base RF gap is simple - no phase correction.
  // The delta time is in seconds
  double delta_time    = 0.;
  syncPart->setTime(syncPart->getTime() + delta_time);

  //wave momentum
  double k             = 2.0 * OrbitConst::PI * RFFreq / OrbitConst::c;
  double ZtoPhi        = k / beta_in;
  double kr            = k / (gamma_in * beta_in);

  //transverse focusing coeff
  double d_rp          = -chargeE0TLsin * k /
                         (2.0 * mass * beta_gap * beta_gap * beta_out *
                          gamma_gap * gamma_gap * gamma_out);
  
  double x, y, r, rp, d_phi;
  double I0, I1;
  for(int i = 0; i < bunch->getSize(); i++)
  {
    x  = bunch->x(i);
    y  = bunch->y(i);
    r  = sqrt(x * x + y * y);
    I0 = bessi0(kr * r);
    I1 = bessi1(kr * r);

    //longitudinal-energy part
    d_phi        = -bunch->z(i) * ZtoPhi;
    bunch->dE(i) =  bunch->dE(i) - chargeE0TLsin * d_phi;

    //transverse focusing
    bunch->xp(i) =  bunch->xp(i) * adbtc_dmp + d_rp * x;
    bunch->yp(i) =  bunch->yp(i) * adbtc_dmp + d_rp * y;
  }
}
