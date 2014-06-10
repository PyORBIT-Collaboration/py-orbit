#include "Dual_Harmonic_Cav.hh"
#include "ParticleMacroSize.hh"

#include <iostream>
#include <cmath>

#include "Bunch.hh"
#include "OrbitConst.hh"

using namespace OrbitUtils;

// Constructor
Dual_Harmonic_Cav::Dual_Harmonic_Cav(double ZtoPhi   ,
                           double RFHNum   ,
                           double RatioRFHNum,
                           double RFVoltage,
                           double RatioVoltage,
                           double RFPhase,
                           double RFPhase2): CppPyWrapper(NULL)
{
  _ZtoPhi    = ZtoPhi;
  _RFHNum    = RFHNum;
  _RatioRFHNum = RatioRFHNum;
  _RFVoltage = RFVoltage;
  _RatioVoltage = RatioVoltage;
  _RFPhase   = RFPhase;
  _RFPhase2   = RFPhase2;

}

// Destructor
Dual_Harmonic_Cav::~Dual_Harmonic_Cav()
{
}

void Dual_Harmonic_Cav::setZtoPhi(double ZtoPhi)
{
  _ZtoPhi = ZtoPhi;
}

double Dual_Harmonic_Cav::getZtoPhi()
{
  return _ZtoPhi;
}


void Dual_Harmonic_Cav::setRatioRFHNum(double RatioRFHNum)
{
  _RatioRFHNum = RatioRFHNum;
}

double Dual_Harmonic_Cav::getRatioRFHNum()
{
  return _RatioRFHNum;
}

void Dual_Harmonic_Cav::setRFHNum(double RFHNum)
{
 _RFHNum = RFHNum;
}

double Dual_Harmonic_Cav::getRFHNum()
{
  return _RFHNum;
}

void Dual_Harmonic_Cav::setRFVoltage(double RFVoltage)
{
  _RFVoltage = RFVoltage;
}

double Dual_Harmonic_Cav::getRFVoltage()
{
  return _RFVoltage;
}

void Dual_Harmonic_Cav::setRatioVoltage(double RatioVoltage)
{
  _RatioVoltage = RatioVoltage;
}

double Dual_Harmonic_Cav::getRatioVoltage()
{
  return _RatioVoltage;
}

void Dual_Harmonic_Cav::setRFPhase(double RFPhase)
{
  _RFPhase = RFPhase;

}

double Dual_Harmonic_Cav::getRFPhase()
{
  return _RFPhase;
}

void Dual_Harmonic_Cav::setRFPhase2(double RFPhase2)
{
  _RFPhase2 = RFPhase2;

}

double Dual_Harmonic_Cav::getRFPhase2()
{
  return _RFPhase2;
}

void Dual_Harmonic_Cav::trackBunch(Bunch* bunch)
{
  double ZtoPhi    = _ZtoPhi;
  double RFHNum    = _RFHNum;
  double RatioRFHNum = _RatioRFHNum;
  double RFVoltage = _RFVoltage;
  double RatioVoltage = _RatioVoltage;
  double RFPhase   = OrbitConst::PI * _RFPhase / 180.0;
  double RFPhase2 = OrbitConst::PI * _RFPhase2 / 180.0;

  double dERF, phase;

  bunch->compress();
  SyncPart* syncPart = bunch->getSyncPart();
  double** arr = bunch->coordArr();

  for(int i = 0; i < bunch->getSize(); i++)
  {
    phase = -ZtoPhi * arr[i][4] * RFHNum ;
	 dERF  = bunch->getCharge()* RFVoltage * (sin(phase) - sin(RFPhase) - RatioVoltage * ( sin(RFPhase + RatioRFHNum * (phase - RFPhase))  - sin(RFPhase2) ));
    arr[i][5] += dERF;
  }
}
