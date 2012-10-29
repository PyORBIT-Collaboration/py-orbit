#include "Harmonic_Cav.hh"
#include "ParticleMacroSize.hh"

#include <iostream>
#include <cmath>

#include "Bunch.hh"
#include "OrbitConst.hh"

using namespace OrbitUtils;

// Constructor
Harmonic_Cav::Harmonic_Cav(double ZtoPhi   ,
                           double dESync   ,
                           double RFHNum   ,
                           double RFVoltage,
                           double RFPhase): CppPyWrapper(NULL)
{
  _ZtoPhi    = ZtoPhi;
  _dESync    = dESync;
  _RFHNum    = RFHNum;
  _RFVoltage = RFVoltage;
  _RFPhase   = RFPhase;
}

// Destructor
Harmonic_Cav::~Harmonic_Cav()
{
}

void Harmonic_Cav::setZtoPhi(double ZtoPhi)
{
  _ZtoPhi = ZtoPhi;
}

double Harmonic_Cav::getZtoPhi()
{
  return _ZtoPhi;
}

void Harmonic_Cav::setdESync(double dESync)
{
  _dESync = dESync;
}

double Harmonic_Cav::getdESync()
{
  return _dESync;
}

void Harmonic_Cav::setRFHNum(double RFHNum)
{
  _RFHNum = RFHNum;
}

double Harmonic_Cav::getRFHNum()
{
  return _RFHNum;
}

void Harmonic_Cav::setRFVoltage(double RFVoltage)
{
  _RFVoltage = RFVoltage;
}

double Harmonic_Cav::getRFVoltage()
{
  return _RFVoltage;
}

void Harmonic_Cav::setRFPhase(double RFPhase)
{
  _RFPhase = RFPhase;
}

double Harmonic_Cav::getRFPhase()
{
  return _RFPhase;
}

void Harmonic_Cav::trackBunch(Bunch* bunch)
{
  double ZtoPhi    = _ZtoPhi;
  double dESync    = _dESync;
  double RFHNum    = _RFHNum;
  double RFVoltage = _RFVoltage;
  double RFPhase   = _RFPhase;

  double dERF, phase;

  bunch->compress();
  SyncPart* syncPart = bunch->getSyncPart();
  double** arr = bunch->coordArr();

  for(int i = 0; i < bunch->getSize(); i++)
  {
    phase = ZtoPhi * arr[i][4];
    dERF  = bunch->getCharge()* RFVoltage * sin(RFHNum * phase + RFPhase);
    arr[i][5] += dERF - dESync;
  }
}
