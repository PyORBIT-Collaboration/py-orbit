#include "Barrier_Cav.hh"
#include "ParticleMacroSize.hh"

#include <iostream>
#include <cmath>

#include "Bunch.hh"
#include "OrbitConst.hh"

using namespace OrbitUtils;

// Constructor
Barrier_Cav::Barrier_Cav(double ZtoPhi,
                         double RFVoltage,
                         double RFPhasep,
                         double RFPhasem,
                         double dRFPhasep,
                         double dRFPhasem): CppPyWrapper(NULL)
{
  _ZtoPhi    = ZtoPhi;
  _RFVoltage = RFVoltage;
  _RFPhasep  = RFPhasep;
  _RFPhasem  = RFPhasem;
  _dRFPhasep = dRFPhasep;
  _dRFPhasem = dRFPhasem;

}

// Destructor
Barrier_Cav::~Barrier_Cav()
{
}

void Barrier_Cav::setZtoPhi(double ZtoPhi)
{
  _ZtoPhi = ZtoPhi;
}

double Barrier_Cav::getZtoPhi()
{
  return _ZtoPhi;
}

void Barrier_Cav::setRFVoltage(double RFVoltage)
{
  _RFVoltage = RFVoltage;
}

double Barrier_Cav::getRFVoltage()
{
  return _RFVoltage;
}

void Barrier_Cav::setRFPhasep(double RFPhasep)
{
  _RFPhasep = RFPhasep;
}

double Barrier_Cav::getRFPhasep()
{
  return _RFPhasep;
}

void Barrier_Cav::setRFPhasem(double RFPhasem)
{
  _RFPhasem = RFPhasem;
}

double Barrier_Cav::getRFPhasem()
{
  return _RFPhasem;
}

void Barrier_Cav::setdRFPhasep(double dRFPhasep)
{
  _dRFPhasep = dRFPhasep;
}

double Barrier_Cav::getdRFPhasep()
{
  return _dRFPhasep;
}

void Barrier_Cav::setdRFPhasem(double dRFPhasem)
{
  _dRFPhasem = dRFPhasem;
}

double Barrier_Cav::getdRFPhasem()
{
  return _dRFPhasem;
}

void Barrier_Cav::trackBunch(Bunch* bunch)
{
  double ZtoPhi    = _ZtoPhi;
  double RFVoltage = _RFVoltage;
  double RFPhasep  = OrbitConst::PI * _RFPhasep / 180.0;
  double RFPhasem  = OrbitConst::PI * _RFPhasem / 180.0;
  double dRFPhasep = OrbitConst::PI * _dRFPhasep / 180.0;
  double dRFPhasem = OrbitConst::PI * _dRFPhasem / 180.0;

  double phase, dERF, arg;

  bunch->compress();
  SyncPart* syncPart = bunch->getSyncPart();
  double** arr = bunch->coordArr();

  for(int i = 0; i < bunch->getSize(); i++)
  {
    phase = -ZtoPhi * arr[i][4];
    if(phase < -OrbitConst::PI) phase += 2.0 * OrbitConst::PI;
    if(phase >  OrbitConst::PI) phase -= 2.0 * OrbitConst::PI;

    dERF = 0.0;

    if((phase > RFPhasep - dRFPhasep) &&
       (phase < RFPhasep + dRFPhasep))
    {
        arg = OrbitConst::PI * (phase - RFPhasep) / (2.0 * dRFPhasep);
        dERF = RFVoltage * cos(arg);
    }

    if((phase > RFPhasem - dRFPhasem) &&
       (phase < RFPhasem + dRFPhasem))
    {
        arg = OrbitConst::PI * (phase - RFPhasem) / (2.0 * dRFPhasem);
        dERF = -RFVoltage * cos(arg);
    }

    arr[i][5] += dERF;
  }
}
