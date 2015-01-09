/////////////////////////////////////////////////////////////////////////////
//
// FILE NAME
//   LSpaceCharge.cc
//
//   05/30/12
//
// DESCRIPTION
//   Calculate the effects of longitudinal space charge and impedance
//   on the bunch (1D) using impedance formulation
//
/////////////////////////////////////////////////////////////////////////////

#include "Grid1D.hh"
#include "BufferStore.hh"
#include "LSpaceChargeCalc.hh"
#include "OrbitConst.hh"
#include <complex>
#include <iostream>
#include <cmath>
#include <cfloat>
#include <complex>

//FFTW library header
#include "fftw3.h"

using namespace OrbitUtils;


LSpaceChargeCalc::LSpaceChargeCalc(double b_a_in, double length_in, int nMacrosMin_in, int useSpaceCharge_in, int nBins_in): CppPyWrapper(NULL)
{
  b_a            = b_a_in;
  length         = length_in;
  nMacrosMin     = nMacrosMin_in;
  useSpaceCharge = useSpaceCharge_in;
  nBins          = nBins_in;
  zGrid          = new Grid1D(nBins, length);

  _fftmagnitude  = new double[nBins / 2];
  _fftphase      = new double[nBins / 2];
  _z             = new double[nBins / 2];
  _chi           = new double[nBins / 2];
  _zImped_n      = new std::complex<double>[nBins / 2];

  for(int n = 0; n < nBins / 2; n++)
  {
    _fftmagnitude[n] = 0.0;
    _fftphase[n]     = 0.0;
    _z[n]            = 0.0;
    _chi[n]          = 0.0;
    _zImped_n[n] = std::complex<double>(0.0, 0.0);
  }

  _in   = (fftw_complex *) fftw_malloc(nBins * sizeof(fftw_complex));
  _out  = (fftw_complex *) fftw_malloc(nBins * sizeof(fftw_complex));
  _plan = fftw_plan_dft_1d(nBins, _in,  _out, FFTW_FORWARD, FFTW_ESTIMATE);
}


LSpaceChargeCalc::~LSpaceChargeCalc()
{
  if(zGrid->getPyWrapper() != NULL)
  {
    Py_DECREF(zGrid->getPyWrapper());
  }
  else
  {
    delete zGrid;
  }
  delete[] _fftmagnitude;
  delete[] _fftphase;
  delete[] _z;
  delete[] _chi;
  delete[] _zImped_n;
  fftw_free(_in);
  fftw_free(_out);
  fftw_destroy_plan(_plan);
}


void LSpaceChargeCalc::assignImpedanceValue(int n, double real, double imag)
{
  _zImped_n[n + 1] = std::complex<double>(real, imag);
}


Grid1D* LSpaceChargeCalc::getLongGrid()
{
  return zGrid;
}


void LSpaceChargeCalc::trackBunch(Bunch* bunch)
{
  int nPartsGlobal = bunch->getSizeGlobal();
  if(nPartsGlobal < nMacrosMin) return;

// Bin the particles

  double zmin, zmax;
  double realPart, imagPart;

  bunchExtremaCalc->getExtremaZ(bunch, zmin, zmax);
  double zextra = (length - (zmax - zmin)) / 2.0;
  zmax += zextra;
  zmin  = zmax - length;
  zGrid->setGridZ(zmin, zmax);
  zGrid->setZero();
  zGrid->binBunchSmoothedByParticle(bunch);
  zGrid->synchronizeMPI(bunch->getMPI_Comm_Local());

// FFT the beam density

  for(int i = 0; i < nBins; i++)
  {
    _in[i][0] = zGrid->getValueOnGrid(i);
    _in[i][1] = 0.;
  }

  fftw_execute(_plan);

// Find the magnitude and phase

  _fftmagnitude[0] = _out[0][0]/(double)nBins;
  _fftphase[0] = 0.;

  for(int n = 1; n < nBins / 2; n++)
  {
    realPart = _out[n][0]/((double)nBins);
    imagPart = _out[n][1]/((double)nBins);
    _fftmagnitude[n] = sqrt(realPart * realPart + imagPart * imagPart);
    _fftphase[n] = atan2(imagPart,realPart);
  }

// Space charge contribution to impedance

  SyncPart* sp = bunch->getSyncPart();

// Set zero if useSpaceCharge = 0

  double zSpaceCharge_n = 0.;

// Otherwise, set positive since space charge is capacitive (Chao convention)

  if(useSpaceCharge != 0)
  {
    double mu_0 = 4.0 * OrbitConst::PI * 1.e-07; // permeability of free space
    double _z_0 = OrbitConst::c * mu_0;

    zSpaceCharge_n = _z_0 * (1.0 + 2.0 * log(b_a)) /
                     (2 * sp->getBeta() * sp->getGamma() * sp->getGamma());
  }

  for(int n = 1; n < nBins / 2; n++)
  {
    realPart = std::real(_zImped_n[n]);
    imagPart = std::imag(_zImped_n[n]) + zSpaceCharge_n;
    _z[n] = n * sqrt(realPart * realPart + imagPart * imagPart);
    _chi[n] = atan2(imagPart, realPart);
  }

// Convert charge to current for a single macroparticle per unit bin length

  double charge2current = bunch->getCharge() * bunch->getMacroSize() *
                          OrbitConst::elementary_charge_MKS * sp->getBeta() *
                          OrbitConst::c / (length / nBins);

// Calculate and add the kick to macroparticles

  double kick, angle, position;

// Convert particle longitudinal coordinate to phi

  double philocal;
  double z;
  double** coords = bunch->coordArr();
  for (int j = 0; j < bunch->getSize(); j++)
  {
    z = bunch->z(j);
    philocal = (z / length) * 2 * OrbitConst::PI;

// Handle cases where the longitudinal coordinate is
// outside of the user-specified length

  if(philocal < -OrbitConst::PI) philocal += 2 * OrbitConst::PI;
  if(philocal >  OrbitConst::PI) philocal -= 2 * OrbitConst::PI;

  double dE = _kick(philocal) * (-1e-9) *
              bunch->getCharge() * charge2current;
  coords[j][5] += dE;
  }
}


///////////////////////////////////////////////////////////////////////////
//
// NAME
//    LSpaceChargeCalc::_kick
//
// DESCRIPTION
//    Calculates the longitudinal space charge and impedance kick
//    to a macroparticle
//
// PARAMETERS
//    angle - the phase angle of the particle (rad)
//
// RETURNS
//    kick
//
///////////////////////////////////////////////////////////////////////////

double LSpaceChargeCalc::_kick(double angle)
{
// n=0 term has no impact (constant in phi)
// f(phi) = _FFTMagnitude(1) + sum (n = 2 -> N / 2) of
//          [2 * _FFTMagnitude(i) * cos(phi * (n - 1) +
//          _FFTPhase(n) + _chi(n))]

  double kick = 0.;
  double cosArg;

  for(int n = 1; n < nBins / 2; n++)
  {
    cosArg = n * (angle + OrbitConst::PI) + _fftphase[n] + _chi[n];
    kick += 2 * _fftmagnitude[n] * _z[n] * cos(cosArg);
  }
  return kick;
}
