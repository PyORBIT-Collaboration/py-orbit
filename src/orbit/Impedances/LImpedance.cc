/////////////////////////////////////////////////////////////////////////////
//
// FILE NAME
//   LImpedance.cc
//
//   12/19/2014
//
// DESCRIPTION
//   Calculate the effects of longitudinal impedance on the bunch
//   Does not include space charge.
//
/////////////////////////////////////////////////////////////////////////////

#include "Grid1D.hh"
#include "BufferStore.hh"
#include "LImpedance.hh"
#include "OrbitConst.hh"
#include <complex>
#include <iostream>
#include <cmath>
#include <cfloat>
#include <complex>

//FFTW library header
#include "fftw3.h"

using namespace OrbitUtils;


LImpedance::LImpedance(double length,
                       int nMacrosMin,
                       int nBins): CppPyWrapper(NULL)
{
  _length        = length;
  _nMacrosMin    = nMacrosMin;
  _nBins         = nBins;
  zGrid          = new Grid1D(_nBins, _length);

  _fftmagnitude  = new double[_nBins / 2];
  _fftphase      = new double[_nBins / 2];
  _z             = new double[_nBins / 2];
  _chi           = new double[_nBins / 2];
  _zImped_n      = new std::complex<double>[_nBins / 2];

  for(int n = 0; n < _nBins / 2; n++)
  {
    _fftmagnitude[n] = 0.0;
    _fftphase[n]     = 0.0;
    _z[n]            = 0.0;
    _chi[n]          = 0.0;
    _zImped_n[n] = std::complex<double>(0.0, 0.0);
  }

  _in   = (fftw_complex *) fftw_malloc(_nBins * sizeof(fftw_complex));
  _out  = (fftw_complex *) fftw_malloc(_nBins * sizeof(fftw_complex));
  _plan = fftw_plan_dft_1d(_nBins, _in,  _out, FFTW_FORWARD, FFTW_ESTIMATE);
}


LImpedance::~LImpedance()
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


void LImpedance::assignImpedanceValue(int n, double real, double imag)
{
  _zImped_n[n + 1] = std::complex<double>(real, imag);
}


Grid1D* LImpedance::getLongGrid()
{
  return zGrid;
}


void LImpedance::trackBunch(Bunch* bunch)
{
  int nPartsGlobal = bunch->getSizeGlobal();
  if(nPartsGlobal < _nMacrosMin) return;

// Bin the particles

  double zmin, zmax;
  double realPart, imagPart;

  bunchExtremaCalc->getExtremaZ(bunch, zmin, zmax);
  double zextra = (_length - (zmax - zmin)) / 2.0;
  zmax += zextra;
  zmin  = zmax - _length;
  zGrid->setGridZ(zmin, zmax);
  zGrid->setZero();
  zGrid->binBunchSmoothedByParticle(bunch);
  zGrid->synchronizeMPI(bunch->getMPI_Comm_Local());

// FFT the beam density

  for(int i = 0; i < _nBins; i++)
  {
    _in[i][0] = zGrid->getValueOnGrid(i);
    _in[i][1] = 0.;
  }

  fftw_execute(_plan);

// Find the magnitude and phase

  _fftmagnitude[0] = _out[0][0]/(double)_nBins;
  _fftphase[0] = 0.;

  for(int n = 1; n < _nBins / 2; n++)
  {
    realPart = _out[n][0]/((double)_nBins);
    imagPart = _out[n][1]/((double)_nBins);
    _fftmagnitude[n] = sqrt(realPart * realPart + imagPart * imagPart);
    _fftphase[n] = atan2(imagPart, realPart);
  }

  SyncPart* sp = bunch->getSyncPart();

  for(int n = 1; n < _nBins / 2; n++)
  {
    realPart = std::real(_zImped_n[n]);
    imagPart = std::imag(_zImped_n[n]);
    _z[n] = n * sqrt(realPart * realPart + imagPart * imagPart);
    _chi[n] = atan2(imagPart, realPart);
  }

// Convert charge to current for a single macroparticle per unit bin length

  double charge2current = bunch->getCharge() * bunch->getMacroSize() *
                          OrbitConst::elementary_charge_MKS * sp->getBeta() *
                          OrbitConst::c / (_length / _nBins);

// Calculate and add the kick to macroparticles

  double kick, angle, position;

// Convert particle longitudinal coordinate to phi

  double philocal;
  double z;
  double phi[bunch->getSize()];
  double** coords = bunch->coordArr();
  for (int j = 0; j < bunch->getSize(); j++)
  {
    z = bunch->z(j);
    philocal = (z / _length) * 2 * OrbitConst::PI;

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
//    LImpedance::_kick
//
// DESCRIPTION
//    Calculates the longitudinal impedance kick to a macroparticle
//
// PARAMETERS
//    angle - the phase angle of the particle (rad)
//
// RETURNS
//    kick
//
///////////////////////////////////////////////////////////////////////////

double LImpedance::_kick(double angle)
{
// n=0 term has no impact (constant in phi)
// f(phi) = _FFTMagnitude(1) + sum (n = 2 -> N / 2) of
//          [2 * _FFTMagnitude(i) * cos(phi * (n - 1) +
//          _FFTPhase(n) + _chi(n))]

  double kick = 0.;
  double cosArg;

  for(int n = 1; n < _nBins / 2; n++)
  {
    cosArg = n * (angle + OrbitConst::PI) + _fftphase[n] + _chi[n];
    kick += 2 * _fftmagnitude[n] * _z[n] * cos(cosArg);
  }
  return kick;
}
