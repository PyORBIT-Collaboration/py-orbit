/////////////////////////////////////////////////////////////////////////////
//
// FILE NAME
//   TImpedance.cc
//
//   12/23/2014
//
// DESCRIPTION
//   Calculate the effects of transverse impedance on the bunch
//
/////////////////////////////////////////////////////////////////////////////

#include "Grid1D.hh"
#include "BufferStore.hh"
#include "TImpedance.hh"
#include "OrbitConst.hh"
#include <complex>
#include <iostream>
#include <cmath>
#include <cfloat>
#include <complex>

//FFTW library header
#include "fftw3.h"

using namespace OrbitUtils;


TImpedance::TImpedance(double length,
                       int nMacrosMin,
                       int nBins,
                       int useX,
                       int useY): CppPyWrapper(NULL)
{
  _length       = length;
  _nMacrosMin   = nMacrosMin;
  _nBins        = nBins;
  _useX         = useX;
  _useY         = useY;
  zGrid         = new Grid1D(_nBins, _length);

  _qX = 0.0;
  _qY = 0.0;
  _charge2TKick = 0.0;
  _Turns        = 0;

  _q       = 0.0;
  _alpha   = 0.0;
  _beta    = 0.0;

  if(_useX)
  {
    _xCentroid      = new double[_nBins];
    _xpCentroid     = new double[_nBins];
    _zXImped_nplus  = new std::complex<double>[_nBins];
    _zXImped_nminus = new std::complex<double>[_nBins];
    for(int n = 0; n < _nBins; n++)
    {
      _xCentroid[n]      = 0.0;
      _xpCentroid[n]     = 0.0;
      _zXImped_nplus[n]  = std::complex<double>(0.0, 0.0);
      _zXImped_nminus[n] = std::complex<double>(0.0, 0.0);
    }
  }
  if(_useY)
  {
    _yCentroid      = new double[_nBins];
    _ypCentroid     = new double[_nBins];
    _zYImped_nplus  = new std::complex<double>[_nBins];
    _zYImped_nminus = new std::complex<double>[_nBins];
    for(int n = 0; n < _nBins; n++)
    {
      _yCentroid[n]      = 0.0;
      _ypCentroid[n]     = 0.0;
      _zYImped_nplus[n]  = std::complex<double>(0.0, 0.0);
      _zYImped_nminus[n] = std::complex<double>(0.0, 0.0);
    }
  }

  _phiCount      = new double[_nBins];
  _FFTResult1    = new std::complex<double>[_nBins];
  _FFTResult2    = new std::complex<double>[_nBins];
  _Centroid      = new double[_nBins];
  _pCentroid     = new double[_nBins];
  _zImped_nplus  = new std::complex<double>[_nBins];
  _zImped_nminus = new std::complex<double>[_nBins];
  for(int n = 0; n < _nBins; n++)
  {
    _phiCount[n]      = 0.0;
    _FFTResult1[n]    = std::complex<double>(0.0, 0.0);
    _FFTResult2[n]    = std::complex<double>(0.0, 0.0);
    _Centroid[n]      = 0.0;
    _pCentroid[n]     = 0.0;
    _zImped_nplus[n]  = std::complex<double>(0.0, 0.0);
    _zImped_nminus[n] = std::complex<double>(0.0, 0.0);
  }

  _in   = (fftw_complex *) fftw_malloc(_nBins * sizeof(fftw_complex));
  _out  = (fftw_complex *) fftw_malloc(_nBins * sizeof(fftw_complex));
  _plan = fftw_plan_dft_1d(_nBins, _in,  _out, FFTW_FORWARD, FFTW_ESTIMATE);
}

TImpedance::~TImpedance()
{
  if(zGrid->getPyWrapper() != NULL)
  {
    Py_DECREF(zGrid->getPyWrapper());
  }
  else
  {
    delete zGrid;
  }
  if(_useX)
  {
    delete[] _xCentroid;
    delete[] _xpCentroid;
    delete[] _zXImped_nplus;
    delete[] _zXImped_nminus;
  }
  if(_useY)
  {
    delete[] _yCentroid;
    delete[] _ypCentroid;
    delete[] _zYImped_nplus;
    delete[] _zYImped_nminus;
  }
  delete[] _FFTResult1;
  delete[] _FFTResult2;  
  delete[] _phiCount;
  delete[] _Centroid;
  delete[] _pCentroid;
  delete[] _zImped_nplus;
  delete[] _zImped_nminus;

  fftw_free(_in);
  fftw_free(_out);
  fftw_destroy_plan(_plan);
}

///////////////////////////////////////////////////////////////////////////
//
// NAME
//   TImpedance::assignLatFuncs
//
// DESCRIPTION
//   Sets the needed lattice functions
//
// PARAMETERS
//   qX     - horizontal lattice tune
//   qY     - vertical   lattice tune
//   alphaX - horizontal Twiss parameter
//   betaX  - horizontal Twiss parameter
//   alphaY - vertical   Twiss parameter
//   betaY  - vertical   Twiss parameter

//
// RETURNS
//   Nothing
//
///////////////////////////////////////////////////////////////////////////

void TImpedance::assignLatFuncs(double qX, double alphaX, double betaX,
                                double qY, double alphaY, double betaY)
{
  _qX     = qX;
  _alphaX = alphaX;
  _betaX  = betaX;
  _qY     = qY;
  _alphaY = alphaY;
  _betaY  = betaY;
}

///////////////////////////////////////////////////////////////////////////
//
// NAME
//   TImpedance::assignImpedanceX
//
// DESCRIPTION
//   Sets the real and imaginary values of the nth modes of
//   the fast and slow contributions to the horizontal impedance
//
// PARAMETERS
//   n
//   realp - real      part of fast contribution
//   imagp - imaginary part of fast contribution
//   realm - real      part of slow contribution
//   imagm - imaginary part of slow contribution
//
// RETURNS
//   Nothing
//
///////////////////////////////////////////////////////////////////////////

void TImpedance::assignImpedanceX(int n,
                                  double realp,
                                  double imagp,
                                  double realm,
                                  double imagm)
{
  _zXImped_nplus[n]  = std::complex<double>(realp, imagp);
  _zXImped_nminus[n] = std::complex<double>(realm, imagm);
}

///////////////////////////////////////////////////////////////////////////
//
// NAME
//   TImpedance::assignImpedanceY
//
// DESCRIPTION
//   Sets the real and imaginary values of the nth modes of
//   the fast and slow contributions to the vertical impedance
//
// PARAMETERS
//   n
//   realp - real      part of fast contribution
//   imagp - imaginary part of fast contribution
//   realm - real      part of slow contribution
//   imagm - imaginary part of slow contribution
//
// RETURNS
//   Nothing
//
///////////////////////////////////////////////////////////////////////////

void TImpedance::assignImpedanceY(int n,
                                  double realp,
                                  double imagp,
                                  double realm,
                                  double imagm)
{
  _zYImped_nplus[n]  = std::complex<double>(realp, imagp);
  _zYImped_nminus[n] = std::complex<double>(realm, imagm);
}

///////////////////////////////////////////////////////////////////////////
//
// NAME
//   TImpedance::trackBunch(Bunch* bunch)
//
// DESCRIPTION
//   Applies the transverse impedance kick to the bunch
//
// PARAMETERS
//   bunch - the bunch to be kicked
//
// RETURNS
//   Nothing
//
///////////////////////////////////////////////////////////////////////////

void TImpedance::trackBunch(Bunch* bunch)
{
  bunch->compress();
  int nPartsGlobal = bunch->getSizeGlobal();
  if(nPartsGlobal < _nMacrosMin) return;

  SyncPart* sp = bunch->getSyncPart();

// coefficient to transform charge harmonics to a transverse kick
// Transverse impedance is not yet included

  _charge2TKick = bunch->getMacroSize() * 1.e-9 *
                  bunch->getCharge() * bunch->getCharge() *
                  OrbitConst::elementary_charge_MKS * OrbitConst::c /
                  (sp->getBeta() * _length *
                  (sp->getEnergy() + sp->getMass()));

// Bin the particles

  int n, i;
  double zmin, zmax;

  bunchExtremaCalc->getExtremaZ(bunch, zmin, zmax);
  double zextra = (_length - (zmax - zmin)) / 2.0;
  zmax += zextra;
  zmin  = zmax - _length;
  zGrid->setGridZ(zmin, zmax);
  zGrid->setZero();
  zGrid->binBunchSmoothedByParticle(bunch);
  for(n = 0; n < _nBins; n++)
  {
    _phiCount[n] = zGrid->getValueOnGrid(n);
  }
  zGrid->binBunchSmoothedMoment(0, bunch, _xCentroid);
  zGrid->binBunchSmoothedMoment(1, bunch, _xpCentroid);
  zGrid->binBunchSmoothedMoment(2, bunch, _yCentroid);
  zGrid->binBunchSmoothedMoment(3, bunch, _ypCentroid);
  zGrid->synchronizeMPI(bunch->getMPI_Comm_Local());

  // FFT of beam and application of kick

  double twopi = 2.0 * OrbitConst::PI;
  double macrophase;
  double** part_coord_arr = bunch->coordArr();

  // Horizontal calculation

  if(_useX)
  {
    _q       = _qX;
    _alpha   = _alphaX;
    _beta    = _betaX;
    for(n = 0; n < _nBins; n++)
    {
      _Centroid[n]      = _xCentroid[n];
      _pCentroid[n]     = _xpCentroid[n];
      _zImped_nplus[n]  = _zXImped_nplus[n];
      _zImped_nminus[n] = _zXImped_nminus[n];
    }

    _prepareToKick();

    for(i = 0; i < bunch->getSize(); i++)
    {
      macrophase = twopi * (part_coord_arr[i][4] - zmin) / _length;
      part_coord_arr[i][1] += _kick(macrophase);
    }
  }
  // end X dimension

  // FFT harmonics for Y

  if(_useY)
  {
    _q       = _qY;
    _alpha   = _alphaY;
    _beta    = _betaY;
    for(n = 0; n < _nBins; n++)
    {
      _Centroid[n]      = _yCentroid[n];
      _pCentroid[n]     = _ypCentroid[n];
      _zImped_nplus[n]  = _zYImped_nplus[n];
      _zImped_nminus[n] = _zYImped_nminus[n];
    }

    _prepareToKick();

    for(i = 0; i < bunch->getSize(); i++)
    {
      macrophase = twopi * (part_coord_arr[i][4] - zmin) / _length;
      part_coord_arr[i][3] += _kick(macrophase);
    }
  }
  // end Y dimension

  _Turns++;
}

///////////////////////////////////////////////////////////////////////////
//
// NAME
//   TImpedance::_prepareToKick()
//
// DESCRIPTION
//   Does FFTs in preparation for kick
//
// PARAMETERS
//   None
//
// RETURNS
//   Nothing
//
///////////////////////////////////////////////////////////////////////////

void TImpedance::_prepareToKick()
{
  int n;
  double twopi = 2.0 * OrbitConst::PI;
  double delta, Phase, dPhase, coeff;

  // calculate the coefficient of cos(omega_betatron * t)

  delta = _Turns * twopi * _q;
  dPhase = twopi * _q / double(_nBins);
  for(n = 0; n < _nBins; n++)
  {
    Phase = n * dPhase + delta;
    coeff = _Centroid[n] * cos(Phase) -
            (_beta  * _pCentroid[n] +
             _alpha * _Centroid[n]) *
            sin(Phase);

    _in[n][0] = coeff * double(_nBins);
    _in[n][1] = 0.;
  }

  // Do the Forward FFT:

  fftw_execute(_plan);
//  fftw(_plan, 1, _in, 1, 0, _out, 1, 0);
  for(n = 0; n < _nBins; n++)
  {
    _FFTResult1[n] = std::complex<double>(_out[n][0], _out[n][1])
                     / double(_nBins);
  }

  // calculate the coefficient of sin(omega_betatron * t)

  for (n = 0; n < _nBins; n++)
  {
    Phase = n * dPhase + delta;
    coeff = _Centroid[n] * sin(Phase) +
            (_beta  * _pCentroid[n] +
             _alpha * _Centroid[n]) *
            cos(Phase);

    _in[n][0] = coeff * double(_nBins);
    _in[n][1] = 0.;
  }

  // Do the Forward FFT:

  fftw_execute(_plan);
//  fftw(_plan, 1, _in, 1, 0, _out, 1, 0);
  for(n = 0; n < _nBins; n++)
  {
    _FFTResult2[n] = std::complex<double>(_out[n][0], _out[n][1])
                     / double(_nBins);
  }
}

///////////////////////////////////////////////////////////////////////////
//
// NAME
//   TImpedance::_kick(macrophase)
//
// DESCRIPTION
//   Returns the transverse kick to a macroparticle
//
// PARAMETERS
//   macrophase - the phase macrophase of the particle (rad)
//
// RETURNS
//   Transverse kick
//
///////////////////////////////////////////////////////////////////////////

double TImpedance::_kick(double macrophase)
{
  double twopi = 2.0 * OrbitConst::PI;
  std::complex<double> sqrtm1 = std::complex<double>(0, 1);;
  double coeff = 0.;
  double dPhasep, dPhasem, delta;

  // delta is betatron phase advance from all previous turns

  delta = _Turns * twopi * _q;

  for(int n = 0; n < _nBins / 2; n++)
  {
    dPhasep = ( _q + n) * (macrophase);
    dPhasem = (-_q + n) * (macrophase);
    coeff += -std::imag((_FFTResult1[n] -
                         sqrtm1 * _FFTResult2[n]) *
                         _zImped_nplus[n] *
                        (cos(dPhasep + delta) +
                         sqrtm1 * sin(dPhasep + delta)) +
                        (_FFTResult1[n] +
                         sqrtm1 * _FFTResult2[n]) *
                        _zImped_nminus[n] *
                        (cos(dPhasem - delta) +
                         sqrtm1 * sin(dPhasem - delta)));
  }

  // now assign negative mode numbers to the upper half modes

  for(int n = _nBins / 2; n < _nBins; n++)
  {
    dPhasep = ( _q + n - _nBins) * (macrophase);
    dPhasem = (-_q + n - _nBins) * (macrophase);
    coeff += -std::imag((_FFTResult1[n] -
                         sqrtm1 * _FFTResult2[n]) *
                         _zImped_nplus[n] *
                        (cos(dPhasep + delta) +
                         sqrtm1 * sin(dPhasep + delta)) +
                        (_FFTResult1[n] +
                         sqrtm1 * _FFTResult2[n]) *
                         _zImped_nminus[n] *
                        (cos(dPhasem - delta) +
                         sqrtm1 * sin(dPhasem - delta)));
  }

  return (-_charge2TKick * coeff / 2.);
  // the - sign above gives instability for positive real Zminus
}

