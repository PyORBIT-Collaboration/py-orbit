// Calculate the transverse impedance effect on the bunch

#ifndef TIMPEDANCE_H
#define TIMPEDANCE_H

//MPI Function Wrappers
#include "orbit_mpi.hh"
#include "wrap_mpi_comm.hh"

#include <cstdlib>
#include <cmath>
#include <complex>

//ORBIT bunch
#include "Bunch.hh"

//pyORBIT utils
#include "CppPyWrapper.hh"
#include "BunchExtremaCalculator.hh"
#include "Grid1D.hh"

//FFTW library header
#include "fftw3.h"

using namespace std;

class TImpedance: public OrbitUtils::CppPyWrapper
{
public:

  /** Constructor */
  TImpedance(double length,
             int nMacrosMin,
             int nBins,
             int useX,
             int useY);

  /** Destructor */
  virtual ~TImpedance();

  /** Assigns the needed lattice functions **/
  void assignLatFuncs(double qX    , double qY,
                      double alphaX, double betaX,
                      double alphaY, double betaY);

  /** Assigns the real and imaginary parts of the
      horizontal impedance for index n **/
  void assignImpedanceX(int n,
                        double realp,
                        double imagp,
                        double realm,
                        double imagm);

  /** Assigns the real and imaginary parts of the
      vertical impedance for index n **/
  void assignImpedanceY(int n,
                        double realp,
                        double imagp,
                        double realm,
                        double imagm);

  /** Calculates the transverse kick to the
      macro-particles in the bunch **/
  void trackBunch(Bunch* bunch);

  /** Does FFTs in preparation for kick **/
  void _prepareToKick();

  /** Routine for calculating the transverse kick to the particle **/
  double _kick(double macrophase);

  /** Returns the 1D grid with a longitudinal density **/
  Grid1D* getLongGrid();

  double _length;
  int _nMacrosMin;
  int _nBins;
  int _useX;
  int _useY;

  double _qX;
  double _qY;
  double _alphaX;
  double _betaX;
  double _alphaY;
  double _betaY;
  double _charge2TKick;
  int _Turns;
  Grid1D* zGrid;
  OrbitUtils::BunchExtremaCalculator* bunchExtremaCalc;

  double* _xCentroid;
  double* _xpCentroid;
  double* _yCentroid;
  double* _ypCentroid;
  double* _phiCount;
  std::complex<double>* _zXImped_nplus;
  std::complex<double>* _zXImped_nminus;
  std::complex<double>* _zYImped_nplus;
  std::complex<double>* _zYImped_nminus;
  std::complex<double>* _FFTResult1;
  std::complex<double>* _FFTResult2;

  double _q;
  double _alpha;
  double _beta;
  double* _Centroid;
  double* _pCentroid;
  std::complex<double>* _zImped_nplus;
  std::complex<double>* _zImped_nminus;

  fftw_complex* _in;
  fftw_complex* _out;
  fftw_plan _plan;
};
//end of TIMPEDANCE_H
#endif

