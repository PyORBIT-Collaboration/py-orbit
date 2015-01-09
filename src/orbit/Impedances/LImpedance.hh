// Calculate the longitudinal impedance effect on the bunch

#ifndef LIMPEDANCE_H
#define LIMPEDANCE_H

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

class LImpedance: public OrbitUtils::CppPyWrapper
{
public:

  /** Constructor */
  LImpedance(double length, int nMacrosMin, int nBins);

  /** Destructor */
  virtual ~LImpedance();

  /** Calculates longitudinal impedance kick to the 
      macro-particles in the bunch */
  void trackBunch(Bunch* bunch);

  /** Assigns the real and imaginary parts of the
      machine impedance for index n**/
  void assignImpedanceValue(int n, double real, double imag);

  /** Routine for calculating the kick to the particle **/
  double _kick(double angle);

  /** Returns the 1D grid with a longitudinal density **/
  Grid1D* getLongGrid();


//private:
  double _length;
  int _nMacrosMin;
  int _nBins;

//protected:
  Grid1D* zGrid;
  OrbitUtils::BunchExtremaCalculator* bunchExtremaCalc;

  //FFT arrays
  double* _fftmagnitude;
  double* _fftphase;
  double* _z;
  double* _chi;

  std::complex<double>* _zImped_n;
  //vector<complex> _zImped_n;

  fftw_plan _plan;
  fftw_complex* _in;
  fftw_complex* _out;
};
//end of LIMPEDANCE_H
#endif

