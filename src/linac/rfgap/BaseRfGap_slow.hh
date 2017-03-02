// This class represents a simple RF gap.
// For this RF gap we know the E0TL parameter only.
// The arrival time is defined by dt = w*z/(beta*c)
// instead of dt = w*z/(beta_synch*c) in BaseRfGap.hh file
// The transformation coefficients are defined for each particle separately,
// so the name of the class has the word slow in the name.

#ifndef BASE_RF_GAP_SLOW_H
#define BASE_RF_GAP_SLOW_H

//MPI Function Wrappers
#include "orbit_mpi.hh"
#include "wrap_mpi_comm.hh"

#include <cstdlib>
#include <cmath>

//ORBIT bunch
#include "Bunch.hh"

//pyORBIT utils
#include "CppPyWrapper.hh"

using namespace std;

/** 
  This class represents a 2D rectangular grid.
*/

class BaseRfGap_slow: public OrbitUtils::CppPyWrapper
{
public:

  /** Constructor for Base RF gap*/
  BaseRfGap_slow();

  /** Destructor */
  virtual ~BaseRfGap_slow();

  /** Tracks the Bunch trough the RF gap. */
  static void trackBunch(Bunch* bunch, double frequency, double E0TL, double phase);

};

#endif