// This class represents a simple RF gap.
// For this RF gap we know the E0TL parameter only.

#ifndef BASE_RF_GAP_H
#define BASE_RF_GAP_H

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

class BaseRfGap: public OrbitUtils::CppPyWrapper
{
public:

  /** Constructor for Base RF gap*/
  BaseRfGap();

  /** Destructor */
  virtual ~BaseRfGap();

  /** Tracks the Bunch trough the RF gap. */
  void trackBunch(Bunch* bunch, double frequency, double E0TL, double phase);

private:

protected:

};

#endif
