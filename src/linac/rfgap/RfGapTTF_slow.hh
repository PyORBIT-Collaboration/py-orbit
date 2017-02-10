/**
   This class represents a Parmila type RF gap. It acts on the coordinates 
   of the particle by using the transit time factors. The model includes  
   non-linearity in transverse direction.
   The transformation coefficients are defined for each particle separately,
   so the name of the class has the word slow in the name.   
*/

#ifndef TTF_RF_GAP_SLOW_H
#define TTF_RF_GAP_SLOW_H

//MPI Function Wrappers
#include "orbit_mpi.hh"
#include "wrap_mpi_comm.hh"


#include <cstdlib>
#include <cmath>

//ORBIT bunch
#include "Bunch.hh"

//pyORBIT utils
#include "CppPyWrapper.hh"
#include "OU_Polynomial.hh"

using namespace std;

/** 
  This class represents a RF gap as a Parmila type gap. 
*/

class RfGapTTF_slow: public OrbitUtils::CppPyWrapper
{
public:
	
	/** Constructor for Parmila's type RF gap with TTF */
  RfGapTTF_slow();
	
  /** Destructor */
  virtual ~RfGapTTF_slow();
	
	/** Tracks the Bunch through the RF gap. */	
	static void trackBunch(Bunch* bunch, double frequency, double E0L, double phase,
	                       OrbitUtils::Polynomial* Tttf,
	                       OrbitUtils::Polynomial* Sttf,
	                       OrbitUtils::Polynomial* Tpttf,
	                       OrbitUtils::Polynomial* Spttf);

};

#endif