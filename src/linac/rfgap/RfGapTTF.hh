/**
   This class represents a Parmila type RF gap. It acts on the coordinates 
   of the particle by using the transit time factors. The model includes  
   non-linearity in transverse direction.
*/

#ifndef TTF_RF_GAP_H
#define TTF_RF_GAP_H

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

class RfGapTTF: public OrbitUtils::CppPyWrapper
{
public:
	
	/** Constructor for Parmila's type RF gap with TTF */
  RfGapTTF();
	
  /** Destructor */
  virtual ~RfGapTTF();
	
	/** Tracks the Bunch through the RF gap. */	
	static void trackBunch(Bunch* bunch, double frequency, double E0L, double phase,
	                       OrbitUtils::Polynomial* Tttf,
	                       OrbitUtils::Polynomial* Sttf,
	                       OrbitUtils::Polynomial* Tpttf,
	                       OrbitUtils::Polynomial* Spttf);

};

#endif
