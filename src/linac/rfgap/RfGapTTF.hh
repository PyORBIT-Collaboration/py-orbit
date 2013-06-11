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
	void trackBunch(Bunch* bunch, double E0, double phase);	
			
	/** 
	    Sets up the gap parameters: T,S, minimal and maximal beta, 
			rf frequency, the gap length,  and the relative amplitude.
	*/
	void setParameters(OrbitUtils::Polynomial* Tttf, OrbitUtils::Polynomial* Tpttf, 
		                 OrbitUtils::Polynomial* Sttf, OrbitUtils::Polynomial* Spttf, 
		                 double beta_min, double beta_max, 
										  double rf_frequency, double gap_length, 
											double relative_amplitude);

  /** Returns the minimal beta. */
	double getBetaMin();
	
  /** Returns the maximal beta. */
	double getBetaMax();
	
	/** Returns T TTF. */
	OrbitUtils::Polynomial* getT_TTF();
	
	/** Returns S TTF. */
	OrbitUtils::Polynomial* getS_TTF();
	
	/** Returns Tp TTF. */
	OrbitUtils::Polynomial* getTp_TTF();
	
	/** Returns Sp TTF. */
	OrbitUtils::Polynomial* getSp_TTF();	
	
  /** Sets the Polynomial to T TTF. */
  void setT_TTF(OrbitUtils::Polynomial* Tttf);
	
  /** Sets the Polynomial to S TTF. */
  void setS_TTF(OrbitUtils::Polynomial* Sttf);
	
  /** Sets the Polynomial to Tp TTF. */
  void setTp_TTF(OrbitUtils::Polynomial* Tpttf);
	
  /** Sets the Polynomial to Sp TTF. */
  void setSp_TTF(OrbitUtils::Polynomial* Spttf);	
	
	/** Returns RF frequency. */
	double getFrequency();
	
	/** Returns the gap length. */
	double getLength();
	
	/** Returns the realtive amplitude. */
	double getRelativeAmplitude();
	
	/** polynomials for T,S TTF as functions of kappa = 2*pi*rf_freq/(c_light*beta) */
	OrbitUtils::Polynomial* Tttf;
	OrbitUtils::Polynomial* Sttf;
	OrbitUtils::Polynomial* Tpttf;
	OrbitUtils::Polynomial* Spttf;	
	
  private:
		
		//the limits where the polynomials are defined
		double beta_min,beta_max;
		
		//this value will be defined from the external code in the method trackBunch(...),
		//but we also need it for auxilary purposes. The units are Hz.
		double rf_frequency;
		
		//gap_length is a gap parameter L to substitute into the E0L expression
		double gap_length;
		
		//the parameter defines the amplitude of the gap relative 
		//to the others gaps in the RF cavity E0 = A*relative_amplitude
		double relative_amplitude;
		
};

#endif
