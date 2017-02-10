/**
   This class represents a RF gap where the transit time factors are calculated
   by using the second order polynomial approximation of the field on the axis.
   The model includes non-linearity in transverse direction.
   The field is defined by three points at positions -dz, 0., and +dz.
   The transformation coefficients are defined for each particle separately,
   so the name of the class has the word slow in the name.
   
*/

#ifndef THREE_POINTS_RF_GAP_SLOW_H
#define THREE_POINTS_RF_GAP_SLOW_H

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
  This class represents a RF gap as a Three Points type gap. 
*/

class RfGapThreePointTTF_slow: public OrbitUtils::CppPyWrapper
{
public:
	
	/** Constructor for Three-Point type RF gap with TTF */
  RfGapThreePointTTF_slow();
	
  /** Destructor */
  virtual ~RfGapThreePointTTF_slow();
	
	/** Tracks the Bunch through the Three-Point RF gap. */	
	static void trackBunch(Bunch* bunch, double dz, double Em, double E0, double Ep, double rf_frequency, double phase);	
	
	/** 
	It calculates the symmetrical TTF for 3-point approximation of the field. 
	This TTF as functions of the kappa variable = 2*pi*f/(c*beta).
  */
	static double Tttf(double dz, double a, double b, double kappa);
	
	/** 
	It calculates the asymmetrical TTF for 3-point approximation of the field. 
	This TTF as functions of the kappa variable = 2*pi*f/(c*beta).
	*/	
	static double Sttf(double dz, double a, double b, double kappa);
	
	/** 
	It calculates the derivative of the symmetrical TTF for 3-point approximation of the field. 
	This TTF as functions of the kappa variable = 2*pi*f/(c*beta).
	*/	
	static double Tpttf(double dz, double a, double b, double kappa);
	
	/** 
	It calculates the derivative of the asymmetrical TTF for 3-point approximation of the field. 
	This TTF as functions of the kappa variable = 2*pi*f/(c*beta).
	*/
	static double Spttf(double dz, double a, double b, double kappa);
	
  private:
		
};

#endif
