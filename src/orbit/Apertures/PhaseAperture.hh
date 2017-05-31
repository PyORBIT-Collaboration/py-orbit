#ifndef PHASE_APERTURE_H
#define PHASE_APERTURE_H
///////////////////////////////////////////////////////////////////////////
//
// NAME
//
//   PhaseAperture::PhaseAperture
//
// Author: A. Shishlo
// Date: 03/29/2017
//
// DESCRIPTION
//   Constructs an aperture for the longitudinal direction.
//
// PARAMETERS
//   The class will remove particles (and put them into lost_bunch) from 
//   the bunch whose longitudinal phase is outside the specified min and 
//   max phases. This class mostly is intended to be used in linacs. 
//   The user have to specify the RF frequency to translate the longitudinal 
//   position from meters to the phase in degrees.
//
///////////////////////////////////////////////////////////////////////////

//pyORBIT utils
#include "CppPyWrapper.hh"
#include "Bunch.hh"

using namespace std;

class PhaseAperture: public OrbitUtils::CppPyWrapper
{
  public:
	
  /**
  The PhaseAperture class constructor. It needs the RF frequency in Hz to translate
  from z coordinate in meters to phase in degrees.
  */	
  PhaseAperture(double frequency);
  
  /**
  The method sets the phase limits particles in the bunch. Phase limits are in degrees.
  */
  void setPhaseLimits(double minPhase, double maxPhase);
  
  /**
  Returns the min limit of the particles' phases. Phase is in degrees.
  */
  double getMinPhase();
  
  /**
  Returns the max limit of the particles' phases. Phase is in degrees.
  */
  double getMaxPhase();
  
  /**
  Sets the position of the Phase Aperture node.
  */
  void setPosition(double position);
  
  /**
  Returns the position of the Phase Aperture node.
  */
  double getPosition();
  
  /**
  Sets the RF frequency of the Phase Aperture node in Hz.
  */
  void setRfFrequency(double frequency);
  
  /**
  Returns the RF frequency of the Phase Aperture node in Hz.
  */
  double getRfFrequency();
  
  /**
  Removes particles with phases outside the phase limits and puts them
  into the lost bunch is it exists.
  */  
  void checkBunch(Bunch* bunch, Bunch* lostbunch);
  
  
  protected:
  	//PhaseAperture parameters
  	double frequency_;
  	double minPhase_;
  	double maxPhase_;
  	double pos_;
  	
};

//end of PHASE_APERTURE_H ifdef
#endif

