#ifndef ENERGY_APERTURE_H
#define ENERGY_APERTURE_H
///////////////////////////////////////////////////////////////////////////
//
// NAME
//
//   EnergyAperture::EnergyAperture
//
// Author: A. Shishlo
// Date: 03/29/2017
//
// DESCRIPTION
//   Constructs an aperture for the energy (relative to synch. part.) variable.
//
// PARAMETERS
//   The class will remove particles (and put them into lost_bunch) from 
//   the bunch whose energy (relative to synch. part.) is outside the specified min and 
//   max values. This class mostly is intended to be used in linacs. 
//   The energy is in GeV.
//
///////////////////////////////////////////////////////////////////////////

//pyORBIT utils
#include "CppPyWrapper.hh"
#include "Bunch.hh"

using namespace std;

class EnergyAperture: public OrbitUtils::CppPyWrapper
{
  public:
	
  /**
  The EnergyAperture class constructor.
  */	
  EnergyAperture();
  
  /**
  The method sets the energy limits in GeV for particles in the bunch.
  */
  void setEnergyLimits(double minEnergy, double maxEnergy);
  
  /**
  Returns the min limit of the particles' energy.
  */
  double getMinEnergy();
  
  /**
  Returns the max limit of the particles' energy.
  */
  double getMaxEnergy();
  
  /**
  Sets the position of the Energy Aperture node.
  */
  void setPosition(double position);
  
  /**
  Returns the position of the Energy Aperture node.
  */
  double getPosition();
  
  /**
  Removes particles with phases outside the energy limits and puts them
  into the lost bunch is it exists.
  */
  void checkBunch(Bunch* bunch, Bunch* lostbunch);
  
  
  protected:
  	//EnergyAperture parameters
  	double minEnergy_;
  	double maxEnergy_;
  	double pos_;
  	
};

//end of ENERGY_APERTURE_H ifdef
#endif

