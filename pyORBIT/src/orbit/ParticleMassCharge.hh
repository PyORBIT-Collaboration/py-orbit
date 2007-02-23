//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//   ParticleMassCharge.hh
//
// AUTHOR
//    M. Perkett
//
// CREATED
//    09/26/2006
//
// DESCRIPTION
//    A subclass of the particle attributes class. This is for use with
//    tracking particles through the injection area of the accumulator
//    ring.  Particles are initially set to mass and charge of proton in
//    init(int particle_index).  Mass is to be stored at index 0 and charge
//    at index 1.  You can either directly change/retrieve the mass and
//    charge through these indeces or simply use one of the following
//    functions: getMass, getCharge, setMass, setCharge.
//
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
//
// INCLUDE FILES
//
///////////////////////////////////////////////////////////////////////////
#ifndef _PARTICLE_MASSCHARGE_H
#define _PARTICLE_MASSCHARGE_H

///////////////////////////////////////////////////////////////////////////
//
// INCLUDE FILES
//
///////////////////////////////////////////////////////////////////////////
#include <string>

///////////////////////////////////////////////////////////////////////////
//
// CLASS NAME
//    ParticleMassCharge
//
///////////////////////////////////////////////////////////////////////////
#include "ParticleAttributes.hh"

class ParticleMassCharge : public ParticleAttributes
{
public:
  //--------------------------------------
  //the public methods of the ParticleMassCharge class
  //--------------------------------------

  ParticleMassCharge(Bunch* bunch);
  ~ParticleMassCharge();

  int getAttSize();
  double getMass(int particle_index);
  double getCharge(int particle_index);
  void setMass(int particle_index, double mass);
  void setCharge(int particle_index, double charge);
  
protected:
  void init(int particle_index);
};

///////////////////////////////////////////////////////////////////////////
//
// END OF FILE
//
///////////////////////////////////////////////////////////////////////////

#endif
