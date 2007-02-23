//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//   ParticleMassCharge.cc
//
// AUTHOR
//    M. Perkett
//
// CREATED
//    09/27/2006
//
// DESCRIPTION
//    
//
//
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
//
// INCLUDE FILES
//
///////////////////////////////////////////////////////////////////////////

#include "Bunch.hh"
#include "OrbitConst.hh"
#include "ParticleMassCharge.hh"
///////////////////////////////////////////////////////////////////////////
//   Constructor and Desctructor
///////////////////////////////////////////////////////////////////////////

ParticleMassCharge::ParticleMassCharge(Bunch* bunch):
  ParticleAttributes(bunch)
{
  cl_name_ = "masscharge";
  attrDescr = "mass_charge";
}


ParticleMassCharge::~ParticleMassCharge()
{
}


int ParticleMassCharge::getAttSize(){
  return 2;
}


void ParticleMassCharge::init(int particle_index)
{
   // initialize to proton charge and mass
   attArr(particle_index)[0] = OrbitConst::mass_proton;           // mass (MeV)
   attArr(particle_index)[1] = OrbitConst::charge_proton*
                               OrbitConst::elementary_charge_CGS; // charge (esu?)
}


double ParticleMassCharge::getMass(int particle_index)
{
   return attValue(particle_index,0);
}


double ParticleMassCharge::getCharge(int particle_index)
{
   return attValue(particle_index,1);
}


void ParticleMassCharge::setMass(int particle_index, double mass)
{
   attArr(particle_index)[0] = mass;
}


void ParticleMassCharge::setCharge(int particle_index, double charge)
{
   attArr(particle_index)[1] = charge;
}

