//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//   ParticleIdNumber.cc
//
// AUTHOR
//    S. Cousineau
//
// CREATED
//    02/13/2014
//
// DESCRIPTION
//    A subclass of a ParticleAttributes class
//    to keep an integer id number parameter for each particle.
//
//
///////////////////////////////////////////////////////////////////////////

#include "Bunch.hh"
#include "ParticleIdNumber.hh"

ParticleIdNumber::ParticleIdNumber(Bunch* bunch):
  ParticleAttributes(bunch,1)
{
  cl_name_ = "ParticleIdNumber";
  attrDescr = "Particle_Id_Number";
}

ParticleIdNumber::~ParticleIdNumber()
{
}

//int& ParticleIdNumber::getIdNumber(int particle_index){
//  return attValue(particle_index,0);//}



