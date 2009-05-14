//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//   ParticleMacroSize.cc
//
// AUTHOR
//    A. Shishlo
//
// CREATED
//    07/19/2005
//
// DESCRIPTION
//    A subclass of a ParticleAttributes class
//    to keep macrosize parameter for each particle.
//
//
///////////////////////////////////////////////////////////////////////////

#include "Bunch.hh"
#include "ParticleMacroSize.hh"

ParticleMacroSize::ParticleMacroSize(Bunch* bunch):
  ParticleAttributes(bunch,1)
{
  cl_name_ = "macrosize";
  attrDescr = "macro_size";
}

ParticleMacroSize::~ParticleMacroSize()
{
}

double& ParticleMacroSize::macrosize(int particle_index){
  return attValue(particle_index,0);
}



