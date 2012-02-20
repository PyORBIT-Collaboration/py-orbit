//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//   LostParticleAttributes.cc
//
// AUTHOR
//   S. Cousineau
//
// CREATED
//    11/08/2011
//
// DESCRIPTION
//    A subclass of a ParticleAttributes class
//
//
///////////////////////////////////////////////////////////////////////////

#include "Bunch.hh"
#include "LostParticleAttributes.hh"

LostParticleAttributes::LostParticleAttributes(Bunch* bunch): 
ParticleAttributes(bunch,1)
{
  cl_name_ = "LostParticleAttributes";
  attrDescr = "Pos (m)";

}

LostParticleAttributes::~LostParticleAttributes()
{
}

double& LostParticleAttributes::getPosition(int particle_index){
	return attValue(particle_index,0);
}
