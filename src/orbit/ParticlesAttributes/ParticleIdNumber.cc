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

/** Returns the Id number for the particle with index. */
int ParticleIdNumber::getIdNumber(int particle_index)
{
	double id = attValue(particle_index,0);
	return int(id);
}

/** Sets the Id number for the particle with index. */
void ParticleIdNumber::setIdNumber(int particle_index, int id)
{
	attValue(particle_index,0) = 1.0*id;
}



