//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//   ParticlePhaseAttributes.cc
//
// AUTHOR
//   S. Cousineau
//
// CREATED
//    05/20/2013
//
// DESCRIPTION
//    A subclass of a ParticlePhaseAttributes class
//
//
///////////////////////////////////////////////////////////////////////////

#include "Bunch.hh"
#include "ParticlePhaseAttributes.hh"

ParticlePhaseAttributes::ParticlePhaseAttributes(Bunch* bunch):
ParticleAttributes(bunch,6)
{
  cl_name_ = "ParticlePhaseAttributes";
  attrDescr = "PhaseVars";
}

ParticlePhaseAttributes::~ParticlePhaseAttributes()
{
}

double& ParticlePhaseAttributes::getLastPhaseX(int particle_index){
	return attValue(particle_index,0);
}

double& ParticlePhaseAttributes::getLastPhaseY(int particle_index){
	return attValue(particle_index,1);
}

double& ParticlePhaseAttributes::getLastTuneX(int particle_index){
	return attValue(particle_index,2);
}

double& ParticlePhaseAttributes::getLastTuneY(int particle_index){
	return attValue(particle_index,3);
}

double& ParticlePhaseAttributes::getLastActionX(int particle_index){
	return attValue(particle_index,4);
}

double& ParticlePhaseAttributes::getLastActionY(int particle_index){
	return attValue(particle_index,5);
}