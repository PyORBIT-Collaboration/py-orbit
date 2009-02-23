//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//   ParticleAttributes.cc
//
// AUTHOR
//    A. Shishlo
//
// CREATED
//    07/14/2005
//
// DESCRIPTION
//    An abstract class for particle attributes. This is a base class for the
//    others classes that keep different attributes, e.g. spin orientation, tunes,
//    quantum amplitudes of exited states etc.
//
//
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
//
// INCLUDE FILES
//
///////////////////////////////////////////////////////////////////////////
#include "Bunch.hh"

#include "ParticleAttributes.hh"

///////////////////////////////////////////////////////////////////////////
//   Constructor and Desctructor
///////////////////////////////////////////////////////////////////////////

ParticleAttributes::ParticleAttributes(Bunch* bunch)
{
  attr_ind_shift_ = 0;
  cl_name_ = "empty";
  bunch_ = bunch;
  attrDescr = "no attributes";
}

ParticleAttributes::~ParticleAttributes()
{
}

const std::string& ParticleAttributes::name(){ return cl_name_;}

const std::string& ParticleAttributes::attrDescription(){ return attrDescr;}

const Bunch* ParticleAttributes::bunch(){ return bunch_;}

double& ParticleAttributes::attValue(int particle_index, int att_index){
  return bunch_->getParticleAttributeVal(particle_index, attr_ind_shift_ + att_index);
}

double* ParticleAttributes::attArr(int particle_index){
  return &bunch_->getParticleAttributeVal(particle_index, attr_ind_shift_);
}

int ParticleAttributes::getAttSize(){
  return 0;
}

int ParticleAttributes::getSize(){
  return bunch_->getSize();
}

int ParticleAttributes::flag(int particle_index){
  return bunch_->flag(particle_index);
}

//-------------------------------------------------------
// Protected members
//-------------------------------------------------------

void ParticleAttributes::init(int particle_index){
  for(int i = 0, n = getAttSize(); i < n; i++){
    attValue(particle_index,i) = 0.0;
  }
}

int ParticleAttributes::getAttrShift(){
  return attr_ind_shift_;
}

void ParticleAttributes::setAttrShift(int attr_ind_shift){
  attr_ind_shift_ = attr_ind_shift;
}
