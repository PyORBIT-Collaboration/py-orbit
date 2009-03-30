//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//   MyAttr.cc
//
// AUTHOR
//    T. Gorlov
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

///////////////////////////////////////////////////////////////////////////
//
// INCLUDE FILES
//
///////////////////////////////////////////////////////////////////////////
#include "Bunch.hh"

#include "AtomPopulations.hh"
///////////////////////////////////////////////////////////////////////////
//   Constructor and Desctructor
///////////////////////////////////////////////////////////////////////////

AtomPopulations::AtomPopulations(Bunch* bunch):
  ParticleAttributes(bunch)
{
  cl_name_ = "Populations";
  attrDescr = "Pops";

}

AtomPopulations::AtomPopulations(Bunch* bunch, int size_in):
  ParticleAttributes(bunch)
{
  cl_name_ = "Populations";
  attrDescr = "Pops";
	size = size_in;
}

AtomPopulations::~AtomPopulations()
{
}



//zero element is phase and the last one is the time of i-th particle in the particle frame
int AtomPopulations::getAttSize(){
  return size;
}

