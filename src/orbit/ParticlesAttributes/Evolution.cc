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

#include "Evolution.hh"
///////////////////////////////////////////////////////////////////////////
//   Constructor and Desctructor
///////////////////////////////////////////////////////////////////////////

Evolution::Evolution(Bunch* bunch):
  ParticleAttributes(bunch,0)
{
  cl_name_ = "Evolution";
  attrDescr = "Evol";

}

Evolution::Evolution(Bunch* bunch, int size_in):
  ParticleAttributes(bunch,size_in)
{
  cl_name_ = "Evolution";
  attrDescr = "Evol";
	size = size_in;
}

Evolution::~Evolution()
{
}



//zero element is phase and the last one is the time of i-th particle in the particle frame
int Evolution::getAttSize(){
  return size;
}

