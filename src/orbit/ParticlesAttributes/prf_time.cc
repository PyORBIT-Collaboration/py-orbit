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

#include "prf_time.hh"
///////////////////////////////////////////////////////////////////////////
//   Constructor and Desctructor
///////////////////////////////////////////////////////////////////////////

prf_time::prf_time(Bunch* bunch):
  ParticleAttributes(bunch)
{
  cl_name_ = "prf_time";
  attrDescr = "time";

}

prf_time::prf_time(Bunch* bunch, int size_in):
  ParticleAttributes(bunch)
{
  cl_name_ = "prf_time";
  attrDescr = "time";
	size = size_in;
}

prf_time::~prf_time()
{
}



//zero element is phase and the last one is the time of i-th particle in the particle frame
int prf_time::getAttSize(){
  return size;
}

