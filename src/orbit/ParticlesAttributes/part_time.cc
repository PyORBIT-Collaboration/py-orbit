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

#include "part_time.hh"
///////////////////////////////////////////////////////////////////////////
//   Constructor and Desctructor
///////////////////////////////////////////////////////////////////////////

part_time::part_time(Bunch* bunch):
  ParticleAttributes(bunch)
{
  cl_name_ = "part_time";
  attrDescr = "time";

}

part_time::part_time(Bunch* bunch, int size_in):
  ParticleAttributes(bunch)
{
  cl_name_ = "part_time";
  attrDescr = "time";
	size = size_in;
}

part_time::~part_time()
{
}



//zero element is phase and the last one is the time of i-th particle in the particle frame
int part_time::getAttSize(){
  return size;
}

