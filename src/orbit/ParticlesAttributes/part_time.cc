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

#include "Bunch.hh"
#include "part_time.hh"

part_time::part_time(Bunch* bunch, int size_in):
  ParticleAttributes(bunch,size_in)
{
  cl_name_ = "part_time";
  attrDescr = "time";
}

part_time::~part_time()
{
}



