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
#include "AtomPopulations.hh"

AtomPopulations::AtomPopulations(Bunch* bunch, int size_in):
  ParticleAttributes(bunch,size_in)
{
  cl_name_ = "Populations";
  attrDescr = "Pops";
}

AtomPopulations::~AtomPopulations()
{
}

