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
#include "Evolution.hh"

Evolution::Evolution(Bunch* bunch, int size_in):
  ParticleAttributes(bunch,size_in)
{
  cl_name_ = "Evolution";
  attrDescr = "Evol";
}

Evolution::~Evolution()
{
}


