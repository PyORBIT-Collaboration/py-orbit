//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//   Evolution.cc
//
// AUTHOR
//    T. Gorlov
//
// CREATED
//    07/19/2005
//
// DESCRIPTION
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


