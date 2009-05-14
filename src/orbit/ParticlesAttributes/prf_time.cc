//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//   prf_time.cc
//
// AUTHOR
//    T. Gorlov
//
// CREATED
//    07/19/2005
//
// DESCRIPTION
//    A subclass of a ParticleAttributes class
//
//
///////////////////////////////////////////////////////////////////////////

#include "Bunch.hh"

#include "prf_time.hh"

prf_time::prf_time(Bunch* bunch, int size_in):
  ParticleAttributes(bunch,size_in)
{
  cl_name_ = "prf_time";
  attrDescr = "time";
}

prf_time::~prf_time()
{
}


