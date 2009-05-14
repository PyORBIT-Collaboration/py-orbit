//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//   pq_coordinates.cc
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
#include "pq_coordinates.hh"

pq_coordinates::pq_coordinates(Bunch* bunch, int size_in):
  ParticleAttributes(bunch,size_in)
{
  cl_name_ = "pq_coords";
  attrDescr = "Coords";
}

pq_coordinates::~pq_coordinates()
{
}

