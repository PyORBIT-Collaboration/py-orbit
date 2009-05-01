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

#include "pq_coordinates.hh"
///////////////////////////////////////////////////////////////////////////
//   Constructor and Desctructor
///////////////////////////////////////////////////////////////////////////

pq_coordinates::pq_coordinates(Bunch* bunch):
  ParticleAttributes(bunch)
{
  cl_name_ = "pq_coords";
  attrDescr = "Coords";

}

pq_coordinates::pq_coordinates(Bunch* bunch, int size_in):
  ParticleAttributes(bunch)
{
  cl_name_ = "pq_coords";
  attrDescr = "Coords";
	size = size_in;
}

pq_coordinates::~pq_coordinates()
{
}



//zero element is phase and the last one is the time of i-th particle in the particle frame
int pq_coordinates::getAttSize(){
  return size;
}

