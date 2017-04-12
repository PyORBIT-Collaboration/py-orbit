//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//   ParticleInitialCoordinates.cc
//
// AUTHOR
//    A. Shishlo
//
// CREATED
//    03/27/2017
//
// DESCRIPTION
//    A subclass of a ParticleAttributes class that keeps the initial
//    coordinates of the particles. This attributes will be used for 
//    different purposes like:
//    1. calculations of the acceptance phase space region for particles that
//       reached a particular point at the lattice. The initial coordinates 
//       should be put into these attributes before tracking.
//    2. Calculation of the transport matrix between the origin (place were 
//       initial coordinates were saved in these attributes) and the another
//       point in the lattice.
//
//
///////////////////////////////////////////////////////////////////////////

#include "Bunch.hh"
#include "ParticleInitialCoordinates.hh"

ParticleInitialCoordinates::ParticleInitialCoordinates(Bunch* bunch):
  ParticleAttributes(bunch,6)
{
  cl_name_ = "ParticleInitialCoordinates";
  attrDescr = "Initial Coordinates";
}

ParticleInitialCoordinates::~ParticleInitialCoordinates()
{
}
