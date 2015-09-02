//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//   ParticlesWithIdFunctions.hh
//
// AUTHOR
//    A. Shishlo
//
// CREATED
//    09/01/2015
//
// DESCRIPTION
//    A set of functions for bunches with the ParticleIdNumber attribute
//
///////////////////////////////////////////////////////////////////////////

#ifndef PARTICLES_WITH_ID_FUNCTIONS_H
#define PARTICLES_WITH_ID_FUNCTIONS_H

#include "ParticlesWithIdFunctions.hh"

//ORBIT bunch
#include "Bunch.hh"

namespace OrbitUtils{
	
	/** A function that will sort bunch according to Id.*/
	void bunch_sort_id(Bunch* bunch);
	
};
///////////////////////////////////////////////////////////////////////////
//
// END OF FILE
//
///////////////////////////////////////////////////////////////////////////

#endif
