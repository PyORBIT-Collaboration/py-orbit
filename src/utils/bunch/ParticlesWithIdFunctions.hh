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

#include "Matrix.hh"

namespace OrbitUtils{
	
	/** A function that will sort bunch according to Id.*/
	void bunch_sort_id(Bunch* bunch);
	
	/** 
	  Calculates transport matrix A: x_out = A*x_in +b. A is a 7x7 matrix. 
		The last column of A is the b vector.
	  It returns 0 if unsuccessful or the size of the statistics otherwise.
	*/
	int transport_mtrx(Bunch* bunch_in, Bunch* bunch_out, Matrix* A_mtr);
	
};
///////////////////////////////////////////////////////////////////////////
//
// END OF FILE
//
///////////////////////////////////////////////////////////////////////////

#endif
