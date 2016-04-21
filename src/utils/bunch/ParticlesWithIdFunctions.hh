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

#include <cmath>

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
	
	/** A function analyzes two bunches assuming the vectors of coordinates 
	    transformation x_out = A*x_in + b, where x_in the initial coordinates,
			and x_out final. The results are matrix A (7x7) with the last column as a vector b.
			We assume that bunch_in and bunch_out are sorted according to Id PartAttr, and
			bunch_out has less particles than bunch_in (some of them could be lost).
			It returns 0 if unsuccessful or the size of the statistics otherwise.
			For statistics it uses the weights according to Gaussian distribution with
			Twiss parameters.
	*/
	int transport_mtrx(Bunch* bunch_in, Bunch* bunch_out, Matrix* A_mtr, int appl_twiss_x, int appl_twiss_y, int appl_twiss_z);	
	
	/** A function analyzes two bunches assuming that they are already
	    sorted and synchronized according to the macro-particles Id. 
	    Coordinates of macro-particles in "in" and "out" bunches will be 
	    multiplied by the same numbers wx*wy*wz where
	    wx = exp(-(x^2+(alphax*x+betax*x')^2)/(2*(betax*emittancex))
	    etc.
	    Alpha, beta, emittance are the Twiss parameters for the corresponding 
	    plane.
	*/
	void apply_twiss_weghts(Bunch* bunch_in, Bunch* bunch_out,int appl_x,int appl_y,int appl_z);	
	
};
///////////////////////////////////////////////////////////////////////////
//
// END OF FILE
//
///////////////////////////////////////////////////////////////////////////

#endif
