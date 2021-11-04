//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//   InitialCoordsAttrFunctions.hh
//
// AUTHOR
//    A. Shishlo
//
// CREATED
//    03/27/2016
//
// DESCRIPTION
//    A set of functions for bunches with the ParticleInitialCoordinates attribute
//
///////////////////////////////////////////////////////////////////////////

#ifndef PARTICLES_INIT_COORDS_FUNCTIONS_H
#define PARTICLES_INIT_COORDS_FUNCTIONS_H

//ORBIT bunch
#include "Bunch.hh"

#include "Matrix.hh"

#include <cmath>

namespace OrbitUtils{
	
	/** 
	  A function that will copy coordinates of the particles to 
	  the 6D initial coordinates Attribute (ParticleInitialCoordinates).
	 */
	void copyCoordsToInitCoordsAttr(Bunch* bunch);
	
	/** 
	  A function that will swap the initial coordinates Attribute 
	  (ParticleInitialCoordinates) and the 6D coordinates of the particles.
	  If the bunch does not have the ParticleInitialCoordinates particles attributes
	  they will be created and full out with 0s, and then they will be swept.
	 */
	void swapInitCoordsAttrAndCoords(Bunch* bunch);
	
	/** 
	  Calculates transport matrix A: x_out = A*x_in +b. A is a 7x7 matrix. 
		The last column of A is the b vector.
	  It returns 0 if unsuccessful or the size of the statistics otherwise.
	  The function uses the initial coordinates attributes as x_in and the usual
	  coordinates as x_out.
	*/
	int transport_mtrx_from_init_coords(Bunch* bunch, Matrix* A_mtr);
	
	/** A function analyzes two bunches assuming the vectors of coordinates 
	    transformation x_out = A*x_in + b, where x_in the initial coordinates,
			and x_out final. The results are matrix A (7x7) with the last column as a vector b.
	    The function uses the initial coordinates attributes as x_in and the usual
	    coordinates as x_out.
			It returns 0 if unsuccessful or the size of the statistics otherwise.
			For statistics it uses the weights according to Gaussian distribution with
			Twiss parameters.
	*/
	int transport_mtrx_from_init_coords(Bunch* bunch, Matrix* A_mtr, int appl_twiss_x, int appl_twiss_y, int appl_twiss_z);	
	
	/** A function analyzes bunch with (ParticleInitialCoordinates) attributes 
	    Macrosizes of macro-particles in the bunch will be 
	    multiplied by the numbers wx0*wy0*wz0*wx1*wy1*wz1 where
	    wx(In or Out) = exp(-(x^2+(alphax*x+betax*x')^2)/(2*(betax*emittancex))
	    etc.
	    Alpha, beta, emittance are the Twiss parameters for the corresponding 
	    plane.
	*/
	void apply_twiss_weights_for_init_coords(Bunch* bunch,int appl_x,int appl_y,int appl_z);	
	
};
///////////////////////////////////////////////////////////////////////////
//
// END OF FILE
//
///////////////////////////////////////////////////////////////////////////

#endif