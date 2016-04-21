//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//   TwissFilteringFunctions.hh
//
// AUTHOR
//    A. Shishlo
//
// CREATED
//    04/20/2016
//
// DESCRIPTION
//    A set of functions for filtering macro-particles from the bunch
//    according to their positions relative to the emittance phase space center 
//
///////////////////////////////////////////////////////////////////////////

#ifndef BUNCH_TWISS_FILTERING_FUNCTIONS_H
#define BUNCH_TWISS_FILTERING_FUNCTIONS_H

#include <cmath>

//ORBIT bunch
#include "Bunch.hh"

namespace OrbitUtils{
	
	/** 
	    Bunch filtering according to the Twiss parameters.
	    bunch_in is the input bunch, and bunch_bad collects
	    the macro-particles that are exluded from the initial bunch.
	    The coefficients are from 0 to infinity and define the limit for 
	    the ratio
	    (gamma*x^2+2*alpha*x*xp+beta*xp^2)/(2*emittance)
	    The function returns the total number of removed particles. 
	*/
	int bunch_twiss_filtering(Bunch* bunch_in,Bunch* bunch_bad, double coeff_x, double coeff_y, double coeff_z);
		
};
///////////////////////////////////////////////////////////////////////////
//
// END OF FILE
//
///////////////////////////////////////////////////////////////////////////

#endif