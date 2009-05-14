
/////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//   MyAttr.hh
//
// AUTHOR
//    T. Gorlov
//
// CREATED
//    07/14/2005
//
// DESCRIPTION
//    A subclass of the particle attributes class. This is container
//    for a macrosize of macro-particles in the bunch.
//
///////////////////////////////////////////////////////////////////////////

#ifndef POPULATIONS_HH_
#define POPULATIONS_HH_

#include <string>

#include "ParticleAttributes.hh"

class AtomPopulations : public ParticleAttributes
{
public:

/** This Attribute describe complex coefficients of Wave functions.
	  * User can specify the number of variables that he wants to reserve.
		*/
	AtomPopulations(Bunch* bunch, int size_in);
	
  ~AtomPopulations();	
	
};

///////////////////////////////////////////////////////////////////////////
//
// END OF FILE
//
///////////////////////////////////////////////////////////////////////////


#endif /*POPULATIONS_HH_*/
