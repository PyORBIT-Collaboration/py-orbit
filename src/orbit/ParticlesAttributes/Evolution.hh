#ifndef EVOLUTION_HH_
#define EVOLUTION_HH_


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

///////////////////////////////////////////////////////////////////////////
//
// INCLUDE FILES
//
///////////////////////////////////////////////////////////////////////////

#include <string>

///////////////////////////////////////////////////////////////////////////
//
// CLASS NAME
//     WaveFunctionAmplitudes
//
///////////////////////////////////////////////////////////////////////////
#include "ParticleAttributes.hh"

class Evolution : public ParticleAttributes
{
public:
  //--------------------------------------
  //the public methods of the ParticleMacroSize class
  //--------------------------------------
	
	/** Constructor. This Attribute describe complex coefficients of Wave functions.
	  * The defailt size is 400. 
		*/
	Evolution(Bunch* bunch);
	
	/** This Attribute describe complex coefficients of Wave functions.
	  * User can specify the number of variables that he wants to reserve.
		*/
	Evolution(Bunch* bunch, int size_in);
	
  ~Evolution();
  
	
  int getAttSize();
	
private:
	int size;
	
	
};

///////////////////////////////////////////////////////////////////////////
//
// END OF FILE
//
///////////////////////////////////////////////////////////////////////////





#endif /*EVOLUTION_HH_*/
