//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//   WaveFunctionAmplitudes.hh
//
// AUTHOR
//    T. Gorlov
//
// CREATED
//    07/14/2005
//
// DESCRIPTION
//    A subclass of the particle attributes class.
//
///////////////////////////////////////////////////////////////////////////
#ifndef WAVE_FUNCTION_AMPLITUDES_H
#define WAVE_FUNCTION_AMPLITUDES_H

#include <string>

#include "ParticleAttributes.hh"

class WaveFunctionAmplitudes : public ParticleAttributes
{
public:
  //--------------------------------------
  //the public methods of the ParticleMacroSize class
  //--------------------------------------

	
	/** This Attribute describe complex coefficients of Wave functions.
	  * User can specify the number of variables that he wants to reserve.
		*/
	WaveFunctionAmplitudes(Bunch* bunch, int size_in);
	
  ~WaveFunctionAmplitudes();
	
};

///////////////////////////////////////////////////////////////////////////
//
// END OF FILE
//
///////////////////////////////////////////////////////////////////////////

#endif   //WAVE_FUNCTION_AMPLITUDES_H definition
