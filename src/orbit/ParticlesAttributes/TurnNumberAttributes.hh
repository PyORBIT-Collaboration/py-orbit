//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//  TurnNumberAttributes.cc
//
// AUTHOR
//   Andrei Shishlo
//
// CREATED
//    03/07/2022
//
// DESCRIPTION
//    A subclass of a ParticleAttributes class
//    
// The main purpose of this class to provide information about the turn
// index in the lost bunch for ring simulations. This attribute should be
// assigned to the tracking main bunch after its instantiating before 
// tracking. During the tracking the collimation nodes will add this 
// attribute to the lost bunch (if it does not exist yet) and will assign
// turn value to the particles lost in this collimator.
//
///////////////////////////////////////////////////////////////////////////

#ifndef TURN_NUMBER_H
#define TURN_NUMBER_H

#include <string>

#include "ParticleAttributes.hh"

class TurnNumberAttributes : public ParticleAttributes
{
public:
  //--------------------------------------
  //the public methods of the TurnNumberAttributes class
  //--------------------------------------

  TurnNumberAttributes(Bunch* bunch);
  ~TurnNumberAttributes();

	/** Returns the turn number for the particle with particle_index in the bunch */
  int getTurn(int particle_index);
	
	/** Sets the turn number for the particle with particle_index in the bunch */
	void setTurn(int particle_index, int turn);	

};

///////////////////////////////////////////////////////////////////////////
//
// END OF FILE
//
///////////////////////////////////////////////////////////////////////////

#endif
