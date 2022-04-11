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

#include "Bunch.hh"
#include "TurnNumberAttributes.hh"

TurnNumberAttributes::TurnNumberAttributes(Bunch* bunch): 
ParticleAttributes(bunch,1)
{
  cl_name_ = "TurnNumber";
  attrDescr = "TurnNumber";
  std::string attr_name_str("TurnNumber");
  //we start counting from the 1st turn
  bunch->setBunchAttribute( attr_name_str, 1);
}

TurnNumberAttributes::~TurnNumberAttributes()
{
}

/**
	Returns the turn number for the particle with particle_index in the bunch.
*/
int TurnNumberAttributes::getTurn(int particle_index){

	return int(attValue(particle_index,0));
}

/**
	Sets the turn number for the particle with particle_index in the bunch.
*/
void TurnNumberAttributes::setTurn(int particle_index, int turn){
	attValue(particle_index,0) = 1.0*turn;
}