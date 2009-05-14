/////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//   prf_time.hh
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

#ifndef PRF_TIME_HH_
#define PRF_TIME_HH_

#include <string>

#include "ParticleAttributes.hh"

class prf_time : public ParticleAttributes
{
public:

	/** This Attribute describe complex coefficients of Wave functions.
	  * User can specify the number of variables that he wants to reserve.
		*/
	prf_time(Bunch* bunch, int size_in);
	
  ~prf_time();

};

#endif /*PRF_TIME_HH_*/
