
//The base class for Apertures. It defines the interface for Aperture
#ifndef FIELDTRACKER_H
#define FIELDTRACKER_H

//pyORBIT utils
#include "CppPyWrapper.hh"
#include "Bunch.hh"

using namespace std;

/**
 The FieldTracker class is used to define how a particle/bunch  propogates through an 
 3-dimensional magnetic field
 */

class FieldTracker: public OrbitUtils::CppPyWrapper
{
public:
	
	/** FieldTracker */
    FieldTracker(); 
    
	/** Routine for transfering particles through a aperture */
	void trackBunch(Bunch* b);
    
protected:
    

};

//end of FieldTracker_H ifdef
#endif
