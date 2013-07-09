/////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//   ParticlePhaseAttributes.hh
//
// AUTHOR
//    S. Cousineau
//
// CREATED
//    05/20/2013
//
// DESCRIPTION
//    A subclass of the particle attributes class. 
//
///////////////////////////////////////////////////////////////////////////
#ifndef PARTICLEPHASEATTRIBUTES_HH_
#define PARTICLEPHASEATTRIBUTES_HH_

///////////////////////////////////////////////////////////////////////////
//
// INCLUDE FILES
//
///////////////////////////////////////////////////////////////////////////
#include "ParticleAttributes.hh"

class ParticlePhaseAttributes : public ParticleAttributes
{
public:
	
	/** This Attribute class contains additional properties of particle phases.
	*/
	ParticlePhaseAttributes(Bunch* bunch);
	~ParticlePhaseAttributes();
	double& getLastPhaseX(int particle_index);
	double& getLastPhaseY(int particle_index);
	double& getLastTuneX(int particle_index);
	double& getLastTuneY(int particle_index);
	double& getLastActionX(int particle_index);
	double& getLastActionY(int particle_index);
};

///////////////////////////////////////////////////////////////////////////
//
// END OF FILE
//
///////////////////////////////////////////////////////////////////////////


#endif
