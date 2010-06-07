//This class repersents a simple RF gap. For this RF gap we know the E0TL parameter only.

#include "BaseRfGap.hh"
#include "ParticleMacroSize.hh"

#include <iostream>

using namespace OrbitUtils;

// Constructor
BaseRfGap::BaseRfGap(): CppPyWrapper(NULL)
{
}

// Destructor
BaseRfGap::~BaseRfGap()
{
}

/** Tracks the Bunch trough the RF gap. */	
void BaseRfGap::trackBunch(Bunch* bunch, double frequency, double E0TL, double phase){
	bunch->compress();
}


