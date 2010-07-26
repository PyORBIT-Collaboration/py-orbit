#ifndef TWISS_ANALYSIS_1D_H
#define TWISS_ANALYSIS_1D_H

#include "StatMoments2D.hh"

//pyORBIT utils
#include "CppPyWrapper.hh"

using namespace std;

/** 
The Moments1D class calculates the arbitrary moments of the (u,up) distribution.
It is used by other classes to calculate Twiss paraemeters etc.
*/

namespace OrbitUtils{ 
	
	class TwissAnalysis1D: public StatMoments2D
	{
	public:
		
		/** Constructor with max order = 2 by default */
		TwissAnalysis1D();
		
		/** Destructor */
		virtual ~TwissAnalysis1D();
		
    /** Returns the emittance */
		double getEmittance();
		
    /** Returns Twiss alpha */
		double getAlpha();
		
    /** Returns Twiss beta */
		double getBeta();
		
    /** Returns Twiss gamma */
		double getGamma();		
		
		/** Returns the rms value of u */ 	
		double getRmsU();
		
		/** Returns the rms value of up */ 	
		double getRmsUP();			
		
	};
	
} //end of OrbitUtils{}

#endif
//endif for TWISS_ANALYSIS_1D_H
