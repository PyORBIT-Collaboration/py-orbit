#ifndef BUNCH_TWISS_ANALYSIS_1D_H
#define BUNCH_TWISS_ANALYSIS_1D_H

#include "TwissAnalysis1D.hh"

//pyORBIT utils
#include "CppPyWrapper.hh"

using namespace std;

/** 
The Moments1D class calculates the arbitrary moments of the (u,up) distribution.
It is used by other classes to calculate Twiss paraemeters etc.
*/

class BunchTwissAnalysis: public OrbitUtils::CppPyWrapper
{
	public:
		
		/** Constructor*/
		BunchTwissAnalysis();
				
		/** Destructor */
		virtual ~BunchTwissAnalysis();
		
};

#endif
//endif for BUNCH_TWISS_ANALYSIS_1D_H
