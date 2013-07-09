#ifndef BUNCH_TUNE_ANALYSIS_H
#define BUNCH_TUNE_ANALYSIS_H

//pyORBIT utils
#include "CppPyWrapper.hh"

#include "Bunch.hh"
#include "BunchTwissAnalysis.hh"

using namespace std;

/** 
  The BunchTuneAnalysis class calculates the particle tunes
*/

class BunchTuneAnalysis: public OrbitUtils::CppPyWrapper
{
	public:
		
		/** Constructor*/
		BunchTuneAnalysis();
				
		/** Destructor */
		virtual ~BunchTuneAnalysis();
		
		/** Performs the Twiss analysis of the bunch */
		void analyzeBunch(Bunch* bunch);
	
		//** Assigns Twiss values at location of calculator */
		void assignTwiss(double bx, double ax, double dx, double dpx, double by, double ay);
		
		/** Returns the average value for coordinate with index ic */
		double getTune(int ic);
				
		
	private:
		//** Twiss */
		double betax;
		double alphax;
		double etax;
		double etapx;
		double betay;
		double alphay;
				
};

#endif
//endif for BUNCH_TUNE_ANALYSIS_H
