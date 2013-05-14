#ifndef BUNCH_TWISS_ANALYSIS_H
#define BUNCH_TWISS_ANALYSIS_H

//pyORBIT utils
#include "CppPyWrapper.hh"

#include "Bunch.hh"

using namespace std;

/** 
  The BunchTwissAnalysis class calculates the average of 6D coordinates and they correlations.
  As results it returns the Twiss parameters for each plane. 
*/

class BunchTwissAnalysis: public OrbitUtils::CppPyWrapper
{
	public:
		
		/** Constructor*/
		BunchTwissAnalysis();
				
		/** Destructor */
		virtual ~BunchTwissAnalysis();
		
		/** Performs the Twiss analysis of the bunch */
		void analyzeBunch(Bunch* bunch);
		
		/** Returns the centered correlation <(x-<x>)*(y-<y>)> = <x*y> - <x>*<y> */
		double getCorrelation(int ic, int jc);
		
		/** Returns the average value for coordinate with index ic */
		double getAverage(int ic);
		
		/** Returns the total number of analysed macroparticles */
		int getGlobalCount();
		
		/** Returns the total macrosize */
		double getGlobalMacrosize();
		
		/** Returns the emittance for index 0,1,2 - x,y,z planes. */
		double getEmittance(int ic);
		
		/** Returns Twiss alpha for index 0,1,2 - x,y,z planes.*/
		double getAlpha(int ic);
		
		/** Returns Twiss beta for index 0,1,2 - x,y,z planes.*/
		double getBeta(int ic);
		
		/** Returns Twiss gamma for index 0,1,2 - x,y,z planes.*/
		double getGamma(int ic);
	
		/** Computes the XY moments of the bunch up to a prescribed order */
		void computeBunchMoments(Bunch* bunch, int order);
	
		/** Returns the XY moment of the beam */
		double getBunchMoment(int i, int j);
		
		
	private:
		
		/** Number of points accounted */
		int count;
		
		/** total macrosize accounted */
		double total_macrosize;
		
		/** array with average values for 6D coordinates */
		double* avg_arr;
		double* avg_arr_MPI;
		
		/** array with correlations between 6D coordinates. 
		    It is a packed 1D array 6*6 = 36. 
				To reach i,j = 0-5 coerrelation use [i+6*j] 
		    It is excessive because <x*y> = <y*x> etc.*/
		double* corr_arr;
		double* corr_arr_MPI;
	
		/** array with XY Moments */
		double** momentXY;
	
		/** order for the moments **/
		double _order;
		
};

#endif
//endif for BUNCH_TWISS_ANALYSIS_H
