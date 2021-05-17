//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//   BunchExtremaCalculator.hh
//
// AUTHOR
//    A. Shishlo
//
// CREATED
//    10/11/2010
//
// DESCRIPTION
//    A class calculates the extrema of the particles coordinates in the bunch.
//
///////////////////////////////////////////////////////////////////////////

#ifndef BUNCH_EXTREMA_CALCULATIONS_H
#define BUNCH_EXTREMA_CALCULATIONS_H

#include "BaseFieldSource.hh"

//ORBIT bunch
#include "Bunch.hh"

namespace OrbitUtils{
	
	/** A class calculates the extrema and averages of the particles coordinates in the bunch.*/
	
	class BunchExtremaCalculator : public CppPyWrapper
	{
		public:
		
			/** Constructor. */
			BunchExtremaCalculator();
			
			/** Destructor */
			~BunchExtremaCalculator();
			
			/** The method calculates the extrema of the particles coordinates in the bunch. */
			void getExtremaXYZ(Bunch* bunch, 
				double& xMin, double& xMax, 
				double& yMin, double& yMax, 
				double& zMin, double& zMax);
			
			/** The method calculates the extrema of the particles coordinates xp, yp, dE in the bunch. */
			void getExtremaXpYpdE(Bunch* bunch, 
			double& xpMin, double& xpMax, 
			double& ypMin, double& ypMax, 
			double& dE_Min, double& dE_Max);
			
			/** The method calculates the z extrema of the particles coordinates in the bunch. */
			void getExtremaZ(Bunch* bunch,  
				double& zMin, double& zMax);
			
			/** The method calculates the sqrt(x**2+y**2) extrema of the particles coordinates in the bunch. */
			void getExtremaR(Bunch* bunch, double& rMax);		
			
	};
};
///////////////////////////////////////////////////////////////////////////
//
// END OF FILE
//
///////////////////////////////////////////////////////////////////////////

#endif
