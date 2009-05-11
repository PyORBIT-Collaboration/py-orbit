//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//   ExternalEffects.hh
//
// CREATED
//    04/18/2008
//
// DESCRIPTION
//    A base class for anything external that acting on the bunch except
//    slow changing magnetic and electric fields.
//
///////////////////////////////////////////////////////////////////////////
#ifndef EXTERNAL_EFFECTS_H
#define EXTERNAL_EFFECTS_H

#include <string>

#include "Bunch.hh"
#include "BaseFieldSource.hh"

#include "CppPyWrapper.hh"

namespace TrackerRK4{
	
	class RungeKuttaTracker;
	
	class ExternalEffects: public OrbitUtils::CppPyWrapper
	{
		//--------------------------------------------------
		// public methods of the ExternalEffects class
		//--------------------------------------------------
	public:
		
		ExternalEffects();
		virtual ~ExternalEffects();
		
		/** It initializes effects. */
		virtual void setupEffects(Bunch* bunch);
		
		/*it memorizes initial coordinates and impulses before rk step*/
		virtual void memorizeInitParams(Bunch* bunch);
		
		/** It finalizes effects. */
		virtual void finalizeEffects(Bunch* bunch);

		/** It applies the external effects to a particle with certain index. */
		virtual void applyEffects(Bunch* bunch, int index, 
	                            double* y_in_vct, double* y_out_vct, 
														  double t, double t_step, 
														  OrbitUtils::BaseFieldSource* fieldSource,
															RungeKuttaTracker* tracker);	
		
		/** It returns the name of the effect to distinguish them later. */
		std::string getName();
		
		/** It sets the name of the effect to distinguish them later. */
		void setName(std::string name);
		

		int getRankSetup();
		int getRankMemorize();
		int getRankApply();
		int getRankFinalize();
		

		void setRankSetup(int i);
		void setRankMemorize(int i);
		void setRankApply(int i);
		void setRankFinalize(int i);
		
	  private:
		   
			std::string name;

			
			int rank_setup;
			int rank_memorize;
			int rank_apply;
			int rank_finalize;
			
			
			
		
	};
	
}; // end of TrackerRK4 name-space


#endif
