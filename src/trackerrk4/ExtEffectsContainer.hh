#ifndef EXTEFFECTSCONTAINER_HH_
#define EXTEFFECTSCONTAINER_HH_



#include "Python.h"

#include "ExternalEffects.hh"
#include <vector>



using namespace TrackerRK4;

namespace LaserStripping{
	
	class  ExtEffectsContainer: public ExternalEffects
	{
	public:
		
		/** Constructor. */
		ExtEffectsContainer();
		
		/** Destructor. */
		~ExtEffectsContainer();
		
		/** Adds the instance of the  ExternalEffects class to the container. */
		void AddEffect(ExternalEffects* eff);
		
		/** It initializes effects. */
		void setupEffects(Bunch* bunch);
		
		/** It finalizes effects. */
		void finalizeEffects(Bunch* bunch);
		
		/** It applies the external effects to a particle with certain index. */
		void applyEffects(Bunch* bunch, int index,
			double* y_in_vct, double* y_out_vct,
			double t, double t_step,
			OrbitUtils::BaseFieldSource* fieldSource,
			RungeKuttaTracker* tracker);
		
		private:
			
			vector<ExternalEffects*>	ref_eff;

	};
};



#endif /*EXTEFFECTSCONTAINER_HH_*/


