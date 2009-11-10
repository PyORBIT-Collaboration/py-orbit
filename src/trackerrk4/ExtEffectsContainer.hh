#ifndef EXTEFFECTSCONTAINER_HH_
#define EXTEFFECTSCONTAINER_HH_



#include "Python.h"

#include "ExternalEffects.hh"
#include <vector>



using namespace TrackerRK4;
using namespace OrbitUtils;


	
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
		
		/*it memorizes initial coordinates and impulses before rk step*/
		void memorizeInitParams(Bunch* bunch);
		
		/** It finalizes effects. */
		void finalizeEffects(Bunch* bunch);
		
		/** It applies the external effects to the bunch as a whole. */
		void applyEffects(Bunch* bunch,
			double t, double t_step,
			BaseFieldSource* fieldSource,
			RungeKuttaTracker* tracker);
		
		/** It applies the external effects to a particle with certain index. */
		void applyEffectsForEach(Bunch* bunch, int index,
			double* y_in_vct, double* y_out_vct,
			double t, double t_step,
			BaseFieldSource* fieldSource,
			RungeKuttaTracker* tracker);
		
		private:
			
			vector<ExternalEffects*>	ref;			
			vector<ExternalEffects*>	ref_setup;
			vector<ExternalEffects*>	ref_memorize;
			vector<ExternalEffects*>	ref_apply;
			vector<ExternalEffects*>	ref_finalize;

	};
	





#endif /*EXTEFFECTSCONTAINER_HH_*/


