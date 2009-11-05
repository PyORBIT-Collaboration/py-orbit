#ifndef FIELDSOURCECONTAINER_HH_
#define FIELDSOURCECONTAINER_HH_




#include "Python.h"

#include "BaseFieldSource.hh"
#include <vector>

namespace OrbitUtils{
	
	class  FieldSourceContainer: public BaseFieldSource
	{
	public:
		
		/** Constructor. */
		FieldSourceContainer();
		
		/** Destructor. */
		~FieldSourceContainer();
		
		/** Adds the instance of the  ExternalEffects class to the container. */
		void AddFieldSource(BaseFieldSource* fs);
		
		/** Adds the instance of the  ExternalEffects class to the container. */
		void getElectricMagneticField(double x, double y, double z, double t, 
						double& E_x, double& E_y, double& E_z,
						double& H_x, double& H_y, double& H_z);
		
		private:

			std::vector<BaseFieldSource*>	ref;

	};
};
	

	










#endif /*FIELDSOURCECONTAINER_HH_*/
