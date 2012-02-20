#ifndef MATERIAL_INTERACTIONS_H
#define MATERIAL_INTERACTIONS_H

//pyORBIT utils
#include "CppPyWrapper.hh"


using namespace std;

/** 
  The MaterialInteractions class contains a set of routines for calculating common interactions of a particle with a
  material, for instance scattering and energy loss in the material. Many of the routines only apply to ~GeV protons. 
  There are several material choices.
*/

class MaterialInteractions 
{
	public:
		
		/** Constructor */
		MaterialInteractions();
				
		/** Destructor */
		virtual ~MaterialInteractions();
		
		/** Routine to generate and apply random, uniformly distributed 2D momentum kicks */
		static void momentumKick(double t, double p, double& dpx, double& dpy);

		/** Routine to apply multiple coulomb scattering kicks following JD Jackson, Chapter 13 */
		static void mcsJackson(double stepsize, double z, double a, double rho, long& idum, double beta, double pfac, double& x, double& y, double& px, double& py);
		
		/** Routine to apply Rutherford scattering following JD Jackson, Chapter 13 */
		static double ruthScattJackson(double stepsize, double z, double a, double rho, long& idum, double beta, int trackit, double pfac, double& thetax, double& thetay);

		/** Routine to generate a random momentum transfer for low energy elastic scattering (<= 0.4 GeV) */
		static double elastic_t(double p, double a, long& idum);

		/** Routine to calculate ionization energy loss. */
		static double ionEnergyLoss(double beta, double z, double a);

	private:

		
};

#endif

