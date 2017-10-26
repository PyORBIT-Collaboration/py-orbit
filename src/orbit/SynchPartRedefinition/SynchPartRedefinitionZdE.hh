#ifndef SYNCH_PARTICLE_REDEFINITION_H
#define SYNCH_PARTICLE_REDEFINITION_H

//pyORBIT utils
#include "CppPyWrapper.hh"

#include "Bunch.hh"

using namespace std;

/** 
  The SynchPartRedefinitionZdE class calculates the average of z and dE coordinates 
  redefines the syncronous particle's energy to move the synch. particle 
  to the center of the bunch's phase space. 
*/

class SynchPartRedefinitionZdE: public OrbitUtils::CppPyWrapper
{
	public:
		
		/** Constructor*/
		SynchPartRedefinitionZdE();
				
		/** Destructor */
		virtual ~SynchPartRedefinitionZdE();
		
		/** Performs the calculation of the z and dE averages of the bunch */
		void analyzeBunch(Bunch* bunch);
				
    /** Move the synch particle energy to the average energy */
    void centerE(Bunch* bunch);	

    /** Shift the synch particle energy */
    void shiftE(Bunch* bunch,double delta_dE);		
		
    /** Move the synch particle's z position to the center of the bunch */
    void centerZ(Bunch* bunch);  

    /** Shift the synch particle's z position */
    void shiftZ(Bunch* bunch,double delta_z);  
    
		/**Returns the average z postion */
		double getAvg_Z();
		
		/** Returns the average dE value */
		double getAvg_dE();
			
	private:
		
		/** array with average values for z and dE coordinates */
		double* z_dE_avg_arr;
		double* z_dE_avg_arr_MPI;
		
		
};

#endif
//endif for SYNCH_PARTICLE_REDEFINITION_H
