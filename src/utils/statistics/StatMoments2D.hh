#ifndef STAT_MOMENTS_2D_H
#define STAT_MOMENTS_2D_H

#include "orbit_mpi.hh"
#include "wrap_mpi_comm.hh"

//pyORBIT utils
#include "CppPyWrapper.hh"

using namespace std;

/** 
The Moments1D class calculates the arbitrary moments of the (u,up) distribution.
It is used by other classes to calculate Twiss paraemeters etc.
*/

namespace OrbitUtils{ 
	
	class StatMoments2D: public CppPyWrapper
	{
	public:
		
		/** Constructor with max order = 2 by default */
		StatMoments2D();
		
		/** Constructor with maximal order parameter */
		StatMoments2D(int maxOrder);
		
		/** Destructor */
		virtual ~StatMoments2D();
		
		/** Sets the maximal order of the moments */   
		void setMaxOrder(int maxOrder);
		
		/** Returns the maximal order of the moments */   
		int getMaxOrder();
		
		/** Initialize all internal arrays to get ready to gather statistical information. */   
		void clean();
		
		/** Takes into account the one point (u,up) */   	
		void account(double u, double up);
		
		/** Returns the statistical moment of the distribution with particular orders for u and up */ 	
		double getStatMoment(int order_u,int order_up);
		
		/** Returns the statistical moment of the distribution with particular order for u*/ 	
		double getStatMomentU(int order_u);
		
		/** Returns the statistical moment of the distribution with particular order for up*/ 	
		double getStatMomentUP(int order_up);
		
		/** Returns the minimal value of u */ 	
		double getMinU();
		
		/** Returns the maximal value of u */ 	
		double getMaxU();
		
		/** Returns the minimal value of up */ 	
		double getMinUP();
		
		/** Returns the maximal value of up */ 	
		double getMaxUP();
		
		/** Returns the number of points in the statistic */ 	
		int getCount();
		
		/** It will synchronize the moments through the MPI communicator */ 	
		int synchronizeMPI(pyORBIT_MPI_Comm* pyComm);
		
    /** Returns the emittance */
		double getEmittance();
		
    /** Returns Twiss alpha */
		double getAlpha();
		
    /** Returns Twiss beta */
		double getBeta();
		
    /** Returns Twiss gamma */
		double getGamma();		
		
		/** Returns the rms value of u */ 	
		double getRmsU();
		
		/** Returns the rms value of up */ 	
		double getRmsUP();				
		
	private:
		
		//make arrays
		void makeArrays();
		
	private:
		
		//max order
		int max_order;
		
		//number of points counted so far
		int count;
		
		//moments array
		double** stat_arr;
		
		double u_max;
		double u_min;
		double up_max;
		double up_min;
		
	};
	
} //end of OrbitUtils{} 

#endif
//endif for STAT_MOMENTS_2D_H
