//The base class for Collimators. It defines the interface for collimator
#ifndef COLLIMATOR_H
#define COLLIMATOR_H

//pyORBIT utils
#include "CppPyWrapper.hh"
#include "Bunch.hh"

using namespace std;

/** 
  The collimator class is used to define how a bunch propogates through a collimator
*/
    
class Collimator: public OrbitUtils::CppPyWrapper
{
public:
	
	/** Constructor */
	Collimator(double length, int ma, 
		              double density_fac, int shape, 
									double a, double b, double c, double d, double angle);

	/** Routine for transfering particles through a collimator */
	void collimateBunch(Bunch* bunch, Bunch* lostbunch);
	
	
private:

	/** check to see if the particle is inside the collimator */
	int checkCollFlag(double x, double y);
	
	/** drift a particle in small steps and check for collimator entry on every step. */
	int driftParticle(int coll_flag, double& zrl, double length, double* coords, SyncPart* syncpart);

	/** drift a particle in small steps and check for collimator entry on every step. */
	double getDirection(double* coords, SyncPart* syncpart);
						
	/** get the fraction offset of the particle momentum from the sync part momentum, 1 +dp/p */
	double getPFactor(double* coords, SyncPart* syncpart);

	/** get the particle beta */
	double getBeta(double* coords, SyncPart* syncpart);
	
	/** get the particle momentum */
	double getP(double* coords, SyncPart* syncpart);
	
	/** check the particle stepsize. Reset it if necessary  */
	void checkStep(double rl, double radlengthfac, double& stepsize, double* coords, SyncPart* syncpart);
	
	/** take a step inside the collimator with MCS and ionization energy loss */
	void takeStep(Bunch* bunch, Bunch* lostbunch, double* coords, SyncPart* syncpart, double z, double a, double density, long& idum, double stepsize, double& zrl, double& rl, int& coll_flag, int ip);
	
	/** delete the particle from the main bunch and add it to the lost particles bunch. */
	void loseParticle(Bunch* bunch, Bunch* lostbunch, int ip, int& nLost, int& coll_flag, double& zrl);

protected:

	//Counters	
	int nHits;
	int nLost;

	//Collimator parameters
	double length_, density_fac_, a_, b_, c_, d_, angle_;
	int ma_, shape_;

};

//end of COLLIMATOR_H ifdef
#endif

