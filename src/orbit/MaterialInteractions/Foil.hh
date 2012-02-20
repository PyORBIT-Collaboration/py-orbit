//The base class for Foils. It defines the interface for collimator
#ifndef FOIL_H
#define FOIL_H

//pyORBIT utils
#include "CppPyWrapper.hh"
#include "Bunch.hh"

using namespace std;

/** 
  The foil class is used to define how a bunch propogates through a foil
*/
    
class Foil: public OrbitUtils::CppPyWrapper
{
public:
	
	/** Constructor */
	Foil(double xmin, double xmax, double ymin, double ymax, double thick);

	/** Routine for transfering particles through a foil with full scattering. Can have particle loss. */
	void traverseFoilFullScatter(Bunch* bunch, Bunch* lostbunch);
	
	/** Routine for transfering particles through a foil with simplified scattering. No particle loss. */
	void traverseFoilSimpleScatter(Bunch* bunch);

	
private:

	/** check to see if the particle is inside the foil */
	int checkFoilFlag(double x, double y);
	
	/** get the direction of the partcile. */
	double getDirection(double* coords, SyncPart* syncpart);
						
	/** get the fraction offset of the particle momentum from the sync part momentum, 1 +dp/p */
	double getPFactor(double* coords, SyncPart* syncpart);

	/** get the particle beta */
	double getBeta(double* coords, SyncPart* syncpart);
	
	/** get the particle momentum */
	double getP(double* coords, SyncPart* syncpart);
	
	/** take a step inside the foil with MCS and ionization energy loss */
	void takeStep(Bunch* bunch, Bunch* lostbunch, double* coords, SyncPart* syncpart, double z, double a, double density, long& idum, double stepsize, double& zrl, double& rl, int& coll_flag, int ip);
	
	/** delete the particle from the main bunch and add it to the lost particles bunch. */
	void loseParticle(Bunch* bunch, Bunch* lostbunch, int ip, int& nLost, int& coll_flag, double& zrl);

protected:

	//Counters	
	int nHits;
	int nLost;

	//Foil parameters
	double xmin_, xmax_, ymin_, ymax_, thick_;
	int ma_;
	double length_;

};

//end of FOIL_H ifdef
#endif

