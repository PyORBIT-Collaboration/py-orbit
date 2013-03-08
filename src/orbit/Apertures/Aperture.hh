//The base class for Apertures. It defines the interface for Aperture
#ifndef APERTURE_H
#define APERTURE_H

//pyORBIT utils
#include "CppPyWrapper.hh"
#include "Bunch.hh"

using namespace std;

/** 
  The aperture class is used to define how a bunch propogates through an aperture
*/
    
class Aperture: public OrbitUtils::CppPyWrapper
{
public:
	
	/** Aperture */
	Aperture(double shape, double a, double b, double c, double d, double angle); //(double shape, double size...)

	/** Routine for transfering particles through a aperture */
	void checkBunch(Bunch* bunch, Bunch* lostbunch);
	
	
private:

	/** check to see if the particle is inside the aperture */
	int checkCollFlag(double x, double y);
	
	/** delete the particle from the main bunch and add it to the lost particles bunch. */
	//void loseParticle(Bunch* bunch, Bunch* lostbunch, int ip, int& nLost, int& coll_flag, double& zrl);

protected:

	//Counters	
	int nLost;

	//Aperture parameters
	//double length_, density_fac_, a_, b_, c_, d_, angle_;
	//int ma_, shape_;

};

//end of APERTURE_H ifdef
#endif

