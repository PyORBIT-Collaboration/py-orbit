#include "Aperture.hh"
#include "SyncPart.hh"
#include "OrbitConst.hh"
#include "Random.hh"

#include <iostream>
#include <cmath>
#include <cfloat>
#include <cstdlib>


// Constructor
///////////////////////////////////////////////////////////////////////////
//
// NAME
//
//   Aperture::Aperture
//
// DESCRIPTION
//   Constructs a aperture
//
// PARAMETERS
// 
//
// RETURNS
//   int.
//
///////////////////////////////////////////////////////////////////////////

Aperture::Aperture(double shape, double a, double b, double c, double d, double angle): CppPyWrapper(NULL)
{

}

void Aperture::checkBunch(Bunch* bunch, Bunch* lostbunch){

	int j = 1, coll_flag = 0, lastArg, trackit;
	double nAvogadro = 6.022045e23;
	
	bunch->compress();
	double m_size = 0.;
	int nParts = bunch->getSize();
	double** part_coord_arr = bunch->coordArr();
	
	
	//Update synchronous particle, compress bunch
	
	bunch->compress();
}

///////////////////////////////////////////////////////////////////////////
//
// NAME
//
//   Aperture::checkCollFlag
//
// DESCRIPTION
//   Checks to see if a particle is located inside a aperture.  Returns
//   1 if the particle is in the aperture, 0 if it isn't.
//
// PARAMETERS
//   x:      x coordinate of particle.
//   y:      y coordinate of particle.
//
// RETURNS
//   int.
//
///////////////////////////////////////////////////////////////////////////

int Aperture::checkCollFlag(double x, double y){

	double xtemp = x, ytemp = y;
	double PI = OrbitConst::PI;
	
	
/*
	if((pow(x, 2) + pow(y, 2)) >= pow(a, 2)){
		return 1;
	
	}
*/	
    return 0; 
}

///////////////////////////////////////////////////////////////////////////
//
// NAME
//
//   Aperture::loseParticle
//
// DESCRIPTION
//   Lose the particle from the alive bunch and add it to the lost bunch
//
// PARAMETERS
//	 coords: particle coordinates
//   bunch: the primary bunch
//	 lostbunch: the lost bunch
//	 nLost:		the number of lost particles in this aperture
//	 coll_flag:	flag for if particle is still in aperture (and alive)
//	 zrt:		remaining length in z the particle has in aperture
//
// RETURNS
//   nothing.
//
///////////////////////////////////////////////////////////////////////////

/*
void Aperture::loseParticle(Bunch* bunch, Bunch* lostbunch, int ip, int& nLost, int& coll_flag, double& zrl){
	double** coords = bunch->coordArr();
	lostbunch->addParticle(coords[ip][0], coords[ip][1], coords[ip][2], coords[ip][3], coords[ip][4], coords[ip][5]);
	
	bunch->deleteParticleFast(ip);
	
	nLost++;
	
}
*/	
	

