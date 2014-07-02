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
//   Constructs an aperture
//
// PARAMETERS
//   Shape can be the intergers 1, 2, or 3. 1 is a circular aperture, 2 is
//   an elipital aperture, and 3 is a rectangular aperture. a is either the
//   radius of the circle or half the length of the aperture in the x 
//   diminsion. b is half the length of the aperture in the y diminsion. 
//   c is the x offset and d is the y offset of the aperture.
//
// RETURNS
//   int.
//
///////////////////////////////////////////////////////////////////////////

Aperture::Aperture(int shape, double a, double b, double c, double d, double pos): CppPyWrapper(NULL)
{
	shape_ = shape;
	a_ = a;
	b_ = b;
	c_ = c;
	d_ = d;
	pos_ = pos;
}

void Aperture::checkBunch(Bunch* bunch, Bunch* lostbunch){
         
  //Removes particle if located outside aperture from main bunch and adds it to the lost bunch.

	int j = 1, coll_flag = 0, lastArg, trackit;
	double a = a_;
	double b = b_;
	double c = c_;
	double d = d_;
	int shape = shape_;
       
	bunch->compress();
	double m_size = 0.;
	int nParts = bunch->getSize();
	double** coord = bunch->coordArr();

	for (int count = 0; count < nParts; count++){
	
	    if(shape == 1)
	      if((pow((coord[count][0]-c), 2) + pow((coord[count][2]-d), 2)) >= pow(a, 2)){
                    lostbunch->addParticle(coord[count][0], coord[count][1], coord[count][2], coord[count][3], coord[count][4], coord[count][5]);
				if (lostbunch->hasParticleAttributes("LostParticleAttributes") > 0) {
					lostbunch->getParticleAttributes("LostParticleAttributes")->attValue(lostbunch->getSize() - 1, 0) = pos_; //position in lattice where particle is lost
					}
                    bunch->deleteParticleFast(count);
	        }

		if(shape == 2)
	      if((pow((coord[count][0]-c), 2)/pow(a,2) + pow((coord[count][2]-d), 2)/pow(b,2)) >= 1){
                    lostbunch->addParticle(coord[count][0], coord[count][1], coord[count][2], coord[count][3], coord[count][4], coord[count][5]);
			  if (lostbunch->hasParticleAttributes("LostParticleAttributes") > 0) {
				  lostbunch->getParticleAttributes("LostParticleAttributes")->attValue(lostbunch->getSize() - 1, 0) = pos_; //position in lattice where particle is lost
			  }
                    bunch->deleteParticleFast(count);
	        }
	
	    if(shape == 3)
	      if((abs((coord[count][0])-c)>=a)||(abs((coord[count][2]-d))>=b)){
	          lostbunch->addParticle(coord[count][0], coord[count][1], coord[count][2], coord[count][3], coord[count][4], coord[count][5]);
			  if (lostbunch->hasParticleAttributes("LostParticleAttributes") > 0) {
				  lostbunch->getParticleAttributes("LostParticleAttributes")->attValue(lostbunch->getSize() - 1, 0) = pos_; //position in lattice where particle is lost
			  }
        	    bunch->deleteParticleFast(count);
	      }
	}
		
	//Update synchronous particle, compress bunch
	bunch->compress();
}

void Aperture::setPosition(double position){
	pos_ = position;
}

	

