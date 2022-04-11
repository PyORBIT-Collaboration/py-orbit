#include "Aperture.hh"
#include "SyncPart.hh"
#include "OrbitConst.hh"
#include "Random.hh"

#include <iostream>
#include <cmath>
#include <cfloat>
#include <cstdlib>

#include "ParticleAttributes.hh"
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
	if(lostbunch != NULL) lostbunch->compress();
	double m_size = 0.;
	int nParts = bunch->getSize();
	double** coord = bunch->coordArr();
	
	
	ParticleAttributes* lostPartAttr = NULL;
	
	ParticleAttributes* partIdNumbAttr = NULL;
	ParticleAttributes* partIdNumbInitAttr = NULL;
	
	ParticleAttributes* partInitCoordsAttr = NULL;
	ParticleAttributes* partInitCoordsInitAttr = NULL;
	
	ParticleAttributes* partMacroAttr = NULL;
	ParticleAttributes* partMacroInitAttr = NULL;	
	
	ParticleAttributes* partTurnNumberAttr = NULL;
	double turn = 0.;
	
	if(lostbunch != NULL) {
		if(lostbunch->hasParticleAttributes("LostParticleAttributes") <= 0){
			std::map<std::string,double> params_dict;
			lostbunch->addParticleAttributes("LostParticleAttributes",params_dict);
		}
		lostPartAttr = lostbunch->getParticleAttributes("LostParticleAttributes");
		
		if(bunch->hasParticleAttributes("ParticleIdNumber") > 0){
			partIdNumbInitAttr = bunch->getParticleAttributes("ParticleIdNumber");
			if(lostbunch->hasParticleAttributes("ParticleIdNumber") <= 0){
				std::map<std::string,double> params_dict;
				lostbunch->addParticleAttributes("ParticleIdNumber",params_dict);
			}	
			partIdNumbAttr = lostbunch->getParticleAttributes("ParticleIdNumber");
		}
		
		if(bunch->hasParticleAttributes("ParticleInitialCoordinates") > 0){
			partInitCoordsInitAttr = bunch->getParticleAttributes("ParticleInitialCoordinates");
			if(lostbunch->hasParticleAttributes("ParticleInitialCoordinates") <= 0){
				std::map<std::string,double> params_dict;
				lostbunch->addParticleAttributes("ParticleInitialCoordinates",params_dict);
			}
			partInitCoordsAttr = lostbunch->getParticleAttributes("ParticleInitialCoordinates");			
		}		
		
		if(bunch->hasParticleAttributes("macrosize") > 0){
			partMacroInitAttr = bunch->getParticleAttributes("macrosize");
			if(lostbunch->hasParticleAttributes("macrosize") <= 0){
				std::map<std::string,double> params_dict;
				lostbunch->addParticleAttributes("macrosize",params_dict);
			}
			partMacroAttr = lostbunch->getParticleAttributes("macrosize");
		}
		
		if (bunch->hasParticleAttributes("TurnNumber") > 0) {
			if (lostbunch->hasParticleAttributes("TurnNumber") <= 0) {
				std::map<std::string,double> part_attr_dict;
				lostbunch->addParticleAttributes("TurnNumber",part_attr_dict);
			}
			std::string attr_name_str("TurnNumber");
			turn = 1.0*bunch->getBunchAttributeInt(attr_name_str);
			partTurnNumberAttr = lostbunch->getParticleAttributes("TurnNumber");
		}
		
		lostbunch->setMacroSize(bunch->getMacroSize());
	}
	
	// shape = 1              ===circular aperture===
	if(shape == 1){
	  for (int count = 0; count < nParts; count++){
	  	if((pow((coord[count][0]-c), 2) + pow((coord[count][2]-d), 2)) >= pow(a, 2)){
	  		if(lostbunch != NULL) {
	  			lostbunch->addParticle(coord[count][0], coord[count][1], coord[count][2], coord[count][3], coord[count][4], coord[count][5]);
	  			//pos_ is a position in lattice where particle is lost
	  			lostPartAttr->attValue(lostbunch->getSize() - 1, 0) = pos_;
	  			if(partIdNumbAttr != NULL){
	  				partIdNumbAttr->attValue(lostbunch->getSize() - 1, 0) = partIdNumbInitAttr->attValue(count,0);
	  			}
	  			if(partInitCoordsAttr != NULL){
	  				for(int j=0; j < 6; ++j){
	  					partInitCoordsAttr->attValue(lostbunch->getSize() - 1, j) = partInitCoordsInitAttr->attValue(count,j);
	  				}
	  			}
	  			if(partMacroAttr != NULL){
	  				partMacroAttr->attValue(lostbunch->getSize() - 1, 0) = partMacroInitAttr->attValue(count,0);
	  			}
					if(partTurnNumberAttr != NULL){
						partTurnNumberAttr->attValue(lostbunch->getSize() - 1, 0) = turn;
					}
	  		}
	  		bunch->deleteParticleFast(count);
	  	}
	  }
	}
	
	// shape = 2              ===elipital aperture===
	if(shape == 2){
	  for (int count = 0; count < nParts; count++){
	  	if((pow((coord[count][0]-c), 2)/pow(a,2) + pow((coord[count][2]-d), 2)/pow(b,2)) >= 1){
	  		if(lostbunch != NULL) {
	  			lostbunch->addParticle(coord[count][0], coord[count][1], coord[count][2], coord[count][3], coord[count][4], coord[count][5]);
	  			//pos_ is a position in lattice where particle is lost
	  			lostPartAttr->attValue(lostbunch->getSize() - 1, 0) = pos_;
	  			if(partIdNumbAttr != NULL){
	  				partIdNumbAttr->attValue(lostbunch->getSize() - 1, 0) = partIdNumbInitAttr->attValue(count,0);
	  			}
	  			if(partInitCoordsAttr != NULL){
	  				for(int j=0; j < 6; ++j){
	  					partInitCoordsAttr->attValue(lostbunch->getSize() - 1, j) = partInitCoordsInitAttr->attValue(count,j);
	  				}
	  			}
	  			if(partMacroAttr != NULL){
	  				partMacroAttr->attValue(lostbunch->getSize() - 1, 0) = partMacroInitAttr->attValue(count,0);
	  			}
					if(partTurnNumberAttr != NULL){
						partTurnNumberAttr->attValue(lostbunch->getSize() - 1, 0) = turn;
					}
	  		}
	  		bunch->deleteParticleFast(count);
	  	}
	  }
	}
	
	// shape = 3              ===rectangular aperture===
	if(shape == 3){
	  for (int count = 0; count < nParts; count++){
	  	if((abs((coord[count][0])-c)>=a)||(abs((coord[count][2]-d))>=b)){
	  		if(lostbunch != NULL) {
	  			lostbunch->addParticle(coord[count][0], coord[count][1], coord[count][2], coord[count][3], coord[count][4], coord[count][5]);
	  			//pos_ is a position in lattice where particle is lost
	  			lostPartAttr->attValue(lostbunch->getSize() - 1, 0) = pos_;
	  			if(partIdNumbAttr != NULL){
	  				partIdNumbAttr->attValue(lostbunch->getSize() - 1, 0) = partIdNumbInitAttr->attValue(count,0);
	  			}
	  			if(partInitCoordsAttr != NULL){
	  				for(int j=0; j < 6; ++j){
	  					partInitCoordsAttr->attValue(lostbunch->getSize() - 1, j) = partInitCoordsInitAttr->attValue(count,j);
	  				}
	  			}
	  			if(partMacroAttr != NULL){
	  				partMacroAttr->attValue(lostbunch->getSize() - 1, 0) = partMacroInitAttr->attValue(count,0);
	  			}
					if(partTurnNumberAttr != NULL){
						partTurnNumberAttr->attValue(lostbunch->getSize() - 1, 0) = turn;
					}
	  		}
	  		bunch->deleteParticleFast(count);
	  	}
	  }
	}
				
	//Update synchronous particle, compress bunch
	bunch->compress();
}

void Aperture::setPosition(double position){
	pos_ = position;
}

	

