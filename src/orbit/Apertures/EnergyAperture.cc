///////////////////////////////////////////////////////////////////////////
//
// NAME
//
//   EnergyAperture::EnergyAperture
//
// Author: A. Shishlo
// Date: 03/29/2017
//
// DESCRIPTION
//   Constructs an aperture for the energy (relative to synch. part.) variable.
//
// PARAMETERS
//   The class will remove particles (and put them into lost_bunch) from 
//   the bunch whose energy (relative to synch. part.) is outside the specified min and 
//   max values. This class mostly is intended to be used in linacs. 
//   The energy is in GeV.
//
///////////////////////////////////////////////////////////////////////////


#include "EnergyAperture.hh"
#include "SyncPart.hh"
#include "OrbitConst.hh"

#include <iostream>
#include <cmath>
#include <cfloat>
#include <cstdlib>

#include "ParticleAttributes.hh"

/**
 The EnergyAperture class constructor.
 */
EnergyAperture::EnergyAperture(): CppPyWrapper(NULL)
{
	minEnergy_ = -1.0e+36;
	maxEnergy_ = +1.0e+36;
	pos_ = 0.;
}

/**
 The method sets the energy limits in GeV for particles in the bunch.
 */
void EnergyAperture::setEnergyLimits(double minEnergy, double maxEnergy){
	minEnergy_ = minEnergy;
	maxEnergy_ = maxEnergy;
}

/**
 Returns the min limit of the particles' energy.
 */
double EnergyAperture::getMinEnergy(){
	return minEnergy_; 
}

/**
 Returns the max limit of the particles' energy.
 */
double EnergyAperture::getMaxEnergy(){
	return maxEnergy_; 
}

/**
 Sets the position of the Energy Aperture node.
 */
void EnergyAperture::setPosition(double position){
	pos_ = position;
}

/**
 Returns the position of the Energy Aperture node.
 */
double EnergyAperture::getPosition(){
	return pos_;
}

/**
 Removes particles with phases outside the energy limits and puts them
 into the lost bunch is it exists.
 */
void EnergyAperture::checkBunch(Bunch* bunch, Bunch* lostbunch){

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
		
		lostbunch->setMacroSize(bunch->getMacroSize());
	}
	
	double dE = 0.;
	for (int count = 0; count < nParts; count++){
		dE = coord[count][5];
		if(dE < minEnergy_ || dE > maxEnergy_){
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
			}
			bunch->deleteParticleFast(count);
		}
	}
	
	//Update synchronous particle, compress bunch
	bunch->compress();
}

	

