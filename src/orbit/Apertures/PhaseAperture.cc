///////////////////////////////////////////////////////////////////////////
//
// NAME
//
//   PhaseAperture::PhaseAperture
//
// Author: A. Shishlo
// Date: 03/29/2017
//
// DESCRIPTION
//   Constructs an aperture for the longitudinal direction.
//
// PARAMETERS
//   The class will remove particles (and put them into lost_bunch) from 
//   the bunch whose longitudinal phase is outside the specified min and 
//   max phases. This class mostly is intended to be used in linacs. 
//   The user have to specify the RF frequency to translate the longitudinal 
//   position from meters to the phase in degrees.
//
///////////////////////////////////////////////////////////////////////////


#include "PhaseAperture.hh"
#include "SyncPart.hh"
#include "OrbitConst.hh"

#include <iostream>
#include <cmath>
#include <cfloat>
#include <cstdlib>

#include "ParticleAttributes.hh"

/**
 The PhaseAperture class constructor. It needs the RF frequency in Hz to translate
 from z coordinate in meters to phase in degrees.
 */
PhaseAperture::PhaseAperture(double frequency): CppPyWrapper(NULL)
{
	frequency_ = frequency;
	minPhase_ = -1.0e+36;
	maxPhase_ = +1.0e+36;
	pos_ = 0.;
}

/**
 The method sets the phase limits particles in the bunch. Phase limits are in degrees.
 */
void PhaseAperture::setPhaseLimits(double minPhase, double maxPhase){
	minPhase_ = minPhase;
	maxPhase_ = maxPhase;
}

/**
 Returns the min limit of the particles' phases. Phase is in degrees.
 */
double PhaseAperture::getMinPhase(){
	return minPhase_; 
}

/**
 Returns the max limit of the particles' phases. Phase is in degrees.
 */
double PhaseAperture::getMaxPhase(){
	return maxPhase_; 
}

/**
 Sets the position of the Phase Aperture node.
 */
void PhaseAperture::setPosition(double position){
	pos_ = position;
}

/**
 Returns the position of the Phase Aperture node.
 */
double PhaseAperture::getPosition(){
	return pos_;
}

/**
 Sets the RF frequency of the Phase Aperture node in Hz.
 */
void PhaseAperture::setRfFrequency(double frequency){
	frequency_ = frequency;
}

/**
 Returns the RF frequency of the Phase Aperture node in Hz.
 */
double PhaseAperture::getRfFrequency(){
	return frequency_;
}

/**
 Removes particles with phases outside the phase limits and puts them
 into the lost bunch is it exists.
 */
void PhaseAperture::checkBunch(Bunch* bunch, Bunch* lostbunch){

	bunch->compress();
	if(lostbunch != NULL) lostbunch->compress();
	double m_size = 0.;
	int nParts = bunch->getSize();
	double** coord = bunch->coordArr();
	
	ParticleAttributes* lostPartAttr = NULL;
	
	ParticleAttributes* partIdNumbAttr = NULL;
	ParticleAttributes* partIdNumbInitAttr = NULL;
	
	ParticleAttributes* partMacroAttr = NULL;
	ParticleAttributes* partMacroInitAttr = NULL;	
	
	ParticleAttributes* partInitCoordsAttr = NULL;
	ParticleAttributes* partInitCoordsInitAttr = NULL;
	
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
		
		if(bunch->hasParticleAttributes("macrosize") > 0){
			partMacroInitAttr = bunch->getParticleAttributes("macrosize");
			if(lostbunch->hasParticleAttributes("macrosize") <= 0){
				std::map<std::string,double> params_dict;
				lostbunch->addParticleAttributes("macrosize",params_dict);
			}
			partMacroAttr = lostbunch->getParticleAttributes("macrosize");
		}	
		
		if(bunch->hasParticleAttributes("ParticleInitialCoordinates") > 0){
			partInitCoordsInitAttr = bunch->getParticleAttributes("ParticleInitialCoordinates");
			if(lostbunch->hasParticleAttributes("ParticleInitialCoordinates") <= 0){
				std::map<std::string,double> params_dict;
				lostbunch->addParticleAttributes("ParticleInitialCoordinates",params_dict);
			}	
			partInitCoordsAttr = lostbunch->getParticleAttributes("ParticleInitialCoordinates");
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

	double beta = bunch->getSyncPart()->getBeta();
	double lambda = OrbitConst::c*beta/frequency_;
	double z_to_phase_coeff = - 360./lambda;
	double z = 0.;
	double phase = 0.;
	
	for (int count = 0; count < nParts; count++){
		z = coord[count][4];
		phase = z*z_to_phase_coeff;
		if(phase < minPhase_ || phase > maxPhase_){
			if(lostbunch != NULL) {
				lostbunch->addParticle(coord[count][0], coord[count][1], coord[count][2], coord[count][3], coord[count][4], coord[count][5]);
				//pos_ is a position in lattice where particle is lost
				lostPartAttr->attValue(lostbunch->getSize() - 1, 0) = pos_;
				if(partIdNumbAttr != NULL){
					partIdNumbAttr->attValue(lostbunch->getSize() - 1, 0) = partIdNumbInitAttr->attValue(count,0);
				}
				if(partMacroAttr != NULL){
					partMacroAttr->attValue(lostbunch->getSize() - 1, 0) = partMacroInitAttr->attValue(count,0);
				}
				if(partInitCoordsInitAttr != NULL){
					for(int init_ind = 0; init_ind < 6; init_ind++){
					  partInitCoordsAttr->attValue(lostbunch->getSize() - 1,init_ind) = partInitCoordsInitAttr->attValue(count,init_ind);
					}
				}
				if(partTurnNumberAttr != NULL){
					partTurnNumberAttr->attValue(lostbunch->getSize() - 1, 0) = turn;
				}				
			}
			bunch->deleteParticleFast(count);
		}
	}
	
	//Update synchronous particle, compress bunch
	bunch->compress();
}

	

