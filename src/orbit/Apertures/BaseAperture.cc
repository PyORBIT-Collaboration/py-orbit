#include "BaseAperture.hh"
#include "ParticleAttributes.hh"

#include <iostream>
#include <cmath>
#include <cfloat>
#include <cstdlib>

///////////////////////////////////////////////////////////////////////////
//
// NAME
//
//   BaseAperture
//
// AUTHOR: 
//   Andrei Shishlo October 2022
//
//   BaseAperture class defines actions with macro particles in a main bunch 
//   and a lost bunch instances with respect of transverse coordinates of 
//   macro particles. It can keep macro particle in the main bunch or remove it
//   after asking the BaseApertureShape class instance about the suitability 
//   of particle's coordinates.
//
///////////////////////////////////////////////////////////////////////////

/** BaseAperture constructor */
BaseAperture::BaseAperture(): CppPyWrapper(NULL)
{
	apertureName = "no_name";
	isActive = 1;
	nLost_ = 0;
	pos_ = 0.;
	apertureShape = NULL;
}

/** BaseAperture decstructor */
BaseAperture::~BaseAperture()
{
	if(apertureShape != NULL){
		Py_XDECREF(apertureShape->getPyWrapper());
	}
}

/** Returns aperture shape */
BaseApertureShape* BaseAperture::getApertureShape(){
	return apertureShape;
}

/** Sets aperture shape */
void BaseAperture::setApertureShape(BaseApertureShape* apertureShapeIn){

	nLost_ = 0;
	
	if(apertureShapeIn == NULL){
		return;
	}
	
	if( ((PyObject*) apertureShapeIn->getPyWrapper()) == NULL){
		ORBIT_MPI_Finalize("BaseAperture class setApertureShape(...): BaseApertureShape Python class needed! Stop.");
	}	
	
	if(apertureShape != NULL){
		if( ((PyObject*) apertureShape->getPyWrapper()) != NULL){
			Py_XDECREF( (PyObject*) apertureShape->getPyWrapper());
		}
	}

	apertureShape = apertureShapeIn;
	Py_INCREF((PyObject*) apertureShape->getPyWrapper());
}

/** 
	Routine for transfering particles through a aperture 
*/
void BaseAperture::checkBunch(Bunch* bunch, Bunch* lostbunch){
	
	nLost_ = 0;
	
	if(isActive != 1){
		return;
	}
	
	if(apertureShape == NULL){
		return;
	}	

	bunch->compress();
	if(lostbunch != NULL) lostbunch->compress();
	double m_size = 0.;
	int nParts = bunch->getSize();
	double** coord = bunch->coordArr();
	
	
	int nPartsGlobal = bunch->getSizeGlobal();
	
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

	//Loop over all particles in the bunch
	for (int count = 0; count < nParts; count++){
		//if particle is not inside the shape we remove it from bunch
		if(apertureShape->inside(bunch,count) != 1){
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
				
	//compress bunch
	bunch->compress();
	
	//Total particle loss across all CPUs in the bunch communicator
	nLost_ = nPartsGlobal - bunch->getSizeGlobal();
}

/** 
	Returns total particle loss across all CPUs in the bunch communicator 
*/
int BaseAperture::getNumberOfLost(){
	return nLost_;
}

/** Returns the aperture name */
string BaseAperture::getName(){
	return apertureName;
}
	
/** Sets the aperture name */
void BaseAperture::setName(string apertureNameIn){
	apertureName = apertureNameIn;
}

/**
	Returns the position of the node in the lattice.
*/
double BaseAperture::getPosition(){
	return pos_;
}

/**
	Sets the position of the node in the lattice.
*/
void BaseAperture::setPosition(double position){
	pos_ = position;
}
	
/**
	Sets the aperture in an active ( 1 )/ not active ( 0 ) state
*/
void BaseAperture::setOnOff(int isActive){
	this->isActive = isActive;
}

/**
	Returns the aperture in an active ( 1 )/ not active ( 0 ) state
*/
int BaseAperture::getOnOff(){
	return isActive;
}