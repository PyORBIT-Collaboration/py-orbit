//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//   ExternalEffects.cc
//
// CREATED
//    04/18/2008
//
// DESCRIPTION
//    A base class for anything external that acting on the bunch except
//    slow changing magnetic and electric fields.
//
///////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <fstream>
#include <cstdlib>

#include "ExternalEffects.hh"
#include "RungeKuttaTracker.hh"

using namespace TrackerRK4;
using namespace OrbitUtils;
		
ExternalEffects::ExternalEffects(){
	name = "empty";
}

ExternalEffects::~ExternalEffects(){
}

/** It initializes effects. */
void ExternalEffects::setupEffects(Bunch* bunch){
}

/*it memorizes initial coordinates and impulses before rk step*/
void ExternalEffects::memorizeInitParams(Bunch* bunch){
}

/** It finalizes effects. */
void ExternalEffects::finalizeEffects(Bunch* bunch){
}

/** It applies the external effects to a particle with certain index. 
    y_in_vct and y_out_vct are double[6] vectors with initial and final
                           		[r,p] coordinates for particular time step.
		t, t_step - initial moment and time step
		fieldSource - electric and magnetic field source 
		tracker - RungeKuttaTracker instance
*/
void ExternalEffects::applyEffects(Bunch* bunch, int index, 
	                                 double* y_in_vct, double* y_out_vct, 
																	 double t, double t_step, 
																	 BaseFieldSource* fieldSource,
																	 RungeKuttaTracker* tracker){
}

/** It returns the name of the effect to distinguish them later. */
std::string ExternalEffects::getName(){
	return name;
}

/** It sets the name of the effect to distinguish them later. */
void ExternalEffects::setName(std::string name){
	this->name = name;
}
