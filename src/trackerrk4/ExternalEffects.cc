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
#include "RungeKuttaTracker.hh"
#include <iostream>
#include <fstream>
#include <cstdlib>

#include "ExternalEffects.hh"


using namespace TrackerRK4;
using namespace OrbitUtils;
		
ExternalEffects::ExternalEffects(){
	name = "empty";

	rank_setup=0;
	rank_memorize=0;
	rank_apply=0;
	rank_finalize=0;
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

/** 
    It applies the external effects to a particle with certain index. 
    y_in_vct and y_out_vct are double[6] vectors with initial and final
                           		[r,p] coordinates for particular time step.
		t, t_step - initial moment and time step
		fieldSource - electric and magnetic field source 
		tracker - RungeKuttaTracker instance
*/
void ExternalEffects::applyEffectsForEach(
	                                Bunch* bunch, int index, 
	                                double* y_in_vct, double* y_out_vct, 
																	 double t, double t_step, 
																	 BaseFieldSource* fieldSource,
																	 RungeKuttaTracker* tracker){
}

/** 
    It applies the external effects to the bunch as a whole. 
		t, t_step - initial moment and time step
		fieldSource - electric and magnetic field source 
		tracker - RungeKuttaTracker instance
*/
void ExternalEffects::applyEffects(Bunch* bunch, 
																	  double t, double t_step, 
																	  BaseFieldSource* fieldSource,
																	  RungeKuttaTracker* tracker){
}


/** It returns the name of the effect to distinguish them later. */
std::string ExternalEffects::getName(){
	return name;
}

int ExternalEffects::getRankSetup()	{
	return rank_setup;
}

int ExternalEffects::getRankMemorize()	{
	return rank_memorize;
}

int ExternalEffects::getRankApply()	{
	return rank_apply;
}

int ExternalEffects::getRankFinalize()	{
	return rank_finalize;
}

void ExternalEffects::setRankSetup(int i)	{
	rank_setup = i;
}

void ExternalEffects::setRankMemorize(int i)	{
	rank_memorize = i;
}

void ExternalEffects::setRankApply(int i)	{
	rank_apply = i;
}

void ExternalEffects::setRankFinalize(int i)	{
	rank_finalize = i;
}

/** It sets the name of the effect to distinguish them later. */
void ExternalEffects::setName(std::string name){
	this->name = name;
}
