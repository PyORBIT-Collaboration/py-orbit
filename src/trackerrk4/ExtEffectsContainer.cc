//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    CppExternalEffects.cc
//
// CREATED
//    06/27/2008
//
// DESCRIPTION
//    The base class for C++ implementation of a external effects
//    during the transport of particles through the external field.
//    It should be sub-classed on Python level and implement
//    setupEffects(Bunch* bunch)
//    finalizeEffects(Bunch* bunch)
//    applyEffects(Bunch* bunch, int index,
//	                            double* y_in_vct, double* y_out_vct,
//														  double t, double t_step,
//														  OrbitUtils::BaseFieldSource* fieldSource)
//    methods.
//    The results of these methods will be available from the c++ level.
//    This is an example of embedding Python in C++ Orbit level.
//
///////////////////////////////////////////////////////////////////////////
#include "orbit_mpi.hh"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>

#include "ExtEffectsContainer.hh"

using namespace LaserStripping;
using namespace OrbitUtils;

ExtEffectsContainer::ExtEffectsContainer(){
}

ExtEffectsContainer::~ExtEffectsContainer(){
	for (int i=0;i<ref_eff.size();i++){
		if(ref_eff[i]->getPyWrapper() == NULL){
			delete ref_eff[i];
		} else {
			Py_XDECREF(ref_eff[i]->getPyWrapper());
		}
	}
}

void ExtEffectsContainer::AddEffect(ExternalEffects* eff)	{
	if(eff->getPyWrapper() != NULL){
		Py_INCREF(eff->getPyWrapper());
	}
	ref_eff.push_back(eff);
}

void ExtEffectsContainer::setupEffects(Bunch* bunch){
	for (int i=0;i<ref_eff.size();i++){
		ref_eff[i]->setupEffects(bunch);
	}
}

void ExtEffectsContainer::finalizeEffects(Bunch* bunch) {
	for (int i=0;i<ref_eff.size();i++){
		ref_eff[i]->finalizeEffects(bunch);
	}
}

void ExtEffectsContainer::applyEffects(Bunch* bunch, int index,
                                       double* y_in_vct, double* y_out_vct,
                                       double t, double t_step,
                                       BaseFieldSource* fieldSource,
                                       RungeKuttaTracker* tracker) {
	for (int i=0;i<ref_eff.size();i++){
		ref_eff[i]->applyEffects(bunch, index, y_in_vct, y_out_vct, t, t_step, fieldSource, tracker);
	}
}




























