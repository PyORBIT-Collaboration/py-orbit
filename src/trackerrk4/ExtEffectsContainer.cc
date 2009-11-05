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




ExtEffectsContainer::ExtEffectsContainer(){
}

ExtEffectsContainer::~ExtEffectsContainer(){
	
	for (int i=0;i<ref.size();i++){ 
		if(ref[i]->getPyWrapper() == NULL){
			delete ref[i];
		} else {
			Py_XDECREF(ref[i]->getPyWrapper());
		}
	}

}

void ExtEffectsContainer::AddEffect(ExternalEffects* eff)	{
	if(eff->getPyWrapper() != NULL){
		Py_INCREF(eff->getPyWrapper());
	}

	ref.push_back(eff);
	
	ref_setup.push_back(eff);
	ref_memorize.push_back(eff);
	ref_apply.push_back(eff);
	ref_finalize.push_back(eff);
}

void ExtEffectsContainer::setupEffects(Bunch* bunch){
	
	ExternalEffects* exchange;
	int min;
	

	for (int i=0; i<ref_setup.size(); i++)	{
		min=ref_setup[i]->getRankSetup();
	for (int j=i+1; j<ref_setup.size(); j++)	{
		
		if(ref_setup[j]->getRankSetup()<min)	{
			min=ref_setup[j]->getRankSetup();		
			exchange = ref_setup[i];
			ref_setup[i] = ref_setup[j];
			ref_setup[j] = exchange;

		}
		
	}
	}
	
	for (int i=0; i<ref_memorize.size(); i++)	{
		min=ref_memorize[i]->getRankMemorize();
	for (int j=i+1; j<ref_memorize.size(); j++)	{
		
		if(ref_memorize[j]->getRankMemorize()<min)	{
			min=ref_memorize[j]->getRankMemorize();		
			exchange = ref_memorize[i];
			ref_memorize[i] = ref_memorize[j];
			ref_memorize[j] = exchange;

		}
		
	}
	}
	
	for (int i=0; i<ref_apply.size(); i++)	{
		min=ref_apply[i]->getRankApply();
	for (int j=i+1; j<ref_apply.size(); j++)	{
		
		if(ref_apply[j]->getRankApply()<min)	{
			min=ref_apply[j]->getRankApply();		
			exchange = ref_apply[i];
			ref_apply[i] = ref_apply[j];
			ref_apply[j] = exchange;

		}
		
	}
	}
	
	for (int i=0; i<ref_finalize.size(); i++)	{
		min=ref_finalize[i]->getRankFinalize();
	for (int j=i+1; j<ref_finalize.size(); j++)	{
		
		if(ref_finalize[j]->getRankFinalize()<min)	{
			min=ref_finalize[j]->getRankFinalize();		
			exchange = ref_finalize[i];
			ref_finalize[i] = ref_finalize[j];
			ref_finalize[j] = exchange;

		}
		
	}
	}
		
	
	
	
	for (int i=0;i<ref_setup.size();i++){
		ref_setup[i]->setupEffects(bunch);
	}



}

void ExtEffectsContainer::memorizeInitParams(Bunch* bunch){	

	for (int i=0;i<ref_memorize.size();i++){
		ref_memorize[i]->memorizeInitParams(bunch);
	}
}

void ExtEffectsContainer::finalizeEffects(Bunch* bunch) {
	for (int i=0;i<ref_finalize.size();i++){
		ref_finalize[i]->finalizeEffects(bunch);
	}
}

void ExtEffectsContainer::applyEffects(Bunch* bunch, int index,
                                       double* y_in_vct, double* y_out_vct,
                                       double t, double t_step,
                                       BaseFieldSource* fieldSource,
                                       RungeKuttaTracker* tracker) {
	for (int i=0;i<ref_apply.size();i++){
		ref_apply[i]->applyEffects(bunch, index, y_in_vct, y_out_vct, t, t_step, fieldSource, tracker);
	}
}




























