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
#include "ExternalEffects.hh"

#include <iostream>
#include <fstream>
#include <cstdlib>

using namespace Tracker3DField;

		
ExternalEffects::ExternalEffects(){
	name = "empty";
}

ExternalEffects::~ExternalEffects(){
}

/** It initializes effects. */
void ExternalEffects::setupEffects(Bunch bunch){
}

/** It finalizes effects. */
void ExternalEffects::finalizeEffects(Bunch bunch){
}

/** It applies the external effects to a particle with certain index. */
void ExternalEffects::applyEffects(Bunch bunch,int index){
}

/** It returns the name of the effect to distinguish them later. */
string ExternalEffects::getName(){
	return name;
}
