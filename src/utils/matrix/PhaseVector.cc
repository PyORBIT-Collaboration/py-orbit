//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//   PhaseVector.cc
//
// AUTHOR
//    A. Shishlo
//
// CREATED
//    03/12/2008
//
// DESCRIPTION
//    A class for a plain double values vector
//
///////////////////////////////////////////////////////////////////////////
#include "orbit_mpi.hh"

#include "PhaseVector.hh"

#include <cstdlib>

using namespace OrbitUtils;

PhaseVector::PhaseVector(int n_in): CppPyWrapper(NULL)
{
	n = n_in;
  v = (double* ) malloc (sizeof(double)*n);
	zero();
}

PhaseVector::PhaseVector(PhaseVector* vIn): CppPyWrapper(NULL)
{
	n = vIn->size();
  v = (double* ) malloc (sizeof(double)*n);
	vIn->copyTo(this);
}

PhaseVector::~PhaseVector()
{
  free(v);
}

void PhaseVector::zero(){
	for(int i = 0; i < n; i++){
			v[i] = 0.;
	}
}

double& PhaseVector::value(int i){ return v[i];}

double* PhaseVector::getArray(){ return v;}

int PhaseVector::size(){ return n;}

int PhaseVector::mult(double val){
	for(int i = 0; i < n; i++){
		  v[i] *= val;
	}
	return 1;
}

int PhaseVector::add(PhaseVector* vIn){
	if(n != vIn->size()){
		ORBIT_MPI_Finalize("PhaseVector: You try to add a vector with wrong size.");
		return 0;
	}
	double* arr = vIn->getArray();
	for(int i = 0; i < n; i++){
		v[i] += arr[i];
	}
	return 1;
}

int PhaseVector::add(double val){
	for(int i = 0; i < n; i++){
		v[i] += val;
	}
	return 1;
}

double PhaseVector::norm(){
	double v2 = 0.;
	for(int i = 0; i < n; i++){
		v2 += v[i]*v[i];
	}
	return v2;
}

double PhaseVector::dot(PhaseVector* vIn){
	if(n != vIn->size()){
		ORBIT_MPI_Finalize("PhaseVector: You try to calculate a scalar product for vectors with wrong size.");
		return 0.;
	}
	double dot = 0.;
	double* arr = vIn->getArray();
	for(int i = 0; i < n; i++){
		dot += v[i]*arr[i];
	}
	return dot;
}

int PhaseVector::copyTo(PhaseVector* vIn){
	if(n != vIn->size()){
	  ORBIT_MPI_Finalize("PhaseVector: You try to copy to a vector with wrong size.");
		return 0;
	}
	double* arr = vIn->getArray();
	for(int i = 0; i < n; i++){
		arr[i] = v[i];
	}
	return 1;
}

