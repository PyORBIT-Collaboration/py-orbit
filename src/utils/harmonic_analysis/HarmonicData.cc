//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    HarmonicData.cc
//
// AUTHOR
//    A. Shishlo
//
// CREATED
//    04/17/2018
//
// DESCRIPTION
//    Keeps the data y(x) for set of x(i), i=1,Ndata values where each x(i) 
//    is a phase between -180 to +180 deg.
//    It also has parameters of the fitting representation of this y(x) data:
//    y_fit(x) = A0 + A1*cos((PI/180.)*(x + phase_shift_1)) + 
//                    A2*cos((PI/180.)*(2*x + parse_shift_2)) +
//                    ...
//    with n_harm the maximal harmonic number in the y_fit function.
//    This class has a method getSumDiff2 that returns
//     sumDiff2 = sum((y(x(i))-y_fit(x(i)))**2, i= 1,Ndata)
//
//    The fast calculation of the sumDiff2 is the main purpose of this 
//    class. This value will be used for a fast preliminary analysis
//    of the phase scan data of RF cavities with one or few RF gaps.
//   
//    This class is used for the harmonics analysis of the scan data for 
//    the RF cavities when we measure the phase response from BPMs.
//
///////////////////////////////////////////////////////////////////////////
#include <cfloat>

#include "orbit_mpi.hh"
#include "HarmonicData.hh"

using namespace OrbitUtils;

HarmonicData::HarmonicData(int order_in, Function* inFunc): CppPyWrapper(NULL)
{  
	if(order_in < 0){
		ORBIT_MPI_Finalize("Orbit Utilites HarmonicData::HarmonicData(order) - order should be >= 0. Stop.");
	}	
	
	// param_arr = [A0,A1,phase1,A2,phase2, ...]
	order = order_in;
	param_arr = new double[2*order+1];
	for (int i = 0; i < (2*order+1); i++){
		param_arr[i] = 0.;
	}
	
	dataFunc = new Function();
	for(int ind = 0; ind < inFunc->getSize(); ind++){
		dataFunc->add(inFunc->x(ind),inFunc->y(ind),inFunc->err(ind));
	}

}

HarmonicData::HarmonicData(HarmonicData* harmonicData): HarmonicData(harmonicData->getOrder(),harmonicData->getDataFunction())
{
}

HarmonicData::~HarmonicData()
{
  delete [] param_arr;
  delete dataFunc;
}

void HarmonicData::setOrder(int order_in)
{
	if(order_in < 0){
		ORBIT_MPI_Finalize("Orbit Utilites HarmonicData::setOrder(order) - order should be >= 0. Stop.");
	}
	
	double* param_new_arr = new double[2*order_in+1];
	for (int i = 0; i < (2*order_in+1); i++){
		param_new_arr[i] = 0.;
	}
	
	int min_order = order_in;
	if(min_order > order){
		min_order = order;
	}
	
	for (int i = 0; i < (2*min_order+1); i++){
		param_new_arr[i] = param_arr[i];
	}	

	delete [] param_arr;
	
	param_arr = param_new_arr;
	order = order_in;
}

int HarmonicData::getOrder()
{
	return order;
}

void HarmonicData::setParam(int index, double val)
{
	if(index < 0 || index > 2*order){
		ORBIT_MPI_Finalize("Orbit Utilites HarmonicData::setParam(index,val) - index out of limits. Stop.");
	}
	param_arr[index] = val;
}

double HarmonicData::getParam(int index)
{
	if(index < 0 || index > 2*order){
		ORBIT_MPI_Finalize("Orbit Utilites HarmonicData::getParam(index) - index out of limits. Stop.");
	}
	return param_arr[index];
}

int HarmonicData::dataSize()
{
	return dataFunc->getSize();
}

double HarmonicData::valueY(int indexX)
{
	return dataFunc->y(indexX);
}

double HarmonicData::valueX(int indexX)
{
	return dataFunc->x(indexX);
}

double HarmonicData::fitValueY(double x)
{
	double sum = param_arr[0];
	double grad_to_rad = PI/180.;
	for(int ind = 0; ind < order; ind++){
		sum += param_arr[2*ind+1]*cos(grad_to_rad*((ind+1)*x +  param_arr[2*ind+2]));
	}
	return sum;
}

void HarmonicData::clean()
{
	for (int i = 0; i < (2*order+1); i++){
		param_arr[i] = 0.;
	}	
	
  dataFunc->clean();
}

double HarmonicData::sumDiff2(){
	double diff = 0.;
	double diff2 = 0.;
	int n_points = dataFunc->getSize();
	for(int ind = 0; ind < n_points; ind++){
		diff = dataFunc->y(ind) - fitValueY(dataFunc->x(ind));
		diff2 += diff*diff;
	}
	return diff2;
}

Function* HarmonicData::getDataFunction(){
	return dataFunc;
}
