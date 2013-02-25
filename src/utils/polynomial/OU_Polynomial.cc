//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    Polynomial.cc
//
// AUTHOR
//    A. Shishlo
//
// CREATED
//    01/31/2013
//
// DESCRIPTION
//    Calculates the polynomial:
//    order = 0   y = coef[0]
//    order = 1   y = coef[0] + coef[1]*x^2
//    etc.
//
///////////////////////////////////////////////////////////////////////////
#include "orbit_mpi.hh"
#include "OU_Polynomial.hh"

using namespace OrbitUtils;

Polynomial::Polynomial(int order_in): CppPyWrapper(NULL)
{  
	if(order_in < 0){
		ORBIT_MPI_Finalize("Orbit Utilites Polynomial::Polynomial(order) - order should be >= 0. Stop.");
	}	
	order = order_in;
	coef_arr = new double[order+1];
	for (int i = 0; i < (order+1); i++){
		coef_arr[i] = 0.;
	}
}

Polynomial::~Polynomial()
{
  delete [] coef_arr;
}

void Polynomial::setOrder(int order_in)
{
	if(order_in < 0){
		ORBIT_MPI_Finalize("Orbit Utilites Polynomial::setOrder(order) - order should be >= 0. Stop.");
	}	
	double* coef_new_arr = new double[order_in+1];
	for (int i = 0; i < (order_in+1); i++){
		coef_new_arr[i] = 0.;
	}	
	int min_order = order_in;
	if(min_order > order){
		min_order = order;
	}
	for (int i = 0; i < (min_order+1); i++){
		coef_new_arr[i] = coef_arr[i];
	}	

	delete [] coef_arr;
	coef_arr = coef_new_arr;
	order = order_in;
}

int Polynomial::getOrder()
{
	return order;
}

void Polynomial::setCoef(int index, double val)
{
	if(index < 0 || index > order){
		ORBIT_MPI_Finalize("Orbit Utilites Polynomial::setCoef(index,val) - index out of limits. Stop.");
	}
	coef_arr[index] = val;
}

double Polynomial::getCoef(int index)
{
	if(index < 0 || index > order){
		ORBIT_MPI_Finalize("Orbit Utilites Polynomial::getCoef(index) - index out of limits. Stop.");
	}
	return coef_arr[index];
}

double Polynomial::value(double x)
{
	double t = 1.0;
	double sum = 0.;
	for (int i = 0; i < (order+1); i++){
		sum = sum + coef_arr[i]*t;
		t = t*x;
	}	
	return sum;
}

double Polynomial::derivative(double x)
{
	double t = 1.0;
	double sum = 0.;
	for (int i = 1; i < (order+1); i++){
		sum = sum + coef_arr[i]*i*t;
		t = t*x;
	}		
	return sum;
}

void Polynomial::derivativeTo(Polynomial* derivP)
{
	if(order == 0){
		derivP->setOrder(0);
		derivP->setCoef(0,0.);
		return;
	}
	derivP->setOrder(order - 1);
	for (int i = 1; i < (order+1); i++){
		derivP->setCoef(i-1,getCoef(i)*i);
	}	
}

void Polynomial::copyTo(Polynomial* p)
{
	p->setOrder(order);
	for (int i = 0; i < (order+1); i++){
		p->setCoef(i,getCoef(i));
	}	
}


