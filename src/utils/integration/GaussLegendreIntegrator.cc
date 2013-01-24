//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    GaussLegendreIntegrator.cc
//
// AUTHOR
//    A. Shishlo
//
// CREATED
//    01/23/2013
//
// DESCRIPTION
//    The integrator for the Gauss-Legendre schema. 
//
///////////////////////////////////////////////////////////////////////////
#include <iomanip>

#include "orbit_mpi.hh"
#include "GaussLegendreIntegrator.hh"
#include "gauss_legendre_points.hh"
#include "OU_Function.hh"
#include "OU_SplineCH.hh"

using namespace OrbitUtils;

GaussLegendreIntegrator::GaussLegendreIntegrator(): CppPyWrapper(NULL)
{  
	n_int_points = 1024;
	x0 = 0.;
	x1 = 1.;
	pw_finc = new Function();
  gauss_legendre_generator(n_int_points,x0,x1,pw_finc);
	n_int_points = pw_finc->getSize();
}

GaussLegendreIntegrator::GaussLegendreIntegrator(int nPoints): CppPyWrapper(NULL)
{  
	n_int_points = nPoints;
	x0 = 0.;
	x1 = 1.;
	pw_finc = new Function();
  gauss_legendre_generator(n_int_points,x0,x1,pw_finc);
	n_int_points = pw_finc->getSize();
}

GaussLegendreIntegrator::GaussLegendreIntegrator(int nPoints, double x_from, double x_to): CppPyWrapper(NULL)
{  
	n_int_points = nPoints;
	x0 = x_from;
	x1 = x_to;
	pw_finc = new Function();
  gauss_legendre_generator(n_int_points,x0,x1,pw_finc);
	n_int_points = pw_finc->getSize();
}

GaussLegendreIntegrator::~GaussLegendreIntegrator()
{
  delete pw_finc;
}

void GaussLegendreIntegrator::setnPoints(int nPoints)
{
	n_int_points = nPoints;
	gauss_legendre_generator(n_int_points,x0,x1,pw_finc);
	n_int_points = pw_finc->getSize();
}

int GaussLegendreIntegrator::getnPoints()
{
  return n_int_points;
}

void GaussLegendreIntegrator::setLimits(double x_from, double x_to)
{
	x0 = x_from;
	x1 = x_to;
	gauss_legendre_generator(n_int_points,x0,x1,pw_finc);
	n_int_points = pw_finc->getSize();
}

Function* GaussLegendreIntegrator::getPointsAndWeightFunc()
{
	return pw_finc;
}

double GaussLegendreIntegrator::integral(Function* func)
{
	double sum = 0.;
	for(int i = 0; i < n_int_points; i++){
		sum += func->getY(pw_finc->x(i))*pw_finc->y(i);
	}
	return sum;
}

double GaussLegendreIntegrator::integral(SplineCH* spline)
{
	double sum = 0.;
	for(int i = 0; i < n_int_points; i++){
		sum += spline->getY(pw_finc->x(i))*pw_finc->y(i);
	}
	return sum;
}

