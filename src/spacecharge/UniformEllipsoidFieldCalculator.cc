/**
  This class caluclates the field of uniformly charged ellipsoid
	*/

#include <iostream>

//MPI Function Wrappers
#include "orbit_mpi.hh"
#include "wrap_mpi_comm.hh"

#include "UniformEllipsoidFieldCalculator.hh"
#include "gauss_legendre_points.hh"

using namespace OrbitUtils;

//macros for max and min
#ifndef max
	#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#ifndef min
	#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif

/** Constructor. There is no parameters */
UniformEllipsoidFieldCalculator::UniformEllipsoidFieldCalculator(): CppPyWrapper(NULL)
{

  uiFunc = new Function();
	intFuncX = new Function();
	intFuncY = new Function();
	intFuncZ = new Function();
	//sets the default integration point number to 128
	setIntegralPointsNumber(64);
	//the number of points
	lambda_function_points = 200;
	//the upper limit of integration lambda_max*10^order_of_integration_limit
	order_of_integration_limit = 3;
	
	//the parameter of ellipses
	a = 1.; b  = 1.; c = 1.;
	x_max = 10.; y_max = 10.; z_max = 10.;
	setEllipsoid(a,b,c,x_max,y_max,z_max);
}

/** Destructor */
UniformEllipsoidFieldCalculator::~UniformEllipsoidFieldCalculator()
{
    delete uiFunc;
		delete intFuncX;
		delete intFuncY;
		delete intFuncZ;
}

	
/** Sets the half-axis of the ellipsoid */
void UniformEllipsoidFieldCalculator::setEllipsoid(double a_in, double b_in, double c_in,
	                                                 double xmax, double ymax, double zmax)
{
	a = a_in; b  = b_in; c = c_in;
	a2 = a*a; b2 = b*b; c2 = c*c;
	x_max = xmax; y_max = ymax; z_max = zmax;
	
	lambda_max = pow(max(max(x_max,y_max),z_max),2);
	lambda_eps = 0.01*lambda_max/(lambda_function_points - 1);
	//calculate integrals from lambda_max to lambda_max*10^order_of_integration_limit	
	double int_valueX0 = 0.;
	double int_valueY0 = 0.;
	double int_valueZ0 = 0.;
	for(int i_order = 0; i_order < order_of_integration_limit; i_order++){
		int_valueX0 += 1.5*this->integralPhi(a2,b2,c2,lambda_max*pow(10.,i_order),lambda_max*pow(10.,i_order+1));
		int_valueY0 += 1.5*this->integralPhi(b2,a2,c2,lambda_max*pow(10.,i_order),lambda_max*pow(10.,i_order+1));
		int_valueZ0 += 1.5*this->integralPhi(c2,b2,a2,lambda_max*pow(10.,i_order),lambda_max*pow(10.,i_order+1));
	}	
	
	//now integrate the part from lambda to lambda_max for different lambda and put values into functions
	double lambda_step = lambda_max/(lambda_function_points - 1);
	double lambda = 0.;	
	double int_valueX = 0.;
	double int_valueY = 0.;
	double int_valueZ = 0.;	
	intFuncX->clean();
	intFuncY->clean();
	intFuncZ->clean();
	for(int iL = 0, nL = (lambda_function_points-1); iL < nL; iL++){
		lambda = lambda_step*iL;
		int_valueX = 1.5*this->integralPhi(a2,b2,c2,lambda,lambda_max);
		int_valueY = 1.5*this->integralPhi(b2,a2,c2,lambda,lambda_max);
		int_valueZ = 1.5*this->integralPhi(c2,b2,a2,lambda,lambda_max);
		intFuncX->add(lambda, int_valueX + int_valueX0);
		intFuncY->add(lambda, int_valueY + int_valueY0);
		intFuncZ->add(lambda, int_valueZ + int_valueZ0);
	}
	lambda = lambda_step*(lambda_function_points-1);
	intFuncX->add(lambda, int_valueX0);
	intFuncY->add(lambda, int_valueY0);
	intFuncZ->add(lambda, int_valueZ0);	
	intFuncX->setConstStep(1);
	intFuncY->setConstStep(1);
	intFuncZ->setConstStep(1);
	if(intFuncX->isStepConst() != 1 || intFuncY->isStepConst() != 1 || intFuncZ->isStepConst() != 1){
		int rank = 0;
		ORBIT_MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		if(rank == 0){
			std::cerr << "UniformEllipsoidFieldCalculator::setEllipsoid(...)" << std::endl 
			<< "The Functions for lambda are not equidistant!!! "<< std::endl 
			<< "Stop."<< std::endl;
		}
		ORBIT_MPI_Finalize();
	}
}

/** Calculates integral for int(1.5*(1/(a^2+s))*1/sqrt((a^2+s)*(b^2+s)*(c^2+s)), over s from start_s to stop_s */
double UniformEllipsoidFieldCalculator::integralPhi(double a_2, double b_2, double c_2, double start_s, double stop_s)
{
	double interval = stop_s - start_s;
	double integral = 0.;
	double a_s;
	double s;
	for(int i = 0, n = uiFunc->getSize(); i < n; i++){
		s = start_s + interval*uiFunc->x(i);
		a_s = a_2 +s;
		integral += uiFunc->y(i)/(a_s*sqrt(a_s*(b_2+s)*(c_2+s)));
	}
	integral *= interval;
	return integral;
}

/** Calculates the field components */
void UniformEllipsoidFieldCalculator::calcField(double x,   double y,   double z, 
	                                                double x2,  double y2,  double z2,
	                                                double& ex, double& ey, double& ez)
{
	double lambda_start = 0.;
	double a2_ls, b2_ls, c2_ls;
	a2_ls = a2+lambda_start; b2_ls = b2+lambda_start; c2_ls = c2+lambda_start; 
	double v_start = x2/a2_ls + y2/b2_ls + z2/c2_ls - 1.0;
	if(v_start <= 0.){
		ex = x*intFuncX->y(0);
		ey = y*intFuncY->y(0);
		ez = z*intFuncZ->y(0);
		return;
	}
	double lambda_stop = lambda_max;
	double v_stop = x2/(a2+lambda_stop) + y2/(b2+lambda_stop) + z2/(c2+lambda_stop) - 1.0;	
	if(v_stop >= 0.){
		ex = 0.; ey = 0.; ez = 0.;
		return;
	}
	double vp_start = x2/pow((a2+lambda_start),2) + y2/pow((b2+lambda_start),2) + z2/pow((c2+lambda_start),2);
	lambda_stop = lambda_start - v_start*(lambda_start - lambda_stop)/(v_start - v_stop);
	lambda_start = lambda_start + v_start/vp_start;
	while(fabs(lambda_start - lambda_stop) > lambda_eps){
		//std::cout << "debug calcFieldlambda_start ="<<lambda_start<<" lambda_stop="<< lambda_stop <<std::endl;
		a2_ls = a2+lambda_start; b2_ls = b2+lambda_start; c2_ls = c2+lambda_start; 
		v_start = x2/a2_ls + y2/b2_ls + z2/c2_ls - 1.0;
		v_stop = x2/(a2+lambda_stop) + y2/(b2+lambda_stop) + z2/(c2+lambda_stop) - 1.0;
		vp_start = x2/(a2_ls*a2_ls) + y2/(b2_ls*b2_ls) + z2/(c2_ls*c2_ls);
	  lambda_stop = lambda_start - v_start*(lambda_start - lambda_stop)/(v_start - v_stop);
	  lambda_start = lambda_start + v_start/vp_start;		
	}
	double lambda = (lambda_start + lambda_stop)/2;
	ex = x*intFuncX->getY(lambda);
	ey = y*intFuncY->getY(lambda);
	ez = z*intFuncZ->getY(lambda);
}

/** Sets the number of integration point in the Gauss-Legendre quadratures */
void UniformEllipsoidFieldCalculator::setIntegralPointsNumber(int nPoints){	
	uiFunc->clean();	
	gauss_legendre_generator(nPoints,0.,1.0,uiFunc);
}
