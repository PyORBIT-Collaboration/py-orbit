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
	intFuncX = new Function();
	intFuncY = new Function();
	intFuncZ = new Function();
	//the number of points
	lambda_function_points = 800;

	//the parameter of ellipses
	a = 1.; b  = 1.; c = 1.;
	x_max = 10.; y_max = 10.; z_max = 10.;
	setEllipsoid(a,b,c,x_max,y_max,z_max);
}

/** Destructor */
UniformEllipsoidFieldCalculator::~UniformEllipsoidFieldCalculator()
{
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
	//now integrate the part from lambda to lambda_max for different lambda and put values into functions
	double lambda_step = lambda_max/(lambda_function_points - 1);
	double lambda = 0.;	
	double int_valueX = 0.;
	double int_valueY = 0.;
	double int_valueZ = 0.;	
	intFuncX->clean();
	intFuncY->clean();
	intFuncZ->clean();
	for(int iL = 0, nL = lambda_function_points; iL < nL; iL++){
		lambda = lambda_step*iL;
		intFuncX->add(lambda, this->integralPhi(a2,b2,c2,lambda));
		intFuncY->add(lambda, this->integralPhi(b2,a2,c2,lambda));
		intFuncZ->add(lambda, this->integralPhi(c2,b2,a2,lambda));
	}
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

/** Calculates integral for int(1.5*(1/(a^2+s))*1/sqrt((a^2+s)*(b^2+s)*(c^2+s)), over s from lambda to infinity */
double UniformEllipsoidFieldCalculator::integralPhi(double a_2, double b_2, double c_2, double lambda)
{
	double x0 = c_2 + lambda;
	double y0 = b_2 + lambda;
	double z0 = a_2 + lambda;
	double x1,y1,z1, mu, nu, X,Y,Z;
	double sum = 0.;
	int nIter = 5;
	for(int i = 0; i < nIter; i++){
		mu = (x0+y0+3*z0)/5;
		X = 1.0 - x0/mu;
		Y = 1.0 - y0/mu;
		Z = 1.0 - z0/mu;
		nu = sqrt(x0*y0) + sqrt(z0*y0) + sqrt(z0*x0);
		x1 = (x0 + nu)*0.25;
		y1 = (y0 + nu)*0.25;
		z1 = (z0 + nu)*0.25;
		sum += 1.0/(pow(4.,i)*(z0+nu)*sqrt(z0));
		x0 = x1; y0 = y1; z0 = z1;
	}
	mu = (x0+y0+3*z0)/5;
	X = 1.0 - x0/mu;
	Y = 1.0 - y0/mu;
	Z = 1.0 - z0/mu;
	double X2,X3,Y2,Y3,Z2,Z3;
  X2 = X*X; X3 = X2*X; Y2 = Y*Y; Y3 = Y2*Y; Z2 = Z*Z; Z3 = Z2*Z;
	double S2 = X2 + Y2 + 3.0*Z2;
	double S3 = X3 + Y3 + 3.0*Z3;
	double S4 = X2*X2 + Y2*Y2 + 3.0*Z2*Z2;
	double S5 = X2*X3 + Y2*Y3 + 3.0*Z2*Z3;
	sum *= 3;
	sum += (1.0/(pow(4.,nIter)*sqrt(mu*mu*mu)))*(1.0+3*S2/7+S3/3+3*S2*S2/22+3*S4/11+3*S2*S3/13+3*S5/13);
	return sum;
}

