/**
  This class calculates the field of uniformly charged ellipsoid by using 
	the symmetric elliptic integral and Carlson formulas for these integrals.
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
	intFuncX0 = new Function();
	intFuncY0 = new Function();
	intFuncZ0 = new Function();
	intFuncX1 = new Function();
	intFuncY1 = new Function();
	intFuncZ1 = new Function();
	intFuncX2 = new Function();
	intFuncY2 = new Function();
	intFuncZ2 = new Function();
	//the number of points
	lambda_function_points0 = 200;
	
	//the number of points
	lambda_function_points1 = 50;
	
	//the number of points
	lambda_function_points2 = 50;

	//Q_total is 1 by dfeault
	Q_total = 1.0;
	
	//the parameter of ellipses
	a = 1.; b  = 1.; c = 1.;
	double r_max = 10.;
	setEllipsoid(a,b,c,r_max);
}

/** Destructor */
UniformEllipsoidFieldCalculator::~UniformEllipsoidFieldCalculator()
{
		delete intFuncX0;
		delete intFuncY0;
		delete intFuncZ0;
		delete intFuncX1;
		delete intFuncY1;
		delete intFuncZ1;
		delete intFuncX2;
		delete intFuncY2;
		delete intFuncZ2;
}

	
/** Sets the half-axis of the ellipsoid */
void UniformEllipsoidFieldCalculator::setEllipsoid(double a_in, double b_in, double c_in, double r_max)
{
	a = a_in; b  = b_in; c = c_in;
	a2 = a*a; b2 = b*b; c2 = c*c;

	lambda_max2 = pow(r_max,2);
	lambda_max0 = pow(2.0*max(max(a,b),c),2);
	lambda_max1 = pow(5.*max(max(a,b),c),2);
	if(lambda_max2 < lambda_max0){
		lambda_max0 = lambda_max2;
		lambda_max1 = lambda_max2;
		intFuncX1->clean();
		intFuncY1->clean();
		intFuncZ1->clean();		
		intFuncX2->clean();
		intFuncY2->clean();
		intFuncZ2->clean();		
	}
	if(lambda_max2 < lambda_max1){
		lambda_max1 = lambda_max2;
		intFuncX2->clean();
		intFuncY2->clean();
		intFuncZ2->clean();			
	}
	//interval from 0 to lambda_max0
	lambda_eps = 0.01*lambda_max0/(lambda_function_points0 - 1);	
	//now integrate the part from lambda to lambda_max for different lambda and put values into functions
	double lambda_step = lambda_max0/(lambda_function_points0 - 1);
	intFuncX0->clean();
	intFuncY0->clean();
	intFuncZ0->clean();
	for(int iL = 0, nL = lambda_function_points0; iL < nL; iL++){
		double lambda = lambda_step*iL;
		intFuncX0->add(lambda, this->integralPhi(a2,b2,c2,lambda));
		intFuncY0->add(lambda, this->integralPhi(b2,a2,c2,lambda));
		intFuncZ0->add(lambda, this->integralPhi(c2,b2,a2,lambda));
	}
	intFuncX0->setConstStep(1);
	intFuncY0->setConstStep(1);
	intFuncZ0->setConstStep(1);
	if(intFuncX0->isStepConst() != 1 || intFuncY0->isStepConst() != 1 || intFuncZ0->isStepConst() != 1){
		int rank = 0;
		ORBIT_MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		if(rank == 0){
			std::cerr << "UniformEllipsoidFieldCalculator::setEllipsoid(...)" << std::endl 
			<< "The Functions for lambda are not equidistant!!! "<< std::endl 
			<< "Stop."<< std::endl;
		}
		ORBIT_MPI_Finalize();
	}
	if(lambda_max1 > lambda_max0){
		lambda_step = (lambda_max1-lambda_max0)/(lambda_function_points1 - 1);
		intFuncX1->clean();
		intFuncY1->clean();
		intFuncZ1->clean();
		for(int iL = 0, nL = lambda_function_points1; iL < nL; iL++){
			double lambda = lambda_max0 + lambda_step*iL;
			intFuncX1->add(lambda, this->integralPhi(a2,b2,c2,lambda)*lambda);
			intFuncY1->add(lambda, this->integralPhi(b2,a2,c2,lambda)*lambda);
			intFuncZ1->add(lambda, this->integralPhi(c2,b2,a2,lambda)*lambda);
		}
		intFuncX1->setConstStep(1);
		intFuncY1->setConstStep(1);
		intFuncZ1->setConstStep(1);
		if(intFuncX1->isStepConst() != 1 || intFuncY1->isStepConst() != 1 || intFuncZ1->isStepConst() != 1){
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
	if(lambda_max2 > lambda_max1){
		lambda_step = (lambda_max2 - lambda_max1)/(lambda_function_points2 - 1);
		intFuncX2->clean();
		intFuncY2->clean();
		intFuncZ2->clean();
		for(int iL = 0, nL = lambda_function_points2; iL < nL; iL++){
			double lambda = lambda_max1 + lambda_step*iL;
			intFuncX2->add(lambda, this->integralPhi(a2,b2,c2,lambda)*lambda);
			intFuncY2->add(lambda, this->integralPhi(b2,a2,c2,lambda)*lambda);
			intFuncZ2->add(lambda, this->integralPhi(c2,b2,a2,lambda)*lambda);
		}
		intFuncX2->setConstStep(1);
		intFuncY2->setConstStep(1);
		intFuncZ2->setConstStep(1);
		if(intFuncX2->isStepConst() != 1 || intFuncY2->isStepConst() != 1 || intFuncZ2->isStepConst() != 1){
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
}

/** Calculates the field components */
void UniformEllipsoidFieldCalculator::calcField(double x,   double y,   double z, 
	double x2,  double y2,  double z2,
	double& ex, double& ey, double& ez)
{
	if((x2/a2+y2/b2+z2/c2) <= 1.){
		ex = Q_total*x*intFuncX0->y(0);
		ey = Q_total*y*intFuncY0->y(0);
		ez = Q_total*z*intFuncZ0->y(0);
		return;
	}
	
	double lambda = this->calcLambda(x,y,z,x2,y2,z2);
	if(lambda < lambda_max0){
		ex = Q_total*x*intFuncX0->getY(lambda);
		ey = Q_total*y*intFuncY0->getY(lambda);
		ez = Q_total*z*intFuncZ0->getY(lambda);
		return;
	} else if(lambda < lambda_max1) {
		ex = Q_total*x*intFuncX1->getY(lambda)/lambda;
		ey = Q_total*y*intFuncY1->getY(lambda)/lambda;
		ez = Q_total*z*intFuncZ1->getY(lambda)/lambda;	
		return;	
	} else if(lambda < lambda_max2) {
		ex = Q_total*x*intFuncX2->getY(lambda)/lambda;
		ey = Q_total*y*intFuncY2->getY(lambda)/lambda;
		ez = Q_total*z*intFuncZ2->getY(lambda)/lambda;	
		return;			
	}
	double r2 = x2+y2+z2;
	double r3 = sqrt(r2)*r2;
	ex = Q_total*x/r3; ey = Q_total*y/r3; ez = Q_total*z/r3;
}

/** Calculates lambda value as a root of eq. x^2/(a^2+s) + y^2/(b^2+s) + z^2/(c^2+s) - 1 = 0 */
double UniformEllipsoidFieldCalculator::calcLambda(double x,   double y,   double z, 
	                                                double x2,  double y2,  double z2)
{
	double lambda_start = 0.;
	double a2_ls, b2_ls, c2_ls;
	a2_ls = a2+lambda_start; b2_ls = b2+lambda_start; c2_ls = c2+lambda_start; 
	double v_start = x2/a2_ls + y2/b2_ls + z2/c2_ls - 1.0;
	double lambda_stop = lambda_max2;
	double v_stop = x2/(a2+lambda_stop) + y2/(b2+lambda_stop) + z2/(c2+lambda_stop) - 1.0;	
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
	return (lambda_start + lambda_stop)/2;
}

/** Calculates integral for int(1.5*(1/(a^2+s))*1/sqrt((a^2+s)*(b^2+s)*(c^2+s)), over s from lambda to infinity */
double UniformEllipsoidFieldCalculator::integralPhi(double a_2, double b_2, double c_2, double lambda)
{
	double x0 = c_2 + lambda;
	double y0 = b_2 + lambda;
	double z0 = a_2 + lambda;
	double x1,y1,z1, mu, nu, X,Y,Z;
	double sum = 0.;
	int nIter = 6;
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

/** Returns the total space charge inside the ellipse. */
double UniformEllipsoidFieldCalculator::getQ()
{
	return Q_total;
}

/** Sets the total space charge inside the ellipse. */
void UniformEllipsoidFieldCalculator::setQ(double Q_in)
{
	Q_total = Q_in;
}


