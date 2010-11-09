#ifndef UNIFORM_SC_ELLIPSOID_HH
#define UNIFORM_SC_ELLIPSOID_HH

#include <cstdlib>
#include <cmath>

//pyORBIT utils
#include "CppPyWrapper.hh"

//Function from OrbitUtils
#include "OU_Function.hh"

using namespace std;

/**
  This class calculates the field of uniformly charged ellipsoid by using 
	the symmetric elliptic integral and Carlson formulas for these integrals.
*/

   
class UniformEllipsoidFieldCalculator: public OrbitUtils::CppPyWrapper
{
  public:
	
  /** Constructor */
  UniformEllipsoidFieldCalculator();
	
  /** Destructor */
  virtual ~UniformEllipsoidFieldCalculator();
	
	/** Sets the half-axis of the ellipsoid and maximal values of radius */
	void setEllipsoid(double a_in, double b_in, double c_in, double r_max);

	/** Calculates the field components */
  void calcField(double x,   double y,   double z, 
	                 double x2,  double y2,  double z2,
	                 double& ex, double& ey, double& ez);
	
	/** Calculates lambda value as a root of eq. x^2/(a^2+s) + y^2/(b^2+s) + z^2/(c^2+s) - 1 = 0 */
	double calcLambda(double x,   double y,   double z, 
	                  double x2,  double y2,  double z2)	;
	
	/** Returns the total space charge inside the ellipse. */
	double getQ();
	
	/** Sets the total space charge inside the ellipse. */
	void setQ(double Q_in);
	
  private:
		
		/** Calculates integral for int(1.5*(1/(a^2+s))*1/sqrt((a^2+s)*(b^2+s)*(c^2+s)), over s from lambda to infinity */
		static double integralPhi(double a_2, double b_2, double c_2, double lambda);
		
	private:

		//total charge Q in the units of the electron cahrge
		double Q_total;
		
		//the half-axis of the ellipsoid
		double a,b,c;
		double a2,b2,c2;
		
		//accuracy for the root of eq. x^2/(a^2+s) + y^2/(b^2+s) + z^2/(c^2+s) - 1 = 0
		double lambda_eps;
		
		//lambda is a value for root of eq. x^2/(a^2+s) + y^2/(b^2+s) + z^2/(c^2+s) - 1 = 0
		//number of points in the function of lambda for interval 0-2.5*max(a,b,c)
		int lambda_function_points0;
		double lambda_max0;
		//number of points in the function of lambda for interval 2.5-5*max(a,b,c)
		int lambda_function_points1;
		double lambda_max1;
		//number of points in the function of lambda for interval 5-lambda_max
		int lambda_function_points2;
		double lambda_max2;	
			
		//values of integral for x,y,z axis lambda for interval 0 to 0-2.0*max(a,b,c)
		// for x int((2/(a^2+s))*1/sqrt((a^2+s)*(b^2+s)*(c^2+s)), over s from lambda to infinity
		// for y int((2/(b^2+s))*1/sqrt((a^2+s)*(b^2+s)*(c^2+s)), over s from lambda to infinity
		// for z int((2/(c^2+s))*1/sqrt((a^2+s)*(b^2+s)*(c^2+s)), over s from lambda to infinity
		OrbitUtils::Function* intFuncX0; 
		OrbitUtils::Function* intFuncY0; 
		OrbitUtils::Function* intFuncZ0; 

		//values of integral for x,y,z axis lambda for interval 2.5-5*max(a,b,c)
		// for x int((2/(a^2+s))*1/sqrt((a^2+s)*(b^2+s)*(c^2+s)), over s from lambda to infinity
		// for y int((2/(b^2+s))*1/sqrt((a^2+s)*(b^2+s)*(c^2+s)), over s from lambda to infinity
		// for z int((2/(c^2+s))*1/sqrt((a^2+s)*(b^2+s)*(c^2+s)), over s from lambda to infinity
		OrbitUtils::Function* intFuncX1; 
		OrbitUtils::Function* intFuncY1; 
		OrbitUtils::Function* intFuncZ1; 
		
		//values of integral for x,y,z axis lambda for interval 5-lambda_max
		// for x int((2/(a^2+s))*1/sqrt((a^2+s)*(b^2+s)*(c^2+s)), over s from lambda to infinity
		// for y int((2/(b^2+s))*1/sqrt((a^2+s)*(b^2+s)*(c^2+s)), over s from lambda to infinity
		// for z int((2/(c^2+s))*1/sqrt((a^2+s)*(b^2+s)*(c^2+s)), over s from lambda to infinity
		OrbitUtils::Function* intFuncX2; 
		OrbitUtils::Function* intFuncY2; 
		OrbitUtils::Function* intFuncZ2; 		
		
};

#endif
