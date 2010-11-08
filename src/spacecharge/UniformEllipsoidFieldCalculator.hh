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
  This class caluclates the field of uniformly charged ellipsoid
*/

   
class UniformEllipsoidFieldCalculator: public OrbitUtils::CppPyWrapper
{
  public:
	
  /** Constructor */
  UniformEllipsoidFieldCalculator();
	
  /** Destructor */
  virtual ~UniformEllipsoidFieldCalculator();
	
	/** Sets the half-axis of the ellipsoid and maximal values of x,y,z */
	void setEllipsoid(double a_in, double b_in, double c_in, double xmax, double ymax, double zmax);

	/** Calculates the field components */
  void calcField(double x,   double y,   double z, 
	                 double x2,  double y2,  double z2,
	                 double& ex, double& ey, double& ez);
  private:
	
		/** Calculates integral for int(1.5*(1/(a^2+s))*1/sqrt((a^2+s)*(b^2+s)*(c^2+s)), over s from lambda to infinity */
		static double integralPhi(double a_2, double b_2, double c_2, double lambda);
		
	private:

		//the half-axis of the ellipsoid
		double a,b,c;
		double a2,b2,c2;
		
		//maximal values for x,y,z
		double x_max,y_max,z_max;
		
		//maximal possible value for root of eq. x^2/(a^2+s) + y^2/(b^2+s) + z^2/(c^2+s) - 1 = 0
		double lambda_max;
		
		//accuracy for the root of eq. x^2/(a^2+s) + y^2/(b^2+s) + z^2/(c^2+s) - 1 = 0
		double lambda_eps;
		
		//number of points in the function of lambda
		int lambda_function_points;
		
		//values of integral for x,y,z axis lambda from 0 to lambda_max = 3*max(x_max^2,y_max^2,z_max^2)
		// for x int((2/(a^2+s))*1/sqrt((a^2+s)*(b^2+s)*(c^2+s)), over s from lambda to infinity
		// for y int((2/(b^2+s))*1/sqrt((a^2+s)*(b^2+s)*(c^2+s)), over s from lambda to infinity
		// for z int((2/(c^2+s))*1/sqrt((a^2+s)*(b^2+s)*(c^2+s)), over s from lambda to infinity
		OrbitUtils::Function* intFuncX; 
		OrbitUtils::Function* intFuncY; 
		OrbitUtils::Function* intFuncZ; 
	
};

#endif
