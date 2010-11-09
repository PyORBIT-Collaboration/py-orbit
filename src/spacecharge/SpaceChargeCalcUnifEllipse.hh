/**
 This class calculates the space charge kicks for bunch. It represent the bunch as the set 
 of uniformly charged ellipses in the center of mass of the bunch system. 
 The space charge kick is transformed later into the lab system.   
*/

#ifndef SC_SPACECHARGE_CALC_UNIFORM_ELLIPSE_HH
#define SC_SPACECHARGE_CALC_UNIFORM_ELLIPSE_HH

//MPI Function Wrappers
#include "orbit_mpi.hh"
#include "wrap_mpi_comm.hh"

#include <cstdlib>
#include <cmath>

//ORBIT bunch
#include "Bunch.hh"

//pyORBIT utils
#include "CppPyWrapper.hh"

#include "UniformEllipsoidFieldCalculator.hh"

using namespace std;

class SpaceChargeCalcUnifEllipse: public OrbitUtils::CppPyWrapper
{
public:
	
	/** Constructor with the "x to y ratio" parameter. */
	SpaceChargeCalcUnifEllipse(int nEllipses_in);

	/** Destructor */
	virtual ~SpaceChargeCalcUnifEllipse();
	
	/** Calculates space charge and applies the transverse and 
	    longitudinal SC kicks to the macro-particles in the bunch. */
	void trackBunch(Bunch* bunch, double length);
	
	/** Analyses the bunch and sets up the ellipsoid filed sources */
  void bunchAnalysis(Bunch* bunch);
	
	/** Calculates the electric filed in the center of the bunch sytem. */
	void calculaField(double x,  double y,  double z, double& ex, double& ey, double& ez)	;

private:
	
protected:

	//number of uniform ellipses 
	int nEllipses;
	
	//parameters of the distribution
	double x_center, y_center, z_center;
	double x2_avg, y2_avg, z2_avg;
 double xMin, xMax, yMin, yMax, zMin, zMax;	

	//sizes of the biggest ellipsoid
	double a_ellips, b_ellips, c_ellips;
	double a2_ellips, b2_ellips, c2_ellips;	
	
	//ellipse calculators
	UniformEllipsoidFieldCalculator** ellipsoidCalc_arr;
	
	//total macrosize in each ellipsoid
	double* macroSizesEll_arr;	
	double* macroSizesEll_MPI_arr;	

};
//end of SC_SPACECHARGE_CALC_UNIFORM_ELLIPSE_HH
#endif
