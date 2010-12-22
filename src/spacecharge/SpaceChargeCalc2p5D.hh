//Calculate the space charge effect of the bunch in the 2.5D 

#ifndef SC_SPACECHARGE_CALC_2P5D_H
#define SC_SPACECHARGE_CALC_2P5D_H

//MPI Function Wrappers
#include "orbit_mpi.hh"
#include "wrap_mpi_comm.hh"

#include <cstdlib>
#include <cmath>

//ORBIT bunch
#include "Bunch.hh"

//pyORBIT utils
#include "CppPyWrapper.hh"
#include "BunchExtremaCalculator.hh"

#include "Grid1D.hh"
#include "Grid2D.hh"
#include "PoissonSolverFFT2D.hh"
#include "BaseBoundary2D.hh"

using namespace std;

class SpaceChargeCalc2p5D: public OrbitUtils::CppPyWrapper
{
public:
	
	/** Constructor */
	SpaceChargeCalc2p5D(int xSize, int ySize, int zSize, double xy_ratio_in);

	SpaceChargeCalc2p5D(int xSize, int ySize, int zSize);
	
	/** Destructor */
	virtual ~SpaceChargeCalc2p5D();
	
	void trackBunch(Bunch* bunch, double length,double pipe_radius);
	virtual void trackBunch(Bunch* bunch, double length, double pipe_radius, BaseBoundary2D* boundary);
	
	/** Returns the 2D rho grid with a transverse density distribution. **/
	Grid2D* getRhoGrid();

	/** Returns the 2D phi grid with a transverse potential. **/
	Grid2D* getPhiGrid();

	/** Returns the 1D grid with a longitudinal density. **/
	Grid1D* getLongGrid();
	
	/** Returns the 1D grid with a derivative of the longitudinal density. **/
	Grid1D* getLongDerivativeGrid();
	
	/** Sets the number of smoothing points to calculate the derivative of the longitudinal density. */
	void setLongAveragingPointsN(int n_points);
	
	/** Returns the number of smoothing points to calculate the derivative of the longitudinal density. */
	int getLongAveragingPointsN();
	
private:
	/** Analyses the bunch and does bining. */
 double bunchAnalysis(Bunch* bunch, double& totalMacrosize, double& x_c, double& y_c, double& a_bunch);	
 
 /** Calculates the derivative of the longitudinal density by using Quadratic Curve Fitting */
 void calculateLongDerivative();
	
protected:
	PoissonSolverFFT2D* poissonSolver;
	Grid2D* rhoGrid;
	Grid2D* phiGrid;
	Grid1D* zGrid;
	Grid1D* zDerivGrid;
	OrbitUtils::BunchExtremaCalculator* bunchExtremaCalc;
	
	double xy_ratio;
	int n_long_avg;
	
	//auxiliary 5x2 array for Quadratic Curve Fitting
	double** S_arr;
};
//end of SC_SPACECHARGE_CALC_2P5D_H
#endif
