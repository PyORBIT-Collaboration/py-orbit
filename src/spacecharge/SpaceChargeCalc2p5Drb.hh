/**
 This class calculates the space charge kicks for bunch using 2.5D approach in transverse direction
 and Rick Baartman's approach (RB) to the longitudinal kicks. There is a hope that will be simular to
 the true 3D discription of space charge. 
*/

#ifndef SC_SPACECHARGE_CALC_2P5D_RB_H
#define SC_SPACECHARGE_CALC_2P5D_RB_H

//MPI Function Wrappers
#include "orbit_mpi.hh"
#include "wrap_mpi_comm.hh"

#include <cstdlib>
#include <cmath>

//ORBIT bunch
#include "Bunch.hh"

//pyORBIT utils
#include "CppPyWrapper.hh"

#include "Grid1D.hh"
#include "Grid2D.hh"
#include "PoissonSolverFFT2D.hh"

using namespace std;

class SpaceChargeCalc2p5Drb: public OrbitUtils::CppPyWrapper
{
public:
	
	/** Constructor with the "x to y ratio" parameter. */
	SpaceChargeCalc2p5Drb(int xSize, int ySize, int zSize, double xy_ratio_in);

	/** Constructor with "ratio" parameter equals 1. */
	SpaceChargeCalc2p5Drb(int xSize, int ySize, int zSize);
	
	/** Destructor */
	virtual ~SpaceChargeCalc2p5Drb();
	
	/** Calculates space charge and applies the transverse and 
	    longitudinal SC kicks to the macro-particles in the bunch. */
	void trackBunch(Bunch* bunch, double length, double pipe_radius);
		
private:
	
	/** Analyses the bunch and does bining. */
 double bunchAnalysis(Bunch* bunch, double& totalMacrosize, double& x_c, double& y_c, double& a_bunch);	
	
protected:
	PoissonSolverFFT2D* poissonSolver;
	Grid2D* rhoGrid;
	Grid2D* phiGrid;
	Grid1D* zGrid;
	
	double xy_ratio;
};
//end of SC_SPACECHARGE_CALC_2P5D_RB_H
#endif
