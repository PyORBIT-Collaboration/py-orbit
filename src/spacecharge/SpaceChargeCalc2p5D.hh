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
	
	/** Calculates space charge and applies the transverse and 
	longitudinal SC kicks to the macro-particles in the bunch. */
	void trackBunch(Bunch* bunch, double length, BaseBoundary2D* boundary);
	
	/** Returns the 2D rho grid with a transverse density distribution. **/
	Grid2D* getRhoGrid();

	/** Returns the 2D phi grid with a transverse potential. **/
	Grid2D* getPhiGrid();

	/** Returns the 1D grid with a longitudinal density. **/
	Grid1D* getLongGrid();	
	
private:
	/** Analyses the bunch and does bining. */
 double bunchAnalysis(Bunch* bunch, double& totalMacrosize, BaseBoundary2D* boundary); 
	
protected:
	PoissonSolverFFT2D* poissonSolver;
	Grid2D* rhoGrid;
	Grid2D* phiGrid;
	Grid1D* zGrid;
	OrbitUtils::BunchExtremaCalculator* bunchExtremaCalc;
	
	double xy_ratio;
};
//end of SC_SPACECHARGE_CALC_2P5D_H
#endif
