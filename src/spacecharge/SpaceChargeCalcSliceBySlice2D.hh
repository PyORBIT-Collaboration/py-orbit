//Calculate the space charge effect of the bunch in 2D slice by slice

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

#include "Grid2D.hh"
#include "Grid3D.hh"
#include "PoissonSolverFFT2D.hh"
#include "BaseBoundary2D.hh"

using namespace std;

class SpaceChargeCalcSliceBySlice2D: public OrbitUtils::CppPyWrapper
{
public:
	
	/** Constructor */
	SpaceChargeCalcSliceBySlice2D(int xSize, int ySize, int zSize, double xy_ratio_in);

	SpaceChargeCalcSliceBySlice2D(int xSize, int ySize, int zSize);
	
	/** Destructor */
	virtual ~SpaceChargeCalcSliceBySlice2D();
	
	/** Calculates space charge and applies the transverse 
	SC kicks to the macro-particles in the bunch. */
	void trackBunch(Bunch* bunch, double length, BaseBoundary2D* boundary);
	
	/** Returns the 3D rho grid with a transverse density distribution. */
	Grid3D* getRhoGrid();

	/** Returns the 3D phi grid with a transverse potential. */
	Grid3D* getPhiGrid();

	
private:
	/** Analyses the bunch and does bining. */
	void bunchAnalysis(Bunch* bunch, double& totalMacrosize, BaseBoundary2D* boundary); 
	
	/** copy a 2D slice of a 3D Grid to another 2D Grid. */
	void copySlice2DtoGrid2D(Grid3D* SourceGrid3D, int iz, Grid2D* TargetGrid2D);
	
	/** copy a 2D Grid to a slice of a 3D Grid. */
	void copyGrid2DtoSlice2D(Grid2D* SourceGrid2D, Grid3D* TargetGrid3D, int iz);
	
	/** Set a slice of a 3D Grid to zero. */
	void setSlice2DZero(Grid3D* TargetGrid3D, int iz);
	
	Grid2D* rhoGrid2D_tmp;
	Grid2D* phiGrid2D_tmp; 
	
protected:
	PoissonSolverFFT2D* poissonSolver;
	Grid3D* rhoGrid3D;	
	Grid3D* phiGrid3D;	
	OrbitUtils::BunchExtremaCalculator* bunchExtremaCalc;
	
	double xy_ratio;
};
//end of SC_SPACECHARGE_CALC_2P5D_H
#endif
