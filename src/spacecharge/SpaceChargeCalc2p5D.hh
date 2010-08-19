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

#include "Grid1D.hh"
#include "Grid2D.hh"
#include "PoissonSolverFFT2D.hh"
#include "BaseBoundary2D.hh"

using namespace std;

class SpaceChargeCalc2p5D: public OrbitUtils::CppPyWrapper
{
public:
	
	/** Constructor */
	SpaceChargeCalc2p5D(int xSize, int ySize, int zSize, double xy_ratio);

	SpaceChargeCalc2p5D(int xSize, int ySize, int zSize,
	             double xMin, double xMax,
	             double yMin, double yMax);
	
	/** Destructor */
	virtual ~SpaceChargeCalc2p5D();
	
	void trackBunch(Bunch* bunch, double length);
	void trackBunch(Bunch* bunch, double length, BaseBoundary2D* boundary);
	double calcMomentumFactor(Bunch* bunch, double length);
	void getBoundaryXY(Bunch* bunch);
	void getBoundaryZ(Bunch* bunch);
	
private:
	//memory allocation and step calculation for dx_ and dy_ 
	void init();
	
protected:
	PoissonSolverFFT2D* poissonSolver;
	Grid2D* rhoGrid;
	Grid2D* phiGrid;
	Grid1D* zGrid;
	
	//Grid size
	int xSize_;
	int ySize_;
	int zSize_;
	
	//Grid limits
	double xMin_,xMax_;
	double yMin_,yMax_;
	double zMin_,zMax_;
	
	double xy_ratio;
};
//end of SC_SPACECHARGE_CALC_2P5D_H
#endif
