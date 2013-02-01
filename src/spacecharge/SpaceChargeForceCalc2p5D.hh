//Calculate the space charge effect of the bunch in the 2.5D 

#ifndef SC_SPACEFORCECHARGE_CALC_2P5D_H
#define SC_SPACECFORCEHARGE_CALC_2P5D_H

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
#include "ForceSolverFFT2D.hh"
#include "BaseBoundary2D.hh"

using namespace std;

class SpaceChargeForceCalc2p5D: public OrbitUtils::CppPyWrapper
{
public:
	
	/** Constructor */
	SpaceChargeForceCalc2p5D(int xSize, int ySize, int zSize, double xy_ratio_in);

	SpaceChargeForceCalc2p5D(int xSize, int ySize, int zSize);
	
	/** Destructor */
	virtual ~SpaceChargeForceCalc2p5D();
	
	/** Calculates space charge and applies the transverse and 
	longitudinal SC kicks to the macro-particles in the bunch. */
	void trackBunch(Bunch* bunch, double length);
	
	/** Returns the 2D rho grid with a transverse density distribution. **/
	Grid2D* getRhoGrid();

	/** Returns the 2D horizontal force grid with a transverse force. **/
	Grid2D* getForceGridX();
	
	/** Returns the 2D vertical force grid with a transverse force. **/
	Grid2D* getForceGridY();

	/** Returns the 1D grid with a longitudinal density. **/
	Grid1D* getLongGrid();	
	
private:
	/** Analyses the bunch and does bining. */
 void bunchAnalysis(Bunch* bunch, double& totalMacrosize); 
	
protected:
	ForceSolverFFT2D* forceSolver;
	Grid2D* rhoGrid;
	Grid2D* phiGrid;
	Grid2D* forceGridX;
	Grid2D* forceGridY;
	Grid1D* zGrid;
	OrbitUtils::BunchExtremaCalculator* bunchExtremaCalc;
	
};
//end of SC_SPACECHARGE_CALC_2P5D_H
#endif
