/**
 This class calculates the space charge kicks for bunch using 3D Poisson Solver.
 The solver implements FFT convolution algorithm to calculate the potential in
 the coordinate system where the bunch is resting. The solver is not parallel in the 
 sense of efficiency, but it is working correctly, and particles can are distributed
 among CPUs (not by the solver).
*/

#ifndef SC_SPACECHARGE_CALC_3D_H
#define SC_SPACECHARGE_CALC_3D_H

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

#include "Grid3D.hh"
#include "PoissonSolverFFT3D.hh"

using namespace std;

class SpaceChargeCalc3D: public OrbitUtils::CppPyWrapper
{
public:
	
	/** Constructor with the 3D grid sizes. */
	SpaceChargeCalc3D(int xSize, int ySize, int zSize);
	
	/** Destructor */
	virtual ~SpaceChargeCalc3D();
	
	/** Calculates space charge and applies 3D kicks to the macro-particles in the bunch. */
	void trackBunch(Bunch* bunch, double length);
		
	/** Returns the 3D rho grid with density distribution. **/
	Grid3D* getRhoGrid();

	/** Returns the 3D phi grid with a potential. **/
	Grid3D* getPhiGrid();
	
	/** Sets the ratio limit for the shape change and Green Function recalculations. */
	void setRatioLimit(double ratio_limit_in);
	
	/** Returns the ratio limit for the shape change and Green Function recalculations. */
	double getRatioLimit();
	
	/** Set number of bunches from both sides for space charge calculations */
	void setNumberOfExternalBunches(int nBunches);

	/** Set frequency of the arrivals of the bunches */
	void setFrequencyOfBunches(double frequency);

	/** Get number of bunches from both sides for space charge calculations */
	int getNumberOfExternalBunches();
	
	/** Get frequency of the arrivals of the bunches */
	double getFrequencyOfBunches();
	
private:
	
	/** Analyses the bunch and does binning. */
	void bunchAnalysis(Bunch* bunch);
 
	/** Analyses the bunch and does binning in the case of accounting for neighboring bunches */
 	void wrappedBunchAnalysis(Bunch* bunch);
	
protected:
	PoissonSolverFFT3D* poissonSolver;
	Grid3D* rhoGrid;
	Grid3D* phiGrid;
	OrbitUtils::BunchExtremaCalculator* bunchExtremaCalc;
	
	double xy_ratio;
	double xz_ratio;
	
	//------------ratio change limit ----------
	//If the shape (x to y and x to z ratios) of 3D region changes more than this
	//limit then we have to change shape and recalculate Green Functions in the Poisson solver	
	double ratio_limit;
	
	//Number of bunches from both sides that should be taken into account.
	//It defines the how many components we will add to Green function.
	//This number will be an even number.
	int nBunches_;
	
	//The frequency of the bunch arrivals in Hz. It defines by the RFQ frequency.
	double frequency_;
};
//end of SC_SPACECHARGE_CALC_3D_H
#endif
