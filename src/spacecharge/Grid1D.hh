#ifndef SC_GRID_1D_H
#define SC_GRID_1D_H

//MPI Function Wrappers
#include "orbit_mpi.hh"
#include "wrap_mpi_comm.hh"

#include <cstdlib>
#include <cmath>

//ORBIT bunch
#include "Bunch.hh"

//pyORBIT utils
#include "CppPyWrapper.hh"

using namespace std;

class Grid1D: public OrbitUtils::CppPyWrapper
{
public:
	/** Constructor with just grid sizes*/
	Grid1D(int zSize);
	
	/** Constructor with grid limits and sizes*/
	Grid1D(int zSize, double zMin, double zMax);
	
	/** Destructor */
	virtual ~Grid1D();
	
	/** Sets the value to the one point*/		
	void setValue(double value,int iZ);
	
	/** Returns the value on grid*/
	double getValueOnGrid(int iZ);
	
	/** Returns the interpolated value*/
	double getValue(double z);
	
	/** Returns the min z in the grid points */
	double getMinZ();
	
	/** Returns the max z in the grid points */
	double getMaxZ();
	
	/** Returns the grid step along x-axis */
	double getStepZ();
	
	/** Sets z-grid */
	void setGridZ(double zMin, double zMax);
	
	double getGridZ(int index);
	
	/** Returns the grid size in z-direction */
	int getSizeZ();
	
	/** synchronizeMPI */
	void synchronizeMPI(pyORBIT_MPI_Comm* comm);
	
	/** Bins the Bunch into the grid. If bunch has a macrosize 
	    particle attribute it will be used.*/	
	void binBunch(Bunch* bunch);
	
	/** Bins the value into the grid */
	void binValue(double macroSize,double z);
	
	/** Calculates gradient at a position (z) */
	void calcGradient(double z,double& ez);
	
	/** Sets all grid points to zero */
	void setZero();
	
private:

  /** Returns the index and the fraction of the grid's point for particular z.
	    The index is a central point in three point interpolation:
	    1 <= ind <= (nBins-2)
	    The fraction will be : 0 <= frac <= 1.0
	*/
	void getIndAndFracZ(double z, int& ind, double& frac);	
	
	
	//memory allocation and step calculation for dx_ and dy_ 
	void init();

protected:
	
	double* arr_;
	
	//Grid size
	int zSize_;
	
	//Grid steps
	double dz_;
	
	//grid limits
	double zMin_,zMax_;
};
#endif
