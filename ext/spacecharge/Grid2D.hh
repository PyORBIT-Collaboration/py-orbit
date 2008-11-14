//This class repersents a 2D rectangular grid 

#ifndef SC_GRID_2D_H
#define SC_GRID_2D_H

//MPI Function Wrappers
#include "orbit_mpi.hh"
#include "wrap_mpi_comm.hh"

#include <cstdlib>
#include <cmath>

//pyORBIT utils
#include "CppPyWrapper.hh"

//local for this module
#include "Boundary2D.hh"

using namespace std;

/** 
  This class repersents a 2D rectangular grid.
*/
    
class Grid2D: public OrbitUtils::CppPyWrapper
{
public:
	
	/** Constructor with just grid sizes*/
  Grid2D(int xBins, int yBins);
	
	/** Constructor with grid limits and sizes*/
	Grid2D(double xMin, double xMax, int xBins,    
	       double yMin, double yMax, int yBins);
	
	/** Constructor with the Boundary2D instance */
  Grid2D(Boundary2D* boundary2D);   
	
  /** Destructor */
  virtual ~Grid2D();
	
	/** Sets the Boundary2D instance and redefine grid parameters */	
	void setBoundary2D(Boundary2D* boundary2D);
	
	/** Returns the Boundary2D instance. It could be NULL. */	
	Boundary2D* getBoundary2D();
	
	/** Sets all grid points to zero */	
	void setZero();
	
  /** Returns the reference to the 2D array */	
	double** getArr();
	
	/** Returns the interpolated value from the 2D grid */	
	double getValue(double x, double y);
	
	/** Sets the value to the one point of the 2D grid  */	
	void setValue(double value, int ix, int iy);
		
	/** Bins the value into the 2D grid */	
	void binValue(double value, double x, double y);	
	
	/** Calculates gradient at a position (x,y) */	
	void calcGradient(double x, double y, double& ex, double& ey);
	
	/** Calculates gradient at a grid point (ix,iy) */	
	void calcGradient(int ix, int iy, double& ex, double& ey);
	
	/** Finds and sets the potential to the external grid */	
	void findPotential(Grid2D* phiGrid2D);
	
  /** Returns the grid size in x-direction */
  int getBinsX(); 

  /** Returns the grid size in y-direction */
  int getBinsY();

  /** Returns the index and the fraction of the grid's point for particular x.
      The index is a central point in three point interpolation:
      1 <= ind <= (nBins-2)
      The fraction will be : 0 <= frac <= 1.0
	*/
  void getIndAndFracX(double x, int& ind, double& frac);
	
  /** Returns the index and the fraction of the grid's point for particular x.
      The index is a central point in three point interpolation:
      1 <= ind <= (nBins-2)
      The fraction will be : 0 <= frac <= 1.0
	*/	
  void getIndAndFracY(double y, int& ind, double& frac);

  /** Returns the grid point x-coordinate for this index. */   
  double getGridX(int index);

  /** Returns the grid point y-coordinate for this index. */   
  double getGridY(int index);

  /** Returns the grid step along x-axis */
  double getStepX();
	
	/** Returns the grid step along y-axis */
  double getStepY();	
	
  /** Sets x-grid */
  void setGridX(double xMin, double xMax, int nPoints);

  /** Sets y-grid */
  void setGridY(double yMin, double yMax, int nPoints);
	
	private:
		void makeGrid(double* grid, double& step, double zMin,double zMax,int n);
	
	protected:
		
		double** arr_;
		
  //Grid size
  int xBins_;
  int yBins_;
	
	//Grid steps
	double dx_;
	double dy_;
	
	//x,y-axis grids 
	double* xGrid_;
  double* yGrid_; 
	
	//boundary
	Boundary2D* boundary2D_;
	
};

#endif
