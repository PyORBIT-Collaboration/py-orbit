//This class repersents a 2D rectangular grid 

#ifndef SC_GRID_2D_H
#define SC_GRID_2D_H

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

/** 
  This class repersents a 2D rectangular grid.
*/
    
class Grid2D: public OrbitUtils::CppPyWrapper
{
public:
	
	/** Constructor with just grid sizes*/
  Grid2D(int xSize, int ySize);
	
	/** Constructor with grid limits and sizes*/
	Grid2D(int xSize, int ySize,
		     double xMin, double xMax,     
	       double yMin, double yMax);

  /** Destructor */
  virtual ~Grid2D();
	
	/** Sets all grid points to zero */	
	void setZero();
	
  /** Returns the reference to the 2D array */	
	double** getArr();
	
	/** Returns the interpolated value from the 2D grid */	
	double getValue(double x, double y);
	
	/** Sets the value to the one point of the 2D grid  */	
	void setValue(double value, int ix, int iy);
		
	/** Bins the Bunch into the 2D grid. If bunch has a macrosize 
	    particle attribute it will be used. 
	*/	
	void binBunch(Bunch* bunch);	
	
	/** Bins the value into the 2D grid */	
	void binValue(double value, double x, double y);	
	
	/** Calculates gradient at a position (x,y) */	
	void calcGradient(double x, double y, double& ex, double& ey);
	
	/** Calculates gradient at a grid point (ix,iy) */	
	void calcGradient(int ix, int iy, double& ex, double& ey);
		
  /** Returns the grid size in x-direction */
  int getSizeX(); 

  /** Returns the grid size in y-direction */
  int getSizeY();
	
	/** Returns 1 if (x,y) is inside the grid region, and 0 otherwise */
	int isInside(double x,double y);

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
	
  /** Returns the maximal value of the grid in x-axis */   
  double getMaxX();

  /** Returns the minimal value of the grid in x-axis */   	
  double getMinX();
	
  /** Returns the maximal value of the grid in y-axis */ 	
  double getMaxY();
	
  /** Returns the minimal value of the grid in y-axis */ 	
  double getMinY();	
	
  /** Sets x-grid */
  void setGridX(double xMin, double xMax);

  /** Sets y-grid */
  void setGridY(double yMin, double yMax);
	
  private:
	
		//memory allocation and step calculation for dx_ and dy_ 
		void init();
		
	protected:
		
		double** arr_;
		
  //Grid size
  int xSize_;
  int ySize_;
	
	//Grid steps
	double dx_;
	double dy_;
	
	//grid limits
	double xMin_,xMax_;
  double yMin_,yMax_; 
	
};

#endif
