//The base class for Poisson Solvers. It defines the interface for solves

#ifndef SC_POISSON_SOLVER_BASE_3D_H
#define SC_POISSON_SOLVER_BASE_3D_H

//MPI Function Wrappers
#include "orbit_mpi.hh"
#include "wrap_mpi_comm.hh"

#include <cstdlib>
#include <cmath>
#include <string>

//pyORBIT utils
#include "CppPyWrapper.hh"

#include "Grid3D.hh"

using namespace std;

/** 
  The PoissonSolver3D class is used to define a boundary geometry
  and to calculate the potential created by charges on the boundary 
	surface.
*/
    
class PoissonSolver3D: public OrbitUtils::CppPyWrapper
{
public:
	
	/** Constructor */
	PoissonSolver3D(int xSize, int ySize, int zSize);
	
	/** Constructor */
	PoissonSolver3D(int xSize, int ySize, int zSize, 
		              double xMin, double xMax, 
									double yMin, double yMax,
									double zMin, double zMax);
	
  /** Destructor */
  virtual ~PoissonSolver3D();

  /** Returns the grid size in x-direction */
  int getSizeX(); 

  /** Returns the grid size in y-direction */
  int getSizeY();
	
  /** Returns the grid size in z-direction */
  int getSizeZ();

  /** Returns the maximal value of the grid in x-axis */   
  double getMaxX();

  /** Returns the minimal value of the grid in x-axis */   	
  double getMinX();
	
  /** Returns the maximal value of the grid in y-axis */ 	
  double getMaxY();
	
  /** Returns the minimal value of the grid in y-axis */ 	
  double getMinY();
	
  /** Returns the maximal value of the grid in z-axis */ 	
  double getMaxZ();
	
  /** Returns the minimal value of the grid in z-axis */ 	
  double getMinZ();

  /** Returns the grid step along x-axis */
  double getStepX();
	
	/** Returns the grid step along y-axis */
  double getStepY();
	
	/** Returns the grid step along z-axis */
  double getStepZ();

  /** Sets x-grid. This method is virtual, because the setting of limits may involve 
	  * some subclass specific actions.     
		*/
  virtual void setGridX(double xMin, double xMax) = 0;

  /** Sets y-grid. This method is virtual, because the setting of limits may involve 
	  * some subclass specific actions.
		*/
  virtual void setGridY(double yMin, double yMax) = 0;	
	
  /** Sets z-grid. This method is virtual, because the setting of limits may involve 
	  * some subclass specific actions.
		*/
  virtual void setGridZ(double zMin, double zMax) = 0;	

  /** Solves the Poisson problem for an external charge distribution and
	    puts results into an external potential grid
	*/
  virtual void findPotential(Grid3D* rhoGrid,Grid3D*  phiGrid) = 0; 
 
protected:
	
	//checks that sizes of Grid3D instances are the same as solver 
	void checkSizes(Grid3D* rhoGrid,Grid3D*  phiGrid); 

	//should be implemented in subclasses. It should be called in both types of constructors.
	virtual void init(int xSize, int ySize, int zSize, double xMin, double xMax, double yMin, double yMax, double zMin, double zMax) = 0;
	
protected:

  //Grid size
  int xSize_;
  int ySize_;
  int zSize_;

   //Grid array and parameters
  double xMin_, xMax_, yMin_, yMax_, zMin_, zMax_;
  double  dx_  ,  dy_,  dz_;

};

//end of SC_POISSON_SOLVER_BASE_3D_H ifdef
#endif

