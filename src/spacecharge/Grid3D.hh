//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    Grid3D.cc
//
// AUTHOR
//    A. Shishlo
//
// CREATED
//    11/05/2010
//
// DESCRIPTION
//    Source code for the class "Grid3D" that is used to operate with 3D arrays
//
///////////////////////////////////////////////////////////////////////////
#ifndef SC_GRID3D_HH
#define SC_GRID3D_HH

#include <iostream> 
#include <cstdlib>

//MPI Function Wrappers
#include "orbit_mpi.hh"
#include "wrap_mpi_comm.hh"

//ORBIT bunch
#include "Bunch.hh"

//pyORBIT utils
#include "CppPyWrapper.hh"

#include "Grid2D.hh"

class Grid3D: public OrbitUtils::CppPyWrapper
{
public:
  //--------------------------------------
  //the public methods of the Grid3D class
  //--------------------------------------
  
  Grid3D(int nX, int nY, int nZ);
	
  virtual ~Grid3D();
	
	/** Returns the reference to the inner 3D array. The array is val[z][x][y].*/
  double*** getArr3D();
	
	/** Returns the reference to one 2D slice of the inner 3D array */
  double**  getSlice2D(int zInd);
  
  /** Returns the reference to Grid2D slice of the inner 3D array */
  Grid2D* getGrid2D(int zInd); 

  /** Returns the grid size in x-direction */
  int getSizeX(); 

  /** Returns the grid size in y-direction */
  int getSizeY();
		
  /** Returns the grid size in z-direction */
  int getSizeZ();	
	
  /** Returns the grid point x-coordinate for this index. */   
  double getGridX(int index);

  /** Returns the grid point y-coordinate for this index. */   
  double getGridY(int index);

  /** Returns the grid point z-coordinate for this index. */   
  double getGridZ(int index);

  /** Returns the grid step along x-axis */
  double getStepX();
	
	/** Returns the grid step along y-axis */
  double getStepY();	

	/** Returns the grid step along z-axis */
  double getStepZ();
	
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
	
  /** Sets the limits for the x-grid */
  void setGridX(double xMin, double xMax);	
	
  /** Sets the limits for the y-grid */
  void setGridY(double yMin, double yMax);	
	
  /** Sets the limits for the z-grid */
  void setGridZ(double zMin, double zMax);	
	
	/** Returns the non-interpolated value on grid*/
	double getValueOnGrid(int ix, int iy, int iz);
	
	/** Sets the value to the one point of the 3D grid  */	
	void setValue(double value, int ix, int iy, int iz);	
	
  /** define all elements as 0.0 */
  void setZero();

  //multiply all elements of Grid3D by constant coefficient
  void multiply(double coeff);
	
  /** Bins the Bunch into the 2D grid. If bunch has a 
	    macrosize particle attribute it will be used. */	
	void binBunch(Bunch* bunch);
	
	/** Bins the value onto grid */
  void binValue(double macroSize, double x, double y, double z);
	
  /** Calculates gradient of Arr3D. gradX = gradient_x(Arr3D), and so on */
  void calcGradient(double x,double& gradX,
	      double y,double& gradY,
		    double z,double& gradZ);
	
  /** Calculates value at the point with coordinates x,y,z */
  double getValue(double x,double y,double z);

  /** returns the sum of all grid points */
  double getSum();

  /** returns the sum of all grid points in the slice with index iZ */
  double getSliceSum(int iZ);

  /** returns the sum of all grid points in the slice with position z */
  double getSliceSum(double z);  

  /** synchronize MPI */
  void synchronizeMPI(pyORBIT_MPI_Comm* pyComm);
	
protected:
  //---------------------------------------
  //the protected methods of the Grid3D class
  //---------------------------------------

  double calcValueOnX(int iX, int iY, int iZ, 
                      double Wxm,double Wx0,double Wxp);

  double calcValueOnY(int iX, int iY, int iZ, 
                      double Wym,double Wy0,double Wyp);

  double calcSheetGradient(int iZ,int iX,int iY,
			   double xm,double x0,double xp,
			   double ym,double y0,double yp);
	
	void getIndAndFracX(double x, int& ind, double& frac);
	
	void getIndAndFracY(double y, int& ind, double& frac);

  void getGridIndAndFrac(double x, int& xInd, double& xFrac,
			 double y, int& yInd, double& yFrac,
			 double z, int& zInd, double& zFrac);
	
protected:
  //---------------------------------------
  //the protected members of the Grid3D class
  //---------------------------------------

	//3D array
  double*** Arr3D;
  
  //Grid2D array
  Grid2D** grid2dArr;

// PRIVATE MEMBERS
//    Arr3D;           holds set of double on each 3D grid points  
//    nZ_,nX_,nY_;     number of bins in x,y,z
//    xMin_,xMax_,yMin_,yMax_,zMin_,zMax_ - the ranges of 3D Grids
//    dx_,dy_,dz_;     the grid steps in x,y,z	
	
  int nZ_,nX_,nY_;
  double xMin_,xMax_,yMin_,yMax_,zMin_,zMax_;
  double dx_, dy_, dz_;

};
#endif

