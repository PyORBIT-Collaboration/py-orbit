#ifndef SC_GRID_1D_H
#define SC_GRID_1D_H

// MPI Function Wrappers
#include "orbit_mpi.hh"
#include "wrap_mpi_comm.hh"

#include <cstdlib>
#include <cmath>

// ORBIT bunch
#include "Bunch.hh"

// pyORBIT utils
#include "CppPyWrapper.hh"

using namespace std;


class Grid1D:public OrbitUtils::CppPyWrapper
{
public:

  /** Constructor with grid size only */
  Grid1D(int zSize);

  /** Constructor with lattice length */
  Grid1D(int zSize, double length);
	
  /** Constructor with grid size and spatial limits */
  Grid1D(int zSize, double zMin, double zMax);

  /** Destructor */
  virtual ~Grid1D();

  /** Sets z-grid */
  void setGridZ(double zMin, double zMax);

  /** Returns the reference to the 1D array */	
  double* getArr();

  /** Returns the min z grid point value */
  double getMinZ();

  /** Returns the max z grid point value */
  double getMaxZ();

  /** Returns the number of grid points */
  int getSizeZ();

  /** Returns the sum of all grid points */
  double getSum();

  /** Returns the grid step size */
  double getStepZ();

  /** Returns the grid point for index */
  double getGridZ(int index);
  
  /** Returns 1 if (z) is inside the grid region, and 0 otherwise */
  int isInside(double z);

  /** Sets arr_ at all grid points to zero */
  void setZero();
  
  /** Multiply all elements of Grid1D by constant coefficient */ 
  void multiply(double coeff);

  /** Sets value of arr_ at one point on the grid */
  void setValue(double value, int iZ);

  /** Returns value of arr_ at one point on the grid */
  double getValueOnGrid(int iZ);

  /** Returns interpolated value of arr_ */
  double getValue(double z);

  /** Returns smoothed interpolated value of arr_ */
  double getValueSmoothed(double z);

  /** Bins the Bunch to the grid incorporating macrosize */
  void binBunch(Bunch* bunch);
  
  /** Bins the Bunch to the grid along any axis */
  void binBunch(Bunch* bunch, int axis_ind); 

  /** Bins the Bunch to the grid using a smoothing algorithm
      and incorporating macrosize */
  void binBunchSmoothed(Bunch* bunch);
  
  /** Bins the Bunch to the grid along any axis using a smoothing algorithm
      and incorporating macrosize */
  void binBunchSmoothed(Bunch* bunch, int axis_ind);  

  /** Bins the Bunch to the grid giving each macroparticle 
      unit weight */
  void binBunchByParticle(Bunch* bunch);
  
  /** Bins the Bunch along the particular coordinate to the grid giving each macro-particle 
      unit weight */
  void binBunchByParticle(Bunch* bunch, int axis_ind);  

  /** Bins the Bunch to the grid using a smoothing algorithm
      and giving each macroparticle unit weight */
  void binBunchSmoothedByParticle(Bunch* bunch);

   /** Bins the Bunch  along the particular coordinate to the grid using a smoothing algorithm
      and giving each macroparticle unit weight */
  void binBunchSmoothedByParticle(Bunch* bunch, int axis_ind); 
  
  /** Bins moment of the Bunch to the grid giving
      each macroparticle unit weight */
  void binBunchMoment(int propindex, Bunch* bunch, double* Moment);

  /** Bins moment of the Bunch to the grid using a smoothing
      algorithm and giving each macroparticle unit weight */
  void binBunchSmoothedMoment(int propindex, Bunch* bunch, double* Moment);

  /** Bins a value to the grid */
  void binValue(double value, double z);

  /** Bins a value to the grid with smoothing */
  void binValueSmoothed(double value, double z);

  /** Bins a moment to the grid */
  void binMoment(double value, double z, double* Moment);

  /** Bins a moment to the grid with smoothing */
  void binMomentSmoothed(double value, double z, double* Moment);

  /** Calculates gradient at a position (z) */
  void calcGradient(double z, double& ez);

  /** Calculates smoothed gradient at a position (z) */
  void calcGradientSmoothed(double z, double& ez);

  /** synchronizeMPI */
  void synchronizeMPI(pyORBIT_MPI_Comm* comm);

private:

  /** Memory allocation and step calculation for dx_ and dy_ */
  void init();

	/** 
			The grid is considered as periodic. The 1st bin is also the next to the last.
			Returns the grid indices and binning/interpolating coefficients for a given z.
			The indices bracket the point of interpolation:
			0 <= iZ0, iZp <= nBins - 1
			The coefficients WZ0 and WZp correspond to iZ0 and iZp 
		*/
  void getIndAndWZ(double z,
                   int& ind0  , int& indp,
                   double& Wz0, double& Wzp);

  /** 
      This is method for interpolation. The grid point responsibility is defined 
      differently for binning and interpolation.   
      Returns the grid index and fractional position for particular z.
      The central index is the central point in smoothed three
      point interpolation:
      0 <= ind <= nBins - 1
      The fraction satisfies -0.5 <= frac <= 0.5 
    */
  void getIndAndWZSmoothed(double z,
                           int& indm   , int& ind0   , int& indp   ,
                           double& Wzm , double& Wz0 , double& Wzp ,
                           double& dWzm, double& dWz0, double& dWzp);
  
  
  /** 
      This is method for binning. The grid point responsibility is defined 
      differently for binning and interpolation.   
      Returns the grid index and fractional position for particular z.
      The central index is the central point in smoothed three
      point interpolation:
      0 <= ind <= nBins - 1
      The fraction satisfies -0.5 <= frac <= 0.5 
    */
  void getBinIndAndWZSmoothed(double z,
                           int& indm   , int& ind0   , int& indp   ,
                           double& Wzm , double& Wz0 , double& Wzp ,
                           double& dWzm, double& dWz0, double& dWzp);  

protected:

  double* arr_;

  /** Grid size */
  int zSize_;

  /** Grid steps */
  double dz_;

  /** grid limits */
  double zMin_, zMax_;
};


#endif
