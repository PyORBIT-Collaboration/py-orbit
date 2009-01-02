//This class is a almost perfect copy of the Boundary class of the
//original ORBIT. 
//Authors: Jeff Holmes, Slava Danilov, John Galambos, Y.Sato, A.Shishlo

#ifndef SC_BOUNDARY_2D_H
#define SC_BOUNDARY_2D_H

//MPI Function Wrappers
#include "orbit_mpi.hh"
#include "wrap_mpi_comm.hh"

#include <cstdlib>
#include <cmath>
#include <string>

//pyORBIT utils
#include "CppPyWrapper.hh"

//#include "rfftw.h"
//#include "fftw.h"
#include "fftw3.h"

using namespace std;

/** 
  The Boundary2D class is used to define a boundary geometry
  and to calculate the potential created by charges on the boundary 
	surface.
*/
    
class Boundary2D: public OrbitUtils::CppPyWrapper
{
public:

  /** Constructor */
  Boundary2D(int xBins   , int yBins   , double xSize, double ySize,
		         int BPoints, const string& BShape, int BModes );

  /** Destructor */
  virtual ~Boundary2D();

  /** The method calculates an impact position on the surface 
	    and a normal vector at the point of entry. 
			The normal vector is directed into the inner volume. 
			The first vector is a position vector and the second one is a 
			normal to the surface vector.
			If the particle did not cross the surface we do not know what to do
			and we will kill it (return TO_BE_KILLED int value).
	*/
  int impactPoint(double x,  double y,  double z,
	                double px, double py, double pz,
									double* r_v,double* n_v);

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

  /** Returns the maximal value of the grid in x-axis */   
  double getGridMaxX();

  /** Returns the minimal value of the grid in x-axis */   	
  double getGridMinX();
	
  /** Returns the maximal value of the grid in y-axis */ 	
  double getGridMaxY();
	
  /** Returns the minimal value of the grid in y-axis */ 	
  double getGridMinY();

  /** Returns the grid step along x-axis */
  double getStepX();
	
	/** Returns the grid step along y-axis */
  double getStepY();

  /** Calculates additional potential phi at one point */
  double calculatePhi(double x, double y);	
	
	/** Returns the pointers to the FFT array of the Green Function values */
	fftw_complex* getOutFFTGreenF();
	
	
  /** Adds potential from the boundary to the external grid */
  void addBoundaryPotential(double** phisc, 
                           int iXmin, int iXmax, 
                           int iYmin, int iYmax);

	/** Adds potential from the boundary to the grid for all points */
	void addBoundaryPotential(double** phisc);
	
  /** Solves the Poisson problem for an external charge distribution and
	    puts results into an external potential grid
	*/
  void findPotential(double** rhosc, double** phisc); 
 

  /** Finds out if a particle is inside the boundary or not */
  int isInside(double x, double y);

  /** Define external grid parameters from inner ones*/
  void defineExtXYgrid(double &xGridMin, double &xGridMax, 
                       double &yGridMin, double &yGridMax,
                       double &dx      , double &dy );

  /** Define external grid's index limits from inner ones */
  void defineExtXYlimits(int &ixMinB, int &ixMaxB, 
                         int &iyMinB, int &iyMaxB);

  /** Returns the x-grid */
  double* getGridX();

  /** Returns the x-grid */
  double* getGridY();

  /** Returns the index of the boundary shape. 
	    1 - Circle, 2 - Ellipse, 3 - Rectangle
	*/
  int getBoundaryShape();

  /** Returns the boundary limit for x-coordinate */
  double getBoundarySizeX();

  /** Returns the boundary limit for y-coordinate */
  double getBoundarySizeY();

  /** Returns the number of boundary points */
  int getNumbBPoints();

  /** Returns the x-coordinate of the boundary point with an index i */
  double getBoundPointX(int i);
	
	/** Returns the y-coordinate of the boundary point with an index i */
  double getBoundPointY(int i);

  ///debug method
  void checkBoundary();

  ///public static members
public:
  const static int IS_INSIDE;
  const static int IS_OUTSIDE;
  const static int TO_BE_KILLED;	

protected:

  ///INITIALIZATION
  void init(int xBins,    int yBins, 
            double xSize, double ySize,
            int BPoints, const string& Shape, int BModes);

  /** Calculates inverse matrix */
  void  _gaussjinv(double **a, int n);

  /** Calculates all of the LSQM functions at one point */
  double* lsq_fuctions(double x, double y);

  /** Defines the FFT of the Green Function */
  void _defineGreenF();

  /** Defines LSQM matrix */
  void _defineLSQMatrix();

  /** Defines LSQM coefficients */
  void _defineLSQcoeff();

  /** Sets array with the interpolation coefficients for boundary points */
  void _setInterpolationCoeff();


protected:

  //Grid size
  int xBins_;
  int yBins_;

  //Twice extended grid size to use convolution method
  int xBins2_;
  int yBins2_; 

  //Green function 
  double** greensF_;
  
  //FFT arrays
  double* in_;
  double* in_res_;
  fftw_complex* out_green_;
  fftw_complex* out_;
  fftw_complex* out_res_;

	fftw_plan planForward_greenF_;
  fftw_plan planForward_;
  fftw_plan planBackward_;

  //boundary points
  //BPoints_ - number of boundary points
  //shape of the boundary  BShape_ = 1 - Circle 2 - Ellipse 3 - Rectangle
  int BPoints_;
  int BShape_;
	int BModes_;

  //Size of the bounding curve from 0 to max [m] or [cm] or [mm] 
  //The center is the point with X,Y coordinates (xSize_/2.,ySize_/2.)
  double xSize_;
  double ySize_; 

  //R_cirle_                - radius of the circle, 
  //BPa_,BPb_               - ellipse parameters, 
  //BPx_length_,BPy_width_  - rectangle parameters
  double R_cirle_,BPa_,BPb_,BPx_length_,BPy_width_;
  double BPrnorm_;
  
  double* BPx_;
  double* BPy_;
  int* iBPx_;
  int* iBPy_;  
  double* theta_;
  double* BPphi_; //external potential at the boundary points   
  
  //Grid array and parameters
  double xGridMin_, xGridMax_, yGridMin_, yGridMax_;
  double  dx_  ,  dy_;

  //Min and max indexes of XY-plane's grid to operate with potential
  int ixMinB_,ixMaxB_,iyMinB_,iyMaxB_;
  
  double* xGrid_;
  double* yGrid_; 

  //LSQM array and matrix 
  //2*BPModes_+1 - number of functions in LSQM 
  int BPModes_; 
  double* func_vector_;
  double* cof_vector_;
  
	int* ipiv_tmp;
	int* indxr_tmp;
	int* indxc_tmp;	
	
  double** LSQ_matrix_;
  double** tmp_matrix_;

  //number of boundary points
  int nBpoints_;

  //array with the interpolation coefficients for boundary points
  //(former Wxm, Wx0, Wxp, Wym, Wy0, Wyp for 9-points scheme)
  double** W_coeff_;

  //PI = 3.141592653589793238462643383279502
  const static double PI;

};

#endif
