//This class is a Poisson Solver based on FFT convolution

#ifndef SC_POISSON_SOLVER_FFT_3D_H
#define SC_POISSON_SOLVER_FFT_3D_H

//MPI Function Wrappers
#include "orbit_mpi.hh"
#include "wrap_mpi_comm.hh"

#include <cstdlib>
#include <cmath>
#include <string>
#include <cfloat>

//pyORBIT utils
#include "CppPyWrapper.hh"

//FFTW library header
#include "fftw3.h"

#include "PoissonSolver3D.hh"

using namespace std;

/** 
  The PoissonSolverFFT3D class calculates electrostatic
  potential of a 3D charge distribution using 3D FFT approach.
*/
    
class PoissonSolverFFT3D: public PoissonSolver3D
{
	public:

		/** Constructor with sizes only*/
		PoissonSolverFFT3D(int xSize, int ySize, int zSize);
		
		/** Constructor wit sizes and limits */
		PoissonSolverFFT3D(int xSize, int ySize, int zSize,
			double xMin, double xMax, 
			double yMin, double yMax,
			double zMin, double zMax);
		
		/** Destructor */
		virtual ~PoissonSolverFFT3D();
		
		/** Set number of bunches from both sides for space charge calculations */
		void setNumberOfExternalBunches(int nBunches);

		/** Set distance between external bunches - period */
		void setSpacingOfExternalBunches(double lambda);	
		
		/** Get number of bunches from both sides for space charge calculations */
		int getNumberOfExternalBunches();
		
		/** Get distance between bunches */
		double getSpacingOfExternalBunches();

	  void setGridX(double xMin, double xMax); 	
	  void setGridY(double yMin, double yMax);
	  void setGridZ(double zMin, double zMax);
		void setGridXYZ(double xMin, double xMax, double yMin, double yMax, double zMin, double zMax);
		
		/** Updates the Green function FFT */
		void updateGreenFunction();
		
		/** Solves the Poisson problem for an external charge distribution and
		puts results into an external potential grid
		*/
		void findPotential(Grid3D* rhoGrid,Grid3D*  phiGrid); 
		
		
	protected:
		
		//initialize the arrays
		void init(int xSize,   int ySize,int zSize, 
			double xMin, double xMax, 
			double yMin, double yMax,
			double zMin, double zMax);
		
		//define green functions table
		void _defineGreenF();
		
	protected:
		
		//Twice extended grid size to use convolution method
		int xSize2_;
		int ySize2_; 
		int zSize2_; 
		
		//Green function 
		double*** greensF_;
		
		//FFT arrays
		double* in_;
		double* in_res_;
		fftw_complex* out_green_;
		fftw_complex* out_;
		fftw_complex* out_res_;
		
		fftw_plan planForward_greenF_;
		fftw_plan planForward_;
		fftw_plan planBackward_;
		
		//Number of bunches from both sides that should be taken into account.
		//It defines the how many components we will add to Green function.
		//This number will be an even number.
		int nBunches_;
		
		//The distance between centers of the neighboring bunches
		double lambda_;
};
//end of SC_POISSON_SOLVER_FFT_3D_H ifdef
#endif

