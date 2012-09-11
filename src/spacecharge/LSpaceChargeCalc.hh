//Calculate the longitudinal space charge effect of the bunch 

#ifndef SC_SPACECHARGE_CALC_L_H
#define SC_SPACECHARGE_CALC_L_H

//MPI Function Wrappers
#include "orbit_mpi.hh"
#include "wrap_mpi_comm.hh"

#include <cstdlib>
#include <cmath>
#include <complex>

//ORBIT bunch
#include "Bunch.hh"

//pyORBIT utils
#include "CppPyWrapper.hh"
#include "BunchExtremaCalculator.hh"
#include "Grid1D.hh"

//FFTW library header
#include "fftw3.h"

using namespace std;

class LSpaceChargeCalc: public OrbitUtils::CppPyWrapper
{
public:
	
	/** Constructor */
	LSpaceChargeCalc(double b_a_in, double length_in, int nMacrosMin_in, int useSpaceCharge_in, int zSize_in);
	
	/** Destructor */
	virtual ~LSpaceChargeCalc();
	
	/** Calculates space charge and applies the transverse and 
	longitudinal SC kicks to the macro-particles in the bunch. */
	void trackBunch(Bunch* bunch);
	
	/** Assigns the real and imaginary peices of the machine impedance for index i**/
	void assignImpedanceValue(int i, double real, double imag);
	
	/** Initialize the machine impedance vector to zero **/
	void initializeImpedance();
	
	/** Routine for calculating the kick to the particle **/
	double _kick(double angle);
	
	/** Returns the 1D grid with a longitudinal density. **/
	Grid1D* getLongGrid();	
	
	
//private:
	double b_a;
	double length;
	int nBins;
	int nMacrosMin;
	int useSpaceCharge;
	
//protected:
	Grid1D* zGrid;
	OrbitUtils::BunchExtremaCalculator* bunchExtremaCalc;
		
	//FFT arrays
	double* _fftmagnitude;
	double* _fftphase;
	double* _z;
	double* _chi;
	
	std::complex<double>* _zImped_n;
	//vector<complex> _zImped_n;
	
	fftw_plan _plan;
	fftw_complex* _in;
	fftw_complex* _out;
	
	//fftw_complex* _zImped_n;
	
};
//end of SC_SPACECHARGE_CALC_2P5D_H
#endif
