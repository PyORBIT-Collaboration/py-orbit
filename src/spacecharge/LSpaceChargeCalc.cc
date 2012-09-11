/////////////////////////////////////////////////////////////////////////////
//
// FILE NAME
//   LSpaceCharge.cc
//
//   05/30/12
//
// DESCRIPTION
//   Calculate the longitudinal space charge effect of the bunch (1D)
//
/////////////////////////////////////////////////////////////////////////////

#include "Grid1D.hh"
#include "BufferStore.hh"
#include "LSpaceChargeCalc.hh"
#include "OrbitConst.hh"
#include <complex>
#include <iostream>
#include <cmath>
#include <cfloat>
#include <complex>
//FFTW library header
#include "fftw3.h"

using namespace OrbitUtils;


LSpaceChargeCalc::LSpaceChargeCalc(double b_a_in, double length_in, int nMacrosMin_in, int useSpaceCharge_in, int nBins_in): CppPyWrapper(NULL)
{
	b_a = b_a_in;
	length = length_in;
	nMacrosMin = nMacrosMin_in;
	useSpaceCharge = useSpaceCharge_in;
	nBins = nBins_in;
	zGrid = new Grid1D(nBins);

	_fftmagnitude = new double[nBins];
	_fftphase = new double[nBins];
	_z = new double[nBins];
	_chi = new double[nBins];
	_zImped_n = new std::complex<double>[nBins];
	//std::complex<double> _zImped_n[nBins];
	//_zImped_n = (fftw_complex *) fftw_malloc(nBins * sizeof(fftw_complex));
	//take this out later and make an explicit assignment from input.
	_in = (fftw_complex *) fftw_malloc(nBins * sizeof(fftw_complex));
	_out = (fftw_complex *) fftw_malloc(nBins * sizeof(fftw_complex));
	_plan = fftw_plan_dft_1d(nBins, _in,  _out, FFTW_FORWARD, FFTW_MEASURE);

}


LSpaceChargeCalc::~LSpaceChargeCalc()
{
	if(zGrid->getPyWrapper() != NULL)
	{
		Py_DECREF(zGrid->getPyWrapper());
	}
	else
	{
		delete zGrid;
	}
	
	delete _zImped_n;
	fftw_free(_in);
	fftw_free(_out);
	fftw_destroy_plan(_plan);
}


void LSpaceChargeCalc::initializeImpedance(){
	for (int n=0; n <= nBins/2 ; n++)
	{
		_zImped_n[n] = std::complex<double>(0.0, 0.0);
	}
}

void LSpaceChargeCalc::assignImpedanceValue(int i, double real, double imag){

 	_zImped_n[i+1] = std::complex<double>(real, imag);
}
									
Grid1D* LSpaceChargeCalc::getLongGrid(){
	return zGrid;
}

void LSpaceChargeCalc::trackBunch(Bunch* bunch)
{
	double zmin, zmax;
	double bunchfactor = 0;
	double zfactor = 0;
	
	int nPartsGlobal = bunch->getSizeGlobal();
	if(nPartsGlobal < 2) return;

	//Bin the particles
	bunchExtremaCalc->getExtremaZ(bunch, zmin, zmax);
	double zextra = (length - (zmax - zmin)) / 2.0;
	zmax += zextra;
	zmin = zmax - length;
	zGrid->setGridZ(zmin, zmax);
	zGrid->setZero();
	zGrid->binBunchSmoothedCount(bunch, length);
	zGrid->synchronizeMPI(bunch->getMPI_Comm_Local());

	if (nPartsGlobal < nMacrosMin) return;

	//FFT the beam
	for(int i=0; i < nBins; i++)
	{
		_in[i][0] = zGrid->getValueOnGrid(i);
		_in[i][1] = 0.;
	}

	fftw_execute(_plan);

	//Find the magnitude and phase
	double realPart, imagPart;
	_fftmagnitude[0] = _out[0][0]/(double)nBins;
	_fftphase[0] = 0.;

	for (int i=1; i< nBins/2; i++)
	{
		realPart = _out[i][0]/((double)nBins);
		imagPart = _out[i][1]/((double)nBins);
		_fftmagnitude[i] = sqrt(realPart*realPart + imagPart*imagPart);
		_fftphase[i] = atan2(imagPart,realPart);
	}

	// Figure space charge part of kick:
	SyncPart* sp = bunch->getSyncPart();

	double zSpaceCharge_n = 0.;  // if useSpaceCharge = 0
	//permeability of free space
	double mu_0 = 4.0 * OrbitConst::PI * 1.e-07;
	double _z_0 = OrbitConst::c * mu_0;
	// Else, set positive since space charge is capacitive (Chao convention)
	if(useSpaceCharge != 0)
		zSpaceCharge_n = _z_0 * (1.0 + 2.0*log(b_a)) /
		  (2* sp->getBeta() * pow(sp->getGamma(),2));
	
	
	for (int n=1; n <= nBins/2; n++)
	{
		_z[n] = n * sqrt( pow(std::real(_zImped_n[n]), 2.0) +
				  pow( (std::imag(_zImped_n[n]) +
				  zSpaceCharge_n), 2.0));
		_chi[n] = atan2((std::imag(_zImped_n[n]) + zSpaceCharge_n),
				(std::real(_zImped_n[n])) );
	}

	//Convert charge to current for a single macroparticle per unit bin length
	double charge2current = bunch->getCharge() * bunch->getMacroSize() *
	  OrbitConst::elementary_charge_MKS * sp->getBeta() * OrbitConst::c /
	  (length/nBins);

	//Calculate and add the kick to macroparticles
	double kick, angle, position;

	//Don't do it if there are not enough particles.
	if(bunch->getSize() < nMacrosMin) return;

	//Convert particle longitudinal coordinate to phi
	double philocal;
	double z;
	double phi[bunch->getSize()];
	double** coords = bunch->coordArr();
	for (int j = 0; j < bunch->getSize(); j++)
	{
		z = bunch->z(j);
		philocal = (z/length)*2*OrbitConst::PI;
		//Handle cases where the longitudinal coordinate is
		//outside of the user specified length;
		if(philocal < -OrbitConst::PI) philocal += 2 * OrbitConst::PI;
		if(philocal > OrbitConst::PI) philocal -= 2 * OrbitConst::PI;

		double dE = _kick(philocal) * (-1e-9) *
		  bunch->getCharge() * charge2current;
		coords[j][5] += dE;
	}
}


///////////////////////////////////////////////////////////////////////////
//
// NAME
//    LSpaceChargeCalc::_kick
//
// DESCRIPTION
//    Returns the longitudinal space charge kick to a macroparticle.
//
// PARAMETERS
//    angle - the phase angle of the particle (rad)
//
// RETURNS
//    Nothing.
//
///////////////////////////////////////////////////////////////////////////

double LSpaceChargeCalc::_kick(double angle)
{
	// n=0 term has no impact (constant with phi)
	// f(phi) = _FFTMagnitude(1) + sum (n=2 -> N/2) of
	//   [2* _FFTMagnitude(i) * cos(phi*(n-1) + _FFTPhase(n) + _chi(n))]
	// Had to convert sign on angle and fft_phase in order to agree with previous SuperCode ORBIT

	double dummy=0.;
	double cosArg;

	for (int n=1; n < nBins/2; n++)
	{
		cosArg = n * (angle + OrbitConst::PI) +
		  _fftphase[n] + _chi[n];
		dummy += 2 * _fftmagnitude[n] * _z[n] * cos(cosArg);
	}
	
	return dummy;
}
