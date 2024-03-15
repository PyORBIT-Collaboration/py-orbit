/////////////////////////////////////////////////////////////////////////////
//
// FILE NAME
//   linac_tracking.hh
//
// AUTHOR
//   Andrei Shishlo, ORNL, shishlo@ornl.gov
//   2017.01.24
//
// DESCRIPTION
//   Define elementary functions for tracking through
//   specific linac accelerator elements without using
//   TEAPOT algorithms. It is necessary for the cases 
//   with a huge energy spread in the bunch.
//
/////////////////////////////////////////////////////////////////////////////

#ifndef LINAC_TRACKING_H
#define LINAC_TRACKING_H

#include "Bunch.hh"

namespace linac_tracking
{
	
	/**
		Drifts a particle bunch in a linac. Length < 0 is not allowed.
	*/	
	void linac_drift(Bunch* bunch, double length);
	
	/**
		Quadrupole element one: linear transport matrix for linac.
		kq = quadrupole field strength [m^(-2)]
	*/		
	void linac_quad1(Bunch* bunch, double length, double kq, int useCharge);
	
	
	/**
		Quadrupole element two: nonlinear piece for linac
	*/
	void linac_quad2(Bunch* bunch, double length);
	
	/**
	     Quadrupole element 3: longitudinal field component of the quad field
	     Bz(z) = x*y*dG(z)/dz
	     This function performs the transverse kicks correction, so the length
	     is just a parameter.
	*/		
	void linac_quad3(Bunch* bunch, double length, double dB_dz);	
	
	
	/**
		Kicker element function: each particle will be kicked with (p_synch/p) coefficient
	*/	
	void kick(Bunch* bunch, double kx, double ky, double kE, int useCharge);

}

#endif  //LINAC_TRACKING_H

