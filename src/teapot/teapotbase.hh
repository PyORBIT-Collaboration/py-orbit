/////////////////////////////////////////////////////////////////////////////
//
// FILE NAME
//   teapotbase.hh
//
// AUTHOR
//   Jeff Holmes, ORNL, jzh@ornl.gov
//   Joshua Abrams, Knox College, jabrams@knox.edu
//   Steven Bunch, University of Tennessee, sbunch2@utk.edu
//
// Modified by Andrei Shishlo
//   12/30/05
//
// Checked by Jeff Holmes
//   02/2012
//
// DESCRIPTION
//   Define elementary functions for different elements
//
/////////////////////////////////////////////////////////////////////////////
#ifndef TEAPOT_BASE_H
#define TEAPOT_BASE_H

#include "Bunch.hh"

namespace teapot_base{

    void init_factorial();
	void delete_factorial();
    
    void rotatexy(Bunch* bunch, double anglexy);
    
    void drifti(Bunch* bunch, int i, double length);
	void drift(Bunch* bunch, double length);

	void multpi(Bunch* bunch, int i, int pole, double kl, int skew);
	void multp(Bunch* bunch, int pole, double kl, int skew);

	void multpfringeIN(Bunch* bunch, int pole, double kl, int skew);
	void multpfringeOUT(Bunch* bunch, int pole, double kl, const int skew);

	void kick(Bunch* bunch, double kx, double ky, double kE);

	void quad1(Bunch* bunch, double length, double kq);
	void quad2(Bunch* bunch, double length);

	void quadfringeIN(Bunch* bunch, double kq);
	void quadfringeOUT(Bunch* bunch, double kq);

	void wedgerotate(Bunch* bunch, double e, int frinout);
	void wedgedrift(Bunch* bunch, double e, int inout);
	void wedgebend(Bunch* bunch, double e, int inout, double rho, int nsteps);

	void bend1(Bunch* bunch, double length, double th);
	void bend2(Bunch* bunch, double length);
	void bend3(Bunch* bunch, double th);
	void bend4(Bunch* bunch, double th);

	void bendfringeIN(Bunch* bunch, double rho);
	void bendfringeOUT(Bunch* bunch, double rho);

	void soln(Bunch* bunch, double length, double B);

	void wedgebendCF(Bunch* bunch, double e, int inout,
                 double rho,
                 int vecnum,
                 std::vector<int>& pole,
                 std::vector<double>& kl,
                 std::vector<int>& skew,
                 int nsteps);

    void RingRF(Bunch* bunch, double ring_length, int harmonic_numb, double voltage, double phase_s);

}

#endif  //TEAPOT_BASE_H
