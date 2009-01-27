#ifndef MATHPOLINOMIAL_HH_
#define MATHPOLINOMIAL_HH_

//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    LorentzTransformationEM.hh
//
// AUTHOR
//    T. Gorlov
//
// CREATED
//    06/28/2008
//
// DESCRIPTION
//    This class provides Lorentz transformations for the electromagnetic field
//    from the laboratory frame to the particle rest frame.
//    mass - mass of the particle in GeV
//    px,py,pz - momentum of the particle the lab frame in GeV/c
//    E_x,E_y,E_z - components of the electric field V/m (parameters are replaced in place) 
//    B_x,B_y,B_z - components of the magnetic field [T] (parameters are replaced in place)   
//
//    OrbitConst::c in [m/sec]
///////////////////////////////////////////////////////////////////////////
#include <complex>
#include <cmath>
typedef std::complex<double>	tcomplex;



class  MathPolinomial
{
public:
	
  static  double 	ReHermite(int n, double x);
  
  static  tcomplex 	ComplexHermite(int n, tcomplex x);
  
  static long int Factorial(int n);
  
};


#endif /*MATHPOLINOMIAL_HH_*/
