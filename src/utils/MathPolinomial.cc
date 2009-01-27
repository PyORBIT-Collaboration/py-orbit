//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    LorentzTransformationEM.cc
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
//
///////////////////////////////////////////////////////////////////////////
#include "MathPolinomial.hh"

double MathPolinomial::ReHermite(int n, double x){
	
	double a,b,ff;
	
	a=1.0;
	b=2.0*x;
	
	if(n==0) ff=a;
	if(n==1) ff=b;
	if(n>=2)	
		for(int i=2; i<=n;i++)	{

			ff=2.0*(x*b-(i-1)*a);
				a=b;
				b=ff;
	}
return ff;
}



tcomplex MathPolinomial::ComplexHermite(int n, tcomplex x){
	
	
	tcomplex a,b,ff;
	
	a=1.0;
	b=2.0*x;
	
	if(n==0) ff=a;
	if(n==1) ff=b;
	if(n>=2)	
		for(int i=2; i<=n;i++)	{

			ff=2.0*(x*b-double(i-1)*a);
				a=b;
				b=ff;
	}
	
return ff;
	
}



long int MathPolinomial::Factorial(int n){
	
	int a=1;
	
	for(int i=1; i<=n; i++)
		a*=i;
	
	return a;
}
