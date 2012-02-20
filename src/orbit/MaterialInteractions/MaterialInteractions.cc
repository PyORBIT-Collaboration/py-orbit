//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//   MaterialInteractions.cc
//
// AUTHOR
//    S. Cousineau
//
// CREATED
//    09/09/2011
//
// DESCRIPTION
//    A class for storing hadron material interaction methods.
//
///////////////////////////////////////////////////////////////////////////
#include "MaterialInteractions.hh"
#include "OrbitConst.hh"
#include "Random.hh"
#include "bessel.hh"

#include <iostream>
#include <cmath>
#include <cfloat>
#include <cstdlib>

/** Constructor */
MaterialInteractions::MaterialInteractions()
{
}

/** Destructor */
MaterialInteractions::~MaterialInteractions()
{
}

///////////////////////////////////////////////////////////////////////////
//
// NAME
//
//   MaterialInteractions::mcsJackson
//
// DESCRIPTION
//   Multiple coulomb scatter a particle through a set distance.
//   Simple formulation from JD Jackson:
//   Classical Electrodynamics, Chapter 13.
//   Units are cgs, esu.
//
// PARAMETERS
//   stepsize: length (thickness) of scattering material {mm}.
//   Z:        atomic number of scattering material.
//   A:        atomic weight of scattering material.
//   rho:      mass density of scattering material.
//   idum:     random number seed.
//   beta:     v/c.
//   pfac:     1+dp/p0.
//   x,y:      horizontal and vertical coordinates {m}.
//   px,py:    horizontal and vertical momenta {radian}.
//
// RETURNS
//   Nothing
//
///////////////////////////////////////////////////////////////////////////

void MaterialInteractions::mcsJackson(double stepsize, double z, double a, double rho, long& idum, double beta, double pfac, double& x, double& y, double& px, double& py){

	//Convert to mm and mrad for this routine.	And density to g/cm3
	x *= 1000.0;
	y *= 1000.0;
	px *= 1000.0;
	py *= 1000.0;
	stepsize *= 1000.0;
	rho /= 1000.0;
	
	double pi = OrbitConst::PI;
	double cvel = OrbitConst::c * 100.0;
	double eesu = OrbitConst::elementary_charge_CGS;
	double hbar = 1.05443e-27;
	double AMU = 1.65979e-24;
	double aBohr = 5.29172e-09;
	double AP = 1.007593;
	double gamma = 1.0 / sqrt(1.0 - beta * beta);
	
	double aMax = 1.4 * aBohr / pow(z, 0.333333);
	double aMin = 1.4e-13 * pow(a, 0.333333);
	
	double N = rho / (a * AMU);
	
	double T = stepsize / 10.0;
	
	double pmom = gamma * beta * cvel * AP * AMU;
	
	double thMin = hbar / (pmom * aMax);
	double thMax = hbar / (pmom * aMin);
	
	double thMin2i = 1.0 / (thMin * thMin);
	double thMax2i = 1.0 / (thMax * thMax);
	double th2iDiff = thMin2i - thMax2i;
	
	double th2s = 2.0 * log(thMax / thMin) / th2iDiff;
	double coeff = 2.0 * z * eesu * eesu / (pmom * beta * cvel);
	
	double sigTot = pi * coeff * coeff * th2iDiff;
	double nColl = N * T * sigTot;
	
	double th2Tot = nColl * th2s;
	
	double probrp = Random::ran1(idum);
	double probxy = 2.0 * pi * Random::ran1(idum);
	
	double angle = sqrt(-th2Tot * log(probrp));
	double anglexMCS = angle * cos(probxy);
	double angleyMCS = angle * sin(probxy);
	
	double xpfac = px / (1000.* pfac);
	double ypfac = py / (1000.* pfac);
	
	double anglex = atan(xpfac) + anglexMCS;
	double angley = atan(ypfac) + angleyMCS;
	
	double tanglex = tan(anglex);
	double tangley = tan(angley);
	
	px = tanglex * (1000.* pfac);
	py = tangley * (1000.* pfac);
	
	double directionfac = sqrt(1.0 + xpfac * xpfac + ypfac * ypfac);
	double zstep = stepsize / directionfac;
	x += zstep * (xpfac + tanglex) / 2.0;
	y += zstep * (ypfac + tangley) / 2.0;

	//Convert back to m and rad.
	x /= 1000.0;
	y /= 1000.0;
	px /= 1000.0;
	py /= 1000.0;
	
}

///////////////////////////////////////////////////////////////////////////
//
// NAME
//
//   MaterialInteractions::getRuthJackson
//
// DESCRIPTION
//   Rutherford Scattering from an atomic nucleus.
//   Simple formulation from JD Jackson:
//   Classical Electrodynamics, Chapter 13.
//   Units are cgs, esu.
//
// PARAMETERS
//   stepsize: length (thickness) of scattering material {mm}.
//   Z:        atomic number of scattering material.
//   A:        atomic weight of scattering material.
//   rho:      mass density of scattering material.
//   idum:     random number seed.
//   beta:     v/c.
//	 trackit: 
//   pfac:     1+dp/p0.
//
//
// RETURNS
//   A double array with thetax and thetay. It also calculates and assigns the 
//   Rutherford cross section rcross. 
//
///////////////////////////////////////////////////////////////////////////

double MaterialInteractions::ruthScattJackson(double stepsize, double z, double a, double rho, long& idum, double beta, int trackit, double pfac, double& thetax, double& thetay){	

	stepsize *= 1000.0; //Convert to mm
	rho /= 1000.0;		//Convert to g/cm3
	
	double theta[2];
	theta[0] = 0.0;
	theta[1] = 0.0;
	
	double pi = OrbitConst::PI;
	double cvel = OrbitConst::c * 100.0;
	double eesu = OrbitConst::elementary_charge_CGS;
	double hbar = 1.05443e-27;
	double AMU = 1.65979e-24;
	double aBohr = 5.29172e-09;
	double AP = 1.007593;
	double gamma = 1.0 / sqrt(1.0 - beta * beta);
	
	double aMax = 1.4 * aBohr / pow(z, 0.333333);
	double aMin = 1.4e-13 * pow(a, 0.333333);
	
	double N = rho / (a * AMU);
	double T = stepsize / 10.0;
	
	double pmom = gamma * beta * cvel * AP * AMU;
	
	double thMin = hbar / (pmom * aMax);
	double thMax = hbar / (pmom * aMin);
	
	double thMin2i = 1.0 / (thMin * thMin);
	double thMax2i = 1.0 / (thMax * thMax);
	double th2iDiff = thMin2i - thMax2i;
	
	double th2s = 2.0 * log(thMax / thMin) / th2iDiff;
	double coeff = 2.0 * z * eesu * eesu / (pmom * beta * cvel);
	
	double sigTot = pi * coeff * coeff * th2iDiff;
	
	double nColl = N * T * sigTot;
	
	double th2Tot = nColl * th2s;
	
	if(thMin < (2.0 * sqrt(th2Tot))) thMin = 2.0 * sqrt(th2Tot);
	
	thMin2i = 1.0 / (thMin * thMin);
	th2iDiff = thMin2i - thMax2i;
	double rcross = 1.e+24 * pi * coeff * coeff * th2iDiff;
	
	if(rcross < 0.0) rcross = 0.0;
	
	thMin *= 1.4;
	thMin2i = 1.0 / (thMin * thMin);
	th2iDiff = thMin2i - thMax2i;
	
	if(trackit != 0)
	{
		double thx = 0.0;
		double thy = 0.0;
		
		if(thMin < thMax)
		{
			double probrp = Random::ran1(idum);
			double probxy = 2.0 * pi * Random::ran1(idum);
			
			double denom2 = probrp * th2iDiff + thMax2i;
			double th = sqrt(1.0 / denom2);
			
			thx = th * cos(probxy);
			thy = th * sin(probxy);
		}
		theta[0] = thx;
		theta[1] = thy;
	}
	
	return rcross;
}  

///////////////////////////////////////////////////////////////////////////
//
// NAME
//
//   MaterialInteractions::momentumKick
//
// DESCRIPTION
//   Generates random, 2D projection of momentum kick.  See Monte 
//   Carlo Methods in PDG, or the K2 code description in the PhD
//   thesis of Nuria Catalan-Lasheras.
//
//
// PARAMETERS
//   t:     magnitude of momentum transfer
//   p:     particle momentum
//	 idum:	seed for random number generator
//
// RETURNS
//   dp: the momentum kick generated in each plane.
//
///////////////////////////////////////////////////////////////////////////


void MaterialInteractions::momentumKick(double t, double p, double& dpx, double& dpy){

	double va, vb, va2, vb2, r2=10., theta;
	long idum = -(unsigned)time(0);
	theta = acos(1 - t/(2*p*p));
	double dp[2];
	dp[0] = 0.0;
	dp[1] = 0.0;	

	while(r2 > 1.)
    {
		va=2.*Random::ran1(idum)-1;
		vb=Random::ran1(idum)-1;
		va2=va*va;
		vb2=vb*vb;
		r2=va2+vb2;
    }
	
	dp[0]=theta * (2.*va*vb)/r2;
	dp[1]=theta * (va2 - vb2)/r2;
	dpx = dp[0];
	dpy = dp[1];

}


///////////////////////////////////////////////////////////////////////////
//
// NAME
//
//   MaterialInteractions::getElast_t
//
// DESCRIPTION
//   Generates random momentum angle for low energy elastic scattering.
//
// PARAMETERS
//   p: particle momentum
//   a: nuclear mass number
//   idum:  random seed.
//
// RETURNS
//   double.
//
///////////////////////////////////////////////////////////////////////////

double MaterialInteractions::elastic_t(double p, double a, long& idum)
{  
	double c = OrbitConst::c;
	double PI = OrbitConst::PI;
	
	int found=0;
	double R, lamda, E, u, sig, cnorm, theta;
	double random, f_x, x=0., angle_cm=0.;
	double M_nuc, E_cm, p_cm, m=.938, h=4.135e-15/1e9;
	double t=0.;
	
	M_nuc = a*.931494;
	E=sqrt(p*p + m*m);
	E_cm = sqrt(m*m + M_nuc*M_nuc + 2*E*M_nuc);             //See PDG.
	p_cm = p*M_nuc/E_cm;
	lamda = h*c/p_cm;
	R = 1.4e-15 * pow(a, 1./3.) + lamda;
	u = 2*R/lamda*sin(PI/2);
	theta = 0.0001;
	
	cnorm = 1./2.*(sqrt(PI)-sqrt(PI)*pow(OrbitUtils::bessj0(u),2)
				   -sqrt(PI)*pow(OrbitUtils::bessj1(u),2))/(sqrt(PI));
	
	random = Random::ran1(idum);
	
	while(theta<=PI)
	{
		x = 2.*R*sin(theta/2.)/lamda;
		f_x=1./2./cnorm*(sqrt(PI)-sqrt(PI)*pow(OrbitUtils::bessj0(x),2)
						 -sqrt(PI)*pow(OrbitUtils::bessj1(x),2))/(sqrt(PI)) - random;
		if( (f_x > -.01) && (f_x < .01) )
		{
			angle_cm = theta;
			theta = PI+1.;
			found=1;
		}
		theta+=1.768e-3;
	}
	
	if(found==0) cout<<"Warning, never found elastic t.\n";
	t=2.*p_cm*p_cm*(1. - cos(angle_cm));
	return t;
}  


///////////////////////////////////////////////////////////////////////////
//
// NAME
//
//   MaterialInteractions::ionEnergyLoss
//
// DESCRIPTION
//   Compute ionization energy loss for a proton.  See Bethe-Bloch
//   equation in the Particle Data Book.
//
// PARAMETERS
//   beta:  relativistic beta.
//   z:     material charge number.
//   a:     atomic number.
//
// RETURNS
//   double. Rate of energy (GeV) loss per stepsize where stepsize is a measure of 
//   mass per unit area (kg/m^2). See PDG handbook. Units are GeV/(kg/m^2)
//
///////////////////////////////////////////////////////////////////////////

double MaterialInteractions::ionEnergyLoss(double beta, double z, double a)
{ 
	double dE, gamma, T_max, m_e=.511, M=938.27, I, arg;
	
	I=z*10.;
	gamma=1/sqrt(1-beta*beta);                  
	T_max=1e6 * (2*m_e*beta*beta*gamma*gamma)/
    (1 + 2*gamma*m_e/M + pow(m_e/M, 2.0));
	arg=2*(m_e*1e6)*beta*beta*gamma*gamma*T_max/(I*I);
	dE=0.307075*z/a/(beta*beta)*(log(arg)/2.0 - beta*beta);
	
	//Unit conversion from MeV/(g/cm^2) to GeV/(kg/m^2)
	dE /= 10000;
	
	return dE;
}  

