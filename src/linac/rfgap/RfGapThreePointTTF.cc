/**
   This class represents a RF gap where the transit time factors are calculated
   by using the second order polynomial approximation of the field on the axis.
   The model includes non-linearity in transverse direction.
   The field is defined by three points at positions -dz, 0., and +dz.
*/

#include <iostream>
#include <cmath>

#include "Bunch.hh"
#include "bessel.hh"
#include "OrbitConst.hh"
#include "RfGapThreePointTTF.hh"

using namespace OrbitUtils;


// Constructor
RfGapThreePointTTF::RfGapThreePointTTF(): CppPyWrapper(NULL)
{
}

// Destructor
RfGapThreePointTTF::~RfGapThreePointTTF()
{
}

/** 
    Tracks the Bunch trough the RF gap as whole using the TTF T and S calculated 
    from polynomial (2nd order) approximation of the axis RF field. It is a thin
    element, and the phase is defined for the center of the gap. The field is defined 
    for 3 points -dz,0, and +dz. Fields are in eV/m, the RF frequency in Hz, and phase
    is in radians. The tracking through the drift -dz and +dz should be done outside of 
    this class.
*/	
void RfGapThreePointTTF::trackBunch(Bunch* bunch, double dz, double Em, double E0, double Ep, double rf_frequency, double phase){
	//energy and mass of particles in GeV and Em,E0, Ep in V/m
	Em = Em/1.0e+9;
	E0 = E0/1.0e+9;
	Ep = Ep/1.0e+9;	
	//field approximation is E(z)= E0*(1+a*z + b*z^2)
	double dz2 = dz*dz;
	double dz3 = dz2*dz;
	double a_param = (Ep-Em)/(2*E0*dz);
	double b_param = (Ep+Em-2*E0)/(2*E0*dz2);
	double E0L = E0*(2*dz+(2.0/3.0)*b_param*dz3);
	bunch->compress();
	SyncPart* syncPart = bunch->getSyncPart();
	double gamma_in = syncPart->getGamma();
	double beta_in = syncPart->getBeta();
	double mass = bunch->getMass();
	double charge = bunch->getCharge();
	double eKin_in = syncPart->getEnergy();
	double kappa_in = 2.0*OrbitConst::PI*rf_frequency/(OrbitConst::c*beta_in);
	double ttf_t = Tttf(dz,a_param,b_param,kappa_in);
	double ttf_s = Sttf(dz,a_param,b_param,kappa_in);
	double delta_eKin = charge*E0L*(ttf_t*cos(phase) - ttf_s*sin(phase));
	//calculate params in the middle of the gap
	syncPart->setMomentum(syncPart->energyToMomentum(eKin_in + delta_eKin/2.0));
	double gamma_gap = syncPart->getGamma();
	double beta_gap = syncPart->getBeta();	
	double kappa_gap = 2.0*OrbitConst::PI*rf_frequency/(OrbitConst::c*beta_gap);
	// T,S,Tp,Sp for kappa = kappa_gap, we assume a small energy spread
  ttf_t = Tttf(dz,a_param,b_param,kappa_gap);
  ttf_s = Sttf(dz,a_param,b_param,kappa_gap);
  double ttf_tp = Tpttf(dz,a_param,b_param,kappa_gap);
  double ttf_sp = Spttf(dz,a_param,b_param,kappa_gap);		
	//the TTF RF gap has the phase correction to simplectic tracking. The delta time in seconds
	double delta_phase = charge*E0L*kappa_gap*(ttf_t*sin(phase) + ttf_s*cos(phase))
	                     /(mass*beta_gap*beta_gap*gamma_gap*gamma_gap*gamma_gap);
	double delta_time = delta_phase/(2.0*OrbitConst::PI*rf_frequency);	
	syncPart->setTime(syncPart->getTime() + delta_time);
	//now move to the end of the gap
	double eKin_out = eKin_in + delta_eKin;
	syncPart->setMomentum(syncPart->energyToMomentum(eKin_out));	
	double gamma_out = syncPart->getGamma();
	double beta_out = syncPart->getBeta();	
	double kappa_out = 2.0*OrbitConst::PI*rf_frequency/(OrbitConst::c*beta_out);
	//(wave momentum)/beta
	double Kr = kappa_gap/gamma_gap;
	double kappa_Kr = kappa_gap/Kr;
	//phase coeff
	double phase_coeff = charge*E0L*kappa_gap/(mass*beta_gap*beta_gap*gamma_gap*gamma_gap*gamma_gap);
  //transverse coeff
	double trans_coeff = charge*E0L/(mass*beta_gap*beta_gap*gamma_gap*gamma_gap);
	double prime_coeff = (beta_in * gamma_in)/(beta_out * gamma_out); 
	double x, y, r, I0,I1, phase_in , phase_out, phase_rf, d_rp;
	double sin_phRf, cos_phRf;
	for(int i = 0, n = bunch->getSize(); i < n; i++){
    x = bunch->x(i);
    y = bunch->y(i);
    r = sqrt(x * x + y * y);
    I0 = bessi0(Kr * r);
    I1 = bessi1(Kr * r);		
		phase_in = bunch->z(i)*kappa_in;
		phase_rf = phase - phase_in;	
		sin_phRf = sin(phase_rf);
		cos_phRf = cos(phase_rf);
		//longitudinal-energy part
		bunch->dE(i) =bunch->dE(i)  + charge*E0L*I0*(ttf_t*cos_phRf - ttf_s*sin_phRf) - delta_eKin;	
		phase_out = phase_in + phase_coeff*(I0*(ttf_tp*sin_phRf + ttf_sp*cos_phRf) +
			                     r*kappa_Kr*I1*(ttf_t*sin_phRf + ttf_s*cos_phRf));
		bunch->z(i) = phase_out/kappa_out;
		//transverse focusing 
		if(r == 0.){
			d_rp = 0.;
		}
		else{
			d_rp = - trans_coeff*I1*(ttf_t*sin_phRf + ttf_s*cos_phRf)/r;
		}
		bunch->xp(i) = bunch->xp(i)*prime_coeff + d_rp*bunch->x(i);
		bunch->yp(i) = bunch->yp(i)*prime_coeff + d_rp*bunch->y(i);		
	}
}	
	
/** 
   It calculates the symmetrical TTF for 3-point approximation of the field. 
   This TTF as functions of the kappa variable = 2*pi*f/(c*beta).
*/
double RfGapThreePointTTF::Tttf(double dz, double a, double b, double kappa)
{ double kappa_dz = kappa*dz;
	double kappa2 = kappa*kappa;
  double ttf = 2*sin(kappa_dz)/kappa *( 1.0 + b*(dz*dz - 2.0/kappa2)) +
               4*b*kappa_dz*cos(kappa_dz)/	(kappa2*kappa);
  ttf = ttf/(2*dz+(2.0/3.0)*b*dz*dz*dz);
  return ttf;
}

/** 
   It calculates the asymmetrical TTF for 3-point approximation of the field. 
   This TTF as functions of the kappa variable = 2*pi*f/(c*beta).
*/
double RfGapThreePointTTF::Sttf(double dz, double a, double b, double kappa)
{ double kappa_dz = kappa*dz;
	double kappa2 = kappa*kappa;
  double ttf = 2*a*(sin(kappa_dz) - kappa_dz*cos(kappa_dz))/kappa2;
  ttf = ttf/(2*dz+(2.0/3.0)*b*dz*dz*dz);
  return ttf;
}


/** 
   It calculates the derivative of the symmetrical TTF for 3-point approximation of the field. 
   This TTF as functions of the kappa variable = 2*pi*f/(c*beta).
*/
double RfGapThreePointTTF::Tpttf(double dz, double a, double b, double kappa)
{ double kappa_dz = kappa*dz;
	double kappa2 = kappa*kappa;
	double ttfp = 2*(dz*cos(kappa_dz)/kappa - sin(kappa_dz)/kappa2)*( 1.0 + b*(kappa_dz*kappa_dz - 2.0)/kappa2)  +
	              8*b*sin(kappa_dz)/(kappa2*kappa2) -
	              4*b*dz*dz*sin(kappa_dz)/kappa2 - 8*b*dz*cos(kappa_dz)/(kappa*kappa2);
  ttfp = ttfp/(2*dz+(2.0/3.0)*b*dz*dz*dz);
  return ttfp;
}

/** 
   It calculates the derivative of the asymmetrical TTF for 3-point approximation of the field. 
   This TTF as functions of the kappa variable = 2*pi*f/(c*beta).
*/
double RfGapThreePointTTF::Spttf(double dz, double a, double b, double kappa)
{ double kappa_dz = kappa*dz;
	double kappa2 = kappa*kappa;
	double ttfp = 2*a*dz*dz*sin(kappa_dz)/kappa - 4*a*(sin(kappa_dz) - kappa_dz*cos(kappa_dz))/(kappa2*kappa);
  ttfp = ttfp/(2*dz+(2.0/3.0)*b*dz*dz*dz);
  return ttfp;
}

