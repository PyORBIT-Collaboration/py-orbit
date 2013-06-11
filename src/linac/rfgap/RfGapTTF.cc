/**
   This class represents a Parmila type RF gap. It acts on the coordinates 
   of the particle by using the transit time factors. The model includes  
   non-linearity in transverse direction.
	 TTFs (T,S,Tp,Sp) are funcftions of the kappa variable = 2*pi*f/(c*beta).
*/

#include <iostream>
#include <cmath>

#include "Bunch.hh"
#include "bessel.hh"
#include "OrbitConst.hh"
#include "RfGapTTF.hh"

using namespace OrbitUtils;


// Constructor
RfGapTTF::RfGapTTF(): CppPyWrapper(NULL)
{

	Tttf = new Polynomial(0);
	Sttf = new Polynomial(0);
	Tpttf = new Polynomial(0);
	Spttf = new Polynomial(0);
	
	beta_min = 0.;
	beta_max = 1.0;
	
	rf_frequency = 402.5e+6;
	
	gap_length = 0.;
	
	relative_amplitude = 1.;
}

// Destructor
RfGapTTF::~RfGapTTF()
{
	if(Tttf->getPyWrapper() != NULL){
		Py_XDECREF(Tttf->getPyWrapper());
	} else {
		delete Tttf;
	}

	if(Sttf->getPyWrapper() != NULL){
		Py_XDECREF(Sttf->getPyWrapper());
	} else {
		delete Sttf;
	}

	if(Tpttf->getPyWrapper() != NULL){
		Py_XDECREF(Tpttf->getPyWrapper());
	} else {
		delete Tpttf;
	}

	if(Spttf->getPyWrapper() != NULL){
		Py_XDECREF(Spttf->getPyWrapper());
	} else {
		delete Spttf;
	}

}

/** Tracks the Bunch trough the RF gap as whole using the TTF T and S. 
    There are a lot of possibilities to modify this algorithm, but for now
		we believe that all changes will give a small difference in results if
		the energy spread in the bunch and the energy gain in the gap are small
		relative to the initial energy.
		A curious user could try to redefine how we calculate the energy gain,
		do we use the same TTF for all particles, recursive procedures etc.
		If you find something interesting, please let us know.
*/	
void RfGapTTF::trackBunch(Bunch* bunch, double E0, double phase){
	//energy and mass of particles in GeV and E0 in V/m
	double E0L = E0*relative_amplitude*gap_length/1.0e+9;
	bunch->compress();
	SyncPart* syncPart = bunch->getSyncPart();
	double gamma_in = syncPart->getGamma();
	double beta_in = syncPart->getBeta();
	double mass = bunch->getMass();
	double charge = bunch->getCharge();
	double eKin_in = syncPart->getEnergy();
	double kappa_in = 2.0*OrbitConst::PI*rf_frequency/(OrbitConst::c*beta_in);
	double delta_eKin = charge*E0L*(Tttf->value(kappa_in)*cos(phase) - Sttf->value(kappa_in)*sin(phase));
	//calculate params in the middle of the gap
	syncPart->setMomentum(syncPart->energyToMomentum(eKin_in + delta_eKin/2.0));
	double gamma_gap = syncPart->getGamma();
	double beta_gap = syncPart->getBeta();	
	double kappa_gap = 2.0*OrbitConst::PI*rf_frequency/(OrbitConst::c*beta_gap);
	// T,S,Tp,Sp for kappa = kappa_gap, we assume a small energy spread
  double ttf_t = Tttf->value(kappa_gap);
  double ttf_s = Sttf->value(kappa_gap);
  double ttf_tp = Tpttf->value(kappa_gap);
  double ttf_sp = Spttf->value(kappa_gap);		
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
Sets up the gap parameters: T,S, minimal and maximal beta, 
rf frequency, the gap length,  and the relative amplitude.
TTFs (T,S,Tp,Sp) are funcftions of the kappa variable = 2*pi*f/(c*beta).
*/
void RfGapTTF::setParameters(Polynomial* Tttf_in, Polynomial* Tpttf_in,
	Polynomial* Sttf_in, Polynomial* Spttf_in,
	double beta_min, double beta_max, 
	double rf_frequency, double gap_length, 
	double relative_amplitude)
{

	setT_TTF(Tttf_in);
	setS_TTF(Sttf_in);
	setTp_TTF(Tpttf_in);
	setSp_TTF(Spttf_in);
	
	this->beta_min = beta_min;
	this->beta_max = beta_max;
	this->rf_frequency = rf_frequency;
	this->gap_length = gap_length;
	this->relative_amplitude = relative_amplitude;
	
}

/** Returns the minimal beta. */
double RfGapTTF::getBetaMin(){
	return beta_min;
}

/** Returns the maximal beta. */
double RfGapTTF::getBetaMax(){
	return beta_max;
}

/** Returns T TTF. */
Polynomial* RfGapTTF::getT_TTF(){
	return Tttf;
}

/** Returns S TTF. */
Polynomial* RfGapTTF::getS_TTF(){
	return Sttf;
}

/** Returns Tp TTF. */
Polynomial* RfGapTTF::getTp_TTF(){
	return Tpttf;
}

/** Returns Sp TTF. */
Polynomial* RfGapTTF::getSp_TTF(){
	return Spttf;
}

/** Sets the Polynomial to T TTF. */
void RfGapTTF::setT_TTF(Polynomial* Tttf_in){
	Tttf_in->copyTo(Tttf);
}

/** Sets the Polynomial to S TTF. */
void RfGapTTF::setS_TTF(Polynomial* Sttf_in){
	Sttf_in->copyTo(Sttf);
}

/** Sets the Polynomial to Tp TTF. */
void RfGapTTF::setTp_TTF(Polynomial* Tpttf_in){
	Tpttf_in->copyTo(Tpttf);
}

/** Sets the Polynomial to Sp TTF. */
void RfGapTTF::setSp_TTF(Polynomial* Spttf_in){
	Spttf_in->copyTo(Spttf);
}

/** Returns RF frequency. */
double RfGapTTF::getFrequency(){
	return rf_frequency;
}

/** Returns the gap length.*/
double RfGapTTF::getLength(){
	return gap_length;
}

/** Returns the realtive amplitude. */
double RfGapTTF::getRelativeAmplitude(){
	return relative_amplitude;
}


