/**
   This class represents a Parmila type RF gap. It acts on the coordinates 
   of the particle by using the transit time factors. The model includes  
   non-linearity in transverse direction.
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
	
	relative_amplitude = 0.;
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

/** Tracks the Bunch trough the RF gap. The formulas are not ready yet. */	
void RfGapTTF::trackBunch(Bunch* bunch, double frequency, double ampl, double E0TL, double phase){
	bunch->compress();
	SyncPart* syncPart = bunch->getSyncPart();
	double gamma = syncPart->getGamma();
	double beta = syncPart->getBeta();
	double mass = bunch->getMass();
	double charge = bunch->getCharge();
	double eKin_in = syncPart->getEnergy();
	double chargeE0TLsin = charge*E0TL*sin(phase);	
	double delta_eKin = charge*E0TL*cos(phase);
	//calculate params in the middle of the gap
	syncPart->setMomentum(syncPart->energyToMomentum(eKin_in + delta_eKin/2.0));
	double gamma_gap = syncPart->getGamma();
	double beta_gap = syncPart->getBeta();	
	//now move to the end of the gap
	double eKin_out = eKin_in + delta_eKin;
	syncPart->setMomentum(syncPart->energyToMomentum(eKin_out));	
	//the base RF gap is simple - no phase correction. The delta time in seconds
	double delta_time = 0.;	
	syncPart->setTime(syncPart->getTime() + delta_time);	
	double gamma_out =	syncPart->getGamma();
	double beta_out = syncPart->getBeta();	
	double prime_coeff = (beta*gamma)/(beta_out*gamma_out);
	//wave momentum
	double k = 2.0*OrbitConst::PI*frequency/OrbitConst::c;
	double phase_time_coeff = k/beta;
	//transverse focusing coeff
	double cappa = - charge*E0TL*k/(2.0*mass*beta_gap*beta_gap*beta_out*gamma_gap*gamma_gap*gamma_out);
	double d_rp = cappa*sin(phase);
	for(int i = 0, n = bunch->getSize(); i < n; i++){
		//longitudinal-energy part
		bunch->dE(i) =bunch->dE(i)  + chargeE0TLsin*phase_time_coeff*bunch->z(i);		
		bunch->z(i) = bunch->z(i)*beta_out/beta;
		//transverse focusing 
		bunch->xp(i) = bunch->xp(i)*prime_coeff + d_rp*bunch->x(i);
		bunch->yp(i) = bunch->yp(i)*prime_coeff + d_rp*bunch->y(i);		
	}
}	
	

/** 
Sets up the gap parameters: T,S, minimal and maximal beta, 
rf frequency, the gap length,  and the relative amplitude.
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

/** Returns the gap length. */
double RfGapTTF::getLength(){
	return gap_length;
}

/** Returns the realtive amplitude. */
double RfGapTTF::getRelativeAmplitude(){
	return relative_amplitude;
}


