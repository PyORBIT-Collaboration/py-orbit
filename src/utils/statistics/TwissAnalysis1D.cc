#include "TwissAnalysis1D.hh"

#include <iostream>
#include <cmath>
#include <cfloat>

/** Constructor */
TwissAnalysis1D::TwissAnalysis1D(): StatMoments2D()
{
}

/** Destructor */
TwissAnalysis1D::~TwissAnalysis1D()
{
}
/** Returns the emittance */
double TwissAnalysis1D::getEmittance()
{
	double x_avg = getStatMomentU(1);
	double xp_avg = getStatMomentUP(1);
	double x2_avg = fabs(getStatMomentU(2) - x_avg*x_avg);
	double xp2_avg = fabs(getStatMomentUP(2) - xp_avg*xp_avg);
	double x_xp_avg = 	getStatMoment(1,1) - x_avg*xp_avg;
  double emitt_rms =  sqrt(x2_avg*xp2_avg - x_xp_avg*x_xp_avg);
	return emitt_rms;
}

/** Returns Twiss alpha */
double TwissAnalysis1D::getAlpha()
{
	double x_avg = getStatMomentU(1);
	double xp_avg = getStatMomentUP(1);
	double x2_avg = fabs(getStatMomentU(2) - x_avg*x_avg);
	double xp2_avg = fabs(getStatMomentUP(2) - xp_avg*xp_avg);
	double x_xp_avg = 	getStatMoment(1,1) - x_avg*xp_avg;
	double emitt2_rms = x2_avg*xp2_avg - x_xp_avg*x_xp_avg;
	if(	emitt2_rms <= 0.) return 0.;
  double emitt_rms =  sqrt(emitt2_rms);
	double alpha = - x_xp_avg/emitt_rms;
	return alpha;
}

/** Returns Twiss beta */
double TwissAnalysis1D::getBeta()
{
	double x_avg = getStatMomentU(1);
	double xp_avg = getStatMomentUP(1);
	double x2_avg = fabs(getStatMomentU(2) - x_avg*x_avg);
	double xp2_avg = fabs(getStatMomentUP(2) - xp_avg*xp_avg);
	double x_xp_avg = 	getStatMoment(1,1) - x_avg*xp_avg;
	double emitt2_rms = x2_avg*xp2_avg - x_xp_avg*x_xp_avg;
	if(	emitt2_rms <= 0.) return 0.;
  double emitt_rms =  sqrt(emitt2_rms);
	double beta = x2_avg/emitt_rms;
	return beta;	
}

/** Returns Twiss gamma */
double TwissAnalysis1D::getGamma()
{	
	double x_avg = getStatMomentU(1);
	double xp_avg = getStatMomentUP(1);
	double x2_avg = fabs(getStatMomentU(2) - x_avg*x_avg);
	double xp2_avg = fabs(getStatMomentUP(2) - xp_avg*xp_avg);
	double x_xp_avg = 	getStatMoment(1,1) - x_avg*xp_avg;
	double emitt2_rms = x2_avg*xp2_avg - x_xp_avg*x_xp_avg;
	if(	emitt2_rms <= 0.) return DBL_MAX;
  double emitt_rms =  sqrt(emitt2_rms);
	double gamma = xp2_avg/emitt_rms;
	return gamma;
}

/** Returns the rms value of u */ 	
double TwissAnalysis1D::getRmsU()
{
	double x_avg = getStatMomentU(1);
	double x2_avg = fabs(getStatMomentU(2) - x_avg*x_avg);
	return x2_avg;
}

/** Returns the rms value of up */ 	
double TwissAnalysis1D::getRmsUP()
{
	double xp_avg = getStatMomentUP(1);
	double xp2_avg = fabs(getStatMomentUP(2) - xp_avg*xp_avg);
	return xp2_avg;
}

