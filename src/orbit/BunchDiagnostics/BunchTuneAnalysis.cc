#include "BunchTuneAnalysis.hh"
#include "SyncPart.hh"
#include "OrbitConst.hh"

#include <iostream>
#include <cmath>
#include <cfloat>
#include <cstdlib>

/** Constructor */
BunchTuneAnalysis::BunchTuneAnalysis(): CppPyWrapper(NULL)
{
	betax = 0;
	alphax = 0;
	etax = 0;
	etapx = 0;
	betay = 0;
	alphay = 0;
}

/** Destructor */
BunchTuneAnalysis::~BunchTuneAnalysis()
{
}

void BunchTuneAnalysis::assignTwiss(double bx, double ax, double dx, double dpx, double by, double ay){
	betax = bx;
	alphax = ax;
	etax = dx;
	etapx = dpx;
	betay = by;
	alphay = ay;
}

/** Performs the Tune analysis of the bunch */
void BunchTuneAnalysis::analyzeBunch(Bunch* bunch){
	
	//initialization
	bunch->compress();
	SyncPart* syncPart = bunch->getSyncPart();
	double beta = syncPart->getBeta();
	double** part_coord_arr = bunch->coordArr();

	if(!bunch->hasParticleAttributes("ParticlePhaseAttributes")){
		cerr<<"adding particle phase information attribute\n";
		std::map<std::string, double> tunemap;
		tunemap.insert(std::make_pair("xLastPhase", 0));
		tunemap.insert(std::make_pair("yLastPhase", 0));
		tunemap.insert(std::make_pair("xLastTune", 0));
		tunemap.insert(std::make_pair("yLastTune", 0));
		tunemap.insert(std::make_pair("xAction", 0));
		tunemap.insert(std::make_pair("yAction", 0));
		bunch->addParticleAttributes("ParticlePhaseAttributes", tunemap);
	}
	
	if(bunch->hasParticleAttributes("ParticlePhaseAttributes")){
		for (int i=0; i < bunch->getSize(); i++)
		{
			double x = part_coord_arr[i][0];
			double xp = part_coord_arr[i][1];
			double y = part_coord_arr[i][2];
			double yp = part_coord_arr[i][3];
			double Etot = syncPart->getEnergy() + syncPart->getMass();
			double dpp = 1/(beta*beta)*part_coord_arr[i][5]/Etot;
			
			double xval = (x - etax * dpp)/sqrt(betax);
			double xpval = (xp - etapx * dpp) * sqrt(betax) + xval * alphax;
			double yval = y / sqrt(betay);
			double ypval = (yp + y * alphay/betay) * sqrt(betay);
			
			double angle = atan2(xpval, xval);
			if(angle < 0.) angle += (2.0*OrbitConst::PI);
			double xPhase = angle;
			double xPhaseOld = bunch->getParticleAttributes("ParticlePhaseAttributes")->attValue(i,0);
			double xTune = (xPhaseOld - xPhase) / (2.0*OrbitConst::PI);
			if(xTune < 0.) xTune += 1.;
			bunch->getParticleAttributes("ParticlePhaseAttributes")->attValue(i, 0) = xPhase;
			bunch->getParticleAttributes("ParticlePhaseAttributes")->attValue(i, 2) = xTune;
			
			angle = atan2(ypval, yval);
			if(angle < 0.) angle += (2.0*OrbitConst::PI);
			double yPhase = angle;
			double yPhaseOld = bunch->getParticleAttributes("ParticlePhaseAttributes")->attValue(i,1);
			double yTune = (yPhaseOld - yPhase) / (2.0*OrbitConst::PI);
			if(yTune < 0.) yTune += 1.;
			bunch->getParticleAttributes("ParticlePhaseAttributes")->attValue(i, 1) = yPhase;
			bunch->getParticleAttributes("ParticlePhaseAttributes")->attValue(i, 3) = yTune;
			
			double xcanonical = x - etax * dpp;
			double ycanonical = y;
			double xpfac = xp - etax * dpp;
			double ypfac = yp;
			double pxcanonical =  xpfac + xcanonical * (alphax/betax);
			double pycanonical =  ypfac + ycanonical * (alphay/betay);
			double xAction = xcanonical  *  xcanonical / betax + pxcanonical * pxcanonical * betax;
			double yAction = ycanonical  *  ycanonical / betay + pycanonical * pycanonical * betay;
	
			bunch->getParticleAttributes("ParticlePhaseAttributes")->attValue(i, 4) = xAction;
			bunch->getParticleAttributes("ParticlePhaseAttributes")->attValue(i, 5) = yAction;
			}
	}

}







