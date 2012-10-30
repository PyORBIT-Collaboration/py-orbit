#include "MaterialInteractions.hh"
#include "Foil.hh"
#include "SyncPart.hh"
#include "cross_sections.hh"
#include "numrecipes.hh"
#include "OrbitConst.hh"
#include "Random.hh"

#include <iostream>
#include <cmath>
#include <cfloat>
#include <cstdlib>


///////////////////////////////////////////////////////////////////////////
//
// NAME
//
//   Foil::Foil
//
// DESCRIPTION
//   Constructs a foil
//
// PARAMETERS
//   xmin: min horizontal foil position
//   xmax: max horizontal foil position
//   ymin: min vertical foil position
//   ymax: max vertical foil position
//
// RETURNS
//   int.
//
///////////////////////////////////////////////////////////////////////////

Foil::Foil(double xmin, double xmax, double ymin, double ymax, double thick): CppPyWrapper(NULL)
{
	xmin_ = xmin;
	xmax_ = xmax;
	ymin_ = ymin;
	ymax_ = ymax;
	thick_ = thick;
	length_ = 0.0;
	ma_ = 0;
}

///////////////////////////////////////////////////////////////////////////
//
// NAME
//
//   Foil::traverseFoilSimpleScatter
//
// DESCRIPTION
//   Material scatters particles through a foil with simplified MCS scatter.
//
// PARAMETERS
//	Bunch - The particle bunch
//  LostBunch 
//
// RETURNS
//   int.
//
///////////////////////////////////////////////////////////////////////////

void Foil::traverseFoilSimpleScatter(Bunch* bunch){

	double BohrRadius=0.52917706e-8;  // hydrogenic Bohr radius in cm
	double hBar = 1.0545887e-27;      // Planck's constant in erg-sec
	double echarge = 4.803242e-10;    // in esu or statcoulombs
	double nAvogadro = 6.022045e23;
	double deg2Rad = 1.74532925e-2;
	double rhofoil = 2.265;
	double muScatter = 1.35;
	double pInj0;
	long idum = (unsigned)time(0);
	idum = -idum;
		
	int foil_flag = 0;
	double length = thick_ / (1.0e3 * OrbitUtils::get_rho(ma_));
	double zrl = 0.0;
	double random1 = 0;
	
    // Momentum in g*cm/sec
	SyncPart* syncPart = bunch->getSyncPart();
    pInj0 = 1.6726e-22 * syncPart->getMass()/OrbitConst::mass_proton * syncPart->getBeta() *
	syncPart->getGamma() * OrbitConst::c;
	
    // Thomas-Fermi atom radius (cm):
	
    double TFRadius = muScatter *  BohrRadius *pow(OrbitUtils::get_z(ma_), -0.33333);
	
    // Minimum scattering angle:
	
    double thetaScatMin =  hBar / (pInj0 * TFRadius);
	
    // Theta max as per Jackson (13.102)
	
    double thetaScatMax = 274.e5 * OrbitConst::mass_electron * OrbitConst::c /
	(pInj0 * pow(OrbitUtils::get_a(ma_), 0.33333));
	
    double pv = 1.e2 * pInj0 * syncPart->getBeta() * OrbitConst::c;
    double term = OrbitUtils::get_z(ma_) * echarge * echarge / pv;
    double sigmacoul = 4*OrbitConst::PI * term * term / (thetaScatMin * thetaScatMin);
	
    // Scattering area per area
	
    double nscatters = nAvogadro * (OrbitUtils::get_rho(ma_)/1000.0) /
	OrbitUtils::get_a(ma_) * length * sigmacoul;
	
    // Mean free path
	
    double lscatter = length/nscatters;
	bunch->compress();
	int nParts = bunch->getSize();
	double** part_coord_arr = bunch->coordArr();
	
	for(int ip = 0; ip < nParts; ip++){
		
		foil_flag = checkFoilFlag(part_coord_arr[ip][0], part_coord_arr[ip][2]);
		
		//If in the foil, tally the hit and start tracking
		if(foil_flag == 1) {
			nHits++;	
			zrl = length;  // distance remaining in foil in cm//
			double thetaX = 0.;
			double thetaY = 0.;
			
			// Generate interaction points until particle exits foil
			
			while (zrl >= 0.0)
			{
				random1 = Random::ran1(idum);
				//cout << "idum "<<idum<<"\n";
				zrl += lscatter * log(random1);
				if(zrl < 0.0) break; // exit foil
				
				// Generate random angles
				
				random1 = Random::ran1(idum);
				double phi = 2*OrbitConst::PI * random1;
				random1 = Random::ran1(idum);
				double theta = thetaScatMin * sqrt(random1 / (1. - random1));
				thetaX += theta * cos(phi);
				thetaY += theta * sin(phi);
				//cout << thetaX <<"\n";
			}
			part_coord_arr[ip][1] += thetaX;
			part_coord_arr[ip][3] += thetaY;
		}
	}

}

	
///////////////////////////////////////////////////////////////////////////
//
// NAME
//
//   Foil::traverseFoilFullScatter
//
// DESCRIPTION
//   Material scatters particles through a foil. MCS and nuclear scattering.
//
// PARAMETERS
//	Bunch - The particle bunch
//  LostBunch 
//
// RETURNS
//   int.
//
///////////////////////////////////////////////////////////////////////////

void Foil::traverseFoilFullScatter(Bunch* bunch, Bunch* lostbunch){

	int j = 1, foil_flag = 0, lastArg, trackit;
	double nAvogadro = 6.022045e23;
	double random, choice, length, dlength, meanfreepath;
	double rl, zrl, stepsize, radlengthfac, directionfac;
	double t, dp_x=0.0, dp_y=0.0, thetax = 0.0, thetay = 0.0, thx = 0.0, thy = 0.0;
	long idum = (unsigned)time(0);
	idum = -idum;
	
	SyncPart* syncPart = bunch->getSyncPart();	
	
	double z = OrbitUtils::get_z(ma_);
	double a = OrbitUtils::get_a(ma_);
	double density = OrbitUtils::get_rho(ma_);
	double b_pN = 14.5 * pow(a, 0.6666667);
	double radlength = OrbitUtils::get_radlength(ma_);

	length = thick_ / (1.0e5 * density);
	length_ = length;
	dlength = length * 1.0e-4;
	
	bunch->compress();
	double m_size = 0.;
	int nParts = bunch->getSize();
	double** part_coord_arr = bunch->coordArr();
	
	for(int ip = 0; ip < nParts; ip++){
	
		int step = 0;
		zrl = length;
		foil_flag = checkFoilFlag(part_coord_arr[ip][0], part_coord_arr[ip][2]);
		
		while(zrl > 0){
			//If in the foil, tally the hit and start tracking
			if(foil_flag == 1) {
				nHits++;	
				directionfac = getDirection(part_coord_arr[ip], syncPart);
				rl = zrl * directionfac;
								
				double beta = Foil::getBeta(part_coord_arr[ip], syncPart);
				double p = Foil::getP(part_coord_arr[ip], syncPart);
				double theta = 0.0136 / (beta * p) / sqrt(radlengthfac);
				double pfac = Foil::getPFactor(part_coord_arr[ip], syncPart);
				double ecross = OrbitUtils::get_elastic_crosssection((syncPart->getEnergy() + part_coord_arr[ip][5]), ma_);
				double icross = OrbitUtils::get_inelastic_crosssection((syncPart->getEnergy() + part_coord_arr[ip][5]), ma_);
				
				if(step == 0){ //If first step, do an iteration with ecross and icross to get first stepsize and first rcross 
					step++;
					double totcross = icross + ecross;
					meanfreepath = (OrbitUtils::get_a(ma_) / (nAvogadro * 1e3) / density / (totcross * 1.0e-28));
					stepsize = -meanfreepath * log(Random::ran1(idum));
				}
				
				double rcross = MaterialInteractions::ruthScattJackson(stepsize, z, a, density, idum, beta, 0, pfac, thetax, thetay);
				double totcross = ecross + icross + rcross;
				meanfreepath = OrbitUtils::get_a(ma_) / ((nAvogadro * 1e3) * density  * (totcross * 1.0e-28));
				stepsize = -meanfreepath * log(Random::ran1(idum));
			
				if(stepsize > rl){ //Take the step but no nuclear scattering event
					stepsize = rl + dlength;
					Foil::takeStep(bunch, lostbunch, part_coord_arr[ip], syncPart, z, a, density, 
								   idum, stepsize, zrl, rl, foil_flag, ip);
					
				}
				if(stepsize <= rl) { //Take the step and allow nuclear scatter
					Foil::takeStep(bunch, lostbunch, part_coord_arr[ip], syncPart, z, a, density, idum, stepsize, zrl, rl, foil_flag, ip);
				
					//If it still exists after MCS and energy loss, nuclear scatter
					if(foil_flag==1 && zrl > 0){
						beta = Foil::getBeta(part_coord_arr[ip], syncPart);
						p = Foil::getP(part_coord_arr[ip], syncPart);
						theta = 0.0136 / (beta * p) / sqrt(radlengthfac);
						pfac = Foil::getPFactor(part_coord_arr[ip], syncPart);
						
						ecross = OrbitUtils::get_elastic_crosssection((syncPart->getEnergy() + part_coord_arr[ip][5]), ma_);
						icross = OrbitUtils::get_inelastic_crosssection((syncPart->getEnergy() + part_coord_arr[ip][5]), ma_);
						rcross = MaterialInteractions::ruthScattJackson(stepsize, z, a, density, idum, beta, 0, pfac, thx, thy);
						
						totcross = ecross + icross + rcross;
						
						double e_frac = ecross/totcross;
						double i_frac = icross/totcross;
						double r_frac = rcross/totcross;
						
						choice = Random::ran1(idum);
						
						// Nuclear Elastic Scattering
						if((choice >= 0.) && (choice <= e_frac))
						{
							if((syncPart->getEnergy() + part_coord_arr[ip][5]) <= 0.4)
							{
								t=MaterialInteractions::elastic_t(p, a, idum);
							}
							if((syncPart->getEnergy() + part_coord_arr[ip][5]) > 0.4)
							{
								t=-log(Random::ran1(idum))/b_pN;
							}

							MaterialInteractions::momentumKick(t, p, dp_x, dp_y);
							part_coord_arr[ip][1] += dp_x * pfac;
							part_coord_arr[ip][3] += dp_y * pfac;
						}
						
						// Rutherford Coulomb scattering
						if((choice > e_frac) && (choice <= (1 - i_frac)))
						{
							rcross = MaterialInteractions::ruthScattJackson(stepsize, z, a, density, idum, beta, 1, pfac, thx, thy);
							
							double xpfac = part_coord_arr[ip][1] / pfac;
							double ypfac = part_coord_arr[ip][3] / pfac;
							
							double anglex = atan(xpfac) + thx;
							double angley = atan(ypfac) + thy;
							
							part_coord_arr[ip][1] = tan(anglex) * pfac;
							part_coord_arr[ip][1] = tan(angley) * pfac;
						}
						
						// Nuclear Inelastic absorption
						if( (choice > (1.-i_frac)) && (choice <= 1.))
						{
							loseParticle(bunch, lostbunch, ip, nLost, foil_flag, zrl);
						}
						
					}
				}
			}
			else{
				//Drift by the foil
				driftParticle(part_coord_arr[ip], syncPart, length);
				zrl = 0;
			}
		}
	}
	
	//Update synchronous particle, compress bunch
	bunch->compress();
	double newtime = syncPart->getTime() + length/( syncPart->getBeta()*OrbitConst::c );
	syncPart->setTime(newtime);
}

///////////////////////////////////////////////////////////////////////////
//
// NAME
//
//   Foil::checkFoilFlag
//
// DESCRIPTION
//   Checks to see if a particle is located inside the foil region.  
//	 Returns 1 if the particle is in the Foil, 0 if it isn't.
//
// PARAMETERS
//   x:      x coordinate of particle.
//   y:      y coordinate of particle.
//
// RETURNS
//   int.
//
///////////////////////////////////////////////////////////////////////////

int Foil::checkFoilFlag(double x, double y){

	double xtemp = x, ytemp = y;
	if((x >= xmin_) && (x <= xmax_) && (y >= ymin_) && (y <= ymax_)){
		return 1;
	}
	else {
		return 0;
	}

}

//////////////////////////////////////////////////////////////////////////
// NAME
//   Drift Particle
//
// DESCRIPTION
//   Drifts a single particle. 
//
// PARAMETERS
//   arr = reference to the particle coordinate array
//   i = particle index
//   length = length of the drift
//
// RETURNS
//   Nothing
//
///////////////////////////////////////////////////////////////////////////

void Foil::driftParticle(double* arr, SyncPart* syncpart, double length)
{
    double KNL, phifac, dp_p;
    double gamma2i = 1.0 / (syncpart->getGamma() * syncpart->getGamma());
    double dp_p_coeff = 1.0 / (syncpart->getMomentum() * syncpart->getBeta());
	
    dp_p = arr[5] * dp_p_coeff;
    KNL  = 1.0 / (1.0 + dp_p);
    arr[0] += KNL * length * arr[1];
    arr[2] += KNL * length * arr[3];
    phifac = (arr[1] * arr[1] + arr[3] * arr[3] + dp_p * dp_p * gamma2i) / 2.0;
    phifac = (phifac * KNL - dp_p * gamma2i) * KNL;
    arr[4] -= length * phifac;
}


///////////////////////////////////////////////////////////////////////////
//
// NAME
//
//  Foil::getDirection
//
// DESCRIPTION
//   Gets the normalized unit vector direction of particle momentum
//
// PARAMETERS
//
//	 coords: particle coordinates
//   syncpart: the relevant synchronous particle
//
// RETURNS
//   int.
//
///////////////////////////////////////////////////////////////////////////

double Foil::getDirection(double* coords, SyncPart* syncpart){
	
	double pfac = Foil::getPFactor(coords, syncpart);
	
	double xpfac = coords[1] / pfac;
	double ypfac = coords[3] / pfac;
	double directionfac = sqrt(1.0 + xpfac * xpfac + ypfac * ypfac);
	
	return directionfac;
	
}

///////////////////////////////////////////////////////////////////////////
//
// NAME
//
//   Foil::getPFac
//
// DESCRIPTION
//   Gets the fractional deviation from sync part momentum (1 +dp/p)
//
// PARAMETERS
//
//	 coords: particle coordinates
//   syncpart: the relevant synchronous particle
//
// RETURNS
//   double.
//
///////////////////////////////////////////////////////////////////////////

double Foil::getPFactor(double* coords, SyncPart* syncpart){
	
	double M = syncpart->getMass();
	double T = syncpart->getEnergy();
	double P_sync = syncpart->energyToMomentum(T);
	double T_part = T + coords[5];
	double E_part = T_part + M;
	double P_part = sqrt(E_part*E_part - M*M);
	double dp = (P_part - P_sync)/P_sync;
	double pfac = 1.0 + dp;
	
	return pfac;
	
}

///////////////////////////////////////////////////////////////////////////
//
// NAME
//
//   Foil::getBeta
//
// DESCRIPTION
//   Gets current beta for particle
//
// PARAMETERS
//
//	 coords: particle coordinates
//   syncpart: the relevant synchronous particle
//
// RETURNS
//   double.
//
///////////////////////////////////////////////////////////////////////////

double Foil::getBeta(double* coords, SyncPart* syncpart){
	
	double T = syncpart->getEnergy();
	double M = syncpart->getMass();
	double T_part = T + coords[5];
	double E_part = T_part + M;
	double beta = sqrt((E_part*E_part - M * M))/E_part;

	return beta;
	
}

///////////////////////////////////////////////////////////////////////////
//
// NAME
//
//   Foil::getP
//
// DESCRIPTION
//   Gets current momentum of particle
//
// PARAMETERS
//
//	 coords: particle coordinates
//   syncpart: the relevant synchronous particle
//
// RETURNS
//   double.
//
///////////////////////////////////////////////////////////////////////////

double Foil::getP(double* coords, SyncPart* syncpart){

	double T = syncpart->getEnergy();
	double M = syncpart->getMass();
	double T_part = T + coords[5];
	double E_part = T_part + M;
	double p_part = sqrt((E_part*E_part - M * M));
	
	return p_part;
	
}

///////////////////////////////////////////////////////////////////////////
//
// NAME
//
//   Foil::takeStep
//
// DESCRIPTION
//   Takes a step in the Foil, allows MCS and ionization energy loss.
//	 Loses the particle if the energy is too low and updates the tracking params.
//
// PARAMETERS
//
//	 coords:	particle coordinates
//   syncpart:	the relevant synchronous particle
//	 z:			z number of material
//	 a:			a number of material
//	 density:	density of material
//	 idum:		a random number seed
//	 stepsize:	the stepsize to be taken
//	 zrl:		remaining Foil length in the z direction
//   rl:		remaining Foil length in the direction of particle momentum
//   foil_flag:	flag for in or out of the Foil.
//
// RETURNS
//   Nothing
//
///////////////////////////////////////////////////////////////////////////
	
void Foil::takeStep(Bunch* bunch, Bunch* lostbunch, double* coords, SyncPart* syncpart, double z, double a, double density, long& idum, double stepsize, double& zrl, double& rl, int& foil_flag, int ip){

	double beta = Foil::getBeta(coords, syncpart);
	double p = Foil::getP(coords, syncpart);
	double pfac = Foil::getPFactor(coords, syncpart);
	
	MaterialInteractions::mcsJackson(stepsize, z, a, density, idum, beta, pfac, coords[0], coords[2], coords[1], coords[3]);
	double dE = MaterialInteractions::ionEnergyLoss(beta, z, a);
	dE = -dE * density * stepsize; //Factors for units m->cm and MeV->GeV
	coords[5] += dE;
	
	if((coords[5] + syncpart->getEnergy()) < 0.02){ 
		Foil::loseParticle(bunch, lostbunch, ip, nLost, foil_flag, zrl);
	}
	else {
		double directionfac = Foil::getDirection(coords, syncpart);
		zrl -= stepsize / directionfac;
		rl = zrl * directionfac;
		foil_flag = checkFoilFlag(coords[0], coords[2]);
	}
}
	

///////////////////////////////////////////////////////////////////////////
//
// NAME
//
//   Foil::loseParticle
//
// DESCRIPTION
//   Lose the particle from the alive bunch and add it to the lost bunch
//
// PARAMETERS
//	 coords: particle coordinates
//   bunch: the primary bunch
//	 lostbunch: the lost bunch
//	 nLost:		the number of lost particles in this Foil
//	 foil_flag:	flag for if particle is still in Foil (and alive)
//	 zrt:		remaining length in z the particle has in Foil
//
// RETURNS
//   nothing.
//
///////////////////////////////////////////////////////////////////////////
	
void Foil::loseParticle(Bunch* bunch, Bunch* lostbunch, int ip, int& nLost, int& foil_flag, double& zrl){

	double** coords = bunch->coordArr();
	lostbunch->addParticle(coords[ip][0], coords[ip][1], coords[ip][2], coords[ip][3], coords[ip][4], coords[ip][5]);
	bunch->deleteParticleFast(ip);
	nLost++;
	foil_flag = 0;
	zrl = -1.;
	
}
	
	

