#include "MaterialInteractions.hh"
#include "Collimator.hh"
#include "SyncPart.hh"
#include "cross_sections.hh"
#include "numrecipes.hh"
#include "OrbitConst.hh"
#include "Random.hh"

#include <iostream>
#include <cmath>
#include <cfloat>
#include <cstdlib>


// Constructor
///////////////////////////////////////////////////////////////////////////
//
// NAME
//
//   Collimator::Collimator
//
// DESCRIPTION
//   Constructs a collimator
//
// PARAMETERS
//   length: length in m
//   ma:	 material number. (1=carbon, 2=aluminum, 3=iron, 4=copper, 5=tantalum, 6=tungstun,
//			 7=platinum, 8=lead, 9 = black absorber)
//   densityfac: density factor (for materials mixed with air or water). 1.0 for pure. 
//   shape:  shape of the collimator: 1=circle, 2=ellipse, 3=one sided
//           flat, 4=two sided flat, 5=rectangular (outside is collimator),
//           6=rectangular (inside is collimator).
//   a:      depending on shape, either (shape = 1) radius, 
//           (shape = 2) semimajor axis, (shape = 3) distance to 
//           flat edge, (shape = 4) minimum edge, (shape=5 or 6) 
//           minimum horizontal edge.
//   b:      depending on shape, either (1) radius, (2) semimajor axis,
//           (3) zero  (4) maximum edge (5) (shape=5 or 6) maximum 
//           horizontal edge.
//   c:      minimum vertical edge (used only in shapes 5 or 6)
//   d:      maximum vertical edge (used only in shapes 5 or 6)
//   angle:  tilt angle of collimator.
///
// RETURNS
//   int.
//
///////////////////////////////////////////////////////////////////////////

Collimator::Collimator(double length, int ma, 
					   double density_fac, int shape, 
					   double a, double b, double c, double d, double angle): CppPyWrapper(NULL)
{
	length_ = length;
	ma_ = ma;
	shape_ = shape;
	density_fac_ = density_fac;
	a_ = a;
	b_ = b;
	c_ = c;
	d_ = d;
	angle_ = angle;
}

void Collimator::collimateBunch(Bunch* bunch, Bunch* lostbunch){

	int j = 1, coll_flag = 0, lastArg, trackit;
	double nAvogadro = 6.022045e23;
	double random, choice, length, dlength, meanfreepath, b_pN;
	double rl, zrl, stepsize, smallstep, radlengthfac, directionfac;
	double t, dp_x=0.0, dp_y=0.0, thetax = 0.0, thetay = 0.0, thx = 0.0, thy = 0.0;
	long idum = (unsigned)time(0);
	idum = -idum;
	
	SyncPart* syncPart = bunch->getSyncPart();	

	length = length_;
	dlength = length * 1.0e-4;
	
	double radlength = OrbitUtils::get_radlength(ma_);
	radlengthfac = radlength / density_fac_;
	smallstep = 0.001 * radlengthfac;
		
	double z = OrbitUtils::get_z(ma_);
	double a = OrbitUtils::get_a(ma_);
	double density = OrbitUtils::get_rho(ma_);
	if(a <= 62.) b_pN = 14.5 * pow(a, 0.6666667);
	if(a >  62.) b_pN = 60.0 * pow(a, 0.3333333);
	
	bunch->compress();
	double m_size = 0.;
	int nParts = bunch->getSize();
	double** part_coord_arr = bunch->coordArr();
	
	for(int ip = 0; ip < nParts; ip++){
		
		int step = 0;
		zrl = length;
		coll_flag = checkCollFlag(part_coord_arr[ip][0], part_coord_arr[ip][2]);
		
		while(zrl > 0){
			//If not inside collimator, drift until the end or entry to the collimator, whichever comes first.
			if(coll_flag == 0 && zrl > 0) coll_flag = driftParticle(coll_flag, zrl, length, part_coord_arr[ip], syncPart);
		
			//If in the collimator, tally the hit and start tracking
			if(coll_flag == 1) {
				nHits++;	
							
				//Check for black absorber.
				if(ma_ == 9) loseParticle(bunch, lostbunch, ip, nLost, coll_flag, zrl);
				
				directionfac = getDirection(part_coord_arr[ip], syncPart);
				rl = zrl * directionfac;
				
				//stepsize = rl;
				
				double beta = Collimator::getBeta(part_coord_arr[ip], syncPart);
				double p = Collimator::getP(part_coord_arr[ip], syncPart);
				double theta = 0.0136 / (beta * p) / sqrt(radlengthfac);
				double pfac = Collimator::getPFactor(part_coord_arr[ip], syncPart);
				double ecross = OrbitUtils::get_elastic_crosssection((syncPart->getEnergy() + part_coord_arr[ip][5]), ma_);
				double icross = OrbitUtils::get_inelastic_crosssection((syncPart->getEnergy() + part_coord_arr[ip][5]), ma_);
				
				if(step == 0){ //If first step, do an iteration with ecross and icross to get first stepsize and first rcross 
					step++;
					double totcross = icross + ecross;
					meanfreepath = (OrbitUtils::get_a(ma_) / (nAvogadro * 1000.0) / (density * density_fac_) / (totcross * 1.0e-28));
					stepsize = -meanfreepath * log(Random::ran1(idum));
				}
				
				double rcross = MaterialInteractions::ruthScattJackson(stepsize, z, a, density, idum, beta, 0, pfac, thetax, thetay);
				double totcross = ecross + icross + rcross;
				meanfreepath = OrbitUtils::get_a(ma_) / ((nAvogadro * 1000.0) * (density * density_fac_) * (totcross * 1.0e-28));
				stepsize = -meanfreepath * log(Random::ran1(idum));
				
				Collimator::checkStep(rl, radlengthfac, stepsize, part_coord_arr[ip], syncPart);
				if(stepsize < smallstep) stepsize = smallstep;
				
				if(stepsize > rl){ //Take the step but no nuclear scattering event
					stepsize = rl + dlength;
					Collimator::checkStep(rl, radlengthfac, stepsize, part_coord_arr[ip], syncPart);
					if(stepsize < smallstep) stepsize = smallstep;
					Collimator::takeStep(bunch, lostbunch, part_coord_arr[ip], syncPart, z, a, density, idum, stepsize, zrl, rl, coll_flag, ip);
					
				}
				
				else{ //Take the step and allow nuclear scatter
					Collimator::takeStep(bunch, lostbunch, part_coord_arr[ip], syncPart, z, a, density, idum, stepsize, zrl, rl, coll_flag, ip);
					
					//If it still exists after MCS and energy loss, nuclear scatter
					if(coll_flag==1 && zrl > 0){
						beta = Collimator::getBeta(part_coord_arr[ip], syncPart);
						p = Collimator::getP(part_coord_arr[ip], syncPart);
						theta = 0.0136 / (beta * p) / sqrt(radlengthfac);
						pfac = Collimator::getPFactor(part_coord_arr[ip], syncPart);
						
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
							loseParticle(bunch, lostbunch, ip, nLost, coll_flag, zrl);
						}
						
					}
				
				}
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
//   Collimator::checkCollFlag
//
// DESCRIPTION
//   Checks to see if a particle is located inside a collimator.  Returns
//   1 if the particle is in the collimator, 0 if it isn't.
//
// PARAMETERS
//   x:      x coordinate of particle.
//   y:      y coordinate of particle.
//
// RETURNS
//   int.
//
///////////////////////////////////////////////////////////////////////////

int Collimator::checkCollFlag(double x, double y){

	double xtemp = x, ytemp = y;
	double PI = OrbitConst::PI;
	double a = a_;
	double b = b_;
	double c = c_;
	double d = d_;
	double angle = angle_;
	int shape = shape_;
	double length = length_;
	
	
	if(shape==1)
    {
		if((pow(x, 2) + pow(y, 2)) >= pow(a, 2)) return 1;
    }
	
	if(shape==2)
    {
		if(angle != 0)
		{
			xtemp =  cos(angle*PI/180.)*x + sin(angle*PI/180.)*y;
			ytemp = -sin(angle*PI/180.)*x + cos(angle*PI/180.)*y;
		}
		if((pow(xtemp/a, 2.0)+pow(ytemp/b, 2.0) >= 1.)) return 1;
    }
	
	if(shape==3)
    {
		if(angle != 0)
		{
			xtemp =  cos(angle*PI/180.)*x + sin(angle*PI/180.)*y;
			ytemp = -sin(angle*PI/180.)*x + cos(angle*PI/180.)*y;
		}
		if(xtemp >= a) return 1;
    }
	
	if(shape==4)
    {
		if(angle != 0)
		{
			xtemp =  cos(angle*PI/180.)*x + sin(angle*PI/180.)*y;
			ytemp = -sin(angle*PI/180.)*x + cos(angle*PI/180.)*y;
		}
		if(xtemp <= a || xtemp >= b) return 1;
    }
	
	if(shape==5)
    {
		if(angle != 0)
		{
			xtemp =  cos(angle*PI/180.)*x + sin(angle*PI/180.)*y;
			ytemp = -sin(angle*PI/180.)*x + cos(angle*PI/180.)*y;
		}
		if(xtemp <=  a ||
		   xtemp >=  b ||
		   ytemp <=  c ||
		   ytemp >=  d ) return 1;
    }
	
	if(shape==6)
    {
		if(angle != 0)
		{
			xtemp =  cos(angle*PI/180.)*x + sin(angle*PI/180.)*y;
			ytemp = -sin(angle*PI/180.)*x + cos(angle*PI/180.)*y;
		}
		
		if(xtemp >=  a &&
		   xtemp <=  b &&
		   ytemp >=  c &&
		   ytemp <=  d ) return 1;
    }
    return 0; 
}

///////////////////////////////////////////////////////////////////////////
//
// NAME
//
//   Collimator::driftParticle
//
// DESCRIPTION
//   Drifts a particle until it either drifts beyond collimator length or enters
//	 collimator aperture. 
//
// PARAMETERS
//   coll_flag: flag for in or out of collimator
//   zrl:	 remaining length
//   length: length of collimator 
//	 coords: particle coordinates
//   syncpart: the relevant synchronous particle
//
// RETURNS
//   int.
//
///////////////////////////////////////////////////////////////////////////
	

int Collimator::driftParticle(int coll_flag, double& zrl, double length, double* coords, SyncPart* syncpart){
	while((coll_flag == 0) && (zrl > 0))
	{
		
		double x = coords[0];
		double xp = coords[1];
		double y = coords[2];
		double yp = coords[3];
		double dlength = length * 1.0e-4;
		double stepsize = 0.001;
		if(stepsize > length / 10.) stepsize = length / 10.;
		if(stepsize > zrl) stepsize = zrl + dlength;
				
		double pfac = Collimator::getPFactor(coords, syncpart);
	
		x += stepsize * xp / pfac;
		y += stepsize * yp / pfac;
		zrl -= stepsize;
		coll_flag = Collimator::checkCollFlag(x, y);
		
		nHits++;
		if(coll_flag == 1) return 1;

	}

	return 0;
}


///////////////////////////////////////////////////////////////////////////
//
// NAME
//
//   Collimator::getDirection
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

double Collimator::getDirection(double* coords, SyncPart* syncpart){
	
	double pfac = Collimator::getPFactor(coords, syncpart);
	
	double xpfac = coords[1] / pfac;
	double ypfac = coords[3] / pfac;
	double directionfac = sqrt(1.0 + xpfac * xpfac + ypfac * ypfac);
	
	return directionfac;
	
}

///////////////////////////////////////////////////////////////////////////
//
// NAME
//
//   Collimator::checkStep
//
// DESCRIPTION
//   Checks to make sure that the randomly generated step size isn't too
//   large (See Catalan-Lasheras, PhD thesis, p 52).  The model is based
//   on the K2 collimation code written by Jean Bernard Jeanneret.
//
// PARAMETERS
//   
//   s:   Length of collimator remaining
//   stepsize: Randomly generated step   
//   angle:  tilt angle of collimator.
//   coords: particle coordinates
//
// RETURNS
//   Nothing
//
///////////////////////////////////////////////////////////////////////////


void Collimator::checkStep(double s, double radlengthfac, double& stepsize, double* coords, SyncPart* syncpart){

	double PI = OrbitConst::PI;
	double beta = Collimator::getBeta(coords, syncpart);
	double p = Collimator::getP(coords, syncpart);
	double theta = 0.0136 / (beta * p) / sqrt(radlengthfac);
	double a = a_;
	double b = b_;
	double c = c_;
	double d = d_;
	double angle = angle_;
	int shape = shape_;
	double length = length_;
	
	double x = coords[0];
	double px = coords[1];
	double y = coords[2];
	double py = coords[3];
	
	float x1=0;
	float x2=length - s;
	int n=50;
	int nb=0;
	float xb1[4];
	float xb2[4];
	float solution[4];
	float smax=0.;
	float xacc=.001;
	int i;
	float r_o, pr_o;
	double arg, theta_o,  x_ellipse, y_ellipse;
	double xtemp=x, ytemp=y;
	
	
	if(shape==1)
    {
		r_o = sqrt(x*x + y*y) - a;
		if(a==0.) r_o=10000.;
    }
	
	
	if(shape==2)
    {
		if(angle != 0)
		{
			xtemp = cos(angle*PI/180.)*x + sin(angle*PI/180.)*y;
			ytemp = -sin(angle*PI/180.)*x + cos(angle*PI/180.)*y;
		}
		arg = ytemp/xtemp;
		if(xtemp>=0) theta_o = atan(arg);
		if(xtemp<0) theta_o = PI+atan(arg);
		y_ellipse = sqrt(pow( a*tan(theta_o) , 2.0)/(1+pow( a*tan(theta_o)/b , 2.0)));
		x_ellipse = sqrt(a*a - pow( a*y_ellipse/b, 2.0));
		r_o = sqrt(xtemp*xtemp + ytemp*ytemp) - 
		sqrt(x_ellipse*x_ellipse + y_ellipse*y_ellipse);
    }	       
	
	
	if(shape==3)
    {
		if(angle != 0)
		{
			xtemp = cos(angle*PI/180.)*x + sin(angle*PI/180.)*y;
			ytemp = -sin(angle*PI/180.)*x + cos(angle*PI/180.)*y;
		}
		r_o = fabs(xtemp - a);
    }
	
	
	if(shape==4)
    {
		if(angle != 0)
		{
			xtemp = cos(angle*PI/180.)*x + sin(angle*PI/180.)*y;
			ytemp = -sin(angle*PI/180.)*x + cos(angle*PI/180.)*y;
		}
		if(xtemp<=a) r_o = fabs(xtemp - a);
		if(xtemp>=b) r_o = fabs(xtemp - b);
    }
	
	
	if(shape==5)
    {
		if(angle != 0)
		{
			xtemp = cos(angle*PI/180.)*x + sin(angle*PI/180.)*y;
			ytemp = -sin(angle*PI/180.)*x + cos(angle*PI/180.)*y;
		}
		
		if(xtemp<=a || xtemp>=b){
			if(xtemp<=a){
				if(ytemp<=c) 
					r_o=sqrt(pow((xtemp - a), 2.0) + pow((ytemp - c), 2.0));
				if(ytemp>=d)
					r_o=sqrt(pow((xtemp - a), 2.0) + pow((ytemp - d), 2.0));
				if(ytemp>=c && ytemp<=d)
					r_o = fabs(xtemp - a);
			}
			if(xtemp>=b){
				if(ytemp<=c)
					r_o=sqrt(pow((xtemp - b), 2.0) + pow((ytemp - c), 2.0));
				if(ytemp>=d)
					r_o=sqrt(pow((xtemp - b), 2.0) + pow((ytemp - d), 2.0));
				if(ytemp>=c && ytemp<=d)
					r_o = fabs(xtemp - b);
			}
		}
		else{
			if(ytemp<=c) r_o = fabs(ytemp - c);
			if(ytemp>=d) r_o = fabs(ytemp - d);
		}
    }
	
	
	if(shape==6)
    {
		if(angle != 0)
		{
			xtemp = cos(angle*PI/180.)*x + sin(angle*PI/180.)*y;
			ytemp = -sin(angle*PI/180.)*x + cos(angle*PI/180.)*y;
		}
		r_o=fabs(xtemp - a);
		if(r_o > fabs(xtemp - b)) r_o = fabs(xtemp - b);
		if(r_o > fabs(ytemp - c)) r_o = fabs(ytemp - c);
		if(r_o > fabs(ytemp - d)) r_o = fabs(ytemp - d);
    }
	
	pr_o = (x*px + y*py)/sqrt(x*x+y*y);

	for(i=0; i<4; i++)
    {
		xb1[i]=0;
		xb2[i]=0;
    }
	
	OrbitUtils::zbrak(OrbitUtils::fstep, x1, x2, n, xb1, xb2, nb, r_o, pr_o, theta);
	
	for(i=1; i<=nb; i++)
    {
		solution[i]=OrbitUtils::rtbis(OrbitUtils::fstep, xb1[i], xb2[i], xacc, r_o, pr_o, theta);
    }
	
	solution[0]=100000.;
	
	for(i=1; i<=nb; i++)
		if(solution[i] < solution[i-1])
		{
			smax=solution[i];
		}
	
	if(smax > 0 && stepsize > smax){
		
	 stepsize=smax;
	 
	 }

}	

///////////////////////////////////////////////////////////////////////////
//
// NAME
//
//   Collimator::getPFac
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

double Collimator::getPFactor(double* coords, SyncPart* syncpart){
	
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
//   Collimator::getBeta
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

double Collimator::getBeta(double* coords, SyncPart* syncpart){
	
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
//   Collimator::getP
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

double Collimator::getP(double* coords, SyncPart* syncpart){

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
//   Collimator::takeStep
//
// DESCRIPTION
//   Takes a step in the collimator, allows MCS and ionization energy loss.
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
//	 zrl:		remaining collimator length in the z direction
//   rl:		remaining collimator length in the direction of particle momentum
//   coll_flag:	flag for in or out of the collimator.
//
// RETURNS
//   Nothing
//
///////////////////////////////////////////////////////////////////////////
	
void Collimator::takeStep(Bunch* bunch, Bunch* lostbunch, double* coords, SyncPart* syncpart, double z, double a, double density, long& idum, double stepsize, double& zrl, double& rl, int& coll_flag, int ip){

	double beta = Collimator::getBeta(coords, syncpart);
	double p = Collimator::getP(coords, syncpart);
	double pfac = Collimator::getPFactor(coords, syncpart);
	
	MaterialInteractions::mcsJackson(stepsize, z, a, density, idum, beta, pfac, coords[0], coords[2], coords[1], coords[3]);
	double dE = MaterialInteractions::ionEnergyLoss(beta, z, a);
	dE = -dE * density * density_fac_ * stepsize; //Factors for units m->cm and MeV->GeV
	coords[5] += dE;

	if((coords[5] + syncpart->getEnergy()) < 0.02){ 
		Collimator::loseParticle(bunch, lostbunch, ip, nLost, coll_flag, zrl);
	}
	else {
		double directionfac = Collimator::getDirection(coords, syncpart);
		zrl -= stepsize / directionfac;
		rl = zrl * directionfac;
		coll_flag = checkCollFlag(coords[0], coords[2]);
	}
}
	

///////////////////////////////////////////////////////////////////////////
//
// NAME
//
//   Collimator::loseParticle
//
// DESCRIPTION
//   Lose the particle from the alive bunch and add it to the lost bunch
//
// PARAMETERS
//	 coords: particle coordinates
//   bunch: the primary bunch
//	 lostbunch: the lost bunch
//	 nLost:		the number of lost particles in this collimator
//	 coll_flag:	flag for if particle is still in collimator (and alive)
//	 zrt:		remaining length in z the particle has in collimator
//
// RETURNS
//   nothing.
//
///////////////////////////////////////////////////////////////////////////
	
void Collimator::loseParticle(Bunch* bunch, Bunch* lostbunch, int ip, int& nLost, int& coll_flag, double& zrl){
	double** coords = bunch->coordArr();
	lostbunch->addParticle(coords[ip][0], coords[ip][1], coords[ip][2], coords[ip][3], coords[ip][4], coords[ip][5]);
	if (lostbunch->hasParticleAttributes("LostParticleAttributes") > 0) {
		lostbunch->getParticleAttributes("LostParticleAttributes")->attValue(lostbunch->getSize() - 1, 0) = length_ - zrl;
	}
	bunch->deleteParticleFast(ip);
	nLost++;
	coll_flag = 0;
	zrl = -1.;
	
}
	
	

