#include "BunchTwissAnalysis.hh"

#include <iostream>
#include <cmath>
#include <cfloat>
#include <cstdlib>

#include "ParticleMacroSize.hh"
#include "SyncPart.hh"

/** Constructor */
BunchTwissAnalysis::BunchTwissAnalysis(): CppPyWrapper(NULL)
{
	
	avg_arr = (double* ) malloc (6*sizeof(double));
	corr_arr = (double* ) malloc (36*sizeof(double));
	
	avg_arr_MPI = (double* ) malloc (6*sizeof(double));
	corr_arr_MPI = (double* ) malloc (36*sizeof(double));
	
	for(int i = 0; i < 6; i++){
		avg_arr[i] = 0.;
		avg_arr_MPI[i] = 0.;		
	}
	
	for(int i = 0; i < 6; i++){
		for(int j = 0; j < 6; j++){
			corr_arr[i+6*j] = 0.;
			corr_arr_MPI[i+6*j] = 0.;
		}	
	}
		
	count = 0;
	_order = 0;
	
	total_macrosize = 0.;
}

/** Destructor */
BunchTwissAnalysis::~BunchTwissAnalysis()
{
	free(avg_arr);
	free(corr_arr);
	free(avg_arr_MPI);
	free(corr_arr_MPI);
}

/** Performs the Twiss analysis of the bunch */		
void BunchTwissAnalysis::analyzeBunch(Bunch* bunch){
	
	//initialization
	for(int i = 0; i < 6; i++){
		avg_arr[i] = 0.;
	}
	
	for(int i = 0; i < 6; i++){
		for(int j = 0; j < 6; j++){
			corr_arr[i+6*j] = 0.;
		}	
	}
	count = 0;
	total_macrosize = 0.;
	
	bunch->compress();
	double m_size = 0.;
	int nParts = bunch->getSize();
	count += nParts;
	double** part_coord_arr = bunch->coordArr();
	int has_msize = bunch->hasParticleAttributes("macrosize");
	if(has_msize > 0){
		ParticleMacroSize* macroSizeAttr = (ParticleMacroSize*) bunch->getParticleAttributes("macrosize");
		double m_size = 0.;
		for(int ip = 0; ip < nParts; ip++){
			m_size = macroSizeAttr->macrosize(ip);
			total_macrosize += m_size;
			for(int i = 0; i < 6; i++){
				avg_arr[i] += m_size*part_coord_arr[ip][i];
			}
			
			for(int i = 0; i < 6; i++){
				for(int j = 0; j < i+1; j++){
					corr_arr[i+6*j] += m_size*part_coord_arr[ip][i]*part_coord_arr[ip][j];
				}	
			}			
		}	
	} else {
		m_size = 1.0;
		for(int ip = 0; ip < nParts; ip++){
			for(int i = 0; i < 6; i++){
				avg_arr[i] += part_coord_arr[ip][i];
			}
			
			for(int i = 0; i < 6; i++){
				for(int j = 0; j < i+1; j++){
					corr_arr[i+6*j] += part_coord_arr[ip][i]*part_coord_arr[ip][j];
				}	
			}	
		}
		total_macrosize += nParts*m_size;	
		for(int i = 0; i < 6; i++){
			avg_arr[i] *= m_size;
		}	
		for(int i = 0; i < 6; i++){
			for(int j = 0; j < i+1; j++){
				corr_arr[i+6*j] *= m_size;
			}	
		}				
	}	
	
	int count_MPI = 0;
	ORBIT_MPI_Allreduce(&count,&count_MPI,1,MPI_INT,MPI_SUM,bunch->getMPI_Comm_Local()->comm);
	count = count_MPI;
	
	double total_macrosize_MPI = 0.;
	ORBIT_MPI_Allreduce(&total_macrosize,&total_macrosize_MPI,1,MPI_DOUBLE,MPI_SUM,bunch->getMPI_Comm_Local()->comm);
	total_macrosize = total_macrosize_MPI;
	
	ORBIT_MPI_Allreduce(avg_arr,avg_arr_MPI,6,MPI_DOUBLE,MPI_SUM,bunch->getMPI_Comm_Local()->comm);
	ORBIT_MPI_Allreduce(corr_arr,corr_arr_MPI,36,MPI_DOUBLE,MPI_SUM,bunch->getMPI_Comm_Local()->comm);
	
	for(int i = 0; i < 6; i++){
		avg_arr[i] = avg_arr_MPI[i]/total_macrosize;
	}
	
	for(int i = 0; i < 6; i++){
		for(int j = 0; j < i+1; j++){
			corr_arr[i+6*j] = corr_arr_MPI[i+6*j]/total_macrosize;
			corr_arr[j+6*i] = corr_arr_MPI[i+6*j]/total_macrosize;
		}	
	}	
	
	SyncPart* syncPart = bunch->getSyncPart();	
	
	bunch_momentum = syncPart->getMomentum();
	bunch_beta = syncPart->getBeta();
	bunch_gamma = syncPart->getGamma();
	bunch_kinenergy = syncPart->getEnergy();
	bunch_mass = syncPart->getMass();
	
}


/** Performs the bunch moments computations */
void BunchTwissAnalysis::computeBunchMoments(Bunch* bunch, int order, int dispersionflag){
	_order = order;
	int i = 0;
	int j = 0;

	bunch->compress();
	
	double dispterm = 0.;
	double m_size = 0.;
	double xAvg = 0.;
	double yAvg = 0.;
	double total_macrosize = 0; //Total macrosize (can different than number of macroparticles if m_size is specified)
	int nParts = bunch->getSize();
	double total_macrosize_MPI = 0.;
	double** part_coord_arr = bunch->coordArr();
	int has_msize = bunch->hasParticleAttributes("macrosize");

	analyzeBunch(bunch);

	if(dispterm > 0){
		if(has_msize > 0){
			ParticleMacroSize* macroSizeAttr = (ParticleMacroSize*) bunch->getParticleAttributes("macrosize");
			double m_size = 0.;
			for(int ip = 0; ip < nParts; ip++){
				m_size = macroSizeAttr->macrosize(ip);
				total_macrosize += m_size;
				if (dispersionflag > 0) {
					dispterm = getDispersion(0) * part_coord_arr[ip][5] / (bunch_kinenergy + bunch_mass) / (bunch_beta*bunch_beta);
				}
				xAvg += m_size*(part_coord_arr[ip][0] - dispterm);
			}
		} else {
			m_size = 1.0;
			for(int ip = 0; ip < nParts; ip++){
				if (dispersionflag > 0) {
					dispterm = getDispersion(0) * part_coord_arr[ip][5] / (bunch_kinenergy + bunch_mass) / (bunch_beta*bunch_beta);
				}
				xAvg += part_coord_arr[ip][0] - dispterm;
			}
			total_macrosize += nParts*m_size;
			xAvg *= m_size;
		}
		
	
		ORBIT_MPI_Allreduce(&total_macrosize,&total_macrosize_MPI,1,MPI_DOUBLE,MPI_SUM,bunch->getMPI_Comm_Local()->comm);
		total_macrosize = total_macrosize_MPI;
	
		double xAvg_MPI = 0;
		ORBIT_MPI_Allreduce(&xAvg,&xAvg_MPI,1,MPI_DOUBLE,MPI_SUM,bunch->getMPI_Comm_Local()->comm);
		xAvg = xAvg_MPI/total_macrosize;
	}
	else{
		xAvg = getAverage(0);
	}
	yAvg = getAverage(2);
	
	double momX [order+1];
	double momY [order+1];
	
	momentXY = new double*[order+1];

	for(i=0; i < order+1; i++)
		momentXY[i] = new double[order+1];
	
	//initialization
	for (int n=0; n < _order+1;n++){
		momX[n] = 0.;
		momY[n] = 0.;
	}
	
	for (int n=0; n<_order+1;n++)
		for (int m=0; m<_order+1;m++)
			momentXY[n][m]=0.;
	
	momX[0]=1.0;
	momY[0]=1.0;
		
	bunch->compress();
	total_macrosize = 0.;
	if(has_msize > 0){
		ParticleMacroSize* macroSizeAttr = (ParticleMacroSize*) bunch->getParticleAttributes("macrosize");
		double m_size = 0.;
		for(int ip = 0; ip < nParts; ip++){
		
			m_size = macroSizeAttr->macrosize(ip);
			total_macrosize += m_size;
			
			if (dispersionflag > 0) {
				dispterm = getDispersion(0) * part_coord_arr[ip][5] / (bunch_kinenergy + bunch_mass) / (bunch_beta*bunch_beta);
			}
			
			for(i = 0; i < _order; i++)
                momX[i+1] = momX[i]*m_size*((part_coord_arr[ip][0] - dispterm) - xAvg);
		
			for(i = 0; i< _order; i++)
                momY[i+1] = momY[i]*m_size*(part_coord_arr[ip][2] - yAvg);
			
			for(j = 0; j<_order; j++)
				for(i=0 ; i< _order+1-j; i++){
					momentXY[i][j] += momX[i]/pow(sqrt(getBeta(0)), double(i)) * momY[j]/pow(sqrt(getBeta(1)), double(j));
					momentXY[i][j] += momX[i] * momY[j];
				}
		}
		
	} else {
		m_size = 1.0;
		for(int ip = 0; ip < nParts; ip++){
			
			if (dispersionflag > 0) {
				dispterm = getDispersion(0) * part_coord_arr[ip][5] / (bunch_kinenergy + bunch_mass) / (bunch_beta*bunch_beta);
			}
			
			for(i = 0; i < _order; i++)
                momX[i+1] = momX[i]*((part_coord_arr[ip][0] - dispterm) - xAvg);
		
			for(i = 0; i< _order; i++)
                momY[i+1] = momY[i]*(part_coord_arr[ip][2] - yAvg);
			
			for(j = 0; j<_order; j++)
				for(i=0 ; i< _order+1-j; i++)
					momentXY[i][j] += momX[i]/pow(sqrt(getBeta(0)), double(i)) * momY[j]/pow(sqrt(getBeta(1)), double(j));
				    momentXY[i][j] += momX[i] * momY[j];
			
		}
		
		total_macrosize += nParts*m_size;
		
	}
	
	//if( nMPIsize_ > 1){
	double* buff_0 = (double *) malloc (sizeof(double)*(_order+1)*(_order+1));
	double* buff_1 = (double *) malloc (sizeof(double)*(_order+1)*(_order+1));
	int count = 0;
	for(j=0; j<_order+1; j++){
		for(i=0 ; i< _order+1-j; i++){
			buff_0[count]= momentXY[i][j];
			count++;
			}
	}
	
	//MPI_Allreduce(buff_0, buff_1, count, MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	
	ORBIT_MPI_Allreduce(buff_0, buff_1, count, MPI_DOUBLE, MPI_SUM, bunch->getMPI_Comm_Local()->comm);
		
	count = 0;
	for(j=0; j<_order+1; j++){
		for(i=0 ; i< _order+1-j; i++){
			momentXY[i][j] = buff_1[count];
			count++;
		}
	}
	free(buff_0);
	free(buff_1);
		
	for(i=0; i< _order+1; i++)
        for(j=0; j< _order+1-i ; j++)
			momentXY[i][j] /= total_macrosize;
	
	momentXY[0][0] = 1.;  // 0th moment
	momentXY[1][0] = xAvg;
	momentXY[0][1] = yAvg;

//return momentXY;
	
}




/** Returns the centered correlation <(x-<x>)*(y-<y>)> = <x*y> - <x>*<y> */
double BunchTwissAnalysis::getCorrelation(int ic, int jc){
	if(ic < 0 || ic > 5 || jc < 0 || jc >5) return 0.;
	return (corr_arr[ic+6*jc] - avg_arr[ic]* avg_arr[jc]);
}

/** Returns the XY moment of the beam */
double BunchTwissAnalysis::getBunchMoment(int i, int j){
	if(i < 0 || i > _order || j < 0 || j > _order) return 0.;
	return momentXY[i][j];
}

/** Returns the average value for coordinate with index ic */
double BunchTwissAnalysis::getAverage(int ic){
	if(ic < 0 || ic > 5 ) return 0.;
	return avg_arr[ic];
}

/** Returns the total number of analysed macroparticles */
int BunchTwissAnalysis::getGlobalCount(){
	return count;
}

/** Returns the total macrosize */
double BunchTwissAnalysis::getGlobalMacrosize(){
	return total_macrosize;
}

/** Returns the emittance for index 0,1,2 - x,y,z planes. */
double BunchTwissAnalysis::getEmittance(int ic)
{	
	// for x and y the pure betatron emittance is computed (subtracting the dispersive contribution)
	if(ic < 0 || ic > 2 ) return 0.;
	double x2_avg = fabs(this->getCorrelation(2*ic,2*ic));
	double xp2_avg = fabs(this->getCorrelation(2*ic+1,2*ic+1));
	double x_xp_avg = this->getCorrelation(2*ic,2*ic+1);
	double x_dE_avg = this->getCorrelation(2*ic,5);
	double xp_dE_avg = this->getCorrelation(2*ic+1,5);	
	double dE2_avg = fabs(this->getCorrelation(5,5));	
	double emitt_rms;
	if(ic==2){
		emitt_rms = sqrt(fabs(x2_avg*xp2_avg - x_xp_avg*x_xp_avg));
	} else {
		emitt_rms = sqrt(fabs( (x2_avg - x_dE_avg*x_dE_avg/dE2_avg) * (xp2_avg - xp_dE_avg*xp_dE_avg/dE2_avg) 
						- (x_xp_avg - x_dE_avg*xp_dE_avg/dE2_avg) * (x_xp_avg - x_dE_avg*xp_dE_avg/dE2_avg) ));
	}						
	return emitt_rms;
}

/** Returns the normalized betatron emittance for index 0,1 - x,y planes. */
double BunchTwissAnalysis::getEmittanceNormalized(int ic)
{
	if(ic < 0 || ic > 1 ) return 0.;
	return this->getEmittance(ic) * bunch_gamma * bunch_beta;
}

/** Returns Twiss alpha (without dispersive part for x,y) for index 0,1,2 - x,y,z planes.*/
double BunchTwissAnalysis::getAlpha(int ic)
{
	if(ic < 0 || ic > 1 ) return 0.;
	double x_xp_avg = this->getCorrelation(2*ic,2*ic+1);
	double x_dE_avg = this->getCorrelation(2*ic, 5);
	double xp_dE_avg = this->getCorrelation(2*ic+1, 5);
	double dE2_avg = fabs(this->getCorrelation(5, 5));
	double alpha;
	if(ic == 2){
		alpha = - x_xp_avg/this->getEmittance(ic);
	} else {
		alpha = -(x_xp_avg - x_dE_avg * xp_dE_avg / dE2_avg) / this->getEmittance(ic);
	}	
	return alpha;
}

/** Returns Twiss beta (without dispersive part for x,y) for index 0,1,2 - x,y,z planes.*/
double BunchTwissAnalysis::getBeta(int ic)
{
	if(ic < 0 || ic > 1 ) return 0.;
	double x2_avg = fabs(this->getCorrelation(2*ic,2*ic));
	double x_dE_avg = this->getCorrelation(2*ic, 5);
	double dE2_avg = fabs(this->getCorrelation(5, 5));
	double beta;
	if(ic == 2){
		beta = x2_avg/this->getEmittance(ic);
	} else {
		beta = (x2_avg - x_dE_avg * x_dE_avg / dE2_avg) / this->getEmittance(ic);
	}
	return beta;	
}

/** Returns Twiss gamma (without dispersive part for x,y) for index 0,1,2 - x,y,z planes.*/
double BunchTwissAnalysis::getGamma(int ic)
{	
	if(ic < 0 || ic > 1 ) return 0.;
	double xp2_avg = fabs(this->getCorrelation(2*ic+1,2*ic+1));
	double xp_dE_avg = this->getCorrelation(2*ic+1, 5);
	double dE2_avg = fabs(this->getCorrelation(5, 5));
	double gamma;
	if(ic == 2){
		gamma = xp2_avg/this->getEmittance(ic);
	} else {
		gamma = (xp2_avg - xp_dE_avg * xp_dE_avg / dE2_avg) / this->getEmittance(ic);
	}
	return gamma;
}

/** Returns Twiss dispersion function for index 0,1 - x,y planes.*/
double BunchTwissAnalysis::getDispersion(int ic)
{	
	if(ic < 0 || ic > 1 ) return 0.;
	double x_dE_avg = this->getCorrelation(2*ic, 5);
	double dE2_avg = fabs(this->getCorrelation(5, 5));
	double dispersion = x_dE_avg/dE2_avg * bunch_momentum * bunch_beta;
	return dispersion;
}

/** Returns Twiss dispersion_prime function for index 0,1 - x,y planes.*/
double BunchTwissAnalysis::getDispersionDerivative(int ic)
{	
	if(ic < 0 || ic > 1 ) return 0.;
	double xp_dE_avg = this->getCorrelation(2*ic+1, 5);
	double dE2_avg = fabs(this->getCorrelation(5, 5));
	double dispersion_prime = xp_dE_avg/dE2_avg * bunch_momentum * bunch_beta;
	return dispersion_prime;
}

/** Returns the effective emittance for index 0,1 - x,y planes. */
double BunchTwissAnalysis::getEffectiveEmittance(int ic)
{
	if(ic < 0 || ic > 1 ) return 0.;
	double x_avg = this->getAverage(2*ic);
	double xp_avg = this->getAverage(2*ic+1);
	double x2_avg = fabs(this->getCorrelation(2*ic,2*ic));
	double xp2_avg = fabs(this->getCorrelation(2*ic+1,2*ic+1));
	double x_xp_avg = this->getCorrelation(2*ic,2*ic+1);
	double emitt_rms =  sqrt(fabs(x2_avg*xp2_avg - x_xp_avg*x_xp_avg));
	return emitt_rms;
}

/** Returns effective Twiss alpha for index 0,1 - x,y planes.*/
double BunchTwissAnalysis::getEffectiveAlpha(int ic)
{
	if(ic < 0 || ic > 1 ) return 0.;
	double x_avg = this->getAverage(2*ic);
	double xp_avg = this->getAverage(2*ic+1);
	double x2_avg = fabs(this->getCorrelation(2*ic,2*ic));
	double xp2_avg = fabs(this->getCorrelation(2*ic+1,2*ic+1));
	double x_xp_avg = this->getCorrelation(2*ic,2*ic+1);
	double emitt2_rms = x2_avg*xp2_avg - x_xp_avg*x_xp_avg;
	if(	emitt2_rms <= 0.) return 0.;
	double emitt_rms =  sqrt(emitt2_rms);
	double alpha = - x_xp_avg/emitt_rms;
	return alpha;
}

/** Returns effective Twiss beta for index 0,1 - x,y planes.*/
double BunchTwissAnalysis::getEffectiveBeta(int ic)
{
	if(ic < 0 || ic > 1 ) return 0.;
	double x_avg = this->getAverage(2*ic);
	double xp_avg = this->getAverage(2*ic+1);
	double x2_avg = fabs(this->getCorrelation(2*ic,2*ic));
	double xp2_avg = fabs(this->getCorrelation(2*ic+1,2*ic+1));
	double x_xp_avg = this->getCorrelation(2*ic,2*ic+1);
	double emitt2_rms = x2_avg*xp2_avg - x_xp_avg*x_xp_avg;
	if(	emitt2_rms <= 0.) return 0.;
	double emitt_rms =  sqrt(emitt2_rms);
	double beta = x2_avg/emitt_rms;
	return beta;	
}

/** Returns effective Twiss gamma for index 0,1 - x,y planes.*/
double BunchTwissAnalysis::getEffectiveGamma(int ic)
{	
	if(ic < 0 || ic > 1 ) return 0.;
	double x_avg = this->getAverage(2*ic);
	double xp_avg = this->getAverage(2*ic+1);
	double x2_avg = fabs(this->getCorrelation(2*ic,2*ic));
	double xp2_avg = fabs(this->getCorrelation(2*ic+1,2*ic+1));
	double x_xp_avg = this->getCorrelation(2*ic,2*ic+1);
	double emitt2_rms = x2_avg*xp2_avg - x_xp_avg*x_xp_avg;
	if(	emitt2_rms <= 0.) return DBL_MAX;
	double emitt_rms =  sqrt(emitt2_rms);
	double gamma = xp2_avg/emitt_rms;
	return gamma;
}


