#include "BunchTwissAnalysis.hh"

#include <iostream>
#include <cmath>
#include <cfloat>
#include <cstdlib>

#include "ParticleMacroSize.hh"

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
}

/** Returns the centered correlation <(x-<x>)*(y-<y>)> = <x*y> - <x>*<y> */
double BunchTwissAnalysis::getCorrelation(int ic, int jc){
	if(ic < 0 || ic > 5 || jc < 0 || jc >5) return 0.;
	return (corr_arr[ic+6*jc] - avg_arr[ic]* avg_arr[jc]);
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
	if(ic < 0 || ic > 2 ) return 0.;
	double x_avg = this->getAverage(2*ic);
	double xp_avg = this->getAverage(2*ic+1);
	double x2_avg = fabs(this->getCorrelation(2*ic,2*ic));
	double xp2_avg = fabs(this->getCorrelation(2*ic+1,2*ic+1));
	double x_xp_avg = this->getCorrelation(2*ic,2*ic+1);
  double emitt_rms =  sqrt(fabs(x2_avg*xp2_avg - x_xp_avg*x_xp_avg));
	return emitt_rms;
}

/** Returns Twiss alpha for index 0,1,2 - x,y,z planes.*/
double BunchTwissAnalysis::getAlpha(int ic)
{
	if(ic < 0 || ic > 2 ) return 0.;
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

/** Returns Twiss beta for index 0,1,2 - x,y,z planes.*/
double BunchTwissAnalysis::getBeta(int ic)
{
	if(ic < 0 || ic > 2 ) return 0.;
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

/** Returns Twiss gamma for index 0,1,2 - x,y,z planes.*/
double BunchTwissAnalysis::getGamma(int ic)
{	
	if(ic < 0 || ic > 2 ) return 0.;
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





