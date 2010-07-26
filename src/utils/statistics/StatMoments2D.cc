#include "StatMoments2D.hh"

#include <iostream>
#include <cmath>
#include <cfloat>

#include "BufferStore.hh"

using namespace OrbitUtils;

/** Constructor with max order = 2 by default */
StatMoments2D::StatMoments2D(): CppPyWrapper(NULL)
{
  max_order = 2;
	stat_arr = NULL;
	makeArrays();	
}

/** Constructor with maximal order parameter */
StatMoments2D::StatMoments2D(int maxOrder): CppPyWrapper(NULL)
{
  max_order = maxOrder;
	if(max_order < 2) max_order = 2;
	stat_arr = NULL;
	makeArrays();
}

/** Destructor */
StatMoments2D::~StatMoments2D()
{
	if(stat_arr != NULL){
			for(int iu = 0, iu_max = max_order + 1; iu < iu_max; iu++){
				delete [] stat_arr[iu];
			}
			delete [] stat_arr;
			stat_arr = NULL;
	}	
}

/** Initialize all internal arrays to get ready to gather statistical information. */   
void StatMoments2D::clean() 
{
	u_max  =  - DBL_MAX;
	u_min  =    DBL_MAX;
	up_max =  - DBL_MAX;
	up_min =    DBL_MAX;
	count = 0;
	for(int iu = 0, iu_max = max_order + 1; iu < iu_max; iu++){
		for(int iup = 0, iup_max = max_order + 1; iup < iup_max; iup++){
			stat_arr[iu][iup] = 0.;
		}
	}	
}

/** make arrays */
void StatMoments2D::makeArrays() 
{
	if(stat_arr != NULL){
			for(int iu = 0, iu_max = max_order + 1; iu < iu_max; iu++){
				delete [] stat_arr[iu];
			}
			delete [] stat_arr;
			stat_arr = NULL;
	}
	
	stat_arr = new double*[max_order + 1];
	for(int iu = 0, iu_max = max_order + 1; iu < iu_max; iu++){
		stat_arr[iu] = new double[max_order + 1];
	}		
	clean();
}
		
/** Sets the maximal order of the moments */   
void StatMoments2D::setMaxOrder(int maxOrder)
{
  max_order = maxOrder;
	if(max_order < 2) max_order = 2;
	makeArrays();	
}
		
/** Returns the maximal order of the moments */   
int StatMoments2D::getMaxOrder()
{
	return max_order;
}

/** Takes into account the one point (u,up) */   	
void StatMoments2D::account(double u, double up)
{
	count += 1;
	
	if(u < u_min) u_min = u;
	if(up < up_min) up_min = up;
	if(u > u_max) u_max = u;
	if(up > up_max) up_max = up;
	
	if(max_order < 3){
		stat_arr[0][0] += 1.;
		stat_arr[0][1] += up;
		stat_arr[0][2] += up*up;
		
		stat_arr[1][0] += u;
		stat_arr[1][1] += u*up;
		stat_arr[1][2] += u*up*up;	
		
		stat_arr[2][0] += u*u;
		stat_arr[2][1] += u*u*up;
		stat_arr[2][2] += u*u*up*up;		
	} else {
		for(int iu = 0, iu_max = max_order + 1; iu < iu_max; iu++){
			for(int iup = 0, iup_max = max_order + 1; iup < iup_max; iup++){
				stat_arr[iu][iup] += pow(u,iu)*pow(up,iup);
			}
		}	
	}
}

/** Returns the statistical moment of the distribution with particular orders for u and up */ 	
double StatMoments2D::getStatMoment(int order_u,int order_up)
{
	if(order_u > max_order || order_up > max_order || count == 0) return 0.;
	return stat_arr[order_u][order_up]/count;
}

/** Returns the statistical moment of the distribution with particular order for u*/ 	
double StatMoments2D::getStatMomentU(int order_u)
{
	if(order_u > max_order || count == 0) return 0.;
	return stat_arr[order_u][0]/count;
}

/** Returns the statistical moment of the distribution with particular order for up*/ 	
double StatMoments2D::getStatMomentUP(int order_up)
{
	if(order_up > max_order || count == 0) return 0.;
	return stat_arr[0][order_up]/count;
}

/** Returns the minimal value of u */ 	
double StatMoments2D::getMinU()
{
	return u_min;
}

/** Returns the maximal value of u */ 	
double StatMoments2D::getMaxU()
{
	return u_max;
}

/** Returns the minimal value of up */ 	
double StatMoments2D::getMinUP()
{
	return up_min;
}

/** Returns the maximal value of up */ 	
double StatMoments2D::getMaxUP()
{
	return up_max;
}

/** Returns the number of points in the statistic */ 	
int StatMoments2D::getCount()
{
	return count;
}

/** It will synchronize the moments through the MPI communicator */ 		
int StatMoments2D::synchronizeMPI(pyORBIT_MPI_Comm* pyComm)
{
	int mpi_size = (max_order+1)*(max_order+1);
	int buff_index0 = -1;
	int buff_index1 = -1;	
  double* inArr  = OrbitUtils::BufferStore::getBufferStore()->getFreeDoubleArr(buff_index0,mpi_size);
  double* outArr = OrbitUtils::BufferStore::getBufferStore()->getFreeDoubleArr(buff_index1,mpi_size);

	int ii = 0;
	for(int iu = 0, iu_max = max_order + 1; iu < iu_max; iu++){
		for(int iup = 0, iup_max = max_order + 1; iup < iup_max; iup++){
			inArr[ii] = stat_arr[iu][iup];
			ii += 1;
		}
	}

	if(pyComm == NULL) {
		ORBIT_MPI_Allreduce(inArr,outArr,mpi_size,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	} else {	
		ORBIT_MPI_Allreduce(inArr,outArr,mpi_size,MPI_DOUBLE,MPI_SUM,pyComm->comm);
	}
		
	ii = 0;
	for(int iu = 0, iu_max = max_order + 1; iu < iu_max; iu++){
		for(int iup = 0, iup_max = max_order + 1; iup < iup_max; iup++){
			stat_arr[iu][iup] = outArr[ii];
			ii += 1;
		}
	}	
	
	int count_MPI = -1;
	if(pyComm == NULL) {
		ORBIT_MPI_Allreduce(&count,&count_MPI,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
	}	else {
		ORBIT_MPI_Allreduce(&count,&count_MPI,1,MPI_INT,MPI_SUM,pyComm->comm);
	}
	count = count_MPI;
	
	OrbitUtils::BufferStore::getBufferStore()->setUnusedDoubleArr(buff_index0);
	OrbitUtils::BufferStore::getBufferStore()->setUnusedDoubleArr(buff_index1);		
}


/** Returns the emittance */
double StatMoments2D::getEmittance()
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
double StatMoments2D::getAlpha()
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
double StatMoments2D::getBeta()
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
double StatMoments2D::getGamma()
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
double StatMoments2D::getRmsU()
{
	double x_avg = getStatMomentU(1);
	double x2_avg = fabs(getStatMomentU(2) - x_avg*x_avg);
	return sqrt(x2_avg);
}

/** Returns the rms value of up */ 	
double StatMoments2D::getRmsUP()
{
	double xp_avg = getStatMomentUP(1);
	double xp2_avg = fabs(getStatMomentUP(2) - xp_avg*xp_avg);
	return sqrt(xp2_avg);
}




