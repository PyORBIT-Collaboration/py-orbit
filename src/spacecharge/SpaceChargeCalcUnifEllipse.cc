/////////////////////////////////////////////////////////////////////////////
//
// FILE NAME
//   SpaceChargeCalcUnifEllipse.cc
//
// AUTHOR
//    A. Shishlo
//
// Created:
//   11/09/10
//
// DESCRIPTION
//  This class calculates the space charge kicks for bunch. It represent the bunch as the set 
//  of uniformly charged ellipses in the center of mass of the bunch system. 
//  The space charge kick is transformed later into the lab system.   
//
/////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <cmath>
#include <cfloat>

#include "SpaceChargeCalcUnifEllipse.hh"
#include "BufferStore.hh"

#include "ParticleMacroSize.hh"

using namespace OrbitUtils;

SpaceChargeCalcUnifEllipse::SpaceChargeCalcUnifEllipse(int nEllipses_in): CppPyWrapper(NULL)
{
	nEllipses = nEllipses_in;
  ellipsoidCalc_arr = new UniformEllipsoidFieldCalculator*[nEllipses];
	for(int ie = 0; ie < nEllipses; ie++){
		ellipsoidCalc_arr[ie] = new UniformEllipsoidFieldCalculator();
	}
	macroSizesEll_arr = (double* ) malloc (sizeof(double)*nEllipses);
	macroSizesEll_MPI_arr = (double* ) malloc (sizeof(double)*nEllipses);
	for(int ie = 0; ie < nEllipses; ie++){
		macroSizesEll_arr[ie] = 0.;
		macroSizesEll_MPI_arr[ie] = 0.;
	}
}


SpaceChargeCalcUnifEllipse::~SpaceChargeCalcUnifEllipse(){
	for(int ie = 0; ie < nEllipses; ie++){
		delete ellipsoidCalc_arr[ie];
	}	
	delete [] ellipsoidCalc_arr;
	
	free(macroSizesEll_arr);
	free(macroSizesEll_MPI_arr);
}


void SpaceChargeCalcUnifEllipse::trackBunch(Bunch* bunch, double length){

	int nPartsGlobal = bunch->getSizeGlobal();
	if(nPartsGlobal < 3) return;

	SyncPart* syncPart = bunch->getSyncPart();	
	double beta = syncPart->getBeta();
	double gamma = syncPart->getGamma();
	
	//analyse the bunch and make the ellipsoid filed sources
	this->bunchAnalysis(bunch);

	double x,y,z,ex,ey,ez;
	for (int i = 0, n = bunch->getSize(); i < n; i++){
		x = bunch->x(i) - x_center;
		y = bunch->y(i) - y_center;
		z = (bunch->z(i) - z_center)*gamma;
		this->calculaField(x,y,z,ex,ey,ez);
		//calculate momentum kicks
	}
}

/** Analyses the bunch and sets up the ellipsoid filed sources */
void SpaceChargeCalcUnifEllipse::bunchAnalysis(Bunch* bunch){

	//average values for x,y,z,x2,y2,z2 and total macrosize
	int buff_index0 = 0;
	int buff_index1 = 0;
	double* coord_avg = BufferStore::getBufferStore()->getFreeDoubleArr(buff_index0,7);
	double* coord_avg_out = BufferStore::getBufferStore()->getFreeDoubleArr(buff_index1,7);
	for (int i = 0; i < 7; i++){
		coord_avg[i] = 0.;
	}	
	
	//caluclate limits and averages
	double** partArr=bunch->coordArr();
	double* coordArr = NULL;
	bunch->compress();
	double** part_coord_arr = bunch->coordArr();
	int has_msize = bunch->hasParticleAttributes("macrosize");
	if(has_msize > 0){
		ParticleMacroSize* macroSizeAttr = (ParticleMacroSize*) bunch->getParticleAttributes("macrosize");
		double m_size = 0.;
		for(int ip = 0, n = bunch->getSize(); ip < n; ip++){
			m_size = macroSizeAttr->macrosize(ip);
			coordArr = partArr[ip];
			coord_avg[0] += m_size*coordArr[0];
			coord_avg[1] += m_size*coordArr[2];
			coord_avg[2] += m_size*coordArr[4];
			coord_avg[3] += m_size*coordArr[0]*coordArr[0];
			coord_avg[4] += m_size*coordArr[2]*coordArr[2];
			coord_avg[5] += m_size*coordArr[4]*coordArr[4];	
			coord_avg[6] += m_size;
		}	
	} else {
		double m_size = bunch->getMacroSize();
		int nParts = bunch->getSize();
		coord_avg[6] = m_size*nParts;	
		for(int ip = 0; ip < nParts; ip++){
			coordArr = partArr[ip];
			coord_avg[0] += coordArr[0];
			coord_avg[1] += coordArr[2];
			coord_avg[2] += coordArr[4];
			coord_avg[3] += coordArr[0]*coordArr[0];
			coord_avg[4] += coordArr[2]*coordArr[2];
			coord_avg[5] += coordArr[4]*coordArr[4];		
		}
		for (int i = 0; i < 6; i++){
			coord_avg[i] *= m_size;
		}	
	}

	//calculates sum over all  CPUs
	ORBIT_MPI_Allreduce(coord_avg,coord_avg_out,7,MPI_DOUBLE,MPI_MAX,bunch->getMPI_Comm_Local()->comm);	
	
	double total_macrosize = coord_avg_out[6];
	
	//calculate the parameters of the biggest ellipse
	x_center = coord_avg_out[0]/total_macrosize;
	y_center = coord_avg_out[1]/total_macrosize;
	z_center = coord_avg_out[2]/total_macrosize;
	x2_avg = fabs(coord_avg_out[3]/total_macrosize - x_center*x_center);
	y2_avg = fabs(coord_avg_out[4]/total_macrosize - y_center*y_center);
	z2_avg = fabs(coord_avg_out[5]/total_macrosize - z_center*z_center);
  a2_ellips = 5.0*x2_avg;
  b2_ellips = 5.0*y2_avg;
  c2_ellips = 5.0*z2_avg;
	a_ellips = sqrt(a2_ellips);
	b_ellips = sqrt(b2_ellips);
	c_ellips = sqrt(c2_ellips);
	
	//free resources
	OrbitUtils::BufferStore::getBufferStore()->setUnusedDoubleArr(buff_index0);
	OrbitUtils::BufferStore::getBufferStore()->setUnusedDoubleArr(buff_index1);	
	
	//check if the beam size is not zero 
  if( x2_avg == 0. || y2_avg == 0.|| z2_avg == 0.){
		int rank = 0;
		ORBIT_MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		if(rank == 0){
			std::cerr << "SpaceChargeCalcUnifEllipse::bunchAnalysis(bunch,...)" << std::endl
         				<< "The bunch coords min and max sizes are wrong! Cannot calculate space charge!" << std::endl
								<<" x2_rms="<< x2_avg << std::endl
								<<" y2_rms="<< y2_avg << std::endl
								<<" z2_rms="<< z2_avg << std::endl
								<< "Stop."<< std::endl;
		}
		ORBIT_MPI_Finalize();
  }
	
	//find the distribution of the macrosizes between nEllipses
	for(int ie = 0; ie < nEllipses; ie++){
		macroSizesEll_arr[ie] = 0.;
	}
	
	double pos = 0.;
	int pos_index = 0;
	if(has_msize > 0){
		ParticleMacroSize* macroSizeAttr = (ParticleMacroSize*) bunch->getParticleAttributes("macrosize");
		double m_size = 0.;
		for(int ip = 0, n = bunch->getSize(); ip < n; ip++){
			m_size = macroSizeAttr->macrosize(ip);
			coordArr = partArr[ip];
			pos = sqrt(coordArr[0]*coordArr[0]/a2_ellips + coordArr[2]*coordArr[2]/b2_ellips + coordArr[4]*coordArr[4]/c2_ellips);
			pos_index = int(pos*nEllipses);
			if(pos_index < 0) pos_index = 0;
			if(pos_index >= nEllipses) pos_index = nEllipses - 1;
			macroSizesEll_arr[pos_index] += m_size;
		}	
	} else {
		double m_size = bunch->getMacroSize();
		int nParts = bunch->getSize();
		for(int ip = 0, n = bunch->getSize(); ip < n; ip++){
			coordArr = partArr[ip];
			pos = sqrt(coordArr[0]*coordArr[0]/a2_ellips + coordArr[2]*coordArr[2]/b2_ellips + coordArr[4]*coordArr[4]/c2_ellips);
			pos_index = int(pos*nEllipses) - 1;
			if(pos_index < 0) pos_index = 0;
			if(pos_index >= nEllipses) pos_index = nEllipses - 1;
			macroSizesEll_arr[pos_index] += m_size;
		}	
	}	
	//calculates sum over all  CPUs
	ORBIT_MPI_Allreduce(macroSizesEll_arr,macroSizesEll_MPI_arr,nEllipses,MPI_DOUBLE,MPI_MAX,bunch->getMPI_Comm_Local()->comm);
	for(int ie = 0; ie < nEllipses; ie++){
		macroSizesEll_arr[ie] = macroSizesEll_MPI_arr[ie];
	}	
	//calculate the relative volume density in each region. This density is a sum of all elipsoids
	for(int ie = 0; ie < nEllipses; ie++){
		macroSizesEll_MPI_arr[ie] /= ((ie+1)*(ie+1)*(ie+1) - ie*ie*ie);
	}
	//calculate the density for each elipsoid
	double rho_sum = 0.;
	for(int ie = (nEllipses-1); ie > 0; ie--){
		macroSizesEll_MPI_arr[ie] -= rho_sum;
		rho_sum += macroSizesEll_MPI_arr[ie];
	}
	
	//now set up the relative total charges in ellipsoids
	double q_sum = 0.;
	for(int ie = 0; ie < nEllipses; ie++){
		macroSizesEll_MPI_arr[ie] = macroSizesEll_MPI_arr[ie]*(ie+1)*(ie+1)*(ie+1);
		q_sum += macroSizesEll_MPI_arr[ie];
	}	
	double q_coeff = total_macrosize/q_sum;
	for(int ie = 0; ie < nEllipses; ie++){
		macroSizesEll_arr[ie] = macroSizesEll_MPI_arr[ie]*q_coeff;
	}
	
	//relativistic factor gamma
	double gamma = bunch->getSyncPart()->getGamma();		
	
	//now we initialize the ellipses filed calculators
	double r_max = a_ellips;
	if(r_max < b_ellips) r_max = b_ellips;
	if(r_max < c_ellips*gamma) r_max = c_ellips*gamma;
	for(int ie = 0; ie < nEllipses; ie++){
		double coeff = (ie+1.)/nEllipses;
		ellipsoidCalc_arr[ie]->setEllipsoid(a_ellips*coeff,b_ellips*coeff,c_ellips*coeff*gamma,10.*r_max*coeff);
		ellipsoidCalc_arr[ie]->setQ(macroSizesEll_arr[ie]);
	}
}

/** Calculates the electric filed in the center of the bunch sytem. */
void SpaceChargeCalcUnifEllipse::calculaField(double x,  double y,  double z, 
	                                            double& ex, double& ey, double& ez)
{
	ex = 0.;  ey = 0.; ez = 0.;
	double x2 = x*x;
	double y2 = y*y;
	double z2 = z*z;
	double ex_l,ey_l,ez_l;
	for(int ie = 0; ie < nEllipses; ie++){
		ellipsoidCalc_arr[ie]->calcField(x,y,z,x2,y2,z2,ex_l,ey_l,ez_l);
		ex += ex_l;
		ey += ey_l;
		ez += ez_l;
	}
}

