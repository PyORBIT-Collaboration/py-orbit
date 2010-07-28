/////////////////////////////////////////////////////////////////////////////
//
// FILE NAME
//   Grid1D.cc
//
//   06/25/10
//
// DESCRIPTION
//   repersents a 1D grid
//
/////////////////////////////////////////////////////////////////////////////

#include "Grid1D.hh"
#include "Bunch.hh"
#include "ParticleMacroSize.hh"
#include "BufferStore.hh"

#include <iostream>

using namespace OrbitUtils;
//zBins - grid size [zBins]
//zSize - geometry parameters [m]

// Constructor
Grid1D::Grid1D(int zSize): CppPyWrapper(NULL){
	zSize_ = zSize;
	zMin_ = -1.0; 
	zMax_ = +1.0; 
	init();
	setZero();
}

Grid1D::Grid1D(int zSize, double zMin, double zMax): CppPyWrapper(NULL){
	zSize_ = zSize;
	zMin_ = zMin;
	zMax_ = zMax;
	init();
	setZero();
}

void Grid1D::init(){
	dz_ = (zMax_ - zMin_)/(zSize_ - 1);
	arr_ = new double[zSize_];
}

// Destructor
Grid1D::~Grid1D(){
	delete [] arr_;
}

/** Sets the value to the one point of the grid  */
void Grid1D::setValue(double value,int iZ){
	arr_[iZ] = value;
}

double Grid1D::getValue(double z){
	int iZ;
	double Wzm, Wz0, Wzp;
	double zFrac, zFrac2;
	getIndAndFracZ(z,iZ,zFrac);
	zFrac2 = zFrac * zFrac;
	
	if(zSize_ > 2){
		Wzm = 0.5 * (0.25 - zFrac + zFrac2);
		Wz0 = 0.75 - zFrac2;
		Wzp = 0.5 * (0.25 + zFrac + zFrac2);
	}
	else if(zSize_ == 2){
		Wzm = 1.0 - zFrac; // for zInd=0
		Wz0 = 0.0;
		Wzp = zFrac;       // for zInd=1
	}
	else if(zSize_ == 1){
		Wzm = 0.0;
		Wz0 = 1.0;
		Wzp = 0.0; 
	}
// Calculate a longitudinal weighting factor(LPosFactor) for Transverse
// space charge calculations: 
// (== local average macroPart density per rad / average density)
	//calc value
	double value= Wzm * arr_[iZ-1] + Wz0 * arr_[iZ] + Wzp * arr_[iZ+1];
	return value;
}

void Grid1D::binBunch(Bunch* bunch){
	bunch->compress();
	
	double** part_coord_arr = bunch->coordArr();
	int has_msize = bunch->hasParticleAttributes("macrosize");
	if(has_msize > 0){
		ParticleMacroSize* macroSizeAttr = (ParticleMacroSize*) bunch->getParticleAttributes("macrosize");
		double m_size = 0.;
		for(int i = 0, n = bunch->getSize(); i < n; i++){
			m_size = macroSizeAttr->macrosize(i);
			binValue(m_size,part_coord_arr[i][4]);
		}	
		return;
	}
	
	double m_size = bunch->getMacroSize();
	for(int i = 0, n = bunch->getSize(); i < n; i++){
		binValue(m_size,part_coord_arr[i][4]);	
	}
}

void Grid1D::binValue(double macroSize,double z){
	int iZ;
	double zFrac,zFrac2;
	double Wzm,Wz0,Wzp;
	Wzm = Wz0 = Wzp = 0.0;
	
	getIndAndFracZ(z,iZ,zFrac);
	zFrac2 = zFrac * zFrac;
	  
	if(zSize_ > 2){
		Wzm = 0.5 * (0.25 - zFrac + zFrac2);
		Wz0 = 0.75 - zFrac2;
		Wzp = 0.5 * (0.25 + zFrac + zFrac2);
	}
	else if(zSize_ == 2){
		Wzm = 1.0 - zFrac; // for zInd=0
		Wz0 = 0.0;
		Wzp = zFrac;       // for zInd=1
	}

	//Add weight of particle
	  if( zSize_ >= 3){
	  	arr_[iZ-1] += Wzm * macroSize;
	  	arr_[iZ] += Wz0 * macroSize;
	  	arr_[iZ+1] += Wzp * macroSize;
	  }
	  else if( zSize_ == 2){
	  	arr_[0] += Wzm * macroSize;
	  	arr_[1] += Wzp * macroSize;	  	  
	  }
	  else if( zSize_ == 1){
	  	arr_[0] += Wz0 * macroSize;
	  }
}

void Grid1D::calcGradient(double z,double& ez){
	double Wzm, Wz0, Wzp;
	double dWzm, dWz0, dWzp;
	int iZ;
	double zFrac,zFrac2;
	Wzm = Wz0 = Wzp = 0.0;
	dWzm = dWz0 = dWzp = 0.0;
	
	getIndAndFracZ(z,iZ,zFrac);
	zFrac2 = zFrac * zFrac; 
	
	if(zSize_ > 2){
		Wzm = 0.5 * (0.25 - zFrac + zFrac2);
		Wz0 = 0.75 - zFrac2;
		Wzp = 0.5 * (0.25 + zFrac + zFrac2);
		dWzm = (-1.0)*  (0.5 - zFrac)/dz_;
		dWz0 = (-1.0)*    2. * zFrac /dz_;
		dWzp = (-1.0)* -(0.5 + zFrac)/dz_;
	}
	else if(zSize_ == 2){
		Wzm = 1.0 - zFrac; // for zInd=0
		Wz0 = 0.0;
		Wzp = zFrac;       // for zInd=1
		dWzm = (-1.0)*  1.0/dz_; // for zInd=0
		dWz0 =  0.0;
		dWzp = (-1.0)* -1.0/dz_; // for zInd=1		
	}
	//calculate gradient
	ez = Wzm * dWzm * arr_[iZ-1] + Wz0 * dWz0 * arr_[iZ] + Wzp * dWzp * arr_[iZ+1];
}



/** Sets all grid points to zero */	
void Grid1D::setZero(){
	for(int i = 0; i < zSize_; i++){
		arr_[i] = 0.;
	}
}
	

void Grid1D::getIndAndFracZ(double z, int& ind, double& frac){
	if(zSize_ > 1){
	    ind  = int ( (z - zMin_)/dz_ + 0.5 );

	    if(zSize_ > 2){
	        //cut off edge for three point interpolation
	        if(ind < 1) ind = 1;
		if(ind > (zSize_-2))
		    ind =  zSize_ - 2;
	    	    frac = (z - (zMin_ + ind * dz_))/dz_;
		}
	    if(zSize_ == 2){
		frac = (z - zMin_)/dz_;
		if(ind < 0) {ind = 0; frac = 0.0;}
		else if(ind > 1) {ind = 1; frac = 1.0;}
		}
	}
	else {dz_=0.0; ind = 0; frac = 0.0;}
}

/** Returns the min z in the grid points */ 
double Grid1D::getMinZ(){return zMin_;};

/** Returns the max z in the grid points */ 
double Grid1D::getMaxZ(){return zMax_;};

/** Sets z-grid */
void Grid1D::setGridZ(double zMin, double zMax){
	zMax_ = zMin;
	zMax_ = zMax;
	dz_ = (zMax_ - zMin_)/(zSize_ -1);
	setZero();
}

double Grid1D::getGridZ(int index){
	return zMin_ + index*dz_;
}

/** Returns the grid size in z-direction */
int Grid1D::getSizeZ(){
	return zSize_;
}

/**synchronizeMPI */
void Grid1D::synchronizeMPI(pyORBIT_MPI_Comm* pyComm){
//MPI
	int size_MPI = zSize_ * zSize_;
	int buff_index0 = 0;
	int buff_index1 = 0;
	double* inArr  = BufferStore::getBufferStore()->getFreeDoubleArr(buff_index0,size_MPI);
	double* outArr = BufferStore::getBufferStore()->getFreeDoubleArr(buff_index1,size_MPI);
	
	for(int i = 0; i < zSize_; i++){
		inArr[i] = arr_[i];
	}
	
	ORBIT_MPI_Allreduce(inArr,outArr,size_MPI,MPI_DOUBLE,MPI_SUM,pyComm->comm);

	for(int i = 0; i < zSize_; i++){
		arr_[i] = outArr[i];	
	}
	
	OrbitUtils::BufferStore::getBufferStore()->setUnusedDoubleArr(buff_index0);
	OrbitUtils::BufferStore::getBufferStore()->setUnusedDoubleArr(buff_index1);	
	
  // ===== MPI end =====
}
