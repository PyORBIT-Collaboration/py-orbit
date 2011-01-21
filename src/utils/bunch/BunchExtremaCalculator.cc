//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//   BunchExtremaCalculator.cc
//
// AUTHOR
//    A. Shishlo
//
// CREATED
//    10/11/2010
//
// DESCRIPTION
//    A class calculates the extrema of the particles coordinates in the bunch.
//
///////////////////////////////////////////////////////////////////////////

#include "orbit_mpi.hh"
#include "BunchExtremaCalculator.hh"
#include "BufferStore.hh"

#include <cstdlib>
#include <iostream>
#include <cmath>
#include <cfloat>

using namespace OrbitUtils;

BunchExtremaCalculator::BunchExtremaCalculator(): CppPyWrapper(NULL)
{
}

BunchExtremaCalculator::~BunchExtremaCalculator()
{
}

/** The method calculates the extrema of the particles coordinates in the bunch. */
void BunchExtremaCalculator::getExtremaXYZ(Bunch* bunch, 
	double& xMin, double& xMax, 
	double& yMin, double& yMax, 
	double& zMin, double& zMax)
{

	int buff_index0 = 0;
	int buff_index1 = 0;
	double* gridLimArr  = BufferStore::getBufferStore()->getFreeDoubleArr(buff_index0,6);
	double* gridLimArr_out = BufferStore::getBufferStore()->getFreeDoubleArr(buff_index1,6);
	for (int i = 0; i < 3; i++){
		gridLimArr[2*i] = DBL_MAX;
		gridLimArr[2*i+1] = -DBL_MAX;
	}
	
	double** partArr=bunch->coordArr();
	double* coordArr = NULL;
	for (int ip = 0, n = bunch->getSize(); ip < n; ip++){
		coordArr = partArr[ip];
		if(coordArr[0] < gridLimArr[0]) gridLimArr[0] = coordArr[0];
		if(coordArr[0] > gridLimArr[1]) gridLimArr[1] = coordArr[0];
		if(coordArr[2] < gridLimArr[2]) gridLimArr[2] = coordArr[2];
		if(coordArr[2] > gridLimArr[3]) gridLimArr[3] = coordArr[2];	
		if(coordArr[4] < gridLimArr[4]) gridLimArr[4] = coordArr[4];
		if(coordArr[4] > gridLimArr[5]) gridLimArr[5] = coordArr[4];		
	}
	
	gridLimArr[0] = - gridLimArr[0];
	gridLimArr[2] = - gridLimArr[2];
	gridLimArr[4] = - gridLimArr[4];

	ORBIT_MPI_Allreduce(gridLimArr,gridLimArr_out,6,MPI_DOUBLE,MPI_MAX,bunch->getMPI_Comm_Local()->comm);

	gridLimArr_out[0] = - gridLimArr_out[0];
	gridLimArr_out[2] = - gridLimArr_out[2];
	gridLimArr_out[4] = - gridLimArr_out[4];
	
  xMin = gridLimArr_out[0];
  xMax = gridLimArr_out[1];	
  yMin = gridLimArr_out[2];
  yMax = gridLimArr_out[3];	
  zMin = gridLimArr_out[4];
  zMax = gridLimArr_out[5];	

	OrbitUtils::BufferStore::getBufferStore()->setUnusedDoubleArr(buff_index0);
	OrbitUtils::BufferStore::getBufferStore()->setUnusedDoubleArr(buff_index1);	
}
	
/** The method calculates the z extrema of the particles coordinates in the bunch. */
void BunchExtremaCalculator::getExtremaZ(Bunch* bunch, 
	double& zMin, double& zMax)
{
	int buff_index0 = 0;
	int buff_index1 = 0;
	double* gridLimArr  = BufferStore::getBufferStore()->getFreeDoubleArr(buff_index0,2);
	double* gridLimArr_out = BufferStore::getBufferStore()->getFreeDoubleArr(buff_index1,2);
	gridLimArr[0] = DBL_MAX;
	gridLimArr[1] = -DBL_MAX;	
	
	double** partArr=bunch->coordArr();
	double* coordArr = NULL;
	for (int ip = 0, n = bunch->getSize(); ip < n; ip++){
		coordArr = partArr[ip];
		if(coordArr[4] < gridLimArr[0]) gridLimArr[0] = coordArr[4];
		if(coordArr[4] > gridLimArr[1]) gridLimArr[1] = coordArr[4];		
	}
	
	gridLimArr[0] = - gridLimArr[0];

	ORBIT_MPI_Allreduce(gridLimArr,gridLimArr_out,2,MPI_DOUBLE,MPI_MAX,bunch->getMPI_Comm_Local()->comm);

	gridLimArr_out[0] = - gridLimArr_out[0];
	
  zMin = gridLimArr_out[0];
  zMax = gridLimArr_out[1];	

	OrbitUtils::BufferStore::getBufferStore()->setUnusedDoubleArr(buff_index0);
	OrbitUtils::BufferStore::getBufferStore()->setUnusedDoubleArr(buff_index1);	
}
