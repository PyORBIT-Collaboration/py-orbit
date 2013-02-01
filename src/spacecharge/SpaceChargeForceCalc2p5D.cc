/////////////////////////////////////////////////////////////////////////////
//
// FILE NAME
//   SpaceChargeForceCalc2p5D.cc
//
//   08/02/10
//
// DESCRIPTION
//   Calculate the space charge effect of the bunch in the 2.5D  
//
/////////////////////////////////////////////////////////////////////////////

#include "Grid1D.hh"
#include "Grid2D.hh"
#include "ForceSolverFFT2D.hh"
#include "SpaceChargeForceCalc2p5D.hh"
#include "BufferStore.hh"

#include <iostream>
#include <cmath>
#include <cfloat>

using namespace OrbitUtils;

SpaceChargeForceCalc2p5D::SpaceChargeForceCalc2p5D(int xSize, int ySize, int zSize): CppPyWrapper(NULL)
{
	forceSolver = new ForceSolverFFT2D(xSize, ySize);
	rhoGrid = new Grid2D(xSize, ySize);
	forceGridX = new Grid2D(xSize, ySize);
	forceGridY = new Grid2D(xSize, ySize);
	zGrid = new Grid1D(zSize);	
	bunchExtremaCalc = new BunchExtremaCalculator();
}

SpaceChargeForceCalc2p5D::~SpaceChargeForceCalc2p5D(){
	delete forceSolver;
	if(rhoGrid->getPyWrapper() != NULL){
		Py_DECREF(rhoGrid->getPyWrapper());
	} else {
		delete rhoGrid;
	}
	if(forceGridX->getPyWrapper() != NULL){
		Py_DECREF(forceGridX->getPyWrapper());
	} else {
		delete forceGridX;
	}
	if(forceGridY->getPyWrapper() != NULL){
		Py_DECREF(forceGridY->getPyWrapper());
	} else {
		delete forceGridY;
	}
	if(zGrid->getPyWrapper() != NULL){
		Py_DECREF(zGrid->getPyWrapper());
	} else {
		delete zGrid;
	}
	delete bunchExtremaCalc;
}

Grid2D* SpaceChargeForceCalc2p5D::getRhoGrid(){
	return rhoGrid;
}

Grid2D* SpaceChargeForceCalc2p5D::getForceGridX(){
	return forceGridX;
}

Grid2D* SpaceChargeForceCalc2p5D::getForceGridY(){
	return forceGridY;
}

Grid1D* SpaceChargeForceCalc2p5D::getLongGrid(){
	return zGrid;
}

void SpaceChargeForceCalc2p5D::trackBunch(Bunch* bunch, double length){

	
	int nPartsGlobal = bunch->getSizeGlobal();
	if(nPartsGlobal < 2) return;	
	
	double totalMacrosize = 0.;
	
	this->bunchAnalysis(bunch, totalMacrosize);
		
	double z_step = zGrid->getStepZ();

	//calculate phiGrid
	forceSolver->findForce(rhoGrid, forceGridX, forceGridY);
	
	SyncPart* syncPart = bunch->getSyncPart();	
	double factor = 2*length*bunch->getClassicalRadius()*pow(bunch->getCharge(),2)/(pow(syncPart->getBeta(),2)*pow(syncPart->getGamma(),3));	
		
	factor = factor/(z_step*totalMacrosize);	
	double Lfactor = 0.;
	double x,y,z,fx,fy;
		
	for (int i = 0, n = bunch->getSize(); i < n; i++){
		x = bunch->x(i);
		y = bunch->y(i);
		z = bunch->z(i);
		
		forceGridX->interpolateBilinear(x,y,fx);
		forceGridY->interpolateBilinear(x,y,fy);
		
		//Lfactor = - zGrid->getValue(z) * factor;
		Lfactor =  zGrid->getValue(z) * factor;
		bunch->xp(i) += fx * Lfactor;
		bunch->yp(i) += fy * Lfactor;
	}
}

void SpaceChargeForceCalc2p5D::bunchAnalysis(Bunch* bunch, double& totalMacrosize){

	double xMin, xMax, yMin, yMax, zMin, zMax;
	
	bunchExtremaCalc->getExtremaXYZ(bunch, xMin, xMax, yMin, yMax, zMin, zMax);
	
	//bunchExtremaCalc->getXY_NRMS(bunch, N, xMin, xMax, yMin, yMax, zMin, zMax);

	//check if the beam size is not zero
	if( xMin >=  xMax || yMin >=  yMax || zMin >=  zMax){
		int rank = 0;
		ORBIT_MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		if(rank == 0){
			std::cerr << "SpaceChargeForceCalc2p5D::bunchAnalysis(bunch,...)" << std::endl
         				<< "The bunch min and max sizes are wrong! Cannot calculate space charge!" << std::endl
								<< "x min ="<< xMin <<" max="<< xMax << std::endl
								<< "y min ="<< yMin <<" max="<< yMax << std::endl
								<< "z min ="<< zMin <<" max="<< zMax << std::endl
								<< "Stop."<< std::endl;
		}
		ORBIT_MPI_Finalize();
	}
	
	//check if the beam size is not zero 
	if(zMin >=  zMax){
		int rank = 0;
		ORBIT_MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		if(rank == 0){
			std::cerr << "SpaceChargeForceCalc2p5D::bunchAnalysis(bunch,...)" << std::endl
         				<< "The bunch min and max sizes are wrong! Cannot calculate space charge!" << std::endl
								<< "z min ="<< zMin <<" max="<< zMax << std::endl
								<< "Stop."<< std::endl;
		}
		ORBIT_MPI_Finalize();
	}

	//set Grids' limits
	rhoGrid->setGridX(xMin,xMax);
	rhoGrid->setGridY(yMin,yMax);
	forceGridX->setGridX(xMin,xMax);
	forceGridX->setGridY(yMin,yMax);
	forceGridY->setGridX(xMin,xMax);
	forceGridY->setGridY(yMin,yMax);
	forceSolver->setGridXY(xMin, xMax, yMin, yMax);
	zGrid->setGridZ(zMin,zMax);

	//sizes of the grids are set up
	//bin rho&z Bunch to the Grid
	rhoGrid->setZero();
	
	rhoGrid->binBunchBilinear(bunch);
	
	zGrid->setZero();
	zGrid->binBunch(bunch);

	rhoGrid->synchronizeMPI(bunch->getMPI_Comm_Local());
	zGrid->synchronizeMPI(bunch->getMPI_Comm_Local());
	
	totalMacrosize = 0.;	
	int nZ = zGrid->getSizeZ();	
	for(int iz = 0; iz < nZ; iz++){
		totalMacrosize += zGrid->getValueOnGrid(iz);		
	}
	
}
