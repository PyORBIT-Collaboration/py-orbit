/////////////////////////////////////////////////////////////////////////////
//
// FILE NAME
//   SpaceChargeCalc2p5D.cc
//
//   08/02/10
//
// DESCRIPTION
//   Calculate the space charge effect of the bunch in the 2.5D  
//
/////////////////////////////////////////////////////////////////////////////

#include "Grid1D.hh"
#include "Grid2D.hh"
#include "PoissonSolverFFT2D.hh"
#include "SpaceChargeCalc2p5D.hh"
#include "BufferStore.hh"

#include <iostream>
#include <cmath>
#include <cfloat>

using namespace OrbitUtils;

SpaceChargeCalc2p5D::SpaceChargeCalc2p5D(int xSize, int ySize, int zSize, double xy_ratio_in): CppPyWrapper(NULL)
{
	xy_ratio = xy_ratio_in;
	poissonSolver = new PoissonSolverFFT2D(xSize, ySize, -xy_ratio, xy_ratio, -1.0, 1.0);
	rhoGrid = new Grid2D(xSize, ySize);
	phiGrid = new Grid2D(xSize, ySize);
	zGrid = new Grid1D(zSize);	
	bunchExtremaCalc = new BunchExtremaCalculator();
}

SpaceChargeCalc2p5D::SpaceChargeCalc2p5D(int xSize, int ySize, int zSize): CppPyWrapper(NULL)
{
	xy_ratio = 1.0;
	poissonSolver = new PoissonSolverFFT2D(xSize, ySize, -xy_ratio, xy_ratio, -1.0, 1.0);
	rhoGrid = new Grid2D(xSize, ySize);
	phiGrid = new Grid2D(xSize, ySize);
	zGrid = new Grid1D(zSize);
	bunchExtremaCalc = new BunchExtremaCalculator();		
}

SpaceChargeCalc2p5D::~SpaceChargeCalc2p5D(){
	delete poissonSolver;
	if(rhoGrid->getPyWrapper() != NULL){
		Py_DECREF(rhoGrid->getPyWrapper());
	} else {
		delete rhoGrid;
	}
	if(phiGrid->getPyWrapper() != NULL){
		Py_DECREF(phiGrid->getPyWrapper());
	} else {
		delete phiGrid;
	}	
	if(zGrid->getPyWrapper() != NULL){
		Py_DECREF(zGrid->getPyWrapper());
	} else {
		delete zGrid;
	}
	delete bunchExtremaCalc;
}

Grid2D* SpaceChargeCalc2p5D::getRhoGrid(){
	return rhoGrid;
}

Grid2D* SpaceChargeCalc2p5D::getPhiGrid(){
	return phiGrid;
}

Grid1D* SpaceChargeCalc2p5D::getLongGrid(){
	return zGrid;
}

void SpaceChargeCalc2p5D::trackBunch(Bunch* bunch, double length, BaseBoundary2D* boundary){

	int nPartsGlobal = bunch->getSizeGlobal();
	if(nPartsGlobal < 2) return;	
	
	double totalMacrosize = 0.;
	this->bunchAnalysis(bunch, totalMacrosize, boundary);
	//std::cerr<<"totalMacrosize="<<totalMacrosize;
	double z_step = zGrid->getStepZ();
	
	//calculate phiGrid
	poissonSolver->findPotential(rhoGrid,phiGrid);
	
	if(boundary != NULL){        
		//update potential with boundary condition		
		boundary->addBoundaryPotential(rhoGrid,phiGrid);
		//std::cerr<<"Boundary ADDED."<<std::endl;		
	}
	
	SyncPart* syncPart = bunch->getSyncPart();	
	double factor = 2*length*bunch->getClassicalRadius()*pow(bunch->getCharge(),2)/(pow(syncPart->getBeta(),2)*pow(syncPart->getGamma(),3));	
	//std::cout<<" debug totalMacrosize="<<totalMacrosize<<" factor="<<factor<<" z_step="<< z_step <<std::endl;	
	
	factor = factor/(z_step*totalMacrosize);	

	double Lfactor = 0.;
	double x,y,z,ex,ey,ez;	
		
	for (int i = 0, n = bunch->getSize(); i < n; i++){
		x = bunch->x(i);
		y = bunch->y(i);
		z = bunch->z(i);	
		
		if(boundary != NULL && boundary->isInside(x,y) == BaseBoundary2D::IS_INSIDE){
			phiGrid->calcGradient(x,y,ex,ey);	
			//std::cout<<" debug ip="<<i<<" x="<<x<<" y="<<y<<" z="<<z<<" ex="<<ex<<" ey="<<ey<<" ez="<<ez<<" rho_z="<< zGrid->getValue(z) <<std::endl;		
			Lfactor = - zGrid->getValue(z) * factor;
			//std::cerr<<" debug zgrid="<<zGrid->getValue(z)<<" lfactor="<<Lfactor;		
			bunch->xp(i) += ex * Lfactor;
			bunch->yp(i) += ey * Lfactor;	
			//std::cerr<<" xp="<<bunch->xp(i)<<" yp="<<bunch->yp(i);
		}
	}
}

void SpaceChargeCalc2p5D::bunchAnalysis(Bunch* bunch, double& totalMacrosize, BaseBoundary2D* boundary){

	double xMin, xMax, yMin, yMax, zMin, zMax;
	
    if(boundary == NULL){
		
	bunchExtremaCalc->getExtremaXYZ(bunch, xMin, xMax, yMin, yMax, zMin, zMax);
	
	//check if the beam size is not zero 
	if( xMin >=  xMax || yMin >=  yMax || zMin >=  zMax){
		int rank = 0;
		ORBIT_MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		if(rank == 0){
			std::cerr << "SpaceChargeCalc2p5D::bunchAnalysis(bunch,...)" << std::endl
         				<< "The bunch min and max sizes are wrong! Cannot calculate space charge!" << std::endl
								<< "x min ="<< xMin <<" max="<< xMax << std::endl
								<< "y min ="<< yMin <<" max="<< yMax << std::endl
								<< "z min ="<< zMin <<" max="<< zMax << std::endl
								<< "Stop."<< std::endl;
		}
		ORBIT_MPI_Finalize();
	}
	
	double xy_ratio_beam = (xMax - xMin)/(yMax - yMin);
	double width, center;
	if(xy_ratio_beam > xy_ratio){
		center = (yMax + yMin)/2.0;
		width = ((yMax - yMin)*(xy_ratio_beam/xy_ratio))/2.0;
		yMin = center - width;
		yMax = center + width;
	} else {
		center = (xMax + xMin)/2.0;
		width = ((xMax - xMin)/(xy_ratio_beam/xy_ratio))/2.0;		
		xMin = center - width;
		xMax = center + width;
	}
    }
    else{   
	xMax = boundary->getMaxX();
	xMin = boundary->getMinX();
	yMax = boundary->getMaxY();
	yMin = boundary->getMinY();
	
	bunchExtremaCalc->getExtremaZ(bunch, zMin, zMax);
	
	//check if the beam size is not zero 
	if(zMin >=  zMax){
		int rank = 0;
		ORBIT_MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		if(rank == 0){
			std::cerr << "SpaceChargeCalc2p5D::bunchAnalysis(bunch,...)" << std::endl
         				<< "The bunch min and max sizes are wrong! Cannot calculate space charge!" << std::endl
								<< "z min ="<< zMin <<" max="<< zMax << std::endl
								<< "Stop."<< std::endl;
		}
		ORBIT_MPI_Finalize();
	}
    }
    
	//set Grids' limits
	rhoGrid->setGridX(xMin,xMax);
	rhoGrid->setGridY(yMin,yMax);	
	
	phiGrid->setGridX(xMin,xMax);
	phiGrid->setGridY(yMin,yMax);
	
	zGrid->setGridZ(zMin,zMax);	
	
	poissonSolver->setGridX(xMin,xMax);
        poissonSolver->setGridY(yMin,yMax);

	//sizes of the grids are set up
	//bin rho&z Bunch to the Grid
	rhoGrid->setZero();
	zGrid->setZero();

	rhoGrid->binBunch(bunch);
	zGrid->binBunch(bunch);

	rhoGrid->synchronizeMPI(bunch->getMPI_Comm_Local());
	zGrid->synchronizeMPI(bunch->getMPI_Comm_Local());
	
	totalMacrosize = 0.;	
	int nZ = zGrid->getSizeZ();	
	for(int iz = 0; iz < nZ; iz++){
		totalMacrosize += zGrid->getValueOnGrid(iz);		
	}
}
