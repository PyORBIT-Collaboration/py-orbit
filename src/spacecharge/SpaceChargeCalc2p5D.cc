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

#include "orbit_mpi.hh"

#include "Grid1D.hh"
#include "Grid2D.hh"
#include "PoissonSolverFFT2D.hh"
#include "SpaceChargeCalc2p5D.hh"
#include "ParticleMacroSize.hh"
#include "ParticleMacroSize.hh"

#include "BaseBoundary2D.hh"

#include <iostream>
#include <cmath>
#include <cfloat>

using namespace OrbitUtils;

// Constructor
SpaceChargeCalc2p5D::SpaceChargeCalc2p5D(int xSize, int ySize, int zSize, double xy_ratio): CppPyWrapper(NULL)
{
	double one_xy_ratio = 1.0/xy_ratio;
	poissonSolver = new PoissonSolverFFT2D(xSize, ySize, -1, 1, -one_xy_ratio, one_xy_ratio);
	rhoGrid = new Grid2D(xSize, ySize);
	phiGrid = new Grid2D(xSize, ySize);
	zGrid = new Grid1D(zSize);

	xSize_ = xSize;
	ySize_ = ySize;
	zSize_ = zSize;
	xMin_ = -1.0;
	xMax_ = +1.0;
	yMin_ = -1.0; 
	yMax_ = +1.0;
	zMin_ = -1.0;
	zMax_ = +1.0;
}
SpaceChargeCalc2p5D::SpaceChargeCalc2p5D(int xSize, int ySize, int zSize,
	             double xMin, double xMax,
	             double yMin, double yMax): CppPyWrapper(NULL)
{
	double xy_ratio = (xMax-xMin)/(yMax-yMin);
	poissonSolver = new PoissonSolverFFT2D(xSize, ySize, xMin, xMax, yMin, yMax);
	rhoGrid = new Grid2D(xSize, ySize, xMin, xMax, yMin, yMax);
	phiGrid = new Grid2D(xSize, ySize, xMin, xMax, yMin, yMax);
	zGrid = new Grid1D(zSize);

	xSize_ = xSize;
	ySize_ = ySize;
	zSize_ = zSize;
	xMin_ = xMin;
	xMax_ = xMax;
	yMin_ = yMin; 
	yMax_ = yMax;	
	zMin_ = -1.0;
	zMax_ = +1.0;
}

// Destructor
SpaceChargeCalc2p5D::~SpaceChargeCalc2p5D(){
	delete poissonSolver;
	delete rhoGrid, phiGrid;
	delete zGrid;
}
/*void SpaceChargeCalc2p5D::init(){
 
}
*/

void SpaceChargeCalc2p5D::trackBunch(Bunch* bunch, BaseBoundary2D* boundary, double length){
	double x,y,z;
	double fx, fy, dx_, dy_, ex, ey;
	double _lambda, factor, Lfactor;
	SyncPart* syncPart = bunch->getSyncPart();
	
	//getBoundary
	if (boundary == NULL){
		getBoundaryXY(bunch);
	} else{
		xMax_ = boundary->getMaxX();
		xMin_ = boundary->getMinX();
		yMax_ = boundary->getMaxY();
		yMin_ = boundary->getMinY();
	}
	getBoundaryZ(bunch);
	
	//setGrid
	rhoGrid->setGridX(xMin_,xMax_);
	rhoGrid->setGridY(yMin_,yMax_);
	phiGrid->setGridX(xMin_,xMax_);
	phiGrid->setGridY(yMin_,yMax_);	
	zGrid->setGridZ(zMin_,zMax_);
	
	//bin rho&z Grid to the bunch
	rhoGrid->binBunch(bunch);
	zGrid->binBunch(bunch);
	
	//calculate phiGrid
	poissonSolver->findPotential(rhoGrid,phiGrid);
	//calculate xp,yp
	_lambda = bunch->getSizeGlobal() / (zMax_ - zMin_);
	
	factor = 2. * _lambda * bunch->getClassicalRadius() / (bunch->getCharge() * pow(syncPart->getBeta(),2) * pow(syncPart->getGamma(),3));        
                
	if (boundary != NULL){ 
		//update potential with boundary condition
		boundary->addBoundaryPotential(rhoGrid,phiGrid);	
	}
	
	//calculate stepX&stepY
	dx_ = phiGrid->getStepX();
	dy_ = phiGrid->getStepY();
	
	for (int i = 0, n = bunch->getSize(); i < n; i++){
		x = bunch->x(i);
		y = bunch->y(i);
		z = bunch->z(i);
		
		phiGrid->calcGradient(x,y,ex,ey);		
		Lfactor = zGrid->getValue(z) * factor;
		
		bunch->xp(i) += ex * dx_ * length * Lfactor;
		bunch->yp(i) += ey * dy_ * length * Lfactor;
		std::cerr<<"xp="<<bunch->xp(i)<<"\n";
		std::cerr<<"yp="<<bunch->yp(i)<<"\n";
	}
}
void SpaceChargeCalc2p5D::getBoundaryXY(Bunch* bunch){
	double xMax_global,yMax_global,xMin_global,yMin_global;
	xMax_ = -DBL_MAX;
	yMax_ = -DBL_MAX;
	xMin_ =  DBL_MAX;
	yMin_ =  DBL_MAX;

	for (int i = 1, n = bunch->getSize(); i < n; i++){
		if(bunch->coordArr()[i][0] > xMax_) xMax_ = bunch->coordArr()[i][0];
		if(bunch->coordArr()[i][2] > yMax_) yMax_ = bunch->coordArr()[i][2];
		if(bunch->coordArr()[i][0] < xMin_) xMin_ = bunch->coordArr()[i][0];
		if(bunch->coordArr()[i][2] < yMin_) yMin_ = bunch->coordArr()[i][2];
	}
	pyORBIT_MPI_Comm* pyComm = bunch->getMPI_Comm_Local();
	if(pyComm == NULL) {
		ORBIT_MPI_Allreduce(&xMax_,&xMax_global,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
		ORBIT_MPI_Allreduce(&yMax_,&yMax_global,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
		ORBIT_MPI_Allreduce(&xMin_,&xMin_global,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
		ORBIT_MPI_Allreduce(&yMin_,&yMin_global,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
		
	} else {
		ORBIT_MPI_Allreduce(&xMax_,&xMax_global,1,MPI_DOUBLE,MPI_MAX,pyComm->comm);
		ORBIT_MPI_Allreduce(&yMax_,&yMax_global,1,MPI_DOUBLE,MPI_MAX,pyComm->comm);
		ORBIT_MPI_Allreduce(&xMin_,&xMin_global,1,MPI_DOUBLE,MPI_MIN,pyComm->comm);
		ORBIT_MPI_Allreduce(&yMin_,&yMin_global,1,MPI_DOUBLE,MPI_MIN,pyComm->comm);
	}	
	xMax_ = xMax_global;
	yMax_ = yMax_global;
	xMin_ = xMin_global;
	yMin_ = yMin_global;
	//clac boundary by xy_ratio
	if((xMax_ - xMin_) > (yMax_ - yMin_)*xy_ratio){
		xMax_ = (xMax_ + xMin_) / 2;
		xMin_ = (3 * xMin_ - xMax_) / 2;
		yMax_ = xMax_ / xy_ratio - (yMax_ - yMin_) / 2;
		yMin_ = xMin_ / xy_ratio - (yMax_ - yMin_) / 2;
	} else{
		xMax_ = yMax_ * xy_ratio - (xMax_ - xMin_) / 2;
		xMin_ = yMin_ * xy_ratio - (xMax_ - xMin_) / 2;
		yMax_ = (yMax_ + yMin_) / 2;
		yMin_ = (3 * yMin_ - yMax_) / 2;
	}
}
void SpaceChargeCalc2p5D::getBoundaryZ(Bunch* bunch){
	double zMax_global,zMin_global;
	zMax_ = -DBL_MAX;
	zMin_ =  DBL_MAX;	
	
	for (int i = 1, n = bunch->getSize(); i < n; i++){
		if(bunch->coordArr()[i][4] > zMax_) zMax_ = bunch->coordArr()[i][4];
		if(bunch->coordArr()[i][4] < zMin_) zMin_ = bunch->coordArr()[i][4];
	}
	pyORBIT_MPI_Comm* pyComm = bunch->getMPI_Comm_Local();
	if(pyComm == NULL) {
		ORBIT_MPI_Allreduce(&zMax_,&zMax_global,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
		ORBIT_MPI_Allreduce(&zMin_,&zMin_global,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
	} else {
		ORBIT_MPI_Allreduce(&zMax_,&zMax_global,1,MPI_DOUBLE,MPI_MAX,pyComm->comm);
		ORBIT_MPI_Allreduce(&zMin_,&zMin_global,1,MPI_DOUBLE,MPI_MIN,pyComm->comm);
	}
	zMax_ = zMax_global;
	zMin_ = zMin_global;
}
