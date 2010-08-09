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

#include "BaseBoundary2D.hh"

#include <iostream>

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
	double Factor;
	
	if (boundary == NULL){
		getBoundaryXY(bunch);
	} else{
		xMax_ = boundary->getMaxX();
		xMin_ = boundary->getMinX();
		yMax_ = boundary->getMaxY();
		yMin_ = boundary->getMinY();
	}
	getBoundaryZ(bunch);
	
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
	
	if (boundary != NULL){ 
		//update potential with boundary condition
		boundary->addBoundaryPotential(rhoGrid,phiGrid);	
	}
	
	//calculate stepX&stepY
	dx_ = phiGrid->getStepX();
	dy_ = phiGrid->getStepY();
	
	//formular **is used here
	for (int i = 0, n = bunch->getSize(); i < n; i++){
		x = bunch->x(i);
		y = bunch->y(i);
		z = bunch->z(i);
		
		phiGrid->calcGradient(x,y,ex,ey);		
		Factor = zGrid->getValue(z);
		
		bunch->xp(i) += ex * dx_ * length * Factor;
		bunch->yp(i) += ey * dy_ * length * Factor;
		std::cerr<<"xp="<<bunch->xp(i)<<"\n";
		std::cerr<<"yp="<<bunch->yp(i)<<"\n";
	}
}
void SpaceChargeCalc2p5D::getBoundaryXY(Bunch* bunch){
	xMax_ = bunch->x(1);
	yMax_ = bunch->y(1);
	xMin_ = bunch->x(1);
	yMin_ = bunch->y(1);

	for (int i = 2, n = bunch->getSize(); i < n; i++){
		if( bunch->x(i) > xMax_) xMax_ = bunch->x(i);
		if( bunch->y(i) > yMax_) yMax_ = bunch->y(i);
		if( bunch->x(i) < xMin_) xMin_ = bunch->x(i);
		if( bunch->y(i) < yMin_) yMin_ = bunch->y(i);
	}
	//clac boundary by xy_ratio 
}
void SpaceChargeCalc2p5D::getBoundaryZ(Bunch* bunch){
	zMin_ = bunch->z(1);
	zMax_ = bunch->z(1);
	for (int i = 2, n = bunch->getSize(); i < n; i++){
		if( bunch->z(i) > zMax_) zMax_ = bunch->z(i);
		if( bunch->z(i) < zMin_) zMin_ = bunch->z(i);
	}
}
