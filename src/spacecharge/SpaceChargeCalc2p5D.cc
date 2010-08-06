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
	poissonSolver = new PoissonSolverFFT2D(xSize, ySize, -1, 1, -1/xy_ratio, 1/xy_ratio);
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
SpaceChargeCalc2p5D::SpaceChargeCalc2p5D(int xSize, int ySize, int zSize, double xy_ratio,
	             double xMin, double xMax,
	             double yMin, double yMax,
	             double zMin, double zMax): CppPyWrapper(NULL)
{
	poissonSolver = new PoissonSolverFFT2D(xSize, ySize, xMin, xMax, xMin/xy_ratio, xMax/xy_ratio);
	rhoGrid = new Grid2D(xSize, ySize, xMin, xMax, yMin, yMax);
	phiGrid = new Grid2D(xSize, ySize, xMin, xMax, yMin, yMax);
	zGrid = new Grid1D(zSize, zMin, zMax);

	xSize_ = xSize;
	ySize_ = ySize;
	zSize_ = zSize;
	xMin_ = xMin;
	xMax_ = xMax;
	yMin_ = yMin; 
	yMax_ = yMax;
	zMin_ = zMin;
	zMax_ = zMax;
}

// Destructor
SpaceChargeCalc2p5D::~SpaceChargeCalc2p5D(){
	delete poissonSolver;
	delete rhoGrid, phiGrid;
	delete zGrid;
}
void SpaceChargeCalc2p5D::init(){
 
}

void SpaceChargeCalc2p5D::trackBunch(Bunch* bunch, BaseBoundary2D* boundary, double length){
	double x,y,z;
	double fx, fy, dx_, dy_, ex, ey;
	double Factor;
	
	if (boundary->getShapeType() > 0){
		xMax_ = boundary->getMaxX();
		xMin_ = boundary->getMinX();
		yMax_ = boundary->getMaxY();
		yMin_ = boundary->getMinY();
	} else{
	}
	zMax_ = zGrid->getMaxZ();
	zMin_ = zGrid->getMinZ();
	
	rhoGrid->setGridX(xMin_,xMax_);
	rhoGrid->setGridY(yMin_,yMax_);
	phiGrid->setGridX(xMin_,xMax_);
	phiGrid->setGridY(yMin_,yMax_);	
	zGrid->setGridZ(zMin_,zMax_);
	
	rhoGrid->binBunch(bunch);
	phiGrid->binBunch(bunch);
	zGrid->binBunch(bunch);
	
	poissonSolver->findPotential(rhoGrid,phiGrid);
	boundary->addBoundaryPotential(rhoGrid,phiGrid);
	
	dx_ = phiGrid->getStepX();
	dy_ = phiGrid->getStepY();
	
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
