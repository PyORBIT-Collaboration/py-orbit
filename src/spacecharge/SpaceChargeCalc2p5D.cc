/////////////////////////////////////////////////////////////////////////////
//
// FILE NAME
//   SpaceChargeCalc2p5D.cc
//
//   06/28/10
//
// DESCRIPTION
//   Calculate the 2.5D space charge 
//
/////////////////////////////////////////////////////////////////////////////

#include "orbit_mpi.hh"

#include "Grid1D.hh"
#include "Grid2D.hh"
#include "BaseBoundary2D.hh"
#include "SpaceChargeCalc2p5D.hh"
#include "ParticleMacroSize.hh"

#include <iostream>

using namespace OrbitUtils;

// Constructor
SpaceChargeCalc2p5D::SpaceChargeCalc2p5D(double lkick): CppPyWrapper(NULL)
//SpaceChargeCalc2p5D::SpaceChargeCalc2p5D(int xSize, int ySize, int zSize): CppPyWrapper(NULL)
{
	/*xSize_ = xSize;
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
	             double yMin, double yMax,
	             double zMin, double zMax):: CppPyWrapper(NULL)
{
	xSize_ = xSize;
	ySize_ = ySize;
	zSize_ = zSize;
	xMin_ = xMin; 
	xMax_ = xMax; 
	yMin_ = yMin; 
	yMax_ = yMax;
	zMin_ = zMin; 
	zMax_ = zMax;*/
}

// Destructor
SpaceChargeCalc2p5D::~SpaceChargeCalc2p5D()
{	
}
void SpaceChargeCalc2p5D::init(){
//MPI
	int rank = 0;
	ORBIT_MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	ORBIT_MPI_Finalize();  
}
/*void SpaceChargeCalc2p5D::calcPotential(Bunch* bunch, Grid2D* phiGrid, Grid1D* zGrid, BaseBoundary* boundary)
{
	BaseBoundary* boundary = NULL;
	if (boundary->getShapeType() > 0)
	
	
}*/
void SpaceChargeCalc2p5D::calcPotential(Bunch* bunch, Grid2D* phiGrid, Grid1D* zGrid)
{
	double x,y,z;
	double fx, fy, dx_, dy_, ex, ey;
	double Factor;
	
	for (int i = 0, n = bunch->getSize(); i < n; i++){
		x = bunch->x(i);
		y = bunch->y(i);
		z = bunch->z(i);
		
		dx_ = phiGrid->getStepX();
		dy_ = phiGrid->getStepY();
		phiGrid->calcGradient(x,y,ex,ey);
		fx = ex * dx_;
		fy = ey * dy_;		
		
		zGrid->binBunch(bunch);
		Factor = zGrid->getValue(z);
		
		bunch->xp(i) += fx * Factor;
		bunch->yp(i) += fy * Factor;
		//bunch->xp(i) += fx * _lkick * Factor;
		//bunch->yp(i) += fy * _lkick * Factor;
	}
}
