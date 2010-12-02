#include "PoissonSolver3D.hh"

#include <iostream>

using namespace OrbitUtils;

// Constructor
PoissonSolver3D::PoissonSolver3D(int xSize, int ySize, int zSize): CppPyWrapper(NULL)
{
	xSize_ = xSize;
	ySize_ = ySize;
	zSize_ = zSize;
	xMin_ = -1.0;
	xMax_ = +1.0;
	yMin_ = -1.0;
	yMax_ = +1.0;
	zMin_ = -1.0;
	zMax_ = +1.0;
	dx_ = (xMax_ - xMin_)/(xSize_ -1);
	dy_ = (yMax_ - yMin_)/(ySize_ -1);	
	dz_ = (zMax_ - zMin_)/(zSize_ -1);	
}

// Constructor
PoissonSolver3D::PoissonSolver3D(int xSize, int ySize, int zSize,
	                               double xMin, double xMax, 
																  double yMin, double yMax,
																  double zMin, double zMax): CppPyWrapper(NULL)
{
	xSize_ = xSize;
	ySize_ = ySize;
	zSize_ = zSize;	
	xMin_ = xMin;
	xMax_ = xMax;
	yMin_ = yMin;
	yMax_ = yMax;
	zMin_ = zMin;
	zMax_ = zMax;	
	dx_ = (xMax_ - xMin_)/(xSize_ -1);
	dy_ = (yMax_ - yMin_)/(ySize_ -1);	
	dz_ = (zMax_ - zMin_)/(zSize_ -1);		
}

// Destructor
PoissonSolver3D::~PoissonSolver3D(){};

/** Returns the number of points in x direction */
int PoissonSolver3D::getSizeX(){return xSize_;}; 

/** Returns the number of points in y direction */
int PoissonSolver3D::getSizeY(){return ySize_;};

/** Returns the number of points in z direction */
int PoissonSolver3D::getSizeZ(){return zSize_;};

/** Returns the max x in the grid points */ 
double PoissonSolver3D::getMaxX(){return xMax_;};

/** Returns the min x in the grid points */ 
double PoissonSolver3D::getMinX(){return xMin_;};

/** Returns the max y in the grid points */ 
double PoissonSolver3D::getMaxY(){return yMax_;};

/** Returns the min y in the grid points */ 
double PoissonSolver3D::getMinY(){return yMin_;};

/** Returns the max z in the grid points */ 
double PoissonSolver3D::getMaxZ(){return zMax_;};

/** Returns the min z in the grid points */ 
double PoissonSolver3D::getMinZ(){return zMin_;};

/** Returns the mesh step in x direction */
double PoissonSolver3D::getStepX(){ return dx_;};

/** Returns the mesh step in y direction */
double PoissonSolver3D::getStepY(){ return dy_;};

/** Returns the mesh step in z direction */
double PoissonSolver3D::getStepZ(){ return dz_;};

void PoissonSolver3D::checkSizes(Grid3D* rhoGrid,Grid3D*  phiGrid){
  if( xSize_ !=  rhoGrid->getSizeX() || ySize_ != rhoGrid->getSizeY() ||
		  xSize_ !=  phiGrid->getSizeX() || ySize_ != phiGrid->getSizeY() ||
		  zSize_ !=  phiGrid->getSizeZ() || zSize_ != rhoGrid->getSizeZ() ||
		  dx_ != rhoGrid->getStepX() || dy_ != rhoGrid->getStepY() ||
			dx_ != phiGrid->getStepX() || dy_ != phiGrid->getStepY() ||
			dz_ != phiGrid->getStepZ() || dz_ != rhoGrid->getStepZ() ||
			xMin_ != rhoGrid->getMinX() || yMin_ != rhoGrid->getMinY() ||
			xMin_ != phiGrid->getMinX() || yMin_ != phiGrid->getMinY() ||
			zMin_ != phiGrid->getMinZ() || zMin_ != rhoGrid->getMinZ()){
		int rank = 0;
		ORBIT_MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		if(rank == 0){
			std::cerr << "PoissonSolver3D (or subclass) checkSizes(...):" 
			<< "The grid sizes are different "<< std::endl 
								<< "number x bins ="<< xSize_ << std::endl
								<< "number y bins ="<< ySize_ << std::endl
								<< "rhoGrid x bins ="<< rhoGrid->getSizeX() <<std::endl
								<< "rhoGrid y bins ="<< rhoGrid->getSizeY() <<std::endl
								<< "rhoGrid z bins ="<< rhoGrid->getSizeZ() <<std::endl
								<< "phiGrid x bins ="<< phiGrid->getSizeX() <<std::endl
								<< "phiGrid y bins ="<< phiGrid->getSizeY() <<std::endl
								<< "phiGrid z bins ="<< phiGrid->getSizeZ() <<std::endl
								<< "dx_  ="<< dx_ <<std::endl
								<< "dy_  ="<< dy_ <<std::endl
								<< "dz_  ="<< dz_ <<std::endl
								<< "rhoGrid dx ="<< rhoGrid->getStepX() <<std::endl
								<< "rhoGrid dy ="<< rhoGrid->getStepY() <<std::endl
								<< "rhoGrid dz ="<< rhoGrid->getStepZ() <<std::endl
								<< "phiGrid dx ="<< phiGrid->getStepX() <<std::endl
								<< "phiGrid dy ="<< phiGrid->getStepY() <<std::endl
								<< "phiGrid dz ="<< phiGrid->getStepZ() <<std::endl
								<< "xMin  ="<< xMin_ <<std::endl
								<< "yMin  ="<< yMin_ <<std::endl
								<< "zMin  ="<< zMin_ <<std::endl
								<< "rhoGrid xMin ="<< rhoGrid->getMinX() <<std::endl
								<< "rhoGrid yMin ="<< rhoGrid->getMinY() <<std::endl
								<< "rhoGrid zMin ="<< rhoGrid->getMinZ() <<std::endl
								<< "phiGrid xMin ="<< phiGrid->getMinX() <<std::endl
								<< "phiGrid yMin ="<< phiGrid->getMinY() <<std::endl
								<< "phiGrid zMin ="<< phiGrid->getMinZ() <<std::endl
								<< "Stop. \n";
		}
		ORBIT_MPI_Finalize();
  }		
}


