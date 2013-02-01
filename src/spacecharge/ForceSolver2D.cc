#include "ForceSolver2D.hh"

#include <iostream>

using namespace OrbitUtils;

// Constructor
ForceSolver2D::ForceSolver2D(int xSize, int ySize): CppPyWrapper(NULL)
{
	xSize_ = xSize;
	ySize_ = ySize;
	xMin_ = -1.0;
	xMax_ = +1.0;
	yMin_ = -1.0;
	yMax_ = +1.0;
	dx_ = (xMax_ - xMin_)/(xSize_ -1);
	dy_ = (yMax_ - yMin_)/(ySize_ -1);	
}

// Destructor
ForceSolver2D::~ForceSolver2D(){};

/** Returns the number of points in x direction */
int ForceSolver2D::getSizeX(){return xSize_;}; 

/** Returns the number of points in y direction */
int ForceSolver2D::getSizeY(){return ySize_;};

/** Returns the max x in the grid points */ 
double ForceSolver2D::getMaxX(){return xMax_;};

/** Returns the min x in the grid points */ 
double ForceSolver2D::getMinX(){return xMin_;};

/** Returns the max y in the grid points */ 
double ForceSolver2D::getMaxY(){return yMax_;};

/** Returns the min y in the grid points */ 
double ForceSolver2D::getMinY(){return yMin_;};

/** Returns the mesh step in x direction */
double ForceSolver2D::getStepX(){ return dx_;};

/** Returns the mesh step in y direction */
double ForceSolver2D::getStepY(){ return dy_;};


void ForceSolver2D::checkSizes(Grid2D* rhoGrid){
  if( xSize_ !=  rhoGrid->getSizeX() || ySize_ != rhoGrid->getSizeY() ||
		  dx_ != rhoGrid->getStepX() || dy_ != rhoGrid->getStepY() ||
			xMin_ != rhoGrid->getMinX() || yMin_ != rhoGrid->getMinY()){
		int rank = 0;
		ORBIT_MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		if(rank == 0){
			std::cerr << "ForceSolver2D (or subclass) checkSizes(...):" 
			<< "The grid sizes are different "<< std::endl 
								<< "number x bins ="<< xSize_ << std::endl
								<< "number y bins ="<< ySize_ << std::endl
								<< "rhoGrid x bins ="<< rhoGrid->getSizeX() <<std::endl
								<< "rhoGrid y bins ="<< rhoGrid->getSizeY() <<std::endl
								<< "dx_  ="<< dx_ <<std::endl
								<< "dy_  ="<< dy_ <<std::endl
								<< "rhoGrid dx ="<< rhoGrid->getStepX() <<std::endl
								<< "rhoGrid dy ="<< rhoGrid->getStepY() <<std::endl
								<< "xMin ="<< xMin_ <<std::endl
								<< "yMin  ="<< yMin_ <<std::endl
								<< "rhoGrid xMin ="<< rhoGrid->getMinX() <<std::endl
								<< "rhoGrid yMin ="<< rhoGrid->getMinY() <<std::endl
								<< "Stop. \n";
		}
		ORBIT_MPI_Finalize();
  }
}


