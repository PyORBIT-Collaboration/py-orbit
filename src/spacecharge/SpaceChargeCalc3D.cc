/////////////////////////////////////////////////////////////////////////////
//
// FILE NAME
//   SpaceChargeCalc3D.cc
//
// AUTHOR
//    A. Shishlo
//
// Created:
//   12/03/10
//
// DESCRIPTION
// This class calculates the space charge kicks for bunch using 3D Poisson Solver.
// The solver implements FFT convolution algorithm to calculate the potential in
// the coordinate system where the bunch is resting. The solver is not parallel in the 
// sense of efficiency, but it is working correctly, and particles can are distributed
// among CPUs (not by the solver).
//
/////////////////////////////////////////////////////////////////////////////

#include "Grid3D.hh"
#include "PoissonSolverFFT3D.hh"
#include "SpaceChargeCalc3D.hh"
#include "BufferStore.hh"
#include "OrbitConst.hh"

#include <iostream>
#include <cmath>
#include <cfloat>

using namespace OrbitUtils;

SpaceChargeCalc3D::SpaceChargeCalc3D(int xSize, int ySize, int zSize): CppPyWrapper(NULL)
{
	xy_ratio = 1.0;
	xz_ratio = 1.0;
	poissonSolver = new PoissonSolverFFT3D(xSize, ySize, zSize, -xy_ratio, xy_ratio, -1./xy_ratio, 1./xy_ratio, -1./xz_ratio, 1./xz_ratio);
	rhoGrid = new Grid3D(xSize, ySize, zSize);
	rhoGrid->setGridX(poissonSolver->getMinX(),poissonSolver->getMaxX());
	rhoGrid->setGridY(poissonSolver->getMinY(),poissonSolver->getMaxY());
	rhoGrid->setGridZ(poissonSolver->getMinZ(),poissonSolver->getMaxZ());
	phiGrid = new Grid3D(xSize, ySize, zSize);
	phiGrid->setGridX(poissonSolver->getMinX(),poissonSolver->getMaxX());
	phiGrid->setGridY(poissonSolver->getMinY(),poissonSolver->getMaxY());
	phiGrid->setGridZ(poissonSolver->getMinZ(),poissonSolver->getMaxZ());	
	bunchExtremaCalc = new BunchExtremaCalculator();
	
	//------------ratio change limit ----------
	//If the shape (x to y and x to z ratios) of 3D region changes more than this
	//limit then we have to change shape and recalculate Green Functions in the Poisson solver
	ratio_limit = 1.2;
	
	//Number of bunches from both sides that should be taken into account for space charge.
	nBunches_ = 0;
	
	// The frequency of the bunch arrivals in Hz. It defines by the RFQ frequency.
	// The non-zero is setup by default to avoid division on zero
	frequency_ = 402.5e+6;
}

SpaceChargeCalc3D::~SpaceChargeCalc3D(){
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
	delete bunchExtremaCalc;
}

Grid3D* SpaceChargeCalc3D::getRhoGrid(){
	return rhoGrid;
}

Grid3D* SpaceChargeCalc3D::getPhiGrid(){
	return phiGrid;
}

void SpaceChargeCalc3D::setNumberOfExternalBunches(int nBunches){
	poissonSolver->setNumberOfExternalBunches(nBunches);
}	

void SpaceChargeCalc3D::setFrequencyOfBunches(double frequency){
	frequency_ = frequency;
}	

int SpaceChargeCalc3D::getNumberOfExternalBunches(){
	return poissonSolver->getNumberOfExternalBunches();
}	

double SpaceChargeCalc3D::getFrequencyOfBunches(){
	return frequency_;
}	

void SpaceChargeCalc3D::trackBunch(Bunch* bunch, double length){

	int nPartsGlobal = bunch->getSizeGlobal();
	if(nPartsGlobal < 2) return;
	
	//calculate max and min of X,Y,Z, bin paricles and set up limits for rhoGrid, phiGrid
	//and bin the particles from the bunch into the rhoGrid
	int nExtBunches = this->getNumberOfExternalBunches();
	if(nExtBunches > 0){
		this->wrappedBunchAnalysis(bunch);
	}
	else{
		this->bunchAnalysis(bunch);
	}
	
	//calculate phiGrid with potential. The z-coordinate is in the center of mass coordinate system
	poissonSolver->findPotential(rhoGrid,phiGrid);
	
	SyncPart* syncPart = bunch->getSyncPart();	
	double gamma = syncPart->getGamma();
	double beta = syncPart->getBeta();
	
	//distance between bunches in the center of mass of the bunch
	double lambda_cm = OrbitConst::c*beta*gamma/frequency_;
	
	double trans_factor =  length*bunch->getClassicalRadius()/(pow(beta,2)*pow(gamma,2));	
	double long_factor =  length*bunch->getClassicalRadius()*bunch->getMass();
	
	double x,y,z,ex,ey,ez;

	double z_center = (phiGrid->getMaxZ() + phiGrid->getMinZ())/2.0;
	
	for (int i = 0, n = bunch->getSize(); i < n; i++){
		x = bunch->x(i);
		y = bunch->y(i);
		z = (bunch->z(i) - z_center)*gamma + z_center;
		
		if(nExtBunches > 0){
			z = remainder(z,lambda_cm);
		}
		
		phiGrid->calcGradient(x,ex,y,ey,z,ez);	
		//std::cout<<"debug ip="<<i<<" x="<<x<<" y="<<y<<" z="<<z<<" ex="<<ex<<" ey="<<ey<<" ez="<<ez<<" rho_z="<< zGrid->getValue(z) <<std::endl;
		//calculate momentum kicks
		bunch->xp(i) += -ex * trans_factor;
		bunch->yp(i) += -ey * trans_factor;
		bunch->dE(i) += -ez * long_factor;
	}
}

void SpaceChargeCalc3D::bunchAnalysis(Bunch* bunch){
	
	double width, center;
	
	double xMin, xMax, yMin, yMax, zMin, zMax;
	
	bunchExtremaCalc->getExtremaXYZ(bunch, xMin, xMax, yMin, yMax, zMin, zMax);
	
	//we are not going to account for neighboring bunches here
	rhoGrid->setLongWrapping(0);
	phiGrid->setLongWrapping(0);
	
	//we will work in the center of the mass of the bunch. 
  //We have to shrink the longitudinal size
	SyncPart* syncPart = bunch->getSyncPart();	
	double gamma = syncPart->getGamma();
	center = (zMax + zMin)/2.0;
	width = (zMax - zMin)*gamma/2.0;		
	zMin = center - width;
	zMax = center + width;	
	
	//check if the beam size is not zero 
  if( xMin >=  xMax || yMin >=  yMax || zMin >=  zMax){
		int rank = 0;
		ORBIT_MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		if(rank == 0){
			std::cerr << "SpaceChargeCalc3D::bunchAnalysis(bunch)" << std::endl
         				<< "The bunch min and max sizes are wrong! Cannot calculate space charge!" << std::endl
								<< "x min ="<< xMin <<" max="<< xMax << std::endl
								<< "y min ="<< yMin <<" max="<< yMax << std::endl
								<< "z min ="<< zMin <<" max="<< zMax << std::endl
								<< "Stop."<< std::endl;
		}
		ORBIT_MPI_Finalize();
  }
	
	double xy_ratio_beam = (xMax - xMin)/(yMax - yMin);
	double xz_ratio_beam = (xMax - xMin)/(zMax - zMin);
	//first we have to decide if there is a need to change 
	//the ratio of x,y,z grid sizes in the Poisson Solver
	int changeRationInfo = -1;
	if(ratio_limit < xy_ratio_beam/xy_ratio || 1.0/ratio_limit > xy_ratio_beam/xy_ratio ||
		 ratio_limit < xz_ratio_beam/xz_ratio || 1.0/ratio_limit > xz_ratio_beam/xz_ratio){
		  changeRationInfo = 1;
	} 
	
	//The change in shape of the bunch is too big, and we have to change 
	//the shape of  the Poisson Solver grid. It will involve a recalculation of 
	//Green functions and will take some time. 
	if(changeRationInfo == 1){
		xy_ratio = xy_ratio_beam;
		xz_ratio = xz_ratio_beam;
		poissonSolver->setGridXYZ(-(xMax - xMin)/2.0,(xMax - xMin)/2.0,-(yMax - yMin)/2.0,(yMax - yMin)/2.0,-(zMax - zMin)/2.0,(zMax - zMin)/2.0);
	}
	
	//now we have to define the sizes of 3D grids for charge density and potential
	//The shape is defined by the Poisson Solver grid, and they should cover all particles
	double scale_coeff_x,scale_coeff_y,scale_coeff_z;
	scale_coeff_x = (xMax - xMin)/(poissonSolver->getMaxX() - poissonSolver->getMinX());
	scale_coeff_y = (yMax - yMin)/(poissonSolver->getMaxY() - poissonSolver->getMinY());	
	scale_coeff_z = (zMax - zMin)/(poissonSolver->getMaxZ() - poissonSolver->getMinZ());	
	double scale_coeff = scale_coeff_x;
	if(scale_coeff < scale_coeff_y) scale_coeff = scale_coeff_y;
	if(scale_coeff < scale_coeff_z) scale_coeff = scale_coeff_z;	
	
	center = (xMax + xMin)/2.0;
	width = (poissonSolver->getMaxX() - poissonSolver->getMinX())*scale_coeff/2.0;		
	xMin = center - width;
	xMax = center + width;
	rhoGrid->setGridX(xMin,xMax);	
	phiGrid->setGridX(xMin,xMax);	

	center = (yMax + yMin)/2.0;
	width = (poissonSolver->getMaxY() - poissonSolver->getMinY())*scale_coeff/2.0;		
	yMin = center - width;
	yMax = center + width;
	rhoGrid->setGridY(yMin,yMax);	
	phiGrid->setGridY(yMin,yMax);	
	
	//for binning we have to use real longitudinal coordinates
	center = (zMax + zMin)/2.0;
	width = (poissonSolver->getMaxZ() - poissonSolver->getMinZ())*scale_coeff/(2.0*gamma);		
	zMin = center - width;
	zMax = center + width;
	rhoGrid->setGridZ(zMin,zMax);	
	phiGrid->setGridZ(zMin,zMax);	
	
	//bin rho&z Bunch to the Grid
	rhoGrid->setZero();
	rhoGrid->binBunch(bunch);
	rhoGrid->synchronizeMPI(bunch->getMPI_Comm_Local());
	
	//after binning we have to move to the the center of mass of the bunch
	center = (zMax + zMin)/2.0;
	width = (zMax - zMin)*gamma/2.0;		
	zMin = center - width;
	zMax = center + width;		
	rhoGrid->setGridZ(zMin,zMax);	
	phiGrid->setGridZ(zMin,zMax);		
}

void SpaceChargeCalc3D::wrappedBunchAnalysis(Bunch* bunch){
	
	double width, center;
	
	double xMin, xMax, yMin, yMax, zMin, zMax;
	
	bunchExtremaCalc->getExtremaXYZ(bunch, xMin, xMax, yMin, yMax, zMin, zMax);
	
	//we will account for neighboring bunches here
	rhoGrid->setLongWrapping(1);
	phiGrid->setLongWrapping(1);
	
	//we will work in the center of the mass of the bunch. 
  //We have to shrink the longitudinal size
	SyncPart* syncPart = bunch->getSyncPart();	
	double gamma = syncPart->getGamma();
	double beta = syncPart->getBeta();
	center = (zMax + zMin)/2.0;
	width = (zMax - zMin)/2.0;		
	zMin = center - width;
	zMax = center + width;
	
	//distance between bunches in the laboratory system
	double lambda = OrbitConst::c*beta/frequency_;
	
	//distance between bunches in the center of mass of the bunch
	double lambda_cm = OrbitConst::c*beta*gamma/frequency_;	
	
	if(fabs(zMax) > lambda/2 || fabs(zMin) > lambda/2){
		// extension should not change the results at all. It just changes the grid size.
		double extension_coeff = 1.0001;
		center = 0.;
		width = extension_coeff*lambda/2;
		zMin = center - width;
		zMax = center + width;		
	}
	
	//check if the beam size is not zero 
  if( xMin >=  xMax || yMin >=  yMax || zMin >=  zMax){
		int rank = 0;
		ORBIT_MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		if(rank == 0){
			std::cerr << "SpaceChargeCalc3D::wrappedBunchAnalysis(bunch)" << std::endl
         				<< "The bunch min and max sizes are wrong! Cannot calculate space charge!" << std::endl
								<< "x min ="<< xMin <<" max="<< xMax << std::endl
								<< "y min ="<< yMin <<" max="<< yMax << std::endl
								<< "z min ="<< zMin <<" max="<< zMax << std::endl
								<< "Stop."<< std::endl;
		}
		ORBIT_MPI_Finalize();
  }
	
	rhoGrid->setGridX(xMin,xMax);	
	phiGrid->setGridX(xMin,xMax);	

	rhoGrid->setGridY(yMin,yMax);	
	phiGrid->setGridY(yMin,yMax);	
	
	rhoGrid->setGridZ(zMin,zMax);	
	phiGrid->setGridZ(zMin,zMax);	
	
	//bin rho&z Bunch to the Grid
	rhoGrid->setZero();
	rhoGrid->binBunch(bunch,lambda);
	rhoGrid->synchronizeMPI(bunch->getMPI_Comm_Local());
	
	//after binning we have to rescale the z-coordinate to the center-of-mass system
	center = (zMax + zMin)/2.0;
	width = (zMax - zMin)*gamma/2.0;		
	zMin = center - width;
	zMax = center + width;		
	rhoGrid->setGridZ(zMin,zMax);	
	phiGrid->setGridZ(zMin,zMax);
	
	//now we set the sizes of the Poisson Solver grids and 
	//calculate Green function FFT inside it
	poissonSolver->setSpacingOfExternalBunches(lambda_cm);
	poissonSolver->setGridXYZ(rhoGrid->getMinX(), rhoGrid->getMaxX(),
			                      rhoGrid->getMinY(), rhoGrid->getMaxY(),
			                      rhoGrid->getMinZ(), rhoGrid->getMaxZ());
}

/** Sets the ratio limit for the shape change and Green Function recalculations. */
void SpaceChargeCalc3D::setRatioLimit(double ratio_limit_in)
{
	ratio_limit = ratio_limit_in;
}

/** Returns the ratio limit for the shape change and Green Function recalculations. */
double SpaceChargeCalc3D::getRatioLimit()
{
	return ratio_limit;
}


