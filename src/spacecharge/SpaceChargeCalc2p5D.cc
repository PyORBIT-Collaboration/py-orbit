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
#include "ParticleMacroSize.hh"

#include "BufferStore.hh"
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
	bunchExtremaCalc = new BunchExtremaCalculator();
	
	xSize_ = xSize;
	ySize_ = ySize;
	zSize_ = zSize;
	xMin_ = -1.0;
	xMax_ = +1.0;
	yMin_ = -one_xy_ratio; 
	yMax_ = one_xy_ratio;
	zMin_ = -1.0;
	zMax_ = +1.0;
	
	xy_ratio_ = xy_ratio;
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
	bunchExtremaCalc = new BunchExtremaCalculator();

	xSize_ = xSize;
	ySize_ = ySize;
	zSize_ = zSize;
	xMin_ = xMin;
	xMax_ = xMax;
	yMin_ = yMin; 
	yMax_ = yMax;	
	zMin_ = -1.0;
	zMax_ = +1.0;

	xy_ratio_ = xy_ratio;
}

// Destructor
SpaceChargeCalc2p5D::~SpaceChargeCalc2p5D(){
	delete poissonSolver;
	delete rhoGrid, phiGrid;
	delete zGrid;
	delete bunchExtremaCalc;
}
//void SpaceChargeCalc2p5D::init(){
 
//}

void SpaceChargeCalc2p5D::trackBunch(Bunch* bunch, double length){
	double x,y,z,_dz;
	double ex,ey;
	double factor, Lfactor;
	
	//getBoundary
	bunchExtremaCalc->getExtremaXYZ(bunch,xMin_,xMax_,yMin_,yMax_,zMin_,zMax_);
	
	bunchAnalysis(bunch);	
  	calcMomentumFactor(bunch,length,factor);
  	
  	//calculate phiGrid
	poissonSolver->findPotential(rhoGrid,phiGrid);
  	
  for (int i = 0, n = bunch->getSize(); i < n; i++){
		x = bunch->x(i);
		y = bunch->y(i);
		z = bunch->z(i);
		
		phiGrid->calcGradient(x,y,ex,ey);
		//std::cerr<<"START: x="<<x<<" y="<<y<<" z="<<z<<"\n";
		Lfactor = zGrid->getValue(z) * factor / bunch->getSizeGlobal();
		//std::cerr<<"Lfactor"<<Lfactor<<"\n";
		bunch->xp(i) += ex * Lfactor;
		bunch->yp(i) += ey * Lfactor;
		//std::cerr<<"xp="<<bunch->xp(i)<<" yp="<<bunch->yp(i)<<"\n";
	}

}
void SpaceChargeCalc2p5D::trackBunch(Bunch* bunch, double length, BaseBoundary2D* boundary){
	double x,y,z,_dz;
	double ex,ey;
	double factor,Lfactor;
	
	//getBoundary	
	bunchExtremaCalc->getExtremaXYZ(bunch,xMin_,xMax_,yMin_,yMax_,zMin_,zMax_);
	
	bunchAnalysis(bunch);
	calcMomentumFactor(bunch,length,factor);
	
	//calculate phiGrid
	poissonSolver->findPotential(rhoGrid,phiGrid);
	
	//update potential with boundary condition
	boundary->addBoundaryPotential(rhoGrid,phiGrid);	

	for (int i = 0, n = bunch->getSize(); i < n; i++){
		x = bunch->x(i);
		y = bunch->y(i);
		z = bunch->z(i);
		
		phiGrid->calcGradient(x,y,ex,ey);
		
		//std::cerr<<"START: x="<<x<<" y="<<y<<" z="<<z<<"\n";
		Lfactor = zGrid->getValue(z) * factor / bunch->getSizeGlobal();
		//std::cerr<<"Lfactor="<<Lfactor<<" ex="<<ex<<" ey="<<ey<<"\n";
		bunch->xp(i) += ex * Lfactor;
		bunch->yp(i) += ey * Lfactor;
		//std::cerr<<"xp="<<bunch->xp(i)<<" yp="<<bunch->yp(i)<<"\n";
	}
}

void SpaceChargeCalc2p5D::calcMomentumFactor(Bunch* bunch, double length, double& factor){
		
	double _dz, _lambda;
	SyncPart* syncPart = bunch->getSyncPart();	
	
	//xp={2rL(lambda)/beta^2*gamma^3}*protential/nSize
	_dz = zGrid->getStepZ();
	factor = zGrid->getValue(0);
	//std::cerr<<"dz="<<_dz<<" ratio="<<xy_ratio_<<" value="<<factor<<"\n";
	_lambda = bunch->getSizeGlobal() / _dz;	
	factor = 2. * length * _lambda * bunch->getClassicalRadius() / (bunch->getCharge() * pow(syncPart->getBeta(),2) * pow(syncPart->getGamma(),3));	
	//std::cerr<<"lambda = "<<_lambda<<" R = "<<bunch->getClassicalRadius()<<" charge = "<<bunch->getCharge()<<" beta^2 = "<<pow(syncPart->getBeta(),2)<<" gamma^3 = "<<pow(syncPart->getGamma(),3)<<"\n";        
	//std::cerr<<"calc factor. factor = "<<factor<<"\n";
}

void SpaceChargeCalc2p5D::bunchAnalysis(Bunch* bunch){	

	//check if the beam size is not zero 
  if( xMin_ >=  xMax_ || yMin_ >=  yMax_ || zMin_ >=  zMax_){
		int rank = 0;
		ORBIT_MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		if(rank == 0){
			std::cerr << "SpaceChargeCalc2p5Drb::bunchAnalysis(bunch,...) \n" 
         				<< "The bunch min and max sizes are wrong! Cannot calculate space charge! \n" 
								<< "x min ="<< xMin_ <<" max="<< xMax_ << std::endl
								<< "y min ="<< yMin_ <<" max="<< yMax_ << std::endl
								<< "z min ="<< zMin_ <<" max="<< zMax_ << std::endl
								<< "Stop."<< std::endl;
		}
		ORBIT_MPI_Finalize();
  }
  
	//clac boundary by xy_ratio
	double delta_ = (xMax_ - xMin_) - (yMax_ - yMin_)*xy_ratio_;;
	if(delta_ > 0){
		yMax_ = yMax_ + 0.5*delta_/xy_ratio_;
		yMin_ = yMin_ - 0.5*delta_/xy_ratio_;
	} else if(delta_ < 0){
		xMax_ = xMax_ + 0.5*delta_;
		xMin_ = xMin_ - 0.5*delta_;
		std::cerr<<" delta="<<delta_<<"\n";
	}
	
	//setGrid
	rhoGrid->setGridX(xMin_,xMax_);
	rhoGrid->setGridY(yMin_,yMax_);	
	
	phiGrid->setGridX(xMin_,xMax_);
	phiGrid->setGridY(yMin_,yMax_);	
	
	zGrid->setGridZ(zMin_,zMax_);
	
	poissonSolver->setGridX(xMin_,xMax_);
	poissonSolver->setGridY(yMin_,yMax_);
	
	std::cerr<<"checkGrid:"<<rhoGrid->getMaxX()<<":"<<rhoGrid->getMinX()<<":"<<rhoGrid->getMaxY()<<":"<<rhoGrid->getMinY()<<":"<<zGrid->getMaxZ()<<":"<<zGrid->getMinZ()<<"\n";
	std::cerr<<"get grid:"<<xMax_<<":"<<xMin_<<":"<<yMax_<<":"<<yMin_<<":"<<zMax_<<":"<<zMin_<<"\n";
	
	//bin rho&z Bunch to the Grid
	rhoGrid->binBunch(bunch);
	zGrid->binBunch(bunch);
}

