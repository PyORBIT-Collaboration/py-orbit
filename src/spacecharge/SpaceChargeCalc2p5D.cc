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
void SpaceChargeCalc2p5D::trackBunch(Bunch* bunch, double length){
	double x,y,z,_dz;
	double ex,ey;
	double factor, Lfactor;
	
	//getBoundary		
	getBoundaryXY(bunch);
	getBoundaryZ(bunch);
	
	factor = calcMomentumFactor(bunch,length);
	
	for (int i = 0, n = bunch->getSize(); i < n; i++){
		x = bunch->x(i);
		y = bunch->y(i);
		z = bunch->z(i);
		
		phiGrid->calcGradient(x,y,ex,ey);		
		//?
		Lfactor = zGrid->getValue(z) * factor / bunch->getSizeGlobal();
	
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
	xMax_ = boundary->getMaxX();
	xMin_ = boundary->getMinX();
	yMax_ = boundary->getMaxY();
	yMin_ = boundary->getMinY();		
	getBoundaryZ(bunch);
	
	factor = calcMomentumFactor(bunch,length);
	
	//update potential with boundary condition
	boundary->addBoundaryPotential(rhoGrid,phiGrid);	

	for (int i = 0, n = bunch->getSize(); i < n; i++){
		x = bunch->x(i);
		y = bunch->y(i);
		z = bunch->z(i);
		
		phiGrid->calcGradient(x,y,ex,ey);		
		//?
		Lfactor = zGrid->getValue(z) * factor / bunch->getSizeGlobal();
	
		bunch->xp(i) += ex * Lfactor;
		bunch->yp(i) += ey * Lfactor;
		//std::cerr<<"xp="<<bunch->xp(i)<<" yp="<<bunch->yp(i)<<"\n";
	}
}

double SpaceChargeCalc2p5D::calcMomentumFactor(Bunch* bunch, double length){
	double _dz, _lambda, factor;
	SyncPart* syncPart = bunch->getSyncPart();	
	
	//setGrid
	rhoGrid->setGridX(xMin_,xMax_);
	rhoGrid->setGridY(yMin_,yMax_);	
	phiGrid->setGridX(xMin_,xMax_);
	phiGrid->setGridY(yMin_,yMax_);	
	zGrid->setGridZ(zMin_,zMax_);	
	poissonSolver->setGridX(xMin_,xMax_);
	poissonSolver->setGridY(yMin_,yMax_);
	
	std::cerr<<"checkGrid:"<<rhoGrid->getMaxX()<<":"<<rhoGrid->getMinX()<<":"<<rhoGrid->getMaxY()<<":"<<rhoGrid->getMinY()<<":"<<zGrid->getMaxZ()<<":"<<zGrid->getMinZ()<<"\n";
	std::cerr<<"final grid:"<<xMax_<<":"<<xMin_<<":"<<yMax_<<":"<<yMin_<<":"<<zMax_<<":"<<zMin_<<"\n";
	
	//bin rho&z Bunch to the Grid
	rhoGrid->binBunch(bunch);
	zGrid->binBunch(bunch);
	
	//calculate phiGrid
	poissonSolver->findPotential(rhoGrid,phiGrid);

	//xp={2rL(lambda)/beta^2*gamma^3}*protential/nSize
	_dz = zGrid->getStepZ();
	_lambda = bunch->getSizeGlobal() / _dz;	
	factor = 2. * length * _lambda * bunch->getClassicalRadius() / (bunch->getCharge() * pow(syncPart->getBeta(),2) * pow(syncPart->getGamma(),3));	
	std::cerr<<"lambda = "<<_lambda<<" R = "<<bunch->getClassicalRadius()<<" charge = "<<bunch->getCharge()<<" beta^2 = "<<pow(syncPart->getBeta(),2)<<" gamma^3 = "<<pow(syncPart->getGamma(),3)<<"\n";        
	std::cerr<<"calc factor. factor = "<<factor<<"\n";	
	return factor;
}

void SpaceChargeCalc2p5D::getBoundaryXY(Bunch* bunch){	
	xMax_ = -DBL_MAX;
	yMax_ = -DBL_MAX;
	xMin_ =  DBL_MAX;
	yMin_ =  DBL_MAX;

	double** partArr=bunch->coordArr();
	for (int i = 0, n = bunch->getSize(); i < n; i++){
		if(partArr[i][0] > xMax_) xMax_ = partArr[i][0];
		if(partArr[i][2] > yMax_) yMax_ = partArr[i][2];
		if(partArr[i][0] < xMin_) xMin_ = partArr[i][0];
		if(partArr[i][2] < yMin_) yMin_ = partArr[i][2];
	}
	std::cerr<<"grid before mpi:"<<xMax_<<":"<<xMin_<<":"<<yMax_<<":"<<yMin_<<"\n";
	pyORBIT_MPI_Comm* pyComm = bunch->getMPI_Comm_Local();
	
	double* gridArr;
	double* gridArr_global;
	gridArr = new double[4];
	gridArr_global = new double[4];
	gridArr[0] = xMax_;
	gridArr[1] = yMax_;
	gridArr[2] = -xMin_;
	gridArr[3] = -yMin_;
	
	for(int i = 0; i < 4; i++){
		gridArr_global[i] = gridArr[i];
	}
	std::cerr<<"grid in arr:"<<gridArr[0]<<":"<<gridArr[2]<<":"<<gridArr[1]<<":"<<gridArr[3]<<"\n";
	
	ORBIT_MPI_Allreduce(gridArr,gridArr_global,4,MPI_DOUBLE,MPI_MAX,pyComm->comm);

	xMax_ = gridArr_global[0];
	yMax_ = gridArr_global[1];
	xMin_ = -gridArr_global[2];
	yMin_ = -gridArr_global[3];
	std::cerr<<"grid after mpi:"<<xMax_<<":"<<xMin_<<":"<<yMax_<<":"<<yMin_<<"\n";
	
	//clac boundary by xy_ratio
	double delta_ = (xMax_ - xMin_) > (yMax_ - yMin_)*xy_ratio;
	if(delta_ > 0){		 
		yMax_ = yMax_ + delta_;
		yMin_ = yMin_ - delta_;
	} else{		
		xMax_ = xMax_ - delta_;
		xMin_ = xMin_ + delta_;
	}
	std::cerr<<"final grid:"<<xMax_<<":"<<xMin_<<":"<<yMax_<<":"<<yMin_<<"\n";
}
void SpaceChargeCalc2p5D::getBoundaryZ(Bunch* bunch){
	double zMax_global,zMin_global;
	zMax_ = -DBL_MAX;
	zMin_ =  DBL_MAX;	
	
	double** partArr=bunch->coordArr();
	for (int i = 0, n = bunch->getSize(); i < n; i++){
		if(partArr[i][4] > zMax_) zMax_ = partArr[i][4];
		if(partArr[i][4] < zMin_) zMin_ = partArr[i][4];
	}
	pyORBIT_MPI_Comm* pyComm = bunch->getMPI_Comm_Local();

	ORBIT_MPI_Allreduce(&zMax_,&zMax_global,1,MPI_DOUBLE,MPI_MAX,pyComm->comm);
	ORBIT_MPI_Allreduce(&zMin_,&zMin_global,1,MPI_DOUBLE,MPI_MIN,pyComm->comm);

	zMax_ = zMax_global;
	zMin_ = zMin_global;	
}
