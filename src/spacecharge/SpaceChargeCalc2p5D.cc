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
	zDerivGrid = new Grid1D(zSize);
	bunchExtremaCalc = new BunchExtremaCalculator();
	//we will use 3 points by default to calculate the longitudinal density derivative
	n_long_avg = 3;
	S_arr = new double*[5];
	for(int i = 0; i < 5; i++){
		S_arr[i] = new double[2];
	}
}

SpaceChargeCalc2p5D::SpaceChargeCalc2p5D(int xSize, int ySize, int zSize): CppPyWrapper(NULL)
{
	xy_ratio = 1.0;
	poissonSolver = new PoissonSolverFFT2D(xSize, ySize, -xy_ratio, xy_ratio, -1.0, 1.0);
	rhoGrid = new Grid2D(xSize, ySize);
	phiGrid = new Grid2D(xSize, ySize);
	zGrid = new Grid1D(zSize);
	zDerivGrid = new Grid1D(zSize);
	bunchExtremaCalc = new BunchExtremaCalculator();	
	//we will use 3 points by default to calculate the longitudinal density derivative
	n_long_avg = 3;
	S_arr = new double*[5];
	for(int i = 0; i < 5; i++){
		S_arr[i] = new double[2];
	}	
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
	if(zDerivGrid->getPyWrapper() != NULL){
		Py_DECREF(zDerivGrid->getPyWrapper());
	} else {
		delete zDerivGrid;
	}	
	delete bunchExtremaCalc;
	for(int i = 0; i < 5; i++){
		delete [] S_arr[i];
	}	
	delete [] S_arr; 
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

Grid1D* SpaceChargeCalc2p5D::getLongDerivativeGrid(){
	return zDerivGrid;
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
	double factor =  2*length*bunch->getClassicalRadius()*pow(bunch->getCharge(),2)/(pow(syncPart->getBeta(),2)*pow(syncPart->getGamma(),3));	
	std::cout<<" debug totalMacrosize="<<totalMacrosize<<" factor="<<factor<<" z_step="<< z_step <<std::endl;	
	
	factor = factor/(z_step*totalMacrosize);	

	double Lfactor = 0.;
	double x,y,z,ex,ey,ez;	
		
	for (int i = 0, n = bunch->getSize(); i < n; i++){
		x = bunch->x(i);
		y = bunch->y(i);
		z = bunch->z(i);		
		
		phiGrid->calcGradient(x,y,ex,ey);		
		//std::cout<<" debug ip="<<i<<" x="<<x<<" y="<<y<<" z="<<z<<" ex="<<ex<<" ey="<<ey<<" ez="<<ez<<" rho_z="<< zGrid->getValue(z) <<std::endl;
		
		Lfactor = - zGrid->getValue(z) * factor;
		//std::cerr<<" debug zgrid="<<zGrid->getValue(z)<<" lfactor="<<Lfactor;
		
		bunch->xp(i) += ex * Lfactor;
		bunch->yp(i) += ey * Lfactor;	
		//std::cerr<<" xp="<<bunch->xp(i)<<" yp="<<bunch->yp(i);
	}
}

double SpaceChargeCalc2p5D::bunchAnalysis(Bunch* bunch, double& totalMacrosize, BaseBoundary2D* boundary){

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
	//clac boundary by xy_ratio
/*	double delta_ = (xMax_ - xMin_) - (yMax_ - yMin_)*xy_ratio_;;
	if(delta_ > 0){
		yMax_ = yMax_ + 0.5*delta_/xy_ratio_;
		yMin_ = yMin_ - 0.5*delta_/xy_ratio_;
	} else if(delta_ < 0){
		xMax_ = xMax_ - 0.5*delta_;
		xMin_ = xMin_ + 0.5*delta_;
		std::cerr<<" delta="<<delta_<<"\n";
	}
*/	
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
	zDerivGrid->setGridZ(zMin,zMax);
	
	poissonSolver->setGridX(xMin,xMax);
        poissonSolver->setGridY(yMin,yMax);

	//sizes of the grids are set up
	//bin rho&z Bunch to the Grid
	rhoGrid->setZero();
	zGrid->setZero();
	zDerivGrid->setZero();
	rhoGrid->binBunch(bunch);
	zGrid->binBunch(bunch);
	rhoGrid->synchronizeMPI(bunch->getMPI_Comm_Local());
	zGrid->synchronizeMPI(bunch->getMPI_Comm_Local());
	
	totalMacrosize = 0.;	
	int nX = rhoGrid->getSizeX();
	int nY = rhoGrid->getSizeY();
	for(int ix = 0; ix < nX; ix++){
		for(int iy = 0; iy < nY; iy++){
			totalMacrosize += rhoGrid->getValueOnGrid(ix,iy);
		}
	}
}

void SpaceChargeCalc2p5D::calculateLongDerivative(){
	int nZ = zGrid->getSizeZ();
	double z_step = zGrid->getStepZ();
	
	if(nZ == 1){
		zDerivGrid->setValue(0.,0);
		return;
	} else if(nZ == 2){
		double zDeriv = 0.;
		zGrid->calcGradient((zGrid->getGridZ(0)+zGrid->getGridZ(1))*0.5,zDeriv);
		zDerivGrid->setValue(zDeriv/z_step,0);
		zDerivGrid->setValue(zDeriv/z_step,1);
		return;
	}
	
	//smoothing the derivative by Quadratic Curve Fitting
	int nAvg = n_long_avg;
	if(nZ < nAvg) nAvg = nZ;
	int iStart, iStop;
	double a,b,c;
	double x,y,z, det;
	double** S = S_arr;
	for(int iz = 0; iz < nZ; iz++){
		z = zGrid->getGridZ(iz);
		
		iStart = iz - nAvg/2;
		if(iStart < 0) iStart = 0;
		if((iStart + nAvg) >= nZ) iStart = nZ - nAvg;
		iStop = iStart + nAvg;
		for(int j = 0; j < 5; j++){
			for(int k = 0; k < 2; k++){
				S[j][k] = 0.;
				for(int i = iStart; i < iStop; i++){
					x = zGrid->getGridZ(i);
					y = zGrid->getValueOnGrid(i);
					S[j][k] += pow(x,j)*pow(y,k);
				}
			}
		}
		det = (S[0][0]*S[2][0]*S[4][0] - S[1][0]*S[1][0]*S[4][0] - S[0][0]*S[3][0]*S[3][0] + 2*S[1][0]*S[2][0]*S[3][0] - S[2][0]*S[2][0]*S[2][0]);
		a = (S[0][1]*S[1][0]*S[3][0] - S[1][1]*S[0][0]*S[3][0] - S[0][1]*S[2][0]*S[2][0]
       + S[1][1]*S[1][0]*S[2][0] + S[2][1]*S[0][0]*S[2][0] - S[2][1]*S[1][0]*S[1][0])/det;
    b = (S[1][1]*S[0][0]*S[4][0] - S[0][1]*S[1][0]*S[4][0] + S[0][1]*S[2][0]*S[3][0]
       - S[2][1]*S[0][0]*S[3][0] - S[1][1]*S[2][0]*S[2][0] + S[2][1]*S[1][0]*S[2][0])/det;
    c = (S[0][1]*S[2][0]*S[4][0] - S[1][1]*S[1][0]*S[4][0] - S[0][1]*S[3][0]*S[3][0]
       + S[1][1]*S[2][0]*S[3][0] + S[2][1]*S[1][0]*S[3][0] - S[2][1]*S[2][0]*S[2][0])/det;
    y = 2*a*z + b;
		zDerivGrid->setValue(y/z_step,iz);
	}
}

/** Sets the number of smoothing points to calculate the derivative of the longitudinal density. */
void SpaceChargeCalc2p5D::setLongAveragingPointsN(int n_points){
	n_long_avg = n_points;
}

/** Returns the number of smoothing points to calculate the derivative of the longitudinal density. */
int SpaceChargeCalc2p5D::getLongAveragingPointsN(){
	return n_long_avg;
}
