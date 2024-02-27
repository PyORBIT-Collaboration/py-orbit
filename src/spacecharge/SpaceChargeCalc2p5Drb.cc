/////////////////////////////////////////////////////////////////////////////
//
// FILE NAME
//   SpaceChargeCalc2p5Drb.cc
//
// AUTHOR
//    A. Shishlo
//
// Created:
//   10/09/10
//
// DESCRIPTION
// This class calculates the space charge kicks for bunch using 2.5D approach in transverse direction
// and Rick Baartman's approach (RB) to the longitudinal kicks. There is a hope that will be simular to
// the true 3D discription of space charge. 
//
/////////////////////////////////////////////////////////////////////////////

#include "Grid1D.hh"
#include "Grid2D.hh"
#include "PoissonSolverFFT2D.hh"
#include "SpaceChargeCalc2p5Drb.hh"
#include "BufferStore.hh"

#include <iostream>
#include <cmath>
#include <cfloat>

using namespace OrbitUtils;

SpaceChargeCalc2p5Drb::SpaceChargeCalc2p5Drb(int xSize, int ySize, int zSize, double xy_ratio_in): CppPyWrapper(NULL)
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

SpaceChargeCalc2p5Drb::SpaceChargeCalc2p5Drb(int xSize, int ySize, int zSize): CppPyWrapper(NULL)
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

SpaceChargeCalc2p5Drb::~SpaceChargeCalc2p5Drb(){
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

Grid2D* SpaceChargeCalc2p5Drb::getRhoGrid(){
	return rhoGrid;
}

Grid2D* SpaceChargeCalc2p5Drb::getPhiGrid(){
	return phiGrid;
}

Grid1D* SpaceChargeCalc2p5Drb::getLongGrid(){
	return zGrid;
}

Grid1D* SpaceChargeCalc2p5Drb::getLongDerivativeGrid(){
	return zDerivGrid;
}

void SpaceChargeCalc2p5Drb::trackBunch(Bunch* bunch, double length, double pipe_radius){

	int nPartsGlobal = bunch->getSizeGlobal();
	if(nPartsGlobal < 2) return;
	
	//calculate max and min of X,Y,Z coordinates, a_bunch**2 = 2*<r^2> for the bunch
	//bin paricles and set up limits for rhoGrid, phiGrid, zGrid 
	double a_bunch = 0.;
	double x_center = 0.;
	double y_center = 0.;	
	double totalMacrosize = 0.;
	this->bunchAnalysis(bunch, totalMacrosize, x_center, y_center, a_bunch);
	double z_step = zGrid->getStepZ();
	
	//calculate phiGrid
	poissonSolver->findPotential(rhoGrid,phiGrid);
	
	SyncPart* syncPart = bunch->getSyncPart();	
	double factor =  2*length*bunch->getClassicalRadius()/(pow(syncPart->getBeta(),2)*pow(syncPart->getGamma(),3));	
	//std::cout<<" debug totalMacrosize="<<totalMacrosize<<" factor="<<factor<<" z_step="<< z_step <<std::endl;	
	

	factor = factor/(z_step*totalMacrosize);

	double Lfactor = 0.;
	double x,y,z,ex,ey,ez, r2;
	
	double long_sc_coeff = 0.;
	double long_sc_factor_in = 1.0+2*log(pipe_radius/a_bunch);
	double long_sc_factor_out = 2*log(pipe_radius);
	double a_bunch_2 = a_bunch*a_bunch;
	double long_sc_factor = - length*bunch->getClassicalRadius()* bunch->getMass()/(pow(syncPart->getGamma(),2));
	//std::cout<<" debug pipe_radius="<<pipe_radius<<" a_bunch="<<a_bunch<<std::endl;	
	//std::cout<<" debug long_sc_factor_in="<<long_sc_factor_in<<" long_sc_factor_out="<<long_sc_factor_out<<std::endl;	
	//std::cout<<" debug long_sc_factor="<<long_sc_factor<<std::endl;
	//std::cout<<" debug z="<<zGrid->getMaxZ()*0.5<<" derivat="<<zDerivGrid->getValue(zGrid->getMaxZ()*0.5)<<std::endl;	
	for (int i = 0, n = bunch->getSize(); i < n; i++){
		x = bunch->x(i);
		y = bunch->y(i);
		z = bunch->z(i);
		r2 = (x - x_center)*(x -x_center)  + (y - y_center)*(y - y_center);
		
		phiGrid->calcGradient(x,y,ex,ey);	
    ez = zDerivGrid->getValue(z);
		//std::cout<<"debug ip="<<i<<" x="<<x<<" y="<<y<<" z="<<z<<" ex="<<ex<<" ey="<<ey<<" ez="<<ez<<" rho_z="<< zGrid->getValue(z) <<std::endl;
		
		Lfactor = - zGrid->getValue(z) * factor;
	
		bunch->xp(i) += ex * Lfactor;
		bunch->yp(i) += ey * Lfactor;
		
		if(r2 <= a_bunch_2){
			long_sc_coeff = long_sc_factor_in - r2/a_bunch_2;
		} else {
			long_sc_coeff = long_sc_factor_out - log(r2);
		}
		bunch->dE(i) += ez*long_sc_factor*long_sc_coeff;
	}
}

void SpaceChargeCalc2p5Drb::bunchAnalysis(Bunch* bunch, double& totalMacrosize, double& x_c, double& y_c, double& a_bunch){

	double xMin, xMax, yMin, yMax, zMin, zMax;
	
	bunchExtremaCalc->getExtremaXYZ(bunch, xMin, xMax, yMin, yMax, zMin, zMax);
	
	//check if the beam size is not zero 
  if( xMin >=  xMax || yMin >=  yMax || zMin >=  zMax){
		int rank = 0;
		ORBIT_MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		if(rank == 0){
			std::cerr << "SpaceChargeCalc2p5Drb::bunchAnalysis(bunch,...)" << std::endl
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

	//set Grids' limits
	rhoGrid->setGridX(xMin,xMax);
	rhoGrid->setGridY(yMin,yMax);	
	
	phiGrid->setGridX(xMin,xMax);
	phiGrid->setGridY(yMin,yMax);
	
	zGrid->setGridZ(zMin,zMax);	
	zDerivGrid->setGridZ(zMin,zMax);	
	
	//sizes of the grids are set up
	//bin rho&z Bunch to the Grid
	rhoGrid->setZero();
	zGrid->setZero();
	zDerivGrid->setZero();
	rhoGrid->binBunch(bunch);
	zGrid->binBunch(bunch);
	rhoGrid->synchronizeMPI(bunch->getMPI_Comm_Local());
	zGrid->synchronizeMPI(bunch->getMPI_Comm_Local());
	
	//calculate the derivative of the longitudinal density that inside the zDerivGrid (Grid1D)
	this->calculateLongDerivative();
	
	//calculate x_avg, x2_avg, y_avg, y2_avg
	double x_avg = 0., x2_avg = 0., y_avg = 0., y2_avg = 0.;
	totalMacrosize = 0.;
	double val = 0., x = 0., y = 0.;
	int nX = rhoGrid->getSizeX();
	int nY = rhoGrid->getSizeY();
	for(int ix = 0; ix < nX; ix++){
		for(int iy = 0; iy < nY; iy++){
			val = rhoGrid->getValueOnGrid(ix,iy);
			x = rhoGrid->getGridX(ix);
			y = rhoGrid->getGridX(iy);
			totalMacrosize += val;
			x_avg += x*val;
			x2_avg += x*x*val;
			y_avg += y*val;
			y2_avg += y*y*val;			
		}
	}
	x_avg /= totalMacrosize;
	x2_avg /= totalMacrosize;
	y_avg /= totalMacrosize;
	y2_avg /= totalMacrosize;	
	x_c = x_avg;
	y_c = y_avg;
	x2_avg = fabs(x2_avg - x_avg*x_avg);
	y2_avg = fabs(y2_avg - y_avg*y_avg);
	a_bunch = sqrt(2*(x2_avg + y2_avg));
}

void SpaceChargeCalc2p5Drb::calculateLongDerivative(){
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
  //The quadratic representation y = a*x^2+b*x+c is built for 
  //x = (z - zGrid->getGridZ(iStart) + z_step)/z_step
  //At the end we transform the derivative back to dy/dz 
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
          x = 1.0*(i - iStart + 1);
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
    
    y = 2*a*(z - zGrid->getGridZ(iStart) + z_step)/z_step + b;
    y = y / z_step;
    zDerivGrid->setValue(y/z_step,iz);
  }
}

/** Sets the number of smoothing points to calculate the derivative of the longitudinal density. */
void SpaceChargeCalc2p5Drb::setLongAveragingPointsN(int n_points){
	n_long_avg = n_points;
}

/** Returns the number of smoothing points to calculate the derivative of the longitudinal density. */
int SpaceChargeCalc2p5Drb::getLongAveragingPointsN(){
	return n_long_avg;
}
