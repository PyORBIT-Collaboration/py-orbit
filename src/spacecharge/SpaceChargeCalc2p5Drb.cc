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
	bunchExtremaCalc = new BunchExtremaCalculator();
}

SpaceChargeCalc2p5Drb::SpaceChargeCalc2p5Drb(int xSize, int ySize, int zSize): CppPyWrapper(NULL)
{
	xy_ratio = 1.0;
	poissonSolver = new PoissonSolverFFT2D(xSize, ySize, -xy_ratio, xy_ratio, -1.0, 1.0);
	rhoGrid = new Grid2D(xSize, ySize);
	phiGrid = new Grid2D(xSize, ySize);
	zGrid = new Grid1D(zSize);
	bunchExtremaCalc = new BunchExtremaCalculator();	
}

SpaceChargeCalc2p5Drb::~SpaceChargeCalc2p5Drb(){
	delete poissonSolver;
	delete rhoGrid, phiGrid;
	delete zGrid;
	delete bunchExtremaCalc;
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

void SpaceChargeCalc2p5Drb::trackBunch(Bunch* bunch, double length, double pipe_radius){

	int nPartsGlobal = bunch->getSizeGlobal();
	if(nPartsGlobal < 2) return;
	
	//calculate max and min of X,Y,Z coordinates, a_bunch**2 = 2*<r^2> for the bunch
	//bin paricles and set up limits for rhoGrid, phiGrid, zGrid 
	double a_bunch = 0.;
	double x_center = 0.;
	double y_center = 0.;	
	double totalMacrosize = 0.;
	bunchAnalysis(bunch, totalMacrosize, x_center, y_center, a_bunch);
	double z_step = zGrid->getStepZ();
	
	//calculate phiGrid
	poissonSolver->findPotential(rhoGrid,phiGrid);
	
	SyncPart* syncPart = bunch->getSyncPart();	
	double factor =  2*length*bunch->getClassicalRadius()*pow(bunch->getCharge(),2)/(pow(syncPart->getBeta(),2)*pow(syncPart->getGamma(),3));	
	factor = factor/(z_step*totalMacrosize);
	
	double Lfactor = 0.;
	double x,y,z,ex,ey,ez, r2;
	
	double long_sc_coeff = 0.;
	double long_sc_factor_in = 1.0+2*log(pipe_radius/a_bunch);
	double long_sc_factor_out = 2*log(pipe_radius);
	double a_bunch_2 = a_bunch*a_bunch;
	double long_sc_factor = - length*bunch->getClassicalRadius()*pow(bunch->getCharge(),2) * bunch->getMass()/(pow(syncPart->getGamma(),2));
	
	for (int i = 0, n = bunch->getSize(); i < n; i++){
		x = bunch->x(i);
		y = bunch->y(i);
		z = bunch->z(i);
		r2 = (x - x_center)*(x -x_center)  + (y - y_center)*(y - y_center);
		
		phiGrid->calcGradient(x,y,ex,ey);	
    zGrid->calcGradient(z,ez);
		//std::cout<<"debug ip="<<i<<" x="<<x<<" y="<<y<<" z="<<z<<" ex="<<ex<<" ey="<<ey<<" ez="<<ez<<" rho_z="<< zGrid->getValue(z) <<std::endl;
		
		Lfactor = zGrid->getValue(z) * factor;
	
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

double SpaceChargeCalc2p5Drb::bunchAnalysis(Bunch* bunch, double& totalMacrosize, double& x_c, double& y_c, double& a_bunch){

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

	//setGrids' limits
	rhoGrid->setGridX(xMin,xMax);
	rhoGrid->setGridY(yMin,yMax);	
	
	phiGrid->setGridX(xMin,xMax);
	phiGrid->setGridY(yMin,yMax);
	
	zGrid->setGridZ(zMin,zMax);	
	
	//sizes of the grids are set up
	//bin rho&z Bunch to the Grid
	rhoGrid->binBunch(bunch);
	zGrid->binBunch(bunch);
	
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


