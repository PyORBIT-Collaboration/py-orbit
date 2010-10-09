/////////////////////////////////////////////////////////////////////////////
//
// FILE NAME
//   SpaceChargeCalc2p5Drb.cc
//
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
#include "ParticleMacroSize.hh"
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
}

SpaceChargeCalc2p5Drb::SpaceChargeCalc2p5Drb(int xSize, int ySize, int zSize): CppPyWrapper(NULL)
{
	xy_ratio = 1.0;
	poissonSolver = new PoissonSolverFFT2D(xSize, ySize, -xy_ratio, xy_ratio, -1.0, 1.0);
	rhoGrid = new Grid2D(xSize, ySize);
	phiGrid = new Grid2D(xSize, ySize);
	zGrid = new Grid1D(zSize);
}

SpaceChargeCalc2p5Drb::~SpaceChargeCalc2p5Drb(){
	delete poissonSolver;
	delete rhoGrid, phiGrid;
	delete zGrid;
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
	double x,y,z,ex,ey,ez;
	for (int i = 0, n = bunch->getSize(); i < n; i++){
		x = bunch->x(i);
		y = bunch->y(i);
		z = bunch->z(i);
		
		phiGrid->calcGradient(x,y,ex,ey);		

		Lfactor = zGrid->getValue(z) * factor;
	
		bunch->xp(i) += ex * Lfactor;
		bunch->yp(i) += ey * Lfactor;
	}
}

double SpaceChargeCalc2p5Drb::bunchAnalysis(Bunch* bunch, double& totalMacrosize, double& x_c, double& y_c, double& a_bunch){

	int buff_index0 = 0;
	int buff_index1 = 0;
	double* gridLimArr  = BufferStore::getBufferStore()->getFreeDoubleArr(buff_index0,6);
	double* gridLimArr_out = BufferStore::getBufferStore()->getFreeDoubleArr(buff_index1,6);
	for (int i = 0; i < 3; i++){
		gridLimArr[2*i] = DBL_MAX;
		gridLimArr[2*i+1] = -DBL_MAX;
	}
	
	double** partArr=bunch->coordArr();
	double* coordArr = NULL;
	for (int ip = 0, n = bunch->getSize(); ip < n; ip++){
		coordArr = coordArr;
		if(coordArr[0] < gridLimArr[0]) gridLimArr[0] = coordArr[0];
		if(coordArr[0] > gridLimArr[1]) gridLimArr[1] = coordArr[0];
		if(coordArr[2] < gridLimArr[2]) gridLimArr[2] = coordArr[2];
		if(coordArr[2] > gridLimArr[3]) gridLimArr[3] = coordArr[2];	
		if(coordArr[4] < gridLimArr[4]) gridLimArr[4] = coordArr[4];
		if(coordArr[4] > gridLimArr[5]) gridLimArr[5] = coordArr[4];		
	}
	
	gridLimArr[0] = - gridLimArr[0];
	gridLimArr[2] = - gridLimArr[2];
	gridLimArr[4] = - gridLimArr[4];
	
	ORBIT_MPI_Allreduce(gridLimArr,gridLimArr_out,6,MPI_DOUBLE,MPI_MAX,bunch->getMPI_Comm_Local()->comm);
	
	gridLimArr[0] = - gridLimArr[0];
	gridLimArr[2] = - gridLimArr[2];
	gridLimArr[4] = - gridLimArr[4];
	
	//check if the beam size is not zero ???
	
	double xy_ratio_beam = (gridLimArr[1] - gridLimArr[0])/(gridLimArr[3] - gridLimArr[2]);
	if(xy_ratio_beam > xy_ratio){
		gridLimArr[3] = (xy_ratio_beam/xy_ratio)*gridLimArr[3];
		gridLimArr[2] = (xy_ratio_beam/xy_ratio)*gridLimArr[2];
	} else {
		gridLimArr[0] = gridLimArr[0]/(xy_ratio_beam/xy_ratio);
		gridLimArr[1] = gridLimArr[1]/(xy_ratio_beam/xy_ratio);
	}

	//setGrid
	rhoGrid->setGridX(gridLimArr[0],gridLimArr[1]);
	rhoGrid->setGridY(gridLimArr[2],gridLimArr[3]);	
	
	phiGrid->setGridX(gridLimArr[0],gridLimArr[1]);
	phiGrid->setGridY(gridLimArr[2],gridLimArr[3]);
	
	zGrid->setGridZ(gridLimArr[4],gridLimArr[5]);	
	
	OrbitUtils::BufferStore::getBufferStore()->setUnusedDoubleArr(buff_index0);
	OrbitUtils::BufferStore::getBufferStore()->setUnusedDoubleArr(buff_index1);	
	
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


