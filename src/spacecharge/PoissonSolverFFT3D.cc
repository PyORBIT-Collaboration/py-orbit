#include "PoissonSolverFFT3D.hh"

#include <iostream>

using namespace OrbitUtils;

// Constructor
PoissonSolverFFT3D::PoissonSolverFFT3D(int xSize, int ySize, int zSize): PoissonSolver3D(xSize,ySize,zSize)
{
  init(xSize,ySize,zSize,xMin_,xMax_,yMin_,yMax_,zMin_,zMax_);
}

PoissonSolverFFT3D::PoissonSolverFFT3D(int xSize, int ySize, int zSize,
	                     double xMin, double xMax, 
											 double  yMin, double yMax,
											 double  zMin, double zMax): PoissonSolver3D(xSize,ySize,zSize,xMin,xMax,yMin,yMax,zMin,zMax)
{
  init(xSize,ySize,zSize,xMin,xMax,yMin,yMax,zMin,zMax);
}


void PoissonSolverFFT3D::init(int xSize, int ySize, int zSize,
                              double xMin, double xMax, 
											        double yMin, double yMax,
															double  zMin, double zMax)
{
  //====== Grid parameters and Green function==============
	//To get double size we need only (xSize -1), (ySize -1), and (zSize -1) additional grid
	//points, but if we get more it will be fine also.
  xSize2_ = 2*xSize;
  ySize2_ = 2*ySize;
  zSize2_ = 2*zSize;

  if( xSize_ < 3 || ySize_ < 3){
		int rank = 0;
		ORBIT_MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		if(rank == 0){
			std::cerr << "PoissonSolverFFT3D::PoissonSolverFFT3D - CONSTRUCTOR \n" 
         				<< "The grid size too small (should be more than 3)! \n" 
								<< "number x bins ="<< xSize_ <<" \n"
								<< "number y bins ="<< ySize_ <<" \n"
								<< "number z bins ="<< zSize_ <<" \n"
								<< "Stop. \n";
		}
		ORBIT_MPI_Finalize();
  }

  greensF_ = new double**[xSize2_];
  for(int ix = 0; ix < xSize2_ ; ix++) {
		greensF_[ix] =  new double* [ySize2_];
		for(int iy = 0; iy < ySize2_ ; iy++) {
			greensF_[ix][iy] =  new double [zSize2_];
		}
  }
	
  in_        = (double *) fftw_malloc(sizeof(double)*xSize2_ * ySize2_ * zSize2_);
  in_res_    = (double *) fftw_malloc(sizeof(double)*xSize2_ * ySize2_ * zSize2_);
  out_green_ = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) *xSize2_ * ySize2_ * (zSize2_/2+1));
  out_       = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) *xSize2_ * ySize2_ * (zSize2_/2+1));
  out_res_   = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) *xSize2_ * ySize2_ * (zSize2_/2+1));

	// FFTW_MEASURE or FFTW_ESTIMATE

  planForward_greenF_ = fftw_plan_dft_r2c_3d(xSize2_ , ySize2_ , zSize2_ , in_,  out_green_, FFTW_ESTIMATE);
  planForward_        = fftw_plan_dft_r2c_3d(xSize2_ , ySize2_ , zSize2_ , in_,  out_,       FFTW_ESTIMATE);
  planBackward_       = fftw_plan_dft_c2r_3d(xSize2_ , ySize2_ , zSize2_ , out_res_, in_res_,FFTW_ESTIMATE);
  
  //define FFT of the Green fuction
  _defineGreenF();
}

// Destructor
PoissonSolverFFT3D::~PoissonSolverFFT3D()
{
	
	//std::cerr<<"debug PoissonSolverFFT3D::~PoissonSolverFFT3D() start! "<<std::endl;
  //delete Green function and FFT input and output arrays

  for(int ix = 0; ix < xSize2_ ; ix++) {
		for(int iy = 0; iy < ySize2_ ; iy++) {	
			delete [] greensF_[ix][iy];
		}   
		delete [] greensF_[ix];
	}
  delete [] greensF_;

  fftw_free(in_);
  fftw_free(in_res_);
  fftw_free(out_green_);
  fftw_free(out_);
  fftw_free(out_res_);
	
  fftw_destroy_plan(planForward_greenF_);
  fftw_destroy_plan(planForward_);
  fftw_destroy_plan(planBackward_);
}


void PoissonSolverFFT3D::setGridX(double xMin, double xMax){
	xMin_ = xMin;
	xMax_ = xMax;
	dx_ = (xMax_ - xMin_)/(xSize_ -1);
	_defineGreenF();
}	

void PoissonSolverFFT3D::setGridY(double yMin, double yMax){
	yMin_ = yMin;
	yMax_ = yMax;
	dy_ = (yMax_ - yMin_)/(ySize_ -1);
	_defineGreenF();
}

void PoissonSolverFFT3D::setGridZ(double zMin, double zMax){
	zMin_ = zMin;
	zMax_ = zMax;
	dz_ = (zMax_ - zMin_)/(zSize_ -1);
	_defineGreenF();
}

void PoissonSolverFFT3D::setGridXYZ(double xMin, double xMax, double yMin, double yMax, double zMin, double zMax){
	xMin_ = xMin;
	xMax_ = xMax;
	dx_ = (xMax_ - xMin_)/(xSize_ -1);	
	yMin_ = yMin;
	yMax_ = yMax;
	dy_ = (yMax_ - yMin_)/(ySize_ -1);
	zMin_ = zMin;
	zMax_ = zMax;
	dz_ = (zMax_ - zMin_)/(zSize_ -1);	
	_defineGreenF();
}

// Defines the FFT of the Green Function: field = lambda/r, potential = - lambda*ln(abs(r))
// Please, keep in mind that the field of point like string 2*lambda*ln(abs(r)) in CGS
void PoissonSolverFFT3D::_defineGreenF()
{
  double rTransY, rTransX, rTransZ, rTot;
  int i, j, k, iY , iX, iZ;
	
	for (iZ = 0; iZ <= zSize2_/2; iZ++)
	{
		rTransZ = iZ * dz_;	
		
		for (iY = 0; iY <= ySize2_/2; iY++)
		{
			rTransY = iY * dy_;
			
			for (iX = 0; iX <= xSize2_/2; iX++)
			{
				rTransX = iX * dx_;	
				rTot = sqrt(rTransX*rTransX + rTransY*rTransY + rTransZ*rTransZ);
				
				if(iX != 0 || iY != 0 || iZ != 0){
					greensF_[iX][iY][iZ] = 1./rTot;
				}
				else{
					greensF_[iX][iY][iZ] = 0.0;
				}
			}
			
			for (iX = xSize2_/2+1; iX < xSize2_; iX++)
			{
				greensF_[iX][iY][iZ] = greensF_[xSize2_-iX][iY][iZ];
			}
		}
		
		for (iY = ySize2_/2+1; iY < ySize2_; iY++)
		{
			for (iX = 0; iX < xSize2_; iX++)
			{
				greensF_[iX][iY][iZ] = greensF_[iX][ySize2_-iY][iZ];
			}
		}
	}
	
	for (iZ = zSize2_/2+1; iZ < zSize2_; iZ++)
	{
		for (iX = 0; iX < xSize2_; iX++)
		{
			for (iY = 0; iY < ySize2_; iY++){
				greensF_[iX][iY][iZ] = greensF_[iX][iY][zSize2_-iZ];
			}
		}
	}		
	//   Calculate the FFT of the Greens Function:
	
	for (i = 0; i < xSize2_; i++)
		for (j = 0; j < ySize2_; j++)
		  for (k = 0; k < zSize2_; k++)
		  {
        in_[k + zSize2_*j + zSize2_*ySize2_*i] = greensF_[i][j][k];
		  }
    
		fftw_execute(planForward_greenF_);
		
		for (i = 0; i < xSize2_; i++)
			for (j = 0; j < ySize2_; j++)
			  for (k = 0; k < zSize2_; k++)
			  {
				  in_[k + zSize2_*j + zSize2_*ySize2_*i] = 0.0;
			  }
			
}

void PoissonSolverFFT3D::findPotential(Grid3D* rhoGrid,Grid3D*  phiGrid)
{
	double shape_diff_limit = 0.0000001;
	//check sizes of the grids
  if( xSize_ !=  rhoGrid->getSizeX() || ySize_ != rhoGrid->getSizeY() ||
		  xSize_ !=  phiGrid->getSizeX() || ySize_ != phiGrid->getSizeY() ||
		  zSize_ !=  phiGrid->getSizeZ() || zSize_ != rhoGrid->getSizeZ() ||
		  fabs(dx_/dy_ - rhoGrid->getStepX()/rhoGrid->getStepY()) > shape_diff_limit ||
			fabs(dx_/dy_ - phiGrid->getStepX()/phiGrid->getStepY()) > shape_diff_limit ||
		  fabs(dz_/dy_ - rhoGrid->getStepZ()/rhoGrid->getStepY()) > shape_diff_limit ||
			fabs(dz_/dy_ - phiGrid->getStepZ()/phiGrid->getStepY()) > shape_diff_limit ||
			phiGrid->getMinX()  != rhoGrid->getMinX()  || phiGrid->getMinY() != rhoGrid->getMinY() ||
			phiGrid->getStepX() != rhoGrid->getStepX() || phiGrid->getStepY() != rhoGrid->getStepY() ||
			phiGrid->getStepZ() != rhoGrid->getStepZ() || phiGrid->getMinZ() != rhoGrid->getMinZ()){
		int rank = 0;
		ORBIT_MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		if(rank == 0){
			std::cerr << "PoissonSolverFFT3D:" 
			<< "The grid sizes or shape are different "<< std::endl 
								<< "number x bins ="<< xSize_ << std::endl
								<< "number y bins ="<< ySize_ << std::endl
								<< "number z bins ="<< zSize_ << std::endl
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
								<< "phiGrid dx ="<< phiGrid->getStepX() <<std::endl
								<< "phiGrid dy ="<< phiGrid->getStepY() <<std::endl
								<< "xMin ="<< xMin_ <<std::endl
								<< "yMin  ="<< yMin_ <<std::endl
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

	double*** rhosc = rhoGrid->getArr3D();
	double*** phisc = phiGrid->getArr3D();
	
	double scale_coeff = dx_/rhoGrid->getStepX();
		
  int i, j, k, index;

  //define the the rho for FFT
  for (i = 0; i < xSize_; i++)
  for (j = 0; j < ySize_; j++)
	for (k = 0; k < zSize_; k++)
  {
    in_[k + zSize2_*(j + ySize2_*i)] = rhosc[k][i][j];
  }

	fftw_execute(planForward_);
	
  //do convolution with the FFT of the Green's function 
	double gr_re = 0.;
	double gr_im = 0.;
  for (i = 0; i < xSize2_; i++)
  for (j = 0; j < ySize2_; j++)
	for (k = 0; k < zSize2_/2+1; k++)
  {
    index = k + (zSize2_/2+1)*(j + ySize2_*i);
		gr_re = out_green_[index][0];
		gr_im = out_green_[index][1];
		//std::cout<<" debug gr_re="<<gr_re<<" gr_im="<<gr_im<<std::endl;
    out_res_[index][0] = out_[index][0]*gr_re - out_[index][1]*gr_im;
    out_res_[index][1] = out_[index][0]*gr_im + out_[index][1]*gr_re;
  }

  //do backward FFT
	fftw_execute(planBackward_);
	
  //set the potential
  double denom =  scale_coeff/ (xSize2_*ySize2_*zSize2_);	
  for (i = 0; i < xSize_; i++)
  for (j = 0; j < ySize_; j++)
	for (k = 0; k < zSize_; k++)
  {
    phisc[k][i][j] = denom * in_res_[k + zSize2_*(j + ySize2_*i)];
  }
}

