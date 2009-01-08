#include "PoissonSolverFFT2D.hh"

#include <iostream>

using namespace OrbitUtils;

// Constructor
PoissonSolverFFT2D::PoissonSolverFFT2D(int xSize, int ySize): PoissonSolver2D(xSize,ySize)
{
  init(xSize,ySize,xMin_,xMax_,yMin_,yMax_);
}

PoissonSolverFFT2D::PoissonSolverFFT2D(int xSize, int ySize, 
	                     double xMin, double xMax, 
											 double yMin, double yMax): PoissonSolver2D(xSize,ySize,xMin,xMax,yMin,yMax)
{
  init(xSize,ySize,xMin,xMax,yMin,yMax);
}


void PoissonSolverFFT2D::init(int xSize, int ySize, 
                              double xMin, double xMax, 
											        double yMin, double yMax)
{
  //====== Grid parameters and Green function==============
	//To get double size we need only (xSize -1) and (ySize -1) additional grid
	//points, but if we get more it will be fine also.
  xSize2_ = 2*xSize;
  ySize2_ = 2*ySize;

  if( xSize_ < 3 || ySize_ < 3){
		int rank = 0;
		ORBIT_MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		if(rank == 0){
			std::cerr << "PoissonSolverFFT2D::PoissonSolverFFT2D - CONSTRUCTOR \n" 
         				<< "The grid size too small (should be more than 3)! \n" 
								<< "number x bins ="<< xSize_ <<" \n"
								<< "number y bins ="<< ySize_ <<" \n"
								<< "Stop. \n";
		}
		ORBIT_MPI_Finalize();
  }

  greensF_ = new double*[xSize2_];
  for(int i = 0; i < xSize2_ ; i++) {
    greensF_[i] =  new double [ySize2_];
  }
	
  in_        = (double *) fftw_malloc(sizeof(double)*xSize2_ * ySize2_);
  in_res_    = (double *) fftw_malloc(sizeof(double)*xSize2_ * ySize2_);
  out_green_ = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) *xSize2_ * (ySize2_/2+1));
  out_       = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) *xSize2_ * (ySize2_/2+1));
  out_res_   = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) *xSize2_ * (ySize2_/2+1));

	// FFTW_MEASURE or FFTW_ESTIMATE

  planForward_greenF_ = fftw_plan_dft_r2c_2d(xSize2_ , ySize2_ , in_,  out_green_, FFTW_ESTIMATE);
  planForward_        = fftw_plan_dft_r2c_2d(xSize2_ , ySize2_ , in_,  out_,       FFTW_ESTIMATE);
  planBackward_       = fftw_plan_dft_c2r_2d(xSize2_ , ySize2_ , out_res_, in_res_,FFTW_ESTIMATE);
  
  //define FFT of the Green fuction
  _defineGreenF();
}

// Destructor
PoissonSolverFFT2D::~PoissonSolverFFT2D()
{
	
	//std::cerr<<"debug PoissonSolverFFT2D::~PoissonSolverFFT2D() start! "<<std::endl;
  //delete Green function and FFT input and output arrays

  for(int i = 0; i < xSize2_ ; i++) {
    delete [] greensF_[i];
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


void PoissonSolverFFT2D::setGridX(double xMin, double xMax){
	xMin_ = xMin;
	xMax_ = xMax;
	dx_ = (xMax_ - xMin_)/(xSize_ -1);
	_defineGreenF();
}	

void PoissonSolverFFT2D::setGridY(double yMin, double yMax){
	yMin_ = yMin;
	yMax_ = yMax;
	dy_ = (yMax_ - yMin_)/(ySize_ -1);
	_defineGreenF();
}

// Defines the FFT of the Green Function
void PoissonSolverFFT2D::_defineGreenF()
{
	
  double rTransY, rTransX, rTot2;
  int i, j, iY , iX;
	
	for (iY = 0; iY <= ySize2_/2; iY++)
	{
		rTransY = iY * dy_;
		
		for (iX = 0; iX <= xSize2_/2; iX++)
		{
			rTransX = iX * dx_;
			rTot2 = rTransX*rTransX + rTransY*rTransY;
			//we can add constant (to get the same numers as in ORBIT)
			//this constant is + log(1000000.0)
			//here in the original ORBIT we deleted this constant 
			if(iX != 0 || iY != 0){
				greensF_[iX][iY] = - log(rTot2)/2;
			}
			else{
				greensF_[iX][iY] = 0.0;
			}
		}
		
		for (iX = xSize2_/2+1; iX < xSize2_; iX++)
		{
			greensF_[iX][iY] = greensF_[xSize2_-iX][iY];
		}
	}
	
	for (iY = ySize2_/2+1; iY < ySize2_; iY++)
	{
		for (iX = 0; iX < xSize2_; iX++)
		{
			greensF_[iX][iY] = greensF_[iX][ySize2_-iY];
		}
	}
	
	//   Calculate the FFT of the Greens Function:
	
	for (i = 0; i < xSize2_; i++)
		for (j = 0; j < ySize2_; j++)
		{
      in_[j + ySize2_*i] = greensF_[i][j];
		}
    
		fftw_execute(planForward_greenF_);
		
		for (i = 0; i < xSize2_; i++)
			for (j = 0; j < ySize2_; j++)
			{
				in_[j + ySize2_*i] = 0.0;
			}
			
}

void PoissonSolverFFT2D::findPotential(Grid2D* rhoGrid,Grid2D*  phiGrid)
{
	
	checkSizes(rhoGrid,phiGrid);
	double** rhosc = rhoGrid->getArr();
	double** phisc = phiGrid->getArr();
		
  int i, j, index;

  //define the the rho for FFT
  for (i = 0; i < xSize_; i++)
  for (j = 0; j < ySize_; j++)
  {
    in_[j + ySize2_*i] = rhosc[i][j];
  }

	fftw_execute(planForward_);
	
  //do convolution with the FFT of the Green's function 
  for (i = 0; i < xSize2_; i++)
  for (j = 0; j < ySize2_/2+1; j++)
  {
    index = j + (ySize2_/2+1)*i;
    out_res_[index][0] = out_[index][0]*out_green_[index][0] - out_[index][1]*out_green_[index][1];
    out_res_[index][1] = out_[index][0]*out_green_[index][1] + out_[index][1]*out_green_[index][0];
  }

  //do backward FFT
	fftw_execute(planBackward_);

  //set the potential
  double denom = 1.0 / (xSize2_*ySize2_);	
  for (i = 0; i < xSize_; i++)
  for (j = 0; j < ySize_; j++)
  {
    phisc[i][j] = denom * in_res_[j + ySize2_*i];
  }
}

