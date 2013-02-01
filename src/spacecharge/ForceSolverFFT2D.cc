#include "ForceSolverFFT2D.hh"

#include <iostream>

using namespace OrbitUtils;

// Constructor
ForceSolverFFT2D::ForceSolverFFT2D(int xSize, int ySize): ForceSolver2D(xSize,ySize)
{
  init(xSize,ySize);
}

void ForceSolverFFT2D::init(int xSize, int ySize)
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
			std::cerr << "ForceSolverFFT2D::ForceSolverFFT2D - CONSTRUCTOR \n" 
         				<< "The grid size too small (should be more than 3)! \n" 
								<< "number x bins ="<< xSize_ <<" \n"
								<< "number y bins ="<< ySize_ <<" \n"
								<< "Stop. \n";
		}
		ORBIT_MPI_Finalize();
	}
	
	greensF_ = new std::complex<double>*[xSize2_];
	
	for(int i = 0; i < xSize2_ ; i++) {
		greensF_[i] =  new std::complex<double>[ySize2_];
	}	
	test_      = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * xSize2_ * ySize2_);
	in_        = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * xSize2_ * ySize2_);
	in_res_    = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * xSize2_ * ySize2_);
	out_green_ = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * xSize2_ * ySize2_);
	out_       = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * xSize2_ * ySize2_);
	out_res_   = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * xSize2_ * ySize2_);

	// FFTW_MEASURE or FFTW_ESTIMATE

	planForward_greenF_ = fftw_plan_dft_2d(xSize2_ , ySize2_ , in_,  out_green_, FFTW_FORWARD, FFTW_MEASURE);
	planForward_        = fftw_plan_dft_2d(xSize2_ , ySize2_ , in_,  out_,       FFTW_FORWARD, FFTW_MEASURE);
	planBackward_       = fftw_plan_dft_2d(xSize2_ , ySize2_ , out_res_, in_res_,FFTW_BACKWARD, FFTW_MEASURE);
  
	//define FFT of the Green fuction
	//_defineGreenF();
}

// Destructor
ForceSolverFFT2D::~ForceSolverFFT2D()
{

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


void ForceSolverFFT2D::setGridX(double xMin, double xMax){
	xMin_ = xMin;
	xMax_ = xMax;
	dx_ = (xMax_ - xMin_)/(xSize_ -1);
	_defineGreenF();
}	

void ForceSolverFFT2D::setGridY(double yMin, double yMax){
	yMin_ = yMin;
	yMax_ = yMax;
	dy_ = (yMax_ - yMin_)/(ySize_ -1);
	_defineGreenF();
}

void ForceSolverFFT2D::setGridXY(double xMin, double xMax, double yMin, double yMax){
	xMin_ = xMin;
	xMax_ = xMax;
	dx_ = (xMax_ - xMin_)/(xSize_ -1);
	yMin_ = yMin;
	yMax_ = yMax;
	dy_ = (yMax_ - yMin_)/(ySize_ -1);
	_defineGreenF();
}

// Defines the FFT of the Green Function: field = lambda/r, potential = - lambda*ln(abs(r))
// Please, keep in mind that the field of point like string 2*lambda*ln(abs(r)) in CGS
void ForceSolverFFT2D::_defineGreenF()
{
	
	double rTransY, rTransX, rTot2;
	int i, j, iY , iX;
	
	for (iY = 0; iY < ySize2_/2; iY++)
	{
		rTransY = iY * dy_;
		
		for (iX = 0; iX < xSize2_/2; iX++)
		{
			rTransX = iX * dx_;
			rTot2 = rTransX*rTransX + rTransY*rTransY;
			greensF_[iX][iY] = complex<double>(rTransX/rTot2, rTransY/rTot2);
		}
		
		greensF_[xSize2_/2][iY] = complex<double>(0,0); //end point
		
		for (iX = xSize2_/2+1; iX < xSize2_; iX++)
		{
			rTransX = (iX - xSize2_) * dx_;
			rTot2 = rTransX*rTransX + rTransY*rTransY;
			greensF_[iX][iY] = complex<double>(rTransX/rTot2, rTransY/rTot2);
		}
	}
	
	for(iX=0; iX < xSize_/2; iX++)   // Null the top row:
	{
		greensF_[iX][ySize2_/2] = complex<double>(0,0);
	}
	
	for (iY = ySize2_/2+1; iY < ySize2_; iY++)  // Bottom rows:
	{
		//rTransY = dy_ * (iY - 1 - ySize2_);
		rTransY = dy_ * (iY - ySize2_);
		for (iX = 0; iX < xSize2_/2; iX++)
		{
			rTransX = iX * dx_;
			rTot2 = rTransX*rTransX + rTransY*rTransY;
			greensF_[iX][iY] = complex<double>(rTransX/rTot2, rTransY/rTot2);
	    }
		
		greensF_[xSize2_/2][iY] = complex<double>(0,0); //end point
				
		for (iX = xSize2_/2+1; iX < xSize2_; iX++)
		{
			rTransX = (iX - xSize2_) * dx_;
			//rTransX = (iX - 1 - xSize2_) * dx_;
			rTot2 = rTransX*rTransX + rTransY*rTransY;
			greensF_[iX][iY] = complex<double>(rTransX/rTot2, rTransY/rTot2);
	    }
	}
	greensF_[0][0] = complex<double>(0,0); //end point
	
	//   Calculate the FFT of the Greens Function:
	
	for (j = 0; j < xSize2_; j++)
		for (i = 0; i < ySize2_; i++)
		{
			int index = j + ySize2_*i;
			in_[index][0] = real(greensF_[i][j]);
			in_[index][1] = imag(greensF_[i][j]);
			test_[index][0] = real(greensF_[i][j]);
			test_[index][1] = imag(greensF_[i][j]);
		}
    
		fftw_execute(planForward_greenF_);
		
	//	out_green_re00_ = out_green_[0][0];
		for (j = 0; j < xSize2_; j++)
			for (i = 0; i < ySize2_; i++)
			{
				int index = j + ySize2_*i;
				in_[index][0] = 0.;
				in_[index][1] = 0.;
			}
			
}

void ForceSolverFFT2D::findForce(Grid2D* rhoGrid, Grid2D* forceGridX, Grid2D* forceGridY)
{
	//check sizes of the grids
	if( xSize_ !=  rhoGrid->getSizeX() || ySize_ != rhoGrid->getSizeY()){
		int rank = 0;
		ORBIT_MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		if(rank == 0){
			std::cerr << "ForceSolverFFT2D:" 
			<< "The grid sizes or shape are different "<< std::endl 
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
	

	//out_green_[0][0] = out_green_re00_ + scale_coeff*(xSize2_*ySize2_);
	double** rhosc = rhoGrid->getArr();
	double** fscx = forceGridX->getArr();
	double** fscy = forceGridY->getArr();

	int i, j, index;

	//define the the rho for FFT
	for (j = 0; j < ySize_; j++)
		for (i = 0; i < xSize_; i++)
		{
			in_[j + ySize2_*i][0] = rhosc[i][j];
			in_[j + ySize2_*i][1] = 0.;
		}

	fftw_execute(planForward_);
	
	//do convolution with the FFT of the Green's function 
	double gr_re = 0.;
	double gr_im = 0.;
	for (j = 0; j < xSize2_; j++)
		for (i = 0; i < ySize2_; i++)
		{
			index = j + ySize2_*i;
			gr_re = out_green_[index][0];
			gr_im = out_green_[index][1];
			
			out_res_[index][0] = out_[index][0]*gr_re - out_[index][1]*gr_im;
			out_res_[index][1] = out_[index][0]*gr_im + out_[index][1]*gr_re;
		}

	//do backward FFT
	fftw_execute(planBackward_);

	double denom = 1.0 / (xSize2_*ySize2_);

	for (int iX = 0; iX < xSize_; iX++)
	{
		for (int iY = 0; iY < ySize_; iY++)
		{
			index = iY + ySize2_ * iX; 
			fscx[iX][iY] = in_res_[index][0] * denom;
			fscy[iX][iY] = in_res_[index][1] * denom;
		}
	}
}

