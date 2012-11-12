#include "BaseBoundary2D.hh"

#include <iostream>
#include <cfloat>

using namespace OrbitUtils;

const int BaseBoundary2D::IS_INSIDE    =  1;
const int BaseBoundary2D::IS_OUTSIDE   = -1;
const int BaseBoundary2D::TO_BE_KILLED =  0;

const double BaseBoundary2D::PI = 3.14159265358979324;


// Constructor
BaseBoundary2D::BaseBoundary2D(int nPoints, int nModes): CppPyWrapper(NULL)
{

	initialized_ = 0;
	
	xMin_ =   DBL_MAX;
	xMax_ = - DBL_MAX;
 	yMin_ =   DBL_MAX;
	yMax_ = - DBL_MAX;
	
	BPrnorm_ = DBL_MAX;
	
	shape_ = "noshape";
	NO_SHAPE = shape_;
	shape_type_ = -1;	
	no_shape_key_ = 1;
	xDim_ = 0.;
	yDim_ = 0.;	
	
	//number of points on the boundary
	nPoints_ = nPoints;
	
	//number of modes in the free-space field
	nModes_ = nModes;

	//The arrays with boundary points 
	bArrX_ = new double[nPoints_];
	bArrY_ = new double[nPoints_];
	for(int i = 0; i < nPoints_; i++){
			bArrX_[i] = DBL_MAX;
			bArrY_[i] = DBL_MAX;
	}
	
	
  if( 2*nModes_+1 >= nPoints_){
		int rank = 0;
		ORBIT_MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		if(rank == 0){
			std::cerr << "BaseBoundary2D::BaseBoundary2D - CONSTRUCTOR " << std::endl
			<< " Boundary class will not work because "<< std::endl
			<< " 2*nModes_+1 = " << 2*nModes_+1 << std::endl
			<< " nPoints_    = " << nPoints_    << std::endl
			<< " Number of points has to be increased."<<std::endl
			<< "Stop."<< std::endl;
		}
		ORBIT_MPI_Finalize();		
  }
	
	bnd_phi_arr_ = new double[nPoints_];
	
  lsq_func_vctr_ = new double[2*nModes_+1];
  lsq_coef_vctr_  = new double[2*nModes_+1];
  for(int i = 0, n =  (2*nModes_+1); i < n; i++) {  
    lsq_coef_vctr_[i] = 0.0;
		lsq_func_vctr_[i] = 0.0;
	}
	
	indxc_tmp= new int[2*nModes_+1];
	indxr_tmp= new int[2*nModes_+1];
	ipiv_tmp = new int[2*nModes_+1];

  tmp_matrix_ =new double*[2*nModes_+1];
  for(int i = 0; i < (2*nModes_+1); i++) {
    tmp_matrix_[i] =  new double [2*nModes_+1];
  } 

  LSQ_matrix_ =new double*[2*nModes_+1];
  for(int i = 0; i < (2*nModes_+1); i++) {
    LSQ_matrix_[i] =  new double [nPoints_];
  }	
  
}

// Destructor
BaseBoundary2D::~BaseBoundary2D()
{
	if(bArrX_ != NULL) delete [] bArrX_;
	if(bArrY_ != NULL) delete [] bArrY_;
	
  //delete LSQM related arrays
	delete [] bnd_phi_arr_;
  delete [] lsq_func_vctr_;
  delete [] lsq_coef_vctr_;

	delete [] ipiv_tmp;
	delete [] indxr_tmp;
	delete [] indxc_tmp;
	
  for(int i = 0; i < (2*nModes_+1); i++) {
    delete [] tmp_matrix_[i];
  }
  delete [] tmp_matrix_;

  for(int i = 0; i < (2*nModes_+1); i++) {
    delete [] LSQ_matrix_[i];
  }
  delete [] LSQ_matrix_;
	
  //std::cerr << "debug Destructor BaseBoundary2D was done.   !!!! \n";
}

/** Returns the maximal value of the grid in x-axis */   
double BaseBoundary2D::getMaxX(){ return xMax_;}

/** Returns the minimal value of the grid in x-axis */   	
double BaseBoundary2D::getMinX(){ return xMin_;}

/** Returns the maximal value of the grid in y-axis */ 	
double BaseBoundary2D::getMaxY(){ return yMax_;}

/** Returns the minimal value of the grid in y-axis */ 	
double BaseBoundary2D::getMinY(){ return yMin_;}

/** Returns the number of boundary points */
int BaseBoundary2D::getNumberOfPoints(){ return nPoints_;}

/** Returns the number of modes in the free-space field */
int BaseBoundary2D::getNumberOfModes(){ return nModes_;}

/** Returns 0 if the boundary is not initialized */
int BaseBoundary2D::isInitialized(){ return initialized_;}

/** Returns IS_INSIDE or IS_OUTSIDE depending on the particle's position */
int BaseBoundary2D::isInside(double x, double y)
{
  int isIns = IS_OUTSIDE;
	if( x > xMin_ && x < xMax_ && y > yMin_ && y < yMax_){
		isIns = IS_INSIDE;
	}
  return isIns;
}

/** Returns the x-coordinate of the boundary point with an index i */
double BaseBoundary2D::getBoundaryPointX(int i){
	if(i < nPoints_ && i >= 0) return bArrX_[i];
	return DBL_MAX;
}

/** Returns the y-coordinate of the boundary point with an index i */
double BaseBoundary2D::getBoundaryPointY(int i){
	if(i < nPoints_ && i >= 0) return bArrY_[i];
	return DBL_MAX;
}

/** Sets the boundary point with index to (x,y) */
void BaseBoundary2D::setBoundaryPoint(int index, double x, double y){
	
	if(no_shape_key_ != 1){
		int rank = 0;
		ORBIT_MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		if(rank == 0){
			std::cerr << "BaseBoundary2D::setBoundaryPoint(...) " << std::endl
			<< " You can not do this! The shape has been defined already! "<< std::endl
			<< " shape = " << shape_ << std::endl
			<< " x size = " << xDim_    << std::endl
			<< " y size = " << yDim_    << std::endl
			<< "Stop."<< std::endl;
		}
		ORBIT_MPI_Finalize();			
	}
	
	if(index < nPoints_ && index >= 0){
		bArrX_[index] = x;
		bArrY_[index] = y;				
	} 
	initialized_ = 0;
}

/** Initializes all arrays related to boundary points  */
void BaseBoundary2D::initializeBPs(){
	xMin_ =   DBL_MAX;
	xMax_ = - DBL_MAX;
 	yMin_ =   DBL_MAX;
	yMax_ = - DBL_MAX;	
	for(int i = 0; i < nPoints_; i++){
		double x = bArrX_[i];
		double y = bArrY_[i];
		if(xMin_ > x) xMin_ = x;
		if(xMax_ < x) xMax_ = x;
		if(yMin_ > y) yMin_ = y;
		if(yMax_ < y) yMax_ = y;		
	}
	
	BPrnorm_ = sqrt((xMin_-xMax_)*(xMin_-xMax_) + (yMin_-yMax_)*(yMin_-yMax_));
	
  for(int i = 0; i < (2*nModes_+1); i++) {
  for(int j = 0; j < (2*nModes_+1); j++) {
    tmp_matrix_[i][j] =  0.0;
  }}

  for (int  iBp = 0; iBp < nPoints_ ; iBp++){
    lsq_fuctions(bArrX_[iBp],bArrY_[iBp]);
     for(int i = 0; i < (2*nModes_+1); i++) {
     for(int j = 0; j < (2*nModes_+1); j++) {
       tmp_matrix_[i][j] += lsq_func_vctr_[i]*lsq_func_vctr_[j];
     }}
  }

  //inverse matrix
  _gaussjinv(tmp_matrix_,2*nModes_+1);

  for(int i = 0; i < (2*nModes_+1); i++) {
  for(int j = 0; j < nPoints_     ; j++) {
    LSQ_matrix_[i][j] = 0.0;
    lsq_fuctions(bArrX_[j],bArrY_[j]);
    for (int k = 0; k < (2*nModes_+1); k++) {
      LSQ_matrix_[i][j] += tmp_matrix_[i][k] * lsq_func_vctr_[k];
    }    
  }} 	
	
	//boundary is initialized
	initialized_ = 1;
}

/** Returns the name of the shape */
string BaseBoundary2D::getShapeName(){
	return shape_;
}

/** Returns the shape index */
int BaseBoundary2D::getShapeType(){
	return shape_type_;
}

void BaseBoundary2D::_gaussjinv(double **a, int n)
{
	//    Taken from Nrecipes and slightly modified.
	//   Get matrix A(nxn) and transform it into A^(-1).
	//   Outputs error if A^(-1) doesn't exist
	
	int *indxc,*indxr,*ipiv;  // indexr and indexc track column permutation
	int i,icol,irow,j,k,l,ll;
	double big,dum,pivinv,temp;
	
	icol = 0;
	irow = 0;
	
	indxc= indxc_tmp;
	indxr= indxr_tmp;
	ipiv = ipiv_tmp;
	
	for (j=0;j<n;j++) ipiv[j]=0;
	for (i=0;i<n;i++) {
		big=0.0;
		// Looking for pivot
		for (j=0;j<n;j++)
			if (ipiv[j] != 1)
			for (k=0;k<n;k++) {
				if (ipiv[k] == 0) {
					if (fabs(a[j][k]) >= big) {
						big=fabs(a[j][k]);
						irow=j;
						icol=k;
					}
				} 
				else if (ipiv[k] > 1) {
					int rank = 0;
					ORBIT_MPI_Comm_rank(MPI_COMM_WORLD, &rank);
					if(rank == 0){
						std::cerr << "Boundary2D::Boundary2D - _gaussjinv"<<std::endl
						<< "Matrix A^(-1) does not_eq exist!"<<std::endl
						<< "Stop 0."<<std::endl;
					}
					ORBIT_MPI_Finalize();					
				};
			}
			++(ipiv[icol]);
			// Pivot found - interchanging rows
			if (irow != icol) {
				for (l=0;l<n;l++) {
					temp = a[irow][l];
					a[irow][l] = a[icol][l];
					a[icol][l] = temp;
				}
			}
			indxr[i]=irow;
			indxc[i]=icol;
			if (a[icol][icol] == 0.0) {
				int rank = 0;
				ORBIT_MPI_Comm_rank(MPI_COMM_WORLD, &rank);
				if(rank == 0){
					std::cerr << "Boundary2D::Boundary2D - _gaussjinv"<<std::endl
					<< "Matrix A^(-1) does not_eq exist!"<<std::endl
					<< "Stop 1."<<std::endl;
				}
				ORBIT_MPI_Finalize();							
			}
			pivinv=1.0/a[icol][icol];
			a[icol][icol]=1.0;
			for (l=0;l<n;l++) a[icol][l] *= pivinv;
			for (ll=0;ll<n;ll++)
			if (ll != icol) {
				dum=a[ll][icol];
				a[ll][icol]=0.0;
				for (l=0;l<n;l++) a[ll][l] -= a[icol][l]*dum;
			}
	}
	for (l=n-1;l>=0;l--) {
		if (indxr[l] != indxc[l])
		for (k=0;k<n;k++){
			temp = a[k][indxr[l]];
			a[k][indxr[l]] = a[k][indxc[l]];
			a[k][indxc[l]] = temp;
		}
	}
}

/** Calculates all of the LSQM functions at one point */
void BaseBoundary2D::lsq_fuctions(double x, double y)
{
  int i = 0;
  double r,r2;
  r2 = x*x + y*y;
  if( r2 == 0. ) {
    lsq_func_vctr_[0] = 1.0; 
     for( i = 1; i < (2*nModes_+1); i++) {
       lsq_func_vctr_[i] = 0.0;
     }
     return;
  }
  
  r = sqrt(r2);
  double sin_f,cos_f, sin_f0,cos_f0, sin_f1,cos_f1;
  double rfac,rj;
  sin_f = y/r;
  cos_f = x/r;  
  sin_f0 = sin_f;
  cos_f0 = cos_f;
  rfac = r / BPrnorm_;
  rj = 1.0;
  lsq_func_vctr_[0] = 1.;
  for( i = 0; i < nModes_; i++) {
    rj *= rfac; 
    lsq_func_vctr_[2*i+1]= rj*cos_f0;
    lsq_func_vctr_[2*i+2]= rj*sin_f0;
    sin_f1 = sin_f*cos_f0 + cos_f*sin_f0;
    cos_f1 = cos_f*cos_f0 - sin_f*sin_f0;
    sin_f0 = sin_f1;
    cos_f0 = cos_f1;
  }
}


/** Adds potential from the boundary to the external grid */
void BaseBoundary2D::addBoundaryPotential(Grid2D* rhoGrid,Grid2D*  phiGrid){

	//check initialization
	if(initialized_ == 0){
		int rank = 0;
		ORBIT_MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		if(rank == 0){
			std::cerr << "BaseBoundary2D (or subclass) addBoundaryPotential(...):" 
			<< "The boundary was not initialized! Call initializeBPs() method!"<< std::endl 
			<< "Stop. \n";
		}
		ORBIT_MPI_Finalize();		
	}
	
	
	//check that the boundary is inside the grids
  if( rhoGrid->getMinX() > xMin_ || rhoGrid->getMaxX() < xMax_ ||
			rhoGrid->getMinY() > yMin_ || rhoGrid->getMaxY() < yMax_){
		int rank = 0;
		ORBIT_MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		if(rank == 0){
			std::cerr << "BaseBoundary2D (or subclass) addBoundaryPotential(...):" 
			<< "The some boundary points are outside the grids:"<< std::endl 
								<< "boundary xMin ="<< xMin_ <<std::endl
								<< "boundary xMax ="<< xMax_ <<std::endl
								<< "boundary yMin ="<< yMin_ <<std::endl
								<< "boundary yMax ="<< yMax_ <<std::endl
								<< "rhoGrid xMin ="<< rhoGrid->getMinX() <<std::endl
								<< "rhoGrid yMin ="<< rhoGrid->getMinY() <<std::endl
								<< "rhoGrid xMax ="<< rhoGrid->getMaxX() <<std::endl
								<< "rhoGrid yMax ="<< rhoGrid->getMaxY() <<std::endl

								<< "Stop. \n";
		}
		ORBIT_MPI_Finalize();
  }		
	
	//check that grids are the same
  if( rhoGrid->getSizeX() != phiGrid->getSizeX() || 
			rhoGrid->getSizeY() != phiGrid->getSizeY() ||
		  rhoGrid->getStepX() != phiGrid->getStepX() ||
			rhoGrid->getStepY() != phiGrid->getStepY() ||
			rhoGrid->getMinX()  != phiGrid->getMinX()  ||
			rhoGrid->getMinY()  != phiGrid->getMinY()){
		int rank = 0;
		ORBIT_MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		if(rank == 0){
			std::cerr << "BaseBoundary2D (or subclass) addBoundaryPotential(...):" 
			<< "The rho and phi grid sizes are different "<< std::endl 
								<< "rhoGrid x bins ="<< rhoGrid->getSizeX() <<std::endl
								<< "rhoGrid y bins ="<< rhoGrid->getSizeY() <<std::endl
								<< "phiGrid x bins ="<< phiGrid->getSizeX() <<std::endl
								<< "phiGrid y bins ="<< phiGrid->getSizeY() <<std::endl
								<< "rhoGrid dx ="<< rhoGrid->getStepX() <<std::endl
								<< "rhoGrid dy ="<< rhoGrid->getStepY() <<std::endl
								<< "phiGrid dx ="<< phiGrid->getStepX() <<std::endl
								<< "phiGrid dy ="<< phiGrid->getStepY() <<std::endl
								<< "rhoGrid xMin ="<< rhoGrid->getMinX() <<std::endl
								<< "rhoGrid yMin ="<< rhoGrid->getMinY() <<std::endl
								<< "phiGrid xMin ="<< phiGrid->getMinX() <<std::endl
								<< "phiGrid yMin ="<< phiGrid->getMinY() <<std::endl
								<< "Stop. \n";
		}
		ORBIT_MPI_Finalize();
  }	
		
	for(int j = 0; j < nPoints_ ; j++) {
		bnd_phi_arr_[j] = phiGrid->getValue(bArrX_[j],bArrY_[j]);
	}
	
  for(int i = 0, n = (2*nModes_+1); i < n; i++) {  
    lsq_coef_vctr_[i] = 0.0;
     for(int j = 0; j < nPoints_ ; j++) {
       lsq_coef_vctr_[i] += LSQ_matrix_[i][j]*bnd_phi_arr_[j];
     } 
		 //std::cout<<"debug i="<<i<<" coef="<<lsq_coef_vctr_[i]<<std::endl;
  }
	
	double phi = 0.;
	double** phi_arr = phiGrid->getArr();
	for(int i = 0, nx = phiGrid->getSizeX(); i < nx; i++){
		for(int j = 0, ny = phiGrid->getSizeY(); j < ny; j++){
			phi = 0.;
			//std::cout<<"debug ================= i="<<i<<"  j="<< j <<std::endl;
			lsq_fuctions(phiGrid->getGridX(i),phiGrid->getGridY(j));
			for(int k = 0, n = (2*nModes_+1); k < n; k++){
				phi += lsq_coef_vctr_[k]*lsq_func_vctr_[k];
				//std::cout<<"debug k="<<k<<" coef="<<lsq_coef_vctr_[k]<<" func="<<lsq_func_vctr_[k]<<std::endl;
			}
			phi_arr[i][j] -=  phi;
			//std::cout<<"debug phi="<<phi<<std::endl;
		}
	}	
}


