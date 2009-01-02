//This class is a almost perfect copy of the Boundary2D class of the
//original ORBIT. 
//Authors: Jeff Holmes, Slava Danilov, John Galambos, Y.Sato, A.Shishlo

#include "Boundary2D.hh"

#include <iostream>

using namespace OrbitUtils;

//xBins, yBins - grid size [xBins]X[yBins]
//xSize, ySize - geometry parameters [m] or [cm] or [mm] 
//BPoints - number of points on the bounary
//BShape - shape of the boundary: Circle - Ellipse - Rectangle
//BModes - number of functions in the LSQ method ( > 10 )

const int Boundary2D::IS_INSIDE    =  1;
const int Boundary2D::IS_OUTSIDE   = -1;
const int Boundary2D::TO_BE_KILLED =  0;
const double Boundary2D::PI = 3.14159265358979324;

// Constructor
Boundary2D::Boundary2D(int xBins,    int yBins, 
			                 double xSize, double ySize,
			                 int BPoints, const string& BShape, int BModes):
            CppPyWrapper(NULL)
{
  init(xBins,yBins,xSize,ySize,BPoints,BShape,BModes);
}

void Boundary2D::init(int xBins, int yBins, 
			       double xSize, double ySize,
			       int BPoints, const string& BShape, 
			       int BModes)
{
  int i;

  xSize_ = xSize; 
  ySize_ = ySize;
 
  //====== Grid parameters and Green function==============

  xBins_ = xBins;
  yBins_ = yBins;

  xBins2_ = 2*xBins;
  yBins2_ = 2*yBins;

  if( xBins_ < 3 || yBins_ < 3){
		int rank = 0;
		ORBIT_MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		if(rank == 0){
			std::cerr << "Boundary2D::Boundary2D - CONSTRUCTOR \n" 
         				<< "The grid size too small (should be more than 3)! \n" 
								<< "number x bins ="<< xBins_ <<" \n"
								<< "number y bins ="<< yBins_ <<" \n"
								<< "Stop. \n";
		}
		ORBIT_MPI_Finalize();
  }

  greensF_ = new double*[xBins2_];
  for( i = 0; i < xBins2_ ; i++) {
    greensF_[i] =  new double [yBins2_];
  }
	
  in_        = (double *) fftw_malloc(sizeof(double)*xBins2_ * yBins2_);
  in_res_    = (double *) fftw_malloc(sizeof(double)*xBins2_ * yBins2_);
  out_green_ = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) *xBins2_ * (yBins2_/2+1));
  out_       = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) *xBins2_ * (yBins2_/2+1));
  out_res_   = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) *xBins2_ * (yBins2_/2+1));

	// FFTW_MEASURE or FFTW_ESTIMATE

  planForward_greenF_ = fftw_plan_dft_r2c_2d(xBins2_ , yBins2_ , in_,  out_green_, FFTW_ESTIMATE);
  planForward_        = fftw_plan_dft_r2c_2d(xBins2_ , yBins2_ , in_,  out_,       FFTW_ESTIMATE);
  planBackward_       = fftw_plan_dft_c2r_2d(xBins2_ , yBins2_ , out_res_, in_res_,FFTW_ESTIMATE);

  //=================================================================
  //================   Border parameters      =======================
  //=================================================================
  BShape_ = -1;
	if( BShape == "Circle") {BShape_ = 1;}
	if( BShape == "Ellipse") {BShape_ = 2;}
	if( BShape == "Rectangle") {BShape_ = 3;}
	
  if( BShape_ < 0){
		int rank = 0;
		ORBIT_MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		if(rank == 0){
			std::cerr << "Boundary2D::Boundary2D - CONSTRUCTOR"<<std::endl 
         			  << "Define the shape!" <<std::endl
                << "Stop."<<std::endl;
		}
		ORBIT_MPI_Finalize();
  }

  BPoints_ = 4* int(BPoints/4);

  BPx_ = new double[BPoints_];
  BPy_ = new double[BPoints_];
  theta_ = new double[BPoints_];
  BPphi_ = new double[BPoints_];


  //BShape_ = 1 - this is the circle =======START=============
  if (BShape_ == 1){
    if(xSize != ySize){
			int rank = 0;
			ORBIT_MPI_Comm_rank(MPI_COMM_WORLD, &rank);
			if(rank == 0){
				std::cerr << "Boundary2D::Boundary2D - CONSTRUCTOR "<< std::endl
				<< "It is not a circle :  xSize != ySize "<< std::endl
				<< "xSize and  ySize = " << xSize <<" "<< ySize << std::endl
				<< "Stop."<< std::endl;
			}
			ORBIT_MPI_Finalize();
    }

    R_cirle_ = xSize/2.0;
    BPrnorm_ = R_cirle_;

    double dtheta =2.0*PI/BPoints_;
    for (i = 0; i < BPoints_; i++){
      theta_[i] = i*dtheta;
      BPx_[i] = R_cirle_*cos(theta_[i]);
      BPy_[i] = R_cirle_*sin(theta_[i]);
    }
  }
  //BShape_ = 1 - this is the circle =======STOP==============

  //BShape_ = 2 - this is the Ellipse ======START=============
  if (BShape_ == 2){
   double ds,denom,resid,th,dtheta;
   int BPP04;
   BPP04 = BPoints_/4;
   BPa_ = xSize/2.0;
   BPb_ = ySize/2.0;
   BPrnorm_ = sqrt(BPa_ * BPb_);
   ds = 2 * PI * BPrnorm_ / BPoints_;
   resid = 1.0;
   while (resid > 1.0e-08 || resid < -1.0e-08){
     denom = BPb_;
     dtheta = ds / denom;
     th = dtheta / 2.;
     theta_[0] = 0.;
      for (i = 0; i < BPP04 ; i++){
        denom = BPb_*BPb_ + (BPa_*BPa_-BPb_*BPb_)*sin(th)*sin(th);
        denom = sqrt(denom);
        dtheta = ds / denom;
        theta_[i+1] = theta_[i] + dtheta;
        th += dtheta;      
      }
    resid = theta_[BPP04] - PI/ 2.; 
    ds *= PI / (2.*theta_[BPP04]); 
   }

   for (i = 0; i < BPP04; i++){
    theta_[2*BPP04-i] = PI - theta_[i];
   }

   for (i = 0; i < 2*BPP04; i++){
    theta_[2* BPP04 + i] = PI + theta_[i];
   }
               
   for (i = 0; i < BPoints_; i++){
     BPx_[i]   = BPa_ * cos(theta_[i]);
     BPy_[i]   = BPb_ * sin(theta_[i]);
     theta_[i] = atan2(BPy_[i],BPx_[i]);
   }
  }
  //BShape_ = 2 - this is the Ellipse ======STOP==============

  //BShape_ = 3 - this is the Rectangle ====START=============
  if (BShape_ == 3){
    double x_min,x_max,y_min,y_max,dx,dy;
    int BPPO4;
    BPx_length_ = xSize;
    BPy_width_  = ySize;
    BPPO4 = BPoints_/4;
    x_max = xSize/2.0;
    x_min = - x_max;   
    y_max = ySize/2.0;
    y_min = - y_max; 
    dx = (x_max - x_min)/BPPO4;
    dy = (y_max - y_min)/BPPO4;
    BPrnorm_ = sqrt((xSize*xSize + ySize*ySize)/4.);
     for (i = 0; i < BPPO4; i++){
       BPx_[i] = x_max;
       BPy_[i] = y_min + i*dy;
       theta_[i] = atan2(BPy_[i],BPx_[i]);
       
       BPx_[BPPO4 + i] = x_max - i*dx;
       BPy_[BPPO4 + i] = y_max;
       theta_[BPPO4 + i] = atan2(BPy_[BPPO4 + i],BPx_[BPPO4 + i]);

       BPx_[2*BPPO4 + i] = x_min;
       BPy_[2*BPPO4 + i] = y_max - i*dy;
       theta_[2*BPPO4 + i] = atan2(BPy_[2*BPPO4 + i],BPx_[2*BPPO4 + i]);

       BPx_[3*BPPO4 + i] = x_min + i*dx;
       BPy_[3*BPPO4 + i] = y_min;
       theta_[3*BPPO4 + i] = atan2(BPy_[3*BPPO4 + i],BPx_[3*BPPO4 + i]);
     }
  }

  //BShape_ = 3 - this is the Rectangle ====STOP==============

  double xBrMin = -xSize_/2.0;
  double xBrMax =  xSize_/2.0;
  double yBrMin = -ySize_/2.0;
  double yBrMax =  ySize_/2.0;

  for (i = 0; i < BPoints_; i++){ 
    if( xBrMin > BPx_[i] ) xBrMin = BPx_[i];
    if( xBrMax < BPx_[i] ) xBrMax = BPx_[i];
    if( yBrMin > BPy_[i] ) yBrMin = BPy_[i];
    if( yBrMax < BPy_[i] ) yBrMax = BPy_[i];
  }

  //extend factor for the grid to use FFT-method for a non-periodic SC-distribution
  //we added 3 factor to increase the size of the grid, so boundary will not touch
  //the grid cells at the edge
  double GridFactorX = 1.0*(xBins_/2 + 3.0)/(xBins_/2);
  double GridFactorY = 1.0*(yBins_/2 + 3.0)/(yBins_/2);

  xGridMin_ = xBrMin*GridFactorX;
  xGridMax_ = xBrMax*GridFactorX;
  yGridMin_ = yBrMin*GridFactorY;
  yGridMax_ = yBrMax*GridFactorY;

  dx_ = (xGridMax_ - xGridMin_)/(xBins_-1);
  dy_ = (yGridMax_ - yGridMin_)/(yBins_-1);

  xGrid_ = new double[xBins_];
  yGrid_ = new double[yBins_];

  for (i = 0; i< xBins_; i++) {xGrid_[i] = xGridMin_ + i * dx_;}
  for (i = 0; i< yBins_; i++) {yGrid_[i] = yGridMin_ + i * dy_;}

  //define the LSQM fuction and matrix 

  BModes_ = BModes;

  if( 2*BModes_+1 >= BPoints_ ){
		int rank = 0;
		ORBIT_MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		if(rank == 0){
			std::cerr << "Boundary2D::Boundary2D - CONSTRUCTOR " << std::endl
			<< " Boundary class will not work because "<< std::endl
			<< " 2*BModes_+1 = " << 2*BModes_+1 << std::endl
			<< " BPoints_    = " << BPoints_    << std::endl
			<< " Number of points has to be increased."<<std::endl
			<< "Stop."<< std::endl;
		}
		ORBIT_MPI_Finalize();		
  }

  func_vector_ = new double[2*BModes_+1];
  cof_vector_  = new double[2*BModes_+1];
  for(int i = 0, n =  (2*BModes_+1); i < n; i++) {  
    cof_vector_[i] = 0.0;
		func_vector_[i] = 0.0;
	}
	
	indxc_tmp= new int[2*BModes_+1];
	indxr_tmp= new int[2*BModes_+1];
	ipiv_tmp = new int[2*BModes_+1];

  tmp_matrix_ =new double*[2*BModes_+1];
  for( i = 0; i < (2*BModes_+1); i++) {
    tmp_matrix_[i] =  new double [2*BModes_+1];
  } 

  LSQ_matrix_ =new double*[2*BModes_+1];
  for( i = 0; i < (2*BModes_+1); i++) {
    LSQ_matrix_[i] =  new double [BPoints_];
  }

  //array with the interpolation coefficients for boundary points
  //there 9-points scheme is used
  W_coeff_ = new double* [BPoints_];
  for( i = 0; i < BPoints_; i++) {  
   W_coeff_[i] = new double [9];
  }

  //indexes of the boundary points
  iBPx_ = new int [BPoints_];
  iBPy_ = new int [BPoints_];
  
  //define FFT of the Green fuction
  _defineGreenF();

  //define tmp_matrix_
  _defineLSQMatrix();

  //define the min and max indexes of XY-plane's grid to operate with potential
  int iX, iY;
  double temp_frac = 0.0;

  ixMinB_ = xBins_;
  iyMinB_ = yBins_;

  ixMaxB_ = 0;
  iyMaxB_ = 0; 
 
  for (i = 0; i < BPoints_; i++){
    getIndAndFracX(BPx_[i],iX,temp_frac);
    getIndAndFracY(BPy_[i],iY,temp_frac);
    iBPx_[i] = iX;
    iBPy_[i] = iY;
    if(ixMinB_ > iX ){ ixMinB_ = iX;}
    if(ixMaxB_ < iX ){ ixMaxB_ = iX;}
    if(iyMinB_ > iY ){ iyMinB_ = iY;}
    if(iyMaxB_ < iY ){ iyMaxB_ = iY;}
  }
  ixMinB_ = ixMinB_ -1;
  ixMaxB_ = ixMaxB_ +1;
  iyMinB_ = iyMinB_ -1;
  iyMaxB_ = iyMaxB_ +1;

  //Sets array with the interpolation coefficients for boundary points
  _setInterpolationCoeff();
}

// Destructor
Boundary2D::~Boundary2D()
{
	
	//std::cerr<<"debug Boundary2D::~Boundary2D() start! "<<std::endl;
   int i;
  //delete Green function and FFT input and output arrays

  for( i = 0; i < xBins2_ ; i++) {
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

  //delete arrays describing boundary
  delete [] BPx_ ;
  delete [] BPy_ ;
  //indexes of the boundary points
  delete [] iBPx_;
  delete [] iBPy_;

  delete [] theta_ ;
  delete [] BPphi_;

  //delete grid
  delete [] xGrid_;
  delete [] yGrid_;

  //delete LSQM related arrays
  delete [] func_vector_;
  delete [] cof_vector_;

	delete [] ipiv_tmp;
	delete [] indxr_tmp;
	delete [] indxc_tmp;
	
  for( i = 0; i < (2*BModes_+1); i++) {
    delete [] tmp_matrix_[i];
  }
  delete [] tmp_matrix_;

  for( i = 0; i < (2*BModes_+1); i++) {
    delete [] LSQ_matrix_[i];
  }
  delete [] LSQ_matrix_;

  //array with the interpolation coefficients for boundary points
  for( i = 0; i < BPoints_; i++) {  
   delete [] W_coeff_[i];
  }
  delete [] W_coeff_;  

  //std::cerr << "debug Destructor Boundary2D was done.   !!!! \n";
}


/** The method calculates an impact position on the surface 
		and a normal vector at the point of entry. 
		The normal vector is directed into the inner volume. 
		The first vector is a position vector and the second one is a 
		normal to the surface vector.
		If the particle did not cross the surface we do not know what to do
		and we will kill it (return TO_BE_KILLED int value).
*/
int Boundary2D::impactPoint(double x,  double y,  double z,
	                          double px, double py, double pz,
														double* r_v,double* n_v)
{
	// NOTICE: in the rectangle case, we could remove corner and boundary plane 
	// conditions if we found out they would not change results so much.	
  double time_impact=-1.0;

  r_v[0]=x;  r_v[1]=y;  r_v[2]=z;
  n_v[0]=0.0;  n_v[1]=0.0;  n_v[2]=0.0;
	
  //we consider the particle outside of the boundary
  int isIns = isInside(x,y);
  if (isIns < 0){
		
    //Circle---------------------------------------------------------------
    if (BShape_ == 1){
      double c1=px*px+py*py, c2=x*px+y*py, c3=x*x+y*y-pow(R_cirle_,2);
      if (c2*c2-c1*c3 <= 0.0){return TO_BE_KILLED;}
      double t1=(c2+sqrt(c2*c2-c1*c3))/c1, t2=(c2-sqrt(c2*c2-c1*c3))/c1;
      if(t1>=0.0 && t2>=0.0){
				time_impact=t2;
      }else{
				//std::cout <<"the particle may have inward momentum \n";
				//std::cout <<"Please check '"<<t1*t2<<"' be positive \n";
				return TO_BE_KILLED;
      }
      double norm = sqrt(x*x+y*y);
      n_v[0] = -x/norm; n_v[1] = -y/norm; n_v[2] = 0.0;
    }
    //-----
		
    //Ellipse--------------------------------------------------------------
    if (BShape_ == 2){
      double c1 = pow(px/BPa_,2) + pow(py/BPb_,2);
      double c2 = x*px/pow(BPa_,2) + y*py/pow(BPb_,2);
      double c3 = pow(x/BPa_,2) + pow(y/BPb_,2) -1;
      if (c2*c2-c1*c3 <= 0.0){return TO_BE_KILLED;}
      double t1=(c2+sqrt(c2*c2-c1*c3))/c1, t2=(c2-sqrt(c2*c2-c1*c3))/c1;
      if(t1>=0.0 && t2>=0.0) {
				time_impact=t2;
      }else{
				//std::cout <<"the particle may have inward momentum \n";
				//std::cout <<"Please check '"<<t1*t2<<"' be positive \n";
				return TO_BE_KILLED;
      }
      double norm = sqrt(pow(BPb_,4)*x*x + pow(BPa_,4)*y*y);
      n_v[0] = -pow(BPb_,2)*x/norm; n_v[1] = -pow(BPa_,2)*y/norm; n_v[2] = 0.0;
    }
    //-----
		
    //Rectangle------------------------------------------------------------
    if (BShape_ == 3){
      double L = BPx_length_,W = BPy_width_; 
			
      double t1,t2,t3,t4, ry1,ry2,rx3,rx4;
      t1= (x - 0.5*L)/px;      ry1=y-t1*py;
      t2= (x + 0.5*L)/px;      ry2=y-t2*py;
      t3= (y - 0.5*W)/py;      rx3=x-t3*px;
      t4= (y + 0.5*W)/py;      rx4=x-t4*px;
			
      //      std::cout<<"debug (t1,t2,t3,t4)= ("
      //	       <<t1<<","<<t2<<","<<t3<<","<<t4<<") \n";
			
      double eps=1.E-10; //to define corners and planes of the boundary 
			
      //the particle is at a corner
      if( fabs(fabs(x)-0.5*L)<eps && fabs(fabs(y)-0.5*W)<eps ){
				return TO_BE_KILLED;
      }
      //the particle is on a boundary plane
      else if( fabs(t1)<eps && fabs(ry1)<0.5*W && px>0.0 ){time_impact =t1;}
      else if( fabs(t2)<eps && fabs(ry2)<0.5*W && px<0.0 ){time_impact =t2;}
      else if( fabs(t3)<eps && fabs(rx3)<0.5*L && py>0.0 ){time_impact =t3;}
      else if( fabs(t4)<eps && fabs(rx4)<0.5*L && py<0.0 ){time_impact =t4;}
      else if( fabs(t1)<eps && fabs(ry1)<0.5*W && px<0.0 ){return IS_INSIDE;}
      else if( fabs(t2)<eps && fabs(ry2)<0.5*W && px>0.0 ){return IS_INSIDE;}
      else if( fabs(t3)<eps && fabs(rx3)<0.5*L && py<0.0 ){return IS_INSIDE;}
      else if( fabs(t4)<eps && fabs(rx4)<0.5*L && py>0.0 ){return IS_INSIDE;}
      else{
				//the particle locates outside boundary
				if(t1>0.0 && fabs(ry1)<0.5*W){time_impact=t1;}
				if(t2>0.0 && fabs(ry2)<0.5*W){
					if(time_impact<=0.0){       time_impact=t2;}
					if(time_impact > t2){       time_impact=t2;}
				}
				if(t3>0.0 && fabs(rx3)<0.5*L){
					if(time_impact<=0.0){       time_impact=t3;} 
					if(time_impact > t3){       time_impact=t3;}
				}
				if(t4>0.0 && fabs(rx4)<0.5*L){ 
					if(time_impact<=0.0){       time_impact=t4;}
					if(time_impact > t4){       time_impact=t4;}
				}	  
				
				//the particle came from a corner
				if( ( fabs(time_impact-t1)<eps && fabs(time_impact-t3)<eps ) ||
					( fabs(time_impact-t1)<eps && fabs(time_impact-t4)<eps ) ||
				( fabs(time_impact-t2)<eps && fabs(time_impact-t3)<eps ) ||
				( fabs(time_impact-t2)<eps && fabs(time_impact-t4)<eps ) ){
				return TO_BE_KILLED;
				}
      }
			
      //      std::cout<<"debug time_impact= ("<<time_impact<<") \n";
			
      if(time_impact<0.0){
				// std::cout<<"the particle has inward momentum \n";
        return TO_BE_KILLED;
      }
			
      if(fabs(time_impact-t1)<eps){n_v[0]=-1.0; n_v[1]= 0.0; n_v[2]= 0.0;}
      if(fabs(time_impact-t2)<eps){n_v[0]= 1.0; n_v[1]= 0.0; n_v[2]= 0.0;}
      if(fabs(time_impact-t3)<eps){n_v[0]= 0.0; n_v[1]=-1.0; n_v[2]= 0.0;}
      if(fabs(time_impact-t4)<eps){n_v[0]= 0.0; n_v[1]= 1.0; n_v[2]= 0.0;}
			
    }
    //-----
		
    //positive "time_impact" means that the particle goes outward 
    r_v[0] = x - time_impact*px; 
    r_v[1] = y - time_impact*py; 
    r_v[2] = z - time_impact*pz;
    return IS_OUTSIDE;
  }
  else{
    return IS_INSIDE;
  }

}

int Boundary2D::getBinsX(){return xBins_;}; 

int Boundary2D::getBinsY(){return yBins_;};

// get the index and the fraction of the grid's point for X coord.
void Boundary2D::getIndAndFracX(double x, int& ind, double& frac)
{
   ind  = int ( (x - xGridMin_)/dx_ + 0.5 );
   if(ind < 1) ind = 1;
   if(ind > (xBins_-2)) ind = xBins_ - 2;
   frac = (x - xGrid_[ind])/dx_; 
}

// get the index and the fraction of the grid's point for Y coord.
void Boundary2D::getIndAndFracY(double y, int& ind, double& frac)
{
   ind  = int ( (y - yGridMin_)/dy_ + 0.5 );
   if(ind < 1) ind = 1;
   if(ind > (yBins_-2)) ind = yBins_ - 2;
   frac = (y - yGrid_[ind])/dy_; 
}

//returns xGridMax_, and so on  
double Boundary2D::getGridMaxX(){return xGridMax_;};
double Boundary2D::getGridMinX(){return xGridMin_;};
double Boundary2D::getGridMaxY(){return yGridMax_;};
double Boundary2D::getGridMinY(){return yGridMin_;};

//returns steps on X and Y axes
double Boundary2D::getStepX(){ return dx_;};
double Boundary2D::getStepY(){ return dy_;};

// Returns the pointers to the FFT array of the Green Function values
fftw_complex* Boundary2D::getOutFFTGreenF()
{
  return out_green_;
}

// Defines the FFT of the Green Function
void Boundary2D::_defineGreenF()
{
	
  double rTransY, rTransX, rTot2;
  int i, j, iY , iX;
	
	for (iY = 0; iY <= yBins2_/2; iY++)
	{
		rTransY = iY * dy_;
		
		for (iX = 0; iX <= xBins2_/2; iX++)
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
		
		for (iX = xBins2_/2+1; iX < xBins2_; iX++)
		{
			greensF_[iX][iY] = greensF_[xBins2_-iX][iY];
		}
	}
	
	for (iY = yBins2_/2+1; iY < yBins2_; iY++)
	{
		for (iX = 0; iX < xBins2_; iX++)
		{
			greensF_[iX][iY] = greensF_[iX][yBins2_-iY];
		}
	}
	
	//   Calculate the FFT of the Greens Function:
	
	for (i = 0; i < xBins2_; i++)
		for (j = 0; j < yBins2_; j++)
		{
      in_[j + yBins2_*i] = greensF_[i][j];
		}
    
		fftw_execute(planForward_greenF_);
		//rfftwnd_one_real_to_complex(planForward_, in_, out_green_);
		
		for (i = 0; i < xBins2_; i++)
			for (j = 0; j < yBins2_; j++)
			{
				in_[j + yBins2_*i] = 0.0;
			}
			
}

// Defines LSQM matrix
void Boundary2D::_defineLSQMatrix()
{

  int i,j;
  for( i = 0; i < (2*BModes_+1); i++) {
  for( j = 0; j < (2*BModes_+1); j++) {
    tmp_matrix_[i][j] =  0.0;
  }}
 
  int iBp;
  for ( iBp = 0; iBp < BPoints_ ; iBp++){
    func_vector_ =  lsq_fuctions( BPx_[iBp],BPy_[iBp]);
     
     for( i = 0; i < (2*BModes_+1); i++) {
     for( j = 0; j < (2*BModes_+1); j++) {
       tmp_matrix_[i][j] += func_vector_[i]*func_vector_[j];
     }}
  }

  //inverse matrix
  _gaussjinv(tmp_matrix_,2*BModes_+1);

  int k;
  for( i = 0; i < (2*BModes_+1); i++) {
  for( j = 0; j < BPoints_     ; j++) {
    LSQ_matrix_[i][j] = 0.0;
    func_vector_ =  lsq_fuctions( BPx_[j],BPy_[j]);
    for ( k = 0; k < (2*BModes_+1); k++) {
      LSQ_matrix_[i][j] += tmp_matrix_[i][k] * func_vector_[k];
    }    
  }} 
}


// Defines LSQM coefficients
void Boundary2D::_defineLSQcoeff()
{
  int i,j;
  for( i = 0; i < (2*BModes_+1); i++) {  
    cof_vector_[i] = 0.0;
     for( j = 0; j < BPoints_ ; j++) {
       cof_vector_[i] += LSQ_matrix_[i][j]*BPphi_[j];
     } 
  }
}

//Defines is a particle into the boundary or not ( returns < 0 - it is not)
int Boundary2D::isInside(double x, double y)
{
  int isIns = IS_OUTSIDE;

  if (BShape_ == 1){
    if( x*x+y*y < R_cirle_*R_cirle_ ) isIns = IS_INSIDE;
  }

  if (BShape_ == 2){
    double xx,yy;
    xx = x/BPa_;
    yy = y/BPb_;
    if( (xx*xx + yy*yy) < 1.0 ) isIns = IS_INSIDE;
  }

  if (BShape_ == 3){
    if( fabs(x/BPx_length_) < 0.5 && fabs(y/BPy_width_) < 0.5 ) isIns = IS_INSIDE;
  }

  return isIns;
}

// Calculates additional potential phi
double Boundary2D::calculatePhi(double x, double y)
{
  //cof_vector_ - should be calculated before
  int i;  
  double phi;
  phi = 0.0;
  func_vector_ =  lsq_fuctions( x, y);
  for( i = 0 ; i < (2*BModes_+1); i++) {
    phi += cof_vector_[i]*func_vector_[i];
  }
  return phi;
}


// Calculates all of the LSQM functions at one point
double* Boundary2D::lsq_fuctions(double x, double y)
{
  int i = 0;
  double r,r2;
  r2 = x*x + y*y;
  if( r2 == 0. ) {
    func_vector_[0] = 1.0; 
     for( i = 1; i < (2*BModes_+1); i++) {
       func_vector_[i] = 0.0;
     }
     return func_vector_;
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
  func_vector_[0] = 1.;
  for( i = 0; i < BModes_; i++) {
    rj *= rfac; 
    func_vector_[2*i+1]= rj*cos_f0;
    func_vector_[2*i+2]= rj*sin_f0;
    sin_f1 = sin_f*cos_f0 + cos_f*sin_f0;
    cos_f1 = cos_f*cos_f0 - sin_f*sin_f0;
    sin_f0 = sin_f1;
    cos_f0 = cos_f1;
  }
  return func_vector_;
}

// Adds potential from the boundary to the grid
void Boundary2D::addBoundaryPotential(double** phisc, 
                                         int iXmin, int iXmax, 
                                         int iYmin, int iYmax)
{
  int iB, iX, iY;
  for (iB = 0; iB < BPoints_; iB++)
  {
 
   iX = iBPx_[iB];
   iY = iBPy_[iB];

    BPphi_[iB] =  W_coeff_[iB][0] * phisc[iX-1][iY-1] +
                  W_coeff_[iB][1] * phisc[iX-1][iY]   +
                  W_coeff_[iB][2] * phisc[iX-1][iY+1] +
                  W_coeff_[iB][3] * phisc[iX]  [iY-1] +
                  W_coeff_[iB][4] * phisc[iX]  [iY]   +
                  W_coeff_[iB][5] * phisc[iX]  [iY+1] +
                  W_coeff_[iB][6] * phisc[iX+1][iY-1] +
                  W_coeff_[iB][7] * phisc[iX+1][iY]   +
                  W_coeff_[iB][8] * phisc[iX+1][iY+1]; 
  }

  //std::cout << "debug Boundary after point's potential calc. \n";

  //define LSQM coefficient
  _defineLSQcoeff();

  //std::cout << "debug Boundary after coeff. calc. \n";

  int ix , iy;
  double x , y, phi;
  for( ix = iXmin; ix <= iXmax; ix++) {
  for( iy = iYmin; iy <= iYmax; iy++) {
    x = xGrid_[ix];
    y = yGrid_[iy];
    phi = calculatePhi(x,y);
    phisc[ix][iy] -= phi;
  }}
  //std::cout << "debug Boundary exit. \n";
}

// Adds potential from the boundary to the grid for all points
void Boundary2D::addBoundaryPotential(double** phisc)
{
	addBoundaryPotential(phisc,0,xBins_ - 1, 0,yBins_ - 1);
}



void Boundary2D::findPotential(double** rhosc, double** phisc)
{
  int i, j, index;

  double denom = 1.0 / (xBins2_*yBins2_);

  //define the the rho for FFT
  for (i = 0; i < xBins_; i++)
  for (j = 0; j < yBins_; j++)
  {
    in_[j + yBins2_*i] = rhosc[i][j];
  }

	fftw_execute(planForward_);
  //rfftwnd_one_real_to_complex(planForward_, in_, out_);

  //do convolution with the FFT of the Green's function 

  for (i = 0; i < xBins2_; i++)
  for (j = 0; j < yBins2_/2+1; j++)
  {
    index = j + (yBins2_/2+1)*i;
    out_res_[index][0] = out_[index][0]*out_green_[index][0] - out_[index][1]*out_green_[index][1];
    out_res_[index][1] = out_[index][0]*out_green_[index][1] + out_[index][1]*out_green_[index][0];
  }

  //do backward FFT
	fftw_execute(planBackward_);
  //rfftwnd_one_complex_to_real(planBackward_, out_res_, in_res_);

  //set the potential
  for (i = 0; i < xBins_; i++)
  for (j = 0; j < yBins_; j++)
  {
    phisc[i][j] = denom * in_res_[j + yBins2_*i];
  }
}



//Sets array with the interpolation coefficients for boundary points
void Boundary2D::_setInterpolationCoeff()
{
  int iB, iX, iY;
  double xFract, yFract, Wxm, Wx0, Wxp, Wym, Wy0, Wyp;
  for (iB = 0; iB < BPoints_; iB++)
  {
    getIndAndFracX(BPx_[iB],iX,xFract);
    getIndAndFracY(BPy_[iB],iY,yFract);

     // TSC interpolation, see Hockney and Eastwood:
      
     Wxm = 0.5 * (0.5 - xFract) * (0.5 - xFract);
     Wx0 = 0.75 - xFract * xFract;
     Wxp = 0.5 * (0.5 + xFract) * (0.5 + xFract);
     Wym = 0.5 * (0.5 - yFract) * (0.5 - yFract);
     Wy0 = 0.75 - yFract * yFract;
     Wyp = 0.5 * (0.5 + yFract) * (0.5 + yFract);
     W_coeff_[iB][0] = Wxm * Wym;
     W_coeff_[iB][1] = Wxm * Wy0;
     W_coeff_[iB][2] = Wxm * Wyp;
     W_coeff_[iB][3] = Wx0 * Wym;
     W_coeff_[iB][4] = Wx0 * Wy0;
     W_coeff_[iB][5] = Wx0 * Wyp;
     W_coeff_[iB][6] = Wxp * Wym;
     W_coeff_[iB][7] = Wxp * Wy0;
     W_coeff_[iB][8] = Wxp * Wyp;
  }
}

//Define external grid's parameters from inner 
void Boundary2D::defineExtXYgrid(double &xGridMin, double &xGridMax, 
                                     double &yGridMin, double &yGridMax,
                                     double &dx      , double &dy)
{
  xGridMin = xGridMin_;
  xGridMax = xGridMax_;
  yGridMin = yGridMin_;
  yGridMax = yGridMax_;
  dx = dx_;
  dy = dy_;
}

//Define external grid's index limits from inner 
void Boundary2D::defineExtXYlimits(int &ixMinB, int &ixMaxB, 
                                       int &iyMinB, int &iyMaxB)
{
  ixMinB = ixMinB_;
  ixMaxB = ixMaxB_;
  iyMinB = iyMinB_;
  iyMaxB = iyMaxB_;
}

//Get X-grid
double* Boundary2D::getGridX(){return xGrid_;};

//Get Y-grid
double* Boundary2D::getGridY(){return yGrid_;};

//returns the index of the boundary shape
int Boundary2D::getBoundaryShape(){return BShape_;};

//Get the bounding curve limit for X-coordinate
double Boundary2D::getBoundarySizeX(){return xSize_;};

//Get the bounding curve limit for Y-coordinate
double Boundary2D::getBoundarySizeY(){return ySize_;};

//Get the number of boundary points
int Boundary2D::getNumbBPoints(){return BPoints_;};

//Get the X and Y coordinates of the boundary points
double Boundary2D::getBoundPointX(int i) {return BPx_[i];};
double Boundary2D::getBoundPointY(int i) {return BPy_[i];};

void Boundary2D::_gaussjinv(double **a, int n)
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
						<< "Stop."<<std::endl;
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
			if (a[icol][icol] == 0.0) exit(1);
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

//debug method
void Boundary2D::checkBoundary()
{

  //define potential at the boundary
  int j;
  for( j = 0; j < BPoints_ ; j++) {
    BPphi_[j] = sin ((BPx_[j]*BPx_[j]*BPy_[j]*BPy_[j]));
  }

  //define LSQM coefficient
  _defineLSQcoeff();

  //check calculatePhi(x,y) - method
  double phi,x,y;

  for( j = 0; j < BPoints_ ; j++) {
    x = BPx_[j];
    y = BPy_[j];
    phi = calculatePhi(x,y);
    std::cout << " x y = \t" << x << " \t" << y << " \t phi0  phi = \t" << BPphi_[j] << "\t " << phi << std::endl;
  }
}
