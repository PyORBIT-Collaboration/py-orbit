#include "ShapedBoundary2D.hh"

#include <iostream>
#include <cfloat>

const int ShapedBoundary2D::IS_INSIDE    =  1;
const int ShapedBoundary2D::IS_OUTSIDE   = -1;
const int ShapedBoundary2D::TO_BE_KILLED =  0;
const double ShapedBoundary2D::PI = 3.14159265358979324;

using namespace OrbitUtils;

// Constructor
ShapedBoundary2D::ShapedBoundary2D(int nPoints, int nModes, string shape, double xDim, double yDim): 
           BaseBoundary2D(4*((int)(nPoints/4)),nModes)
{
	nPoints = getNumberOfPoints();
	
	shape_ = shape;
	xDim_ = xDim;
	yDim_ = yDim;
	
	r_circle_ = 0.;
  a_ellipse_ = 0.;
	b_ellipse_ = 0.;
	a_rect_ = 0.;
	b_rect_ = 0.;	
	
	shape_type_ = -1;
	if( shape_ == "Circle") {shape_type_ = 1; r_circle_ = xDim/2; }
	if( shape_ == "Ellipse") {shape_type_ = 2; a_ellipse_ = xDim/2; b_ellipse_ = yDim/2;}
	if( shape_ == "Rectangle") {shape_type_ = 3; a_rect_ = xDim; b_rect_ = yDim; }	

	if(shape_type_ < 0){
		int rank = 0;
		ORBIT_MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		if(rank == 0){
			std::cerr << "ShapedBoundary2D::ShapedBoundary2D - CONSTRUCTOR " << std::endl
			<< " Boundary class will not work because the shape parameter should be "<< std::endl
			<< " Circle or Ellipse or Rectangle"<< std::endl
			<< "Stop."<< std::endl;
		}
		ORBIT_MPI_Finalize();			
	}
	
	if(shape_type_ == 1 && xDim_ != yDim_){
		int rank = 0;
		ORBIT_MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		if(rank == 0){
			std::cerr << "ShapedBoundary2D::ShapedBoundary2D - CONSTRUCTOR " << std::endl
			<< " The Circle should have the same x and y sizes. "<< std::endl
			<< " xDimension = "<< xDim_<< std::endl
			<< " yDimension = "<< yDim_<< std::endl
			<< "Stop."<< std::endl;
		}
		ORBIT_MPI_Finalize();			
	}
	
	//after initializing we have to put this key to 0
	no_shape_key_ = 1;
	
	//initialize circle 
	if(shape_type_ == 1){
		for(int i = 0; i < nPoints; i++){
			double x = r_circle_*sin(2*i*PI/nPoints);
			double y = r_circle_*cos(2*i*PI/nPoints);
			setBoundaryPoint(i,x,y);
		}
	} 
	
	//initialize ellipse
	if(shape_type_ == 2){
		double ds,denom,resid,th,dtheta;
		ds = 2 * PI * sqrt(a_ellipse_*b_ellipse_)/nPoints;
		resid = 1.0;
		double* theta_arr = new double[nPoints];
		while (resid > 1.0e-08 || resid < -1.0e-08){
			denom = b_ellipse_;
			dtheta = ds / denom;
			th = dtheta / 2.;
			theta_arr[0] = 0.;
      for (int i = 0; i < nPoints/4 ; i++){
        denom = b_ellipse_*b_ellipse_ + (a_ellipse_*a_ellipse_-b_ellipse_*b_ellipse_)*sin(th)*sin(th);
        denom = sqrt(denom);
        dtheta = ds / denom;
        theta_arr[i+1] = theta_arr[i] + dtheta;
        th += dtheta;      
      }
			resid = theta_arr[nPoints/4] - PI/2.; 
			ds *= PI / (2.*theta_arr[nPoints/4]); 
		}
		
		for (int i = 0; i < nPoints/4; i++){
			theta_arr[nPoints/2-i] = PI - theta_arr[i];
		}
		
		for (int i = 0; i < nPoints/2; i++){
			theta_arr[nPoints/2 + i] = PI + theta_arr[i];
		}
		
		for (int i = 0; i < nPoints; i++){
			double x = a_ellipse_ * cos(theta_arr[i]);
			double y = b_ellipse_ * sin(theta_arr[i]);
			setBoundaryPoint(i,x,y);		 
		}	
		delete [] theta_arr;
	}
	
	//initialize rectangle
	if(shape_type_ == 3){
		double dx = a_rect_/(nPoints/4);
		double dy = b_rect_/(nPoints/4);
		for (int i = 0; i < nPoints/4; i++){
			double x = (-a_rect_/2.0) + i*dx;
			double y = b_rect_/2.0;		
			setBoundaryPoint(i,x,y);	
		}
		for (int i = 0; i < nPoints/4; i++){
			double x = a_rect_/2.0;
			double y = b_rect_/2.0 - i*dy;		
			setBoundaryPoint(i,x,y);	
		}
		for (int i = 0; i < nPoints/4; i++){
			double x = a_rect_/2.0 - i*dx;
			double y = -b_rect_/2.0;		
			setBoundaryPoint(i,x,y);	
		}
		for (int i = 0; i < nPoints/4; i++){
			double x = -a_rect_/2.0;
			double y = (-b_rect_/2.0) + i*dy;		
			setBoundaryPoint(i,x,y);	
		}
	}
	
	//end of initialization
	initializeBPs();
	no_shape_key_ = 0;
}

// Destructor
ShapedBoundary2D::~ShapedBoundary2D()
{
}


