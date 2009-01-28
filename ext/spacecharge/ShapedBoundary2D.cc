#include "ShapedBoundary2D.hh"

#include <iostream>
#include <cfloat>

const int ShapedBoundary2D::IS_INSIDE    =  1;
const int ShapedBoundary2D::IS_OUTSIDE   = -1;
const int ShapedBoundary2D::TO_BE_KILLED =  0;

const double ShapedBoundary2D::PI = 3.14159265358979324;

using namespace OrbitUtils;

/** Constructor */
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

/** Destructor */
ShapedBoundary2D::~ShapedBoundary2D()
{
}

/** The method calculates an impact position on the surface 
and a normal vector at the point of entry. 
The normal vector is directed into the inner volume. 
The first vector is a position vector and the second one is a 
normal to the surface vector.
If the particle did not cross the surface we do not know what to do
and we will kill it (return TO_BE_KILLED int value).
*/
int ShapedBoundary2D::impactPoint(double x,  double y,  double z,
                                  double px, double py, double pz,
	                                double* r_v,double* n_v)
{
	// NOTICE: in the rectangle case, we could remove corner and boundary plane 
	// conditions if we found out they would not change results so much.	
  double time_impact=-1.0;
	
  r_v[0]=x;    r_v[1]=y;    r_v[2]=z;
  n_v[0]=0.0;  n_v[1]=0.0;  n_v[2]=0.0;
	
  //we consider the particle outside of the boundary
  int isIns = isInside(x,y);
  if (isIns == IS_OUTSIDE){
		
    //Circle-------------------------------------------------
    if (shape_type_ == 1){
      double c1=px*px+py*py, c2=x*px+y*py, c3=x*x+y*y-r_circle_*r_circle_;
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
    //-------------------------------------------------------
		
    //Ellipse------------------------------------------------
    if (shape_type_ == 2){
      double c1 = pow(px/a_ellipse_,2) + pow(py/b_ellipse_,2);
      double c2 = x*px/pow(a_ellipse_,2) + y*py/pow(b_ellipse_,2);
      double c3 = pow(x/a_ellipse_,2) + pow(y/b_ellipse_,2) -1;
      if (c2*c2-c1*c3 <= 0.0){return TO_BE_KILLED;}
      double t1=(c2+sqrt(c2*c2-c1*c3))/c1, t2=(c2-sqrt(c2*c2-c1*c3))/c1;
      if(t1>=0.0 && t2>=0.0) {
				time_impact=t2;
      }else{
				//std::cout <<"the particle may have inward momentum \n";
				//std::cout <<"Please check '"<<t1*t2<<"' be positive \n";
				return TO_BE_KILLED;
      }
      double norm = sqrt(pow(b_ellipse_,4)*x*x + pow(a_ellipse_,4)*y*y);
      n_v[0] = -pow(b_ellipse_,2)*x/norm; n_v[1] = -pow(a_ellipse_,2)*y/norm; n_v[2] = 0.0;
    }
    //-------------------------------------------------------
		
    //Rectangle----------------------------------------------
    if (shape_type_ == 3){
      double L = a_rect_,W = b_rect_; 
			
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
    //-------------------------------------------------------
		
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

/** Returns IS_INSIDE or IS_OUTSIDE depending on the particle's position */
int ShapedBoundary2D::isInside(double x, double y)
{
  int isIns = IS_OUTSIDE;

  if (shape_type_ == 1){
    if( x*x+y*y < r_circle_*r_circle_ ) isIns = IS_INSIDE;
  }

  if (shape_type_ == 2){
    double xx,yy;
    xx = x/a_ellipse_;
    yy = y/b_ellipse_;
    if( (xx*xx + yy*yy) < 1.0 ) isIns = IS_INSIDE;
  }

  if (shape_type_ == 3){
    if( fabs(x/a_rect_) < 0.5 && fabs(y/b_rect_) < 0.5 ) isIns = IS_INSIDE;
  }

  return isIns;
}


