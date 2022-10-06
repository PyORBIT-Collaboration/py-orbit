#include "orbit_mpi.hh"

#include "ConvexApertureShape.hh"

#include <iostream>
#include <cmath>
#include <cfloat>
#include <cstdlib>

///////////////////////////////////////////////////////////////////////////
//
// NAME
//
//   ConvexApertureShape
//
// AUTHOR: 
//   Andrei Shishlo October 2022
//
//   ConvexApertureShape is an implementation of BaseApertureShape class
//   and a collection of (x,y) points describing convex shape. User should 
//   define these points.
///////////////////////////////////////////////////////////////////////////

/** ConvexApertureShape constructor */
ConvexApertureShape::ConvexApertureShape(): BaseApertureShape()
{
		shapeName = "convex";
		typeName = "convex";
}

/** ConvexApertureShape decstructor */
ConvexApertureShape::~ConvexApertureShape()
{
}

/** Return 1 if the particular macro-particle is inside this shape */
int ConvexApertureShape::inside(Bunch* bunch, int count){
	
	double** coord = bunch->coordArr();
	
	double x = coord[count][0] - x_center;
	double y = coord[count][2] - y_center;
	
	int nPoints = convexX.size();
	for(int ind = 0; ind < (nPoints-1); ind++){
		if(this->checkOnePoint(ind,ind+1,x,y) != 1){
			return 0;
		}
	}
	if(nPoints == 2) return 1;
	return this->checkOnePoint(nPoints-1,0,x,y);
}

/** Adds the new aperture shape to the collection */ 
void ConvexApertureShape::addPoint(double x, double y)
{
	convexX.push_back(x);
	convexY.push_back(y);
}

/** Returns vector of x-coordinates of points */ 
std::vector<double> ConvexApertureShape::getPointsX()
{
	return convexX;
}

/** Returns vector of y-coordinates of points */ 
std::vector<double> ConvexApertureShape::getPointsY()
{
	return convexY;
}

/** Removes all (x,y) points */ 
void ConvexApertureShape::removeAllPoints()
{
	convexX.clear();
	convexY.clear();
}

/** Checks that we have a convex shape*/ 
int ConvexApertureShape::checkAllPoints()
{
	int nPoints = convexX.size();
	if(nPoints == 0 || nPoints == 1){
		ORBIT_MPI_Finalize("ConvexApertureShape.checkAllPoints() - number of convex shape points = 0 or 1. Stop.");
		return 0;
	}
	if(nPoints == 2){
		return 1;
	}

	
	for(int ind = 0; ind < (nPoints-2); ind++){
		if(this->checkOnePoint(ind,ind+1,ind+2) != 1){
			ORBIT_MPI_Finalize("ConvexApertureShape.checkAllPoints() - not a convex shape. Stop.");
			return 0;				
		}
	}
	
	if(this->checkOnePoint(nPoints - 2, nPoints - 1, 0) != 1){
		ORBIT_MPI_Finalize("ConvexApertureShape.checkAllPoints() - not a convex shape. Stop.");
		return 0;
	}
	
	if(this->checkOnePoint(nPoints - 1, 0, nPoints - 2) != 1){
		ORBIT_MPI_Finalize("ConvexApertureShape.checkAllPoints() - not a convex shape. Stop.");
		return 0;
	}
	return 1;
}

/** Checks that a point with index i that it is at clockwise half-plain relative to 0->1 line */ 
int ConvexApertureShape::checkOnePoint(int i0, int i1, int i)
{
	double x = convexX[i];
	double y = convexY[i];
	return this->checkOnePoint(i0,i1,x,y);
}

/** Checks that a point (x,y) that it is at clockwise half-plain relative to 0->1 line */ 
int ConvexApertureShape::checkOnePoint(int i0, int i1, double x, double y)
{
	//double res = ((convexX[i1] - convexX[i0])*(y - convexY[i0]) -(x - convexX[i0])*(convexY[i1] - convexY[i0]));
	//std::cout<<"debug x,y="<< x <<"  "<< y <<" i0="<< i0 <<" i1="<< i1 <<"  res="<< res <<std::endl;
	if( ((convexX[i1] - convexX[i0])*(y - convexY[i0]) - 
		  (x - convexX[i0])*(convexY[i1] - convexY[i0])) < 0.){
		return 1;
	}
	return 0;
}

