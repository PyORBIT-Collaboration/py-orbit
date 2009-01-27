#ifndef SC_SHAPED_BOUNDARY_2D_H
#define SC_SHAPED_BOUNDARY_2D_H

#include <string>

#include "Grid2D.hh"
#include "BaseBoundary2D.hh"

using namespace std;

/** 
The ShapedBoundary2D class defines a boundary geometry for three cases:
circle, ellipse, and rectangular shape
and calculates the potential created by charges on the boundary 
surface.
*/

class ShapedBoundary2D: public BaseBoundary2D
{
	public:
	
		/** Constructor */
		ShapedBoundary2D(int nPoints, int nModes, string shape, double xDim, double yDim);
		
		/** Destructor */
		virtual ~ShapedBoundary2D();
	
		///public static members
		public:
			const static int IS_INSIDE;
			const static int IS_OUTSIDE;
			const static int TO_BE_KILLED;
			
			const static double PI;	
			
	protected:
		
		double r_circle_;
		double a_ellipse_,b_ellipse_;
		double a_rect_,b_rect_;
};

#endif
//endif for SC_SHAPED_BOUNDARY_2D_H
