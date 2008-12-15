//This class repersents a 2D rectangular grid 

#include "Grid2D.hh"

#include <iostream>

using namespace OrbitUtils;

// Constructor
Grid2D::Grid2D(int xBins,    int yBins): CppPyWrapper(NULL)
{
	boundary2D_ = NULL;
	xBins_ = xBins;
	yBins_ = yBins;
	dx_ = 1.0;
	dy_ = 1.0;
	xGrid_ = new double[xBins_];
	yGrid_ = new double[yBins_];	
	makeGrid(xGrid_,dx_,0.,1.0*(xBins_ - 1), xBins_);
	makeGrid(yGrid_,dy_,0.,1.0*(yBins_ - 1), yBins_);	
	arr_ = new double*[xBins_];
	for(int i = 0; i < xBins_; i++){
		arr_[i] = new double[yBins_];
	}
	setZero();
}

Grid2D::Grid2D(double xMin, double xMax, int xBins,    
	             double yMin, double yMax, int yBins): CppPyWrapper(NULL)
{
	boundary2D_ = NULL;
	xBins_ = xBins;
	yBins_ = yBins;
	xGrid_ = new double[xBins_];
	yGrid_ = new double[yBins_];	
	makeGrid(xGrid_,dx_,xMin,xMax, xBins_);
	makeGrid(yGrid_,dy_,yMin,yMax, yBins_);	
	arr_ = new double*[xBins_];
	for(int i = 0; i < xBins_; i++){
		arr_[i] = new double[yBins_];
	}
	setZero();
}

Grid2D::Grid2D(Boundary2D* boundary2D): CppPyWrapper(NULL)
{
	boundary2D_ = boundary2D;
	xBins_ = boundary2D_->getBinsX();
	yBins_ = boundary2D_->getBinsX();
	xGrid_ = new double[xBins_];
	yGrid_ = new double[yBins_];
	makeGrid(xGrid_, dx_, boundary2D_->getGridMinX(),boundary2D_-> getGridMaxX(), xBins_);
	makeGrid(yGrid_, dy_, boundary2D_->getGridMinY(),boundary2D_-> getGridMaxY(), yBins_);
	arr_ = new double*[xBins_];
	for(int i = 0; i < xBins_; i++){
		arr_[i] = new double[yBins_];
	}
	setZero();
}

// Destructor
Grid2D::~Grid2D()
{
	delete[] xGrid_;
	delete[] yGrid_;
	for(int i = 0; i < xBins_; i++){
		delete [] arr_[i];
	}
	delete [] arr_;
}

void Grid2D::makeGrid(double* grid, double& step, double zMin,double zMax,int n){
	grid[0] = zMin;
	step = (zMax - zMin)/(n - 1);
	for(int i = 1; i < n; i++){
		grid[i] =  grid[0] + ((zMax - zMin)*i)/(n - 1);
	}
	grid[n-1] = zMax;
}

Boundary2D* Grid2D::getBoundary2D(){
	return boundary2D_;
}

void Grid2D::setBoundary2D(Boundary2D* boundary2D){
	boundary2D_ = boundary2D;
	int xBins = boundary2D_->getBinsX();
	int yBins = boundary2D_->getBinsX();
	double xMin = boundary2D_->getGridMinX();
	double xMax = boundary2D_->getGridMaxX();
	double yMin = boundary2D_->getGridMinY();
	double yMax = boundary2D_->getGridMaxY();
	if(	xBins != xBins_ || yBins != yBins_ ||
		xMin  != xGrid_[0] || xMax != xGrid_[xBins_-1] ||
		yMin  != yGrid_[0] || yMax != yGrid_[yBins_-1]){
	  if(xBins != xBins_ || yBins != yBins_){
			delete[] xGrid_;
			delete[] yGrid_;
			for(int i = 0; i < xBins_; i++){
				delete [] arr_[i];
			}
			delete [] arr_;
			xBins_ = boundary2D_->getBinsX();
			yBins_ = boundary2D_->getBinsX();		
			arr_ = new double*[xBins_];
			for(int i = 0; i < xBins_; i++){
				arr_[i] = new double[yBins_];
			}			
			xGrid_ = new double[xBins_];
			yGrid_ = new double[yBins_];
		}
		makeGrid(xGrid_, dx_, boundary2D_->getGridMinX(),boundary2D_-> getGridMaxX(), xBins_);
		makeGrid(yGrid_, dy_, boundary2D_->getGridMinY(),boundary2D_-> getGridMaxY(), yBins_);
	  setZero();		
	}	
}
	
/** Sets the value to the one point of the 2D grid  */	
void Grid2D::setValue(double value, int ix, int iy){
	arr_[ix][iy] = value;
}

/** Returns the interpolated value from the 2D grid */	
double Grid2D::getValue(double x, double y){
	int iX, iY;
	double Wxm, Wx0, Wxp, Wym, Wy0, Wyp;
	double xFract,  yFract;
	double xFract2, yFract2;
	getIndAndFracX(x,iX,xFract);
	getIndAndFracY(y,iY,yFract);
	xFract2 = xFract * xFract;
	yFract2 = yFract * yFract;
	Wxm = 0.5 * (0.25 - xFract + xFract2);
	Wx0 = 0.75 - xFract2;
	Wxp = 0.5 * (0.25 + xFract + xFract2);
	Wym = 0.5 * (0.25 - yFract + yFract2);
	Wy0 = 0.75 - yFract2;
	Wyp = 0.5 * (0.25 + yFract + yFract2);
	double value =
	Wxm * Wym * arr_[iX-1][iY-1] +
	Wxm * Wy0 * arr_[iX-1][iY]   +
	Wxm * Wyp * arr_[iX-1][iY+1] +
	Wx0 * Wym * arr_[iX]  [iY-1] +
	Wx0 * Wy0 * arr_[iX]  [iY]   +
	Wx0 * Wyp * arr_[iX]  [iY+1] +
	Wxp * Wym * arr_[iX+1][iY-1] +
	Wxp * Wy0 * arr_[iX+1][iY]   +
	Wxp * Wyp * arr_[iX+1][iY+1]; 
	return value;
}	

/** Bins the value into the 2D grid */	
void Grid2D::binValue(double value, double x, double y){
	int iX, iY;
	double Wxm, Wx0, Wxp, Wym, Wy0, Wyp;
	double xFract,  yFract;
	double xFract2, yFract2;
	getIndAndFracX(x,iX,xFract);
	getIndAndFracY(y,iY,yFract);
	xFract2 = xFract * xFract;
	yFract2 = yFract * yFract;
	Wxm = 0.5 * (0.25 - xFract + xFract2);
	Wx0 = 0.75 - xFract2;
	Wxp = 0.5 * (0.25 + xFract + xFract2);
	Wym = 0.5 * (0.25 - yFract + yFract2);
	Wy0 = 0.75 - yFract2;
	Wyp = 0.5 * (0.25 + yFract + yFract2);
	arr_[iX-1][iY-1] += Wxm * Wym * value;      
  arr_[iX-1][iY]   += Wxm * Wy0 * value;   
  arr_[iX-1][iY+1] += Wxm * Wyp * value; 
  arr_[iX]  [iY-1] += Wx0 * Wym * value; 
  arr_[iX]  [iY]   += Wx0 * Wy0 * value; 
  arr_[iX]  [iY+1] += Wx0 * Wyp * value; 
  arr_[iX+1][iY-1] += Wxp * Wym * value; 
  arr_[iX+1][iY]   += Wxp * Wy0 * value; 
  arr_[iX+1][iY+1] += Wxp * Wyp * value;             
}

/** Calculates gradient at a position (x,y) */	
void Grid2D::calcGradient(double x, double y, double& ex, double& ey){
	double dWxm,dWx0,dWxp,dWym,dWy0,dWyp;
	double Wxm, Wx0, Wxp, Wym, Wy0, Wyp;
	int iX, iY;
	double xFract,  yFract;
	double xFract2, yFract2;
	getIndAndFracX(x,iX,xFract);
	getIndAndFracY(y,iY,yFract);	
	xFract2 = xFract * xFract;
	yFract2 = yFract * yFract;
	Wxm = 0.5 * (0.25 - xFract + xFract2);
	Wx0 = 0.75 - xFract2;
	Wxp = 0.5 * (0.25 + xFract + xFract2);
	Wym = 0.5 * (0.25 - yFract + yFract2);
	Wy0 = 0.75 - yFract2;
	Wyp = 0.5 * (0.25 + yFract + yFract2);	
  dWxm = (-1.0)*(0.5 - xFract); 
	dWx0 = (-1.0)*(-2.) * xFract; 
	dWxp = (-1.0)*(-(0.5 + xFract)); 
	dWym = (-1.0)*(0.5 - yFract); 
	dWy0 = (-1.0)*(+2.) * yFract; 
	dWyp = (-1.0)*(-(0.5 + yFract));
  ex =
	dWxm * Wym * arr_[iX-1][iY-1] +
	dWxm * Wy0 * arr_[iX-1][iY]   +
	dWxm * Wyp * arr_[iX-1][iY+1] +
	dWx0 * Wym * arr_[iX]  [iY-1] +
	dWx0 * Wy0 * arr_[iX]  [iY]   +
	dWx0 * Wyp * arr_[iX]  [iY+1] +
	dWxp * Wym * arr_[iX+1][iY-1] +
	dWxp * Wy0 * arr_[iX+1][iY]   +
	dWxp * Wyp * arr_[iX+1][iY+1]; 
	ex = ex / dx_;
  ey =
	Wxm * dWym * arr_[iX-1][iY-1] +
	Wxm * dWy0 * arr_[iX-1][iY]   +
	Wxm * dWyp * arr_[iX-1][iY+1] +
	Wx0 * dWym * arr_[iX]  [iY-1] +
	Wx0 * dWy0 * arr_[iX]  [iY]   +
	Wx0 * dWyp * arr_[iX]  [iY+1] +
	Wxp * dWym * arr_[iX+1][iY-1] +
	Wxp * dWy0 * arr_[iX+1][iY]   +
	Wxp * dWyp * arr_[iX+1][iY+1];
	ey = ey / dy_;
}

/** Calculates gradient at a grid point (ix,iy) */	
void Grid2D::calcGradient(int iX, int iY, double& ex, double& ey){
	if(iX != 0 && iX != (xBins_ - 1) && iY != 0 && iY != (yBins_ - 1)){
		ex =
		(- 0.125) * arr_[iX-1][iY-1] +
		(- 0.75 )* arr_[iX-1][iY]   +
		(- 0.125) * arr_[iX-1][iY+1] +
		0.125 * arr_[iX+1][iY-1] +
		0.75  * arr_[iX+1][iY]   +
		0.125 * arr_[iX+1][iY+1]; 	
		ex = 0.5 * ex / dx_;
		ey =
		(-0.125) * arr_[iX-1][iY-1] +
		( 0.125) * arr_[iX-1][iY+1] +
		(-0.75 ) * arr_[iX]  [iY-1] +
		( 0.75 ) * arr_[iX]  [iY+1] +
		(-0.125) * arr_[iX+1][iY-1] +
		( 0.125) * arr_[iX+1][iY+1];
		ey = 0.5 * ey / dy_;
		return;
	} 
	double x = xGrid_[iX];
	double y = yGrid_[iY];
	calcGradient(x,y,ex,ey);
}

/** Finds and sets the potential to the external grid */	
void Grid2D::findPotential(Grid2D* phiGrid2D){
	if(getBoundary2D() == NULL || getBoundary2D() != phiGrid2D->getBoundary2D()){
		int rank = 0;
		ORBIT_MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		if(rank == 0){
			std::cerr << " Grid2D::findPotential(Grid2D* phiGrid2D) "<< std::endl
			<< "local boundar !=  phiGrid2D boundary "<< std::endl
			<< "this->getBoundary2D() = " << getBoundary2D() << std::endl
			<< "phiGrid2D->getBoundary2D() = " << phiGrid2D->getBoundary2D() << std::endl
			<< "Stop."<< std::endl;
		}
		ORBIT_MPI_Finalize();
	}
	getBoundary2D()->findPotential(getArr(),phiGrid2D->getArr());
}

/** Adds the boundary induced potential to itself. It should be the potential Grid2D. */	
void Grid2D::addBoundaryPotential(){
	if(getBoundary2D() == NULL){
		int rank = 0;
		ORBIT_MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		if(rank == 0){
			std::cerr << " Grid2D::addBoundaryPotential() "<< std::endl
			<< "No local boundary!"<< std::endl
			<< "this->getBoundary2D() = " << getBoundary2D() << std::endl
			<< "Stop."<< std::endl;
		}
		ORBIT_MPI_Finalize();
	}
	getBoundary2D()->addBoundaryPotential(getArr());
}

/** Sets all grid points to zero */	
void Grid2D::setZero(){
	for(int i = 0; i < xBins_; i++){
		for(int j = 0; j < yBins_; j++){
			arr_[i][j] = 0.;
		}
	}
}

/** Returns the reference to the 2D array */	
double** Grid2D::getArr(){
	return arr_;
}
	
/** Returns the grid size in x-direction */
int Grid2D::getBinsX(){
	return xBins_;
}

/** Returns the grid size in y-direction */
int Grid2D::getBinsY(){
	return yBins_;
}

void Grid2D::getIndAndFracX(double x, int& ind, double& frac){
   ind  = int ( (x - xGrid_[0])/dx_ + 0.5 );
   if(ind < 1) ind = 1;
   if(ind > (xBins_-2)) ind = xBins_ - 2;
   frac = (x - xGrid_[ind])/dx_; 	
}

void Grid2D::getIndAndFracY(double y, int& ind, double& frac){
   ind  = int ( (y - yGrid_[0])/dy_ + 0.5 );
   if(ind < 1) ind = 1;
   if(ind > (yBins_-2)) ind = yBins_ - 2;
   frac = (y - yGrid_[ind])/dy_; 	
}	

/** Returns the grid point x-coordinate for this index. */   
double Grid2D::getGridX(int index){
	return xGrid_[index];
}

/** Returns the grid point y-coordinate for this index. */   
double Grid2D::getGridY(int index){
	return yGrid_[index];
}

/** Returns the grid step along x-axis */
double Grid2D::getStepX(){
	return dx_;
}

/** Returns the grid step along y-axis */
double Grid2D::getStepY(){
	return dy_;
}

/** Sets x-grid */
void Grid2D::setGridX(double xMin, double xMax, int nPoints){
	if(nPoints != xBins_ || xMin != xGrid_[0] || xMax != xGrid_[xBins_-1]){
		if(nPoints != xBins_){
			delete[] xGrid_;
			for(int i = 0; i < xBins_; i++){
				delete [] arr_[i];
			}
			delete [] arr_;
			xBins_ = nPoints;
			xGrid_ = new double[xBins_];
			arr_ = new double*[xBins_];
			for(int i = 0; i < xBins_; i++){
				arr_[i] = new double[yBins_];
			}
		}
		makeGrid(xGrid_, dx_, xMin, xMax, xBins_);		
		setZero();	
	}
}

/** Sets y-grid */
void Grid2D::setGridY(double yMin, double yMax, int nPoints){
	if(nPoints != yBins_ || yMin != yGrid_[0] || yMax != yGrid_[yBins_-1]){
		if(nPoints != yBins_){
			delete[] yGrid_;
			for(int i = 0; i < xBins_; i++){
				delete [] arr_[i];
			}
			delete [] arr_;
			yBins_ = nPoints;
			yGrid_ = new double[yBins_];
			arr_ = new double*[xBins_];
			for(int i = 0; i < xBins_; i++){
				arr_[i] = new double[yBins_];
			}
		}
		makeGrid(yGrid_, dy_, yMin, yMax, yBins_);
		setZero();	
	}		
}
	
