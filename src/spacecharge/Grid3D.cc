//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    Grid3D.cc
//
// AUTHOR
//    A. Shishlo
//
// CREATED
//    11/05/2010
//
// DESCRIPTION
//    Source code for the class "Grid3D" that is used to operate with 3D arrays
//
///////////////////////////////////////////////////////////////////////////
#include <iostream>

#include "ParticleMacroSize.hh"
#include "BufferStore.hh"

#include "Grid3D.hh"

using namespace OrbitUtils;

/** Constructor */
Grid3D::Grid3D(int nX, int nY, int nZ): CppPyWrapper(NULL)
{
  nX_ = nX;
  nY_ = nY;
  nZ_ = nZ;
	
  dx_ = 0.;
  dy_ = 0.;
  dz_ = 0.;
	
  xMin_=0.0; xMax_=0.0;
  yMin_=0.0; yMax_=0.0;
  zMin_=0.0; zMax_=0.0;


  //Allocate memory for the 3D distribution
  Arr3D = new double**[nZ_];
  for(int iz=0 ; iz < nZ_; iz++){
    Arr3D[iz] = new double*[nX_];
    for(int ix=0 ; ix < nX_; ix++){
      Arr3D[iz][ix] = new double[nY_];
      for(int iy=0 ; iy < nY_; iy++){
	     Arr3D[iz][ix][iy] = 0.0;
      }
    }
  }
}

Grid3D::~Grid3D()
{
  for(int iz=0 ; iz < nZ_; iz++){
    for(int ix=0 ; ix < nX_; ix++){
      delete [] Arr3D[iz][ix];
    }   
    delete [] Arr3D[iz];
  }  
  delete [] Arr3D;
}

/** Returns the reference to the inner 3D array */
double*** Grid3D::getArr3D(){return Arr3D;}

/** Returns the reference to one 2D slice of the inner 3D array */
double** Grid3D::getSlice2D(int zInd){return Arr3D[zInd];}


/** Returns the grid size in x-direction */
int Grid3D::getSizeX(){
	return nX_;
}

/** Returns the grid size in y-direction */
int Grid3D::getSizeY(){
	return nY_;
}

/** Returns the grid size in z-direction */
int Grid3D::getSizeZ(){
	return nZ_;
}

/** Returns the grid point x-coordinate for this index. */   
double Grid3D::getGridX(int index){
	return xMin_ + index*dx_;
}

/** Returns the grid point y-coordinate for this index. */   
double Grid3D::getGridY(int index){
	return yMin_ + index*dy_;
}

/** Returns the grid point z-coordinate for this index. */   
double Grid3D::getGridZ(int index){
	return zMin_ + index*dz_;
}

/** Returns the grid step along x-axis */
double Grid3D::getStepX(){
	return dx_;
}

/** Returns the grid step along y-axis */
double Grid3D::getStepY(){
	return dy_;
}

/** Returns the grid step along z-axis */
double Grid3D::getStepZ(){
	return dz_;
}

/** Returns the max x in the grid points */ 
double Grid3D::getMaxX(){return xMax_;};

/** Returns the min x in the grid points */ 
double Grid3D::getMinX(){return xMin_;};

/** Returns the max y in the grid points */ 
double Grid3D::getMaxY(){return yMax_;};

/** Returns the min y in the grid points */ 
double Grid3D::getMinY(){return yMin_;};

/** Returns the max z in the grid points */ 
double Grid3D::getMaxZ(){return zMax_;};

/** Returns the min z in the grid points */ 
double Grid3D::getMinZ(){return zMin_;};

/** Sets the limits for the x-grid */
void Grid3D::setGridX(double xMin, double xMax){
	xMin_ = xMin;
	xMax_ = xMax;
	dx_ = (xMax_ - xMin_)/(nX_ -1);
	setZero();
}

/** Sets the limits for the y-grid */
void Grid3D::setGridY(double yMin, double yMax){
	yMin_ = yMin;
	yMax_ = yMax;
	dy_ = (yMax_ - yMin_)/(nY_ -1);		
	setZero();
}

/** Sets the limits for the z-grid */
void Grid3D::setGridZ(double zMin, double zMax)
{
  zMin_ = zMin;
  zMax_ = zMax;
  if( nZ_ > 1){
   dz_ = (zMax_-zMin_)/( nZ_ -1);
  }
  else{
   dz_ = zMax_-zMin_;
  }
	setZero();	
}

/** Set all array values to zero */
void Grid3D::setZero()
{
  for(int k = 0; k < nZ_; k++){
    for(int i = 0; i < nX_; i++){
      for(int j = 0; j < nY_; j++){
	     Arr3D[k][i][j] = 0.;
      }
    }
  }
}

/** Multiply all elements of Grid3D by constant coefficient */
void Grid3D::multiply(double coeff)
{
  for(int k = 0; k < nZ_; k++){
    for(int i = 0; i < nX_; i++){
      for(int j = 0; j < nY_; j++){
	      Arr3D[k][i][j] = coeff*Arr3D[k][i][j];
      }
    }
  }
}

void Grid3D::getIndAndFracX(double x, int& ind, double& frac){
   ind  = int ( (x - xMin_)/dx_ + 0.5 );
   if(ind < 1) ind = 1;
   if(ind > (nX_-2)) ind = nX_ - 2;
   frac = (x - (xMin_ + ind*dx_))/dx_; 	
}

void Grid3D::getIndAndFracY(double y, int& ind, double& frac){
   ind  = int ( (y - yMin_)/dy_ + 0.5 );
   if(ind < 1) ind = 1;
   if(ind > (nY_-2)) ind = nY_ - 2;
   frac = (y - (yMin_ + ind*dy_))/dy_; 	
}	

void Grid3D::getGridIndAndFrac(double x, int& xInd, double& xFrac,
	double y, int& yInd, double& yFrac,
	double z, int& zInd, double& zFrac)
{
  this->getIndAndFracX( x, xInd, xFrac);
  this->getIndAndFracY( y, yInd, yFrac);
	
  //z direction
  if( nZ_ > 1){
    zInd  = int ( (z - zMin_)/dz_ + 0.5 );
		
    if( nZ_ > 2){
      //cut off edge for three point interpolation
      if(zInd <= 0) zInd = 1;
      if(zInd >=  nZ_ -1) zInd =  nZ_  - 2;
      zFrac = (z - (zMin_ + zInd * dz_))/dz_;
    }
    if( nZ_ == 2){
      zFrac = (z - zMin_)/dz_;
      if(zInd < 0) {
				zInd = 0;
				zFrac = 0.0;
      }
      else if(zInd > 1) {
				zInd = 1;
				zFrac = 1.0;
      }
    }
  }
  else{
    dz_=0.0; zInd = 0;zFrac = 0.0;
  }
}

/** Sets the value to the one point of the 3D grid  */	
void Grid3D::setValue(double value, int ix, int iy, int iz){
	Arr3D[iz][ix][iy] = value;
}

/** Returns the value on grid*/
double Grid3D::getValueOnGrid(int ix, int iy, int iz){
	return Arr3D[iz][ix][iy];
}

/** Bins the Bunch into the 2D grid. If bunch has a macrosize particle attribute it will be used. */	
void Grid3D::binBunch(Bunch* bunch){
	bunch->compress();
	double** part_coord_arr = bunch->coordArr();
	int has_msize = bunch->hasParticleAttributes("macrosize");
	if(has_msize > 0){
		ParticleMacroSize* macroSizeAttr = (ParticleMacroSize*) bunch->getParticleAttributes("macrosize");
		double m_size = 0.;
		for(int i = 0, n = bunch->getSize(); i < n; i++){
			m_size = macroSizeAttr->macrosize(i);
			this->binValue(m_size,part_coord_arr[i][0],part_coord_arr[i][2],part_coord_arr[i][4]);
		}	
		return;
	}
	double m_size = bunch->getMacroSize();
	int nParts = bunch->getSize();
	for(int i = 0; i < nParts; i++){
		this->binValue(m_size,part_coord_arr[i][0],part_coord_arr[i][2],part_coord_arr[i][4]);	
	}
  this->synchronizeMPI(bunch->getMPI_Comm_Local());
}


/** Bins the value onto grid */
void Grid3D::binValue(double macroSize, double x, double y, double z)
{
	if(x < xMin_ || x > xMax_ || y < yMin_ || y > yMax_ || z < zMin_ || z > zMax_) return;	
  int iX, iY, iZ;
  double xFrac, yFrac, zFrac;
  getGridIndAndFrac(x, iX, xFrac, y, iY, yFrac, z, iZ, zFrac);

  //Calculate interpolation weight
  double Wxm,Wx0,Wxp,Wym,Wy0,Wyp;
  double Wzm,Wz0,Wzp;
  Wzm = Wz0 = Wzp = 0.0;
 
  Wxm = 0.5 * (0.5 - xFrac) * (0.5 - xFrac);
  Wx0 = 0.75 - xFrac * xFrac;
  Wxp = 0.5 * (0.5 + xFrac) * (0.5 + xFrac);
  Wym = 0.5 * (0.5 - yFrac) * (0.5 - yFrac);
  Wy0 = 0.75 - yFrac * yFrac;
  Wyp = 0.5 * (0.5 + yFrac) * (0.5 + yFrac);

  if( nZ_ >= 3){
    Wzm = 0.5 * (0.5 - zFrac) * (0.5 - zFrac);
    Wz0 = 0.75 - zFrac * zFrac;
    Wzp = 0.5 * (0.5 + zFrac) * (0.5 + zFrac);
  }
  if( nZ_ == 2){
    Wzm = 1.0 - zFrac; // for zInd=0
    Wz0 = 0.0; 
    Wzp = zFrac;       // for zInd=1
  }

  //Add weight of particle to Arr3D

  double tmp;

  if( nZ_ >= 3){
    tmp = Wym * Wzm *macroSize;
    Arr3D[iZ-1][iX-1][iY-1] += Wxm * tmp;
    Arr3D[iZ-1][iX  ][iY-1] += Wx0 * tmp;
    Arr3D[iZ-1][iX+1][iY-1] += Wxp * tmp;
    tmp = Wy0 * Wzm *macroSize;
    Arr3D[iZ-1][iX-1][iY  ] += Wxm * tmp;
    Arr3D[iZ-1][iX  ][iY  ] += Wx0 * tmp;
    Arr3D[iZ-1][iX+1][iY  ] += Wxp * tmp;
    tmp = Wyp * Wzm *macroSize;
    Arr3D[iZ-1][iX-1][iY+1] += Wxm * tmp;
    Arr3D[iZ-1][iX  ][iY+1] += Wx0 * tmp;
    Arr3D[iZ-1][iX+1][iY+1] += Wxp * tmp;
    tmp = Wym * Wz0 *macroSize;
    Arr3D[iZ  ][iX-1][iY-1] += Wxm * tmp;
    Arr3D[iZ  ][iX  ][iY-1] += Wx0 * tmp;
    Arr3D[iZ  ][iX+1][iY-1] += Wxp * tmp;
    tmp = Wy0 * Wz0 *macroSize;
    Arr3D[iZ  ][iX-1][iY  ] += Wxm * tmp;
    Arr3D[iZ  ][iX  ][iY  ] += Wx0 * tmp;
    Arr3D[iZ  ][iX+1][iY  ] += Wxp * tmp;
    tmp = Wyp * Wz0 *macroSize;
    Arr3D[iZ  ][iX-1][iY+1] += Wxm * tmp;
    Arr3D[iZ  ][iX  ][iY+1] += Wx0 * tmp;
    Arr3D[iZ  ][iX+1][iY+1] += Wxp * tmp;
    tmp = Wym * Wzp *macroSize;
    Arr3D[iZ+1][iX-1][iY-1] += Wxm * tmp;
    Arr3D[iZ+1][iX  ][iY-1] += Wx0 * tmp;
    Arr3D[iZ+1][iX+1][iY-1] += Wxp * tmp;
    tmp = Wy0 * Wzp *macroSize;
    Arr3D[iZ+1][iX-1][iY  ] += Wxm * tmp;
    Arr3D[iZ+1][iX  ][iY  ] += Wx0 * tmp;
    Arr3D[iZ+1][iX+1][iY  ] += Wxp * tmp;
    tmp = Wyp * Wzp *macroSize;
    Arr3D[iZ+1][iX-1][iY+1] += Wxm * tmp;
    Arr3D[iZ+1][iX  ][iY+1] += Wx0 * tmp;
    Arr3D[iZ+1][iX+1][iY+1] += Wxp * tmp;  
  }
  if( nZ_ == 2){
    tmp = Wym * Wzm *macroSize;
    Arr3D[0][iX-1][iY-1] += Wxm * tmp;
    Arr3D[0][iX  ][iY-1] += Wx0 * tmp;
    Arr3D[0][iX+1][iY-1] += Wxp * tmp;
    tmp = Wy0 * Wzm *macroSize;
    Arr3D[0][iX-1][iY  ] += Wxm * tmp;
    Arr3D[0][iX  ][iY  ] += Wx0 * tmp;
    Arr3D[0][iX+1][iY  ] += Wxp * tmp;
    tmp = Wyp * Wzm *macroSize;
    Arr3D[0][iX-1][iY+1] += Wxm * tmp;
    Arr3D[0][iX  ][iY+1] += Wx0 * tmp;
    Arr3D[0][iX+1][iY+1] += Wxp * tmp;
    tmp = Wym * Wzp *macroSize;
    Arr3D[1][iX-1][iY-1] += Wxm * tmp;
    Arr3D[1][iX  ][iY-1] += Wx0 * tmp;
    Arr3D[1][iX+1][iY-1] += Wxp * tmp;
    tmp = Wy0 * Wzp *macroSize;
    Arr3D[1][iX-1][iY  ] += Wxm * tmp;
    Arr3D[1][iX  ][iY  ] += Wx0 * tmp;
    Arr3D[1][iX+1][iY  ] += Wxp * tmp;
    tmp = Wyp * Wzp *macroSize;
    Arr3D[1][iX-1][iY+1] += Wxm * tmp;
    Arr3D[1][iX  ][iY+1] += Wx0 * tmp;
    Arr3D[1][iX+1][iY+1] += Wxp * tmp;  
  }
  if( nZ_ == 1){
    tmp = Wym * macroSize;
    Arr3D[0][iX-1][iY-1] += Wxm * tmp;
    Arr3D[0][iX  ][iY-1] += Wx0 * tmp;
    Arr3D[0][iX+1][iY-1] += Wxp * tmp;
    tmp = Wy0 * macroSize;
    Arr3D[0][iX-1][iY  ] += Wxm * tmp;
    Arr3D[0][iX  ][iY  ] += Wx0 * tmp;
    Arr3D[0][iX+1][iY  ] += Wxp * tmp;
    tmp = Wyp * macroSize;
    Arr3D[0][iX-1][iY+1] += Wxm * tmp;
    Arr3D[0][iX  ][iY+1] += Wx0 * tmp;
    Arr3D[0][iX+1][iY+1] += Wxp * tmp;
  }
}

/** Calculates gradient of Arr3D. gradX = gradient_x(Arr3D), and so on */
void Grid3D::calcGradient(double x,double& gradX,
			  double y,double& gradY,
			  double z,double& gradZ)
{
  int iX, iY, iZ;
  double xFrac, yFrac, zFrac;
  getGridIndAndFrac(x, iX, xFrac, y, iY, yFrac, z, iZ, zFrac); 

  //coeff. of Arr3D
  double Wxm,Wx0,Wxp,Wym,Wy0,Wyp;
  double Wzm, Wz0, Wzp;
  Wzm = Wz0 = Wzp = 0.0;

  Wxm = 0.5 * (0.5 - xFrac) * (0.5 - xFrac);
  Wx0 = 0.75 - xFrac * xFrac;
  Wxp = 0.5 * (0.5 + xFrac) * (0.5 + xFrac);
  Wym = 0.5 * (0.5 - yFrac) * (0.5 - yFrac);
  Wy0 = 0.75 - yFrac * yFrac;
  Wyp = 0.5 * (0.5 + yFrac) * (0.5 + yFrac);
  Wzm = 0.0; Wz0 = 0.0; Wzp = 0.0;
  if( nZ_ >= 3){
    Wzm = 0.5 * (0.5 - zFrac) * (0.5 - zFrac);
    Wz0 = 0.75 - zFrac * zFrac;
    Wzp = 0.5 * (0.5 + zFrac) * (0.5 + zFrac);
  }
  else if( nZ_ == 2){
    Wzm = 1.0 - zFrac; // for zInd=0
    Wz0 = 0.0; 
    Wzp = zFrac;       // for zInd=1
  }

  //coeff. for grad(Arr3D) through the Gradient weights, dWs
  double dWxm,dWx0,dWxp,dWym,dWy0,dWyp;
  double dWzm, dWz0, dWzp; 
  dWzm = dWz0 = dWzp = 0.0;

  dWxm = (-1.0)*  (0.5 - xFrac)/dx_;
  dWx0 = (-1.0)*    2. * xFrac /dx_;
  dWxp = (-1.0)* -(0.5 + xFrac)/dx_;
  dWym = (-1.0)*  (0.5 - yFrac)/dy_;
  dWy0 = (-1.0)*    2. * yFrac /dy_;
  dWyp = (-1.0)* -(0.5 + yFrac)/dy_;
  if( nZ_ >= 3){
    dWzm = (-1.0)*  (0.5 - zFrac)/dz_;
    dWz0 = (-1.0)*    2. * zFrac /dz_;
    dWzp = (-1.0)* -(0.5 + zFrac)/dz_;
  }
  else if( nZ_ == 2){
    dWzm = (-1.0)*  1.0/dz_; // for zInd=0
    dWz0 =  0.0; 
    dWzp = (-1.0)* -1.0/dz_; // for zInd=1   
  }

  //calculate gradient
  if( nZ_ >= 3){
    gradX = 
      calcSheetGradient(iZ-1,iX,iY,dWxm,dWx0,dWxp,Wym,Wy0,Wyp) *Wzm +
      calcSheetGradient(iZ  ,iX,iY,dWxm,dWx0,dWxp,Wym,Wy0,Wyp) *Wz0 +
      calcSheetGradient(iZ+1,iX,iY,dWxm,dWx0,dWxp,Wym,Wy0,Wyp) *Wzp;
    gradY = 
      calcSheetGradient(iZ-1,iX,iY,Wxm,Wx0,Wxp,dWym,dWy0,dWyp) *Wzm +
      calcSheetGradient(iZ  ,iX,iY,Wxm,Wx0,Wxp,dWym,dWy0,dWyp) *Wz0 +
      calcSheetGradient(iZ+1,iX,iY,Wxm,Wx0,Wxp,dWym,dWy0,dWyp) *Wzp;
    gradZ = 
      calcSheetGradient(iZ-1,iX,iY,Wxm,Wx0,Wxp,Wym,Wy0,Wyp) *dWzm +
      calcSheetGradient(iZ  ,iX,iY,Wxm,Wx0,Wxp,Wym,Wy0,Wyp) *dWz0 +
      calcSheetGradient(iZ+1,iX,iY,Wxm,Wx0,Wxp,Wym,Wy0,Wyp) *dWzp;
  } 
  if( nZ_ == 2){
    gradX = 
      calcSheetGradient(0,iX,iY,dWxm,dWx0,dWxp,Wym,Wy0,Wyp) *Wzm +
      calcSheetGradient(1,iX,iY,dWxm,dWx0,dWxp,Wym,Wy0,Wyp) *Wzp;
    gradY = 
      calcSheetGradient(0,iX,iY,Wxm,Wx0,Wxp,dWym,dWy0,dWyp) *Wzm +
      calcSheetGradient(1,iX,iY,Wxm,Wx0,Wxp,dWym,dWy0,dWyp) *Wzp;
    gradZ = 
      calcSheetGradient(0,iX,iY,Wxm,Wx0,Wxp,Wym,Wy0,Wyp) *dWzm +
      calcSheetGradient(1,iX,iY,Wxm,Wx0,Wxp,Wym,Wy0,Wyp) *dWzp;
  }
  if( nZ_ == 1){
    gradX = 
      calcSheetGradient(0,iX,iY,dWxm,dWx0,dWxp,Wym,Wy0,Wyp);
    gradY = 
      calcSheetGradient(0,iX,iY,Wxm,Wx0,Wxp,dWym,dWy0,dWyp);
    gradZ = 0.0;
  }

}

/** Calculates value at the point with coordinates x,y,z */
double Grid3D::getValue(double x,double y,double z)
{
  int iX,iY,iZ;
  double fracX,fracY,fracZ;
  getGridIndAndFrac( x, iX, fracX, y, iY, fracY, z, iZ, fracZ);

  double value = 0.0; 

  double Wxm, Wx0, Wxp, Wym, Wy0, Wyp;
  double Wzm = 0.0, Wz0 = 0.0, Wzp = 0.0;

  Wxm = 0.5 - fracX;
  Wxm *=0.5*Wxm; 
  Wxp = 0.5 + fracX;
  Wxp *=0.5*Wxp;
  Wx0 = 0.75 - fracX*fracX;

  Wym = 0.5 - fracY;
  Wym *=0.5*Wym; 
  Wyp = 0.5 + fracY;
  Wyp *=0.5*Wyp;
  Wy0 = 0.75 - fracY*fracY;

  if( nZ_ > 2 ){
    Wzm = 0.5 - fracZ;
    Wzm *=0.5*Wzm; 
    Wzp = 0.5 + fracZ;
    Wzp *=0.5*Wzp;
    Wz0 = 0.75 - fracZ*fracZ;
  } 
  else if( nZ_ == 2 ){
    Wzm = 1.0 - fracZ;
    Wz0 = 0.0;
    Wzp = fracZ;
  }
  else if( nZ_ == 1 ){
    Wzm = 0.0;
    Wz0 = 1.0;
    Wzp = 0.0;    
  }

  double Vm, V0, Vp;
  double Vzm, Vz0, Vzp;

  if( nZ_ > 2 ){
    Vm = calcValueOnX(iX,iY-1,iZ-1,Wxm,Wx0,Wxp);
    V0 = calcValueOnX(iX,iY  ,iZ-1,Wxm,Wx0,Wxp);
    Vp = calcValueOnX(iX,iY+1,iZ-1,Wxm,Wx0,Wxp);
    Vzm = Wym*Vm + Wy0*V0 + Wyp*Vp;

    Vm = calcValueOnX(iX,iY-1,iZ,Wxm,Wx0,Wxp);
    V0 = calcValueOnX(iX,iY  ,iZ,Wxm,Wx0,Wxp);
    Vp = calcValueOnX(iX,iY+1,iZ,Wxm,Wx0,Wxp);
    Vz0 = Wym*Vm + Wy0*V0 + Wyp*Vp;

    Vm = calcValueOnX(iX,iY-1,iZ+1,Wxm,Wx0,Wxp);
    V0 = calcValueOnX(iX,iY  ,iZ+1,Wxm,Wx0,Wxp);
    Vp = calcValueOnX(iX,iY+1,iZ+1,Wxm,Wx0,Wxp);
    Vzp = Wym*Vm + Wy0*V0 + Wyp*Vp;

    value = Wzm*Vzm + Wz0*Vz0 + Wzp*Vzp;   
  }
  else if(nZ_ == 2){
    Vm = calcValueOnX(iX,iY-1,0,Wxm,Wx0,Wxp);
    V0 = calcValueOnX(iX,iY  ,0,Wxm,Wx0,Wxp);
    Vp = calcValueOnX(iX,iY+1,0,Wxm,Wx0,Wxp);
    Vzm = Wym*Vm + Wy0*V0 + Wyp*Vp; 

    Vm = calcValueOnX(iX,iY-1,1,Wxm,Wx0,Wxp);
    V0 = calcValueOnX(iX,iY  ,1,Wxm,Wx0,Wxp);
    Vp = calcValueOnX(iX,iY+1,1,Wxm,Wx0,Wxp);
    Vzp = Wym*Vm + Wy0*V0 + Wyp*Vp; 

    value = Wzm*Vzm + Wzp*Vzp; 
  }
  else if(nZ_ == 1){
    Vm = calcValueOnX(iX,iY-1,0,Wxm,Wx0,Wxp);
    V0 = calcValueOnX(iX,iY  ,0,Wxm,Wx0,Wxp);
    Vp = calcValueOnX(iX,iY+1,0,Wxm,Wx0,Wxp);

    value = Wym*Vm + Wy0*V0 + Wyp*Vp; 
  }

  return value;
}

/** returns the sum of all grid points */
double Grid3D::getSum()
{
  double sum = 0.;
  for(int iZ = 0; iZ < nZ_; iZ++){
    for(int i = 0; i < nX_; i++){
      for(int j = 0; j < nY_; j++){
	      sum += Arr3D[iZ][i][j];
      }
    }
  }
  return sum;
}

/** returns the sum of all grid points in the slice with index iZ */
double Grid3D::getSliceSum(int iZ)
{
  double sum = 0.;
  if( iZ <   0 ) return sum;
  if( iZ > (nZ_-1)) return sum;
  for(int i = 0; i < nX_; i++){
    for(int j = 0; j < nY_; j++){
      sum += Arr3D[iZ][i][j];
    }
  }
  return sum;
}

/** returns the sum of all grid points in the slice with  position z */
double Grid3D::getSliceSum(double z)
{
  double sum = 0.;
  if(z > getMaxZ()) return sum;
  if(z < getMinZ()) return sum;
  int zInd = 0;
  double zFrac = 0.;
  if(nZ_ > 2){
    double sumM,sum0,sumP;
    zInd  = int ( (z - zMin_)/dz_ + 0.5 );
    //cut off edge for three point interpolation
    if(zInd <= 0) zInd = 1;
    if(zInd >=  nZ_ -1) zInd =  nZ_  - 2;
    zFrac = (z - (zMin_ + zInd * dz_))/dz_;
    sumM = getSliceSum(zInd - 1);
    sum0 = getSliceSum(zInd);
    sumP = getSliceSum(zInd + 1);
    sum = 0.5*(0.5 - zFrac)*(0.5 - zFrac)*sumM + 
          (0.75 - zFrac*zFrac)*sum0 +
          0.5*(0.5 + zFrac)*(0.5 + zFrac)*sumP;
    return sum;
  }
  if(nZ_ == 2){
   zFrac = (z - zMin_)/dz_;
   sum = (1.0 - zFrac)*getSliceSum(0) + zFrac*getSliceSum(1);
   return sum;
  }
  if(nZ_ == 1){
    sum = getSliceSum(0);
    return sum;
  }
  return sum;
}

/** Calculates interpolated value along x-axis for fixed y and z */
double Grid3D::calcValueOnX(int iX, int iY, int iZ, double Wxm,double Wx0,double Wxp)
{
  return ( Wxm*Arr3D[iZ][iX-1][iY]+Wx0*Arr3D[iZ][iX][iY]+Wxp*Arr3D[iZ][iX+1][iY]);
}

/** Calculates interpolated value along y-axis for fixed x and z */
double Grid3D::calcValueOnY(int iX, int iY, int iZ, double Wym,double Wy0,double Wyp)
{
  return ( Wym*Arr3D[iZ][iX][iY-1]+Wy0*Arr3D[iZ][iX][iY]+Wyp*Arr3D[iZ][iX][iY+1]);
}

/** Calculates Gradient from each z-sheet without interpolating in z */
double Grid3D::calcSheetGradient(int iZ,int iX,int iY,
				 double xm,double x0,double xp,
				 double ym,double y0,double yp)		     
{
  double sheetGradient = 
    xm * ym * Arr3D[iZ][iX-1][iY-1] +
    x0 * ym * Arr3D[iZ][iX  ][iY-1] +
    xp * ym * Arr3D[iZ][iX+1][iY-1] +
    xm * y0 * Arr3D[iZ][iX-1][iY  ] +
    x0 * y0 * Arr3D[iZ][iX  ][iY  ] +
    xp * y0 * Arr3D[iZ][iX+1][iY  ] +
    xm * yp * Arr3D[iZ][iX-1][iY+1] +
    x0 * yp * Arr3D[iZ][iX  ][iY+1] +
    xp * yp * Arr3D[iZ][iX+1][iY+1];
  return sheetGradient;
}

/**synchronize MPI */
void Grid3D::synchronizeMPI(pyORBIT_MPI_Comm* pyComm){
  // ====== MPI  start ========
	int size_MPI = nX_ * nY_*nZ_;
	int buff_index0 = 0;
	int buff_index1 = 0;
	double* inArr  = BufferStore::getBufferStore()->getFreeDoubleArr(buff_index0,size_MPI);
	double* outArr = BufferStore::getBufferStore()->getFreeDoubleArr(buff_index1,size_MPI);
	
	int count = 0;
	for(int iz = 0; iz < nZ_; iz++){
		for(int ix = 0; ix < nX_; ix++){
			for(int iy = 0; iy < nY_; iy++){
				inArr[count] = Arr3D[iz][ix][iy];
				count += 1;
			}
		}
	}
	
	if(pyComm == NULL) {
		ORBIT_MPI_Allreduce(inArr,outArr,size_MPI,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	} else {
		ORBIT_MPI_Allreduce(inArr,outArr,size_MPI,MPI_DOUBLE,MPI_SUM,pyComm->comm);
	}
	
	count = 0;
	for(int iz = 0; iz < nZ_; iz++){	
		for(int ix = 0; ix < nX_; ix++){
			for(int iy = 0; iy < nY_; iy++){
				Arr3D[iz][ix][iy] = outArr[count];
				count += 1;
			}
		}
	}	

	OrbitUtils::BufferStore::getBufferStore()->setUnusedDoubleArr(buff_index0);
	OrbitUtils::BufferStore::getBufferStore()->setUnusedDoubleArr(buff_index1);	
	
  // ===== MPI end =====
}
