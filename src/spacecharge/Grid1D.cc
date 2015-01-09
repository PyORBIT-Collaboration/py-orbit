/////////////////////////////////////////////////////////////////////////////
//
// FILE NAME
//   Grid1D.cc
//
//   06/25/10
//
// DESCRIPTION
//   Provides a 1D grid and binning routines on that grid
//
/////////////////////////////////////////////////////////////////////////////

#include "Grid1D.hh"
#include "Bunch.hh"
#include "ParticleMacroSize.hh"
#include "BufferStore.hh"

#include <iostream>

using namespace OrbitUtils;


/** Constructor with grid size only */
Grid1D::Grid1D(int zSize):CppPyWrapper(NULL)
{
  zSize_ = zSize;
  zMin_  = -1.0;
  zMax_  = +1.0;
  length_ = 0.0;
  init();
  setZero();
}

/** Constructor with grid size and lattice length */
Grid1D::Grid1D(int zSize, double length):CppPyWrapper(NULL)
{
	zSize_ = zSize;
	length_ = length;
	zMin_  = -1.0;
	zMax_  = +1.0;
	init();
	setZero();
}

/** Constructor with grid size and spatial limits */
Grid1D::Grid1D(int zSize, double zMin, double zMax):CppPyWrapper(NULL)
{
  zSize_ = zSize;
  zMin_  = zMin;
  zMax_  = zMax;
  init();
  setZero();
}


/** Destructor */
Grid1D::~Grid1D()
{
  delete [] arr_;
}


/** Memory allocation and step calculation for dx_ and dy_ */
void Grid1D::init()
{
  if(zSize_ < 1)
  {
    int rank = 0;
    ORBIT_MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(rank == 0)
    {
      std::cerr << "Grid1D::Grid1D - CONSTRUCTOR" << std::endl
                << "The grid size is too small (should be more than 0)!"
                << std::endl << "number of bins = " << zSize_ << std::endl
                << "Stop." << std::endl;
    }
    ORBIT_MPI_Finalize();
  }
  dz_  = (zMax_ - zMin_) / zSize_;
  arr_ = new double[zSize_];
}


/** Sets z-grid */
void Grid1D::setGridZ(double zMin, double zMax)
{
  zMin_ = zMin;
  zMax_ = zMax;
  dz_ = (zMax_ - zMin_) / zSize_;
}


/** Returns the min z grid point value */
double Grid1D::getMinZ()
{
  return zMin_;
}


/** Returns the max z grid point value */
double Grid1D::getMaxZ()
{
  return zMax_;
}


/** Returns the number of grid points */
int Grid1D::getSizeZ()
{
  return zSize_;
}


/** Returns the grid step size */
double Grid1D::getStepZ()
{
  return dz_;
}


/** Returns the grid point for index */
double Grid1D::getGridZ(int index)
{
  return zMin_ + index * dz_;
}


/** Sets arr_ at all grid points to zero */
void Grid1D::setZero()
{
  for(int i = 0; i < zSize_; i++)
  {
    arr_[i] = 0.;
  }
}


/** Sets value of arr_ at one point on the grid */
void Grid1D::setValue(double value, int iZ)
{
  arr_[iZ] = value;
}


/** Returns value of arr_ at one point on the grid */
double Grid1D::getValueOnGrid(int iZ)
{
  return arr_[iZ];
}


/** Returns interpolated value of arr_ */
double Grid1D::getValue(double z)
{
  int iZ0, iZp;
  double WZ0, WZp;
  double arrval;

  getIndAndWZ(z, iZ0, iZp, WZ0, WZp);
  if(zSize_ > 1)
  {
    arrval = WZ0 * arr_[iZ0] + WZp * arr_[iZp];
  }
  else if(zSize_ == 1)
  {
    arrval = arr_[iZ0];
  }
  return arrval;
}


/** Returns smoothed interpolated value of arr_ */
double Grid1D::getValueSmoothed(double z)
{
  double  WZm,  WZ0,  WZp;
  double dWZm, dWZ0, dWZp;
  double arrval;
  int iZm, iZ0, iZp;
  getIndAndWZSmoothed(z, iZm,  iZ0,  iZp,
                         WZm,  WZ0,  WZp,
                        dWZm, dWZ0, dWZp);
  if(zSize_ > 2)
  {
    arrval = WZm * arr_[iZm] + WZ0 * arr_[iZ0] + WZp * arr_[iZp];
  }
  else if(zSize_ == 2)
  {
    arrval = WZ0 * arr_[iZ0] + WZp * arr_[iZp];
  }
  else if(zSize_ == 1)
  {
    arrval = arr_[iZ0];
  }
  return arrval;
}


/** Bins the Bunch to the grid incorporating macrosize */
void Grid1D::binBunch(Bunch* bunch)
{
  double m_size;
  bunch->compress();
  double** part_coord_arr = bunch->coordArr();
  int has_msize = bunch->hasParticleAttributes("macrosize");
  if(has_msize > 0)
  {
    ParticleMacroSize* macroSizeAttr =
            (ParticleMacroSize*) bunch->getParticleAttributes("macrosize");
    m_size = 0.;
    for(int i = 0; i < bunch->getSize(); i++)
    {
      m_size = macroSizeAttr->macrosize(i);
      binValue(m_size, part_coord_arr[i][4]);
    }
  }
  else
  {
    m_size = bunch->getMacroSize();
    for(int i = 0; i < bunch->getSize(); i++)
    {
      binValue(m_size, part_coord_arr[i][4]);
    }
  }
}


/** Bins the Bunch to the grid using a smoothing algorithm
    and incorporating macrosize */
void Grid1D::binBunchSmoothed(Bunch* bunch)
{
  double m_size;
  bunch->compress();
  double** part_coord_arr = bunch->coordArr();
  int has_msize = bunch->hasParticleAttributes("macrosize");
  if(has_msize > 0)
  {
    ParticleMacroSize* macroSizeAttr =
            (ParticleMacroSize*) bunch->getParticleAttributes("macrosize");
    m_size = 0.;
    for(int i = 0; i < bunch->getSize(); i++)
    {
      m_size = macroSizeAttr->macrosize(i);
      binValueSmoothed(m_size, part_coord_arr[i][4]);
    }
  }
  else
  {
    m_size = bunch->getMacroSize();
    for(int i = 0; i < bunch->getSize(); i++)
    {
      binValueSmoothed(m_size, part_coord_arr[i][4]);
    }
  }
}


/** Bins the Bunch to the grid giving each macroparticle 
    unit weight */
void Grid1D::binBunchByParticle(Bunch* bunch)
{
  bunch->compress();
  double** part_coord_arr = bunch->coordArr();
  for(int i = 0; i < bunch->getSize(); i++)
  {
    binValue(1.0, part_coord_arr[i][4]);
  }
}


/** Bins the Bunch to the grid using a smoothing algorithm
    and giving each macroparticle unit weight */
void Grid1D::binBunchSmoothedByParticle(Bunch* bunch)
{
  bunch->compress();
  double** part_coord_arr = bunch->coordArr();
  for(int i = 0; i < bunch->getSize(); i++)
  {
    binValueSmoothed(1.0, part_coord_arr[i][4]);
  }
}


/** Bins property of the Bunch to the grid giving
    each macroparticle unit weight */
void Grid1D::binBunchMoment(int propindex, Bunch* bunch, double* Moment)
{
  double property;
  for(int i = 0; i < zSize_; i++)
  {
    Moment[i] = 0.0;
  }
  bunch->compress();
  double** part_coord_arr = bunch->coordArr();
  for(int i = 0; i < bunch->getSize(); i++)
  {
    property = part_coord_arr[i][propindex];
    binMoment(property, part_coord_arr[i][4], Moment);
  }
}


/** Bins moment of the Bunch to the grid using a smoothing
    algorithm and giving each macroparticle unit weight */
void Grid1D::binBunchSmoothedMoment(int propindex, Bunch* bunch, double* Moment)
{
  double property;
  for(int i = 0; i < zSize_; i++)
  {
    Moment[i] = 0.0;
  }
  bunch->compress();
  double** part_coord_arr = bunch->coordArr();
  for(int i = 0; i < bunch->getSize(); i++)
  {
    property = part_coord_arr[i][propindex];
    binMomentSmoothed(property, part_coord_arr[i][4], Moment);
  }
}


/** Bins a value to the grid */
void Grid1D::binValue(double value, double z)
{
  if(z < zMin_ || z > zMax_ ) return;
  double WZ0, WZp;
  int iZ0, iZp;
  getIndAndWZ(z, iZ0, iZp, WZ0, WZp);
  if(zSize_ > 1)
  {
    arr_[iZ0] += WZ0 * value;
    arr_[iZp] += WZp * value;
  }
  else if(zSize_ == 1)
  {
    arr_[iZ0] += value;
  }
}


/** Bins a value to the grid with smoothing */
void Grid1D::binValueSmoothed(double value, double z)
{
  if(z < zMin_ || z > zMax_ ) return;
  double  WZm,  WZ0,  WZp;
  double dWZm, dWZ0, dWZp;
  int iZm, iZ0, iZp;
  getIndAndWZSmoothed(z, iZm,  iZ0,  iZp,
                         WZm,  WZ0,  WZp,
                        dWZm, dWZ0, dWZp);
  if(zSize_ > 2)
  {
    arr_[iZm] += WZm * value;
    arr_[iZ0] += WZ0 * value;
    arr_[iZp] += WZp * value;
  }
  else if(zSize_ == 2)
  {
    arr_[iZ0] += WZ0 * value;
    arr_[iZp] += WZp * value;
  }
  else if(zSize_ == 1)
  {
    arr_[iZ0] += value;
  }
}


/** Bins a moment to the grid */
void Grid1D::binMoment(double value, double z, double* Moment)
{
  if(z < zMin_ || z > zMax_ ) return;
  double WZ0, WZp;
  int iZ0, iZp;
  getIndAndWZ(z, iZ0, iZp, WZ0, WZp);
  if(zSize_ > 1)
  {
    Moment[iZ0] += WZ0 * value;
    Moment[iZp] += WZp * value;
  }
  else if(zSize_ == 1)
  {
    Moment[iZ0] += value;
  }
}


/** Bins a moment to the grid with smoothing */
void Grid1D::binMomentSmoothed(double value, double z, double* Moment)
{
  if(z < zMin_ || z > zMax_ ) return;
  double  WZm,  WZ0,  WZp;
  double dWZm, dWZ0, dWZp;
  int iZm, iZ0, iZp;
  getIndAndWZSmoothed(z, iZm,  iZ0,  iZp,
                         WZm,  WZ0,  WZp,
                        dWZm, dWZ0, dWZp);
  if(zSize_ > 2)
  {
    Moment[iZm] += WZm * value;
    Moment[iZ0] += WZ0 * value;
    Moment[iZp] += WZp * value;
  }
  else if(zSize_ == 2)
  {
    Moment[iZ0] += WZ0 * value;
    Moment[iZp] += WZp * value;
  }
  else if(zSize_ == 1)
  {
    Moment[iZ0] += value;
  }
}


/** Calculates gradient at a position (z) */
void Grid1D::calcGradient(double z, double& ez)
{
  double WZ0, WZp;
  int iZ0, iZp;
  getIndAndWZ(z, iZ0, iZp, WZ0, WZp);
  if(zSize_ > 1)
  {
    ez = (arr_[iZp] - arr_[iZ0]) / dz_;
  }
  else if(zSize_ == 1)
  {
    ez = 0.;
  }
}


/** Calculates smoothed gradient at a position (z) */
void Grid1D::calcGradientSmoothed(double z, double& ez)
{
  double  WZm,  WZ0,  WZp;
  double dWZm, dWZ0, dWZp;
  int iZm, iZ0, iZp;
  getIndAndWZSmoothed(z, iZm,  iZ0,  iZp,
                         WZm,  WZ0,  WZp,
                        dWZm, dWZ0, dWZp);
  if(zSize_ > 2)
  {
    ez = dWZm * arr_[iZm] + dWZ0 * arr_[iZ0] + dWZp * arr_[iZp];
  }
  else if(zSize_ == 2)
  {
    ez = (arr_[iZp] - arr_[iZ0]) / dz_;
  }
  else if(zSize_ == 1)
  {
    ez = 0.;
  }
}


void Grid1D::getIndAndFracZ(double z, int& ind, double& frac)
{
  if(zSize_ > 1)
  {
    double zgloc = (z - zMin_) / dz_;
    ind = int(zgloc);
    frac = zgloc - double(ind);

    // Keep indices in bounds
    if(ind < 0) ind = zSize_ - 1;
    if(ind > (zSize_ - 1)) ind =  0;
  }
  else
  {
    ind = 0;
    frac = 1.0;
  }
}


void Grid1D::getIndAndFracZSmoothed(double z, int& ind, double& frac)
{
  if(zSize_ > 2)
  {
    double zgloc = (z - zMin_) / dz_;
    ind  = int(zgloc + 0.5);
    frac = zgloc - double(ind);

    // Keep indices in bounds
    if(ind < 0) ind = zSize_ - 1;
    if(ind > (zSize_ - 1)) ind =  0;
  }
  else if(zSize_ == 2)
  {
    ind = 0;
    frac = (z - zMin_) / dz_;
  }
  else if(zSize_ == 1)
  {
    ind = 0;
    frac = 1.0;
  }
}


/** Returns the grid indices and interpolation coefficients for a given z.
    The indices bracket the point of interpolation:
    0 <= iZ0, iZp <= nBins - 1
    The coefficients WZ0 and WZp correspond to iZ0 and iZp */
void Grid1D::getIndAndWZ(double z,
                         int& iZ0   , int& iZp,
                         double& WZ0, double& WZp)
{
  if(zSize_ > 1)
  {
    double zgloc = (z - zMin_) / dz_;
    iZ0 = int(zgloc);
    WZp = zgloc - double(iZ0);
    WZ0 = 1.0 - WZp;

    // Keep indices in bounds
    if(iZ0 < 0) iZ0 = zSize_ - 1;
    if(iZ0 > (zSize_ - 1)) iZ0 =  0;
    iZp = iZ0 + 1;
    if(iZ0 == (zSize_ - 1)) iZp = 0;
  }
  else
  {
    iZ0 = 0;
    iZp = 0;
    WZ0  = 1.0;
    WZp  = 0.0;
  }
}


/** Returns the grid index and fractional position for particular z.
    The central index is the central point in smoothed three
    point interpolation:
    0 <= iZm, iZ0, iZp <= nBins - 1
    The coefficients WZm, WZ0, and WZp
    correspond to iZm, iZ0, and iZp */
void Grid1D::getIndAndWZSmoothed(double z,
                                 int& iZm    , int& iZ0    , int& iZp    ,
                                 double& WZm , double& WZ0 ,  double& WZp,
                                 double& dWZm, double& dWZ0, double& dWZp)
{
  if(zSize_ > 2)
  {
	double zgloc = 0;
    //double zgloc = (z - zMin_) / dz_;
	if(length_ == 0)
    {
		std::cerr<<"Warning!: Lattice length should not be zero for smoothed binning.\n";
		zgloc = (z - zMin_) / dz_;
	}
	else
    {
		zgloc = (z + length_ / 2.0) / dz_;
	}
	  
    iZ0 = int(zgloc + 0.5);
    double frac = zgloc - double(iZ0);
    WZm  = 0.5 * (0.5 - frac) * (0.5 - frac);
    WZ0  = 0.75 - frac * frac;
    WZp  = 0.5 * (0.5 + frac) * (0.5 + frac);
    dWZm = (frac - 0.5) / dz_;
    dWZ0 = -2.0 * frac / dz_;
    dWZp = (frac + 0.5) / dz_;

    // Keep indices in bounds
    if(iZ0 < 0) iZ0 = zSize_ - 1;
    if(iZ0 > (zSize_ - 1)) iZ0 =  0;
    iZm = iZ0 - 1;
    if(iZ0 == 0) iZm = zSize_ - 1;
    iZp = iZ0 + 1;
    if(iZ0 == (zSize_ - 1)) iZp = 0;
  }
  else if(zSize_ == 2)
  {
    iZ0 = 0;
    iZm = 1;
    iZp = 1;
    WZm = 0.0;
    WZp = (z - zMin_) / dz_;
    WZ0 = 1.0 - WZp;
    dWZm =  0.0;
    dWZ0 = -1.0 / dz_;
    dWZp =  1.0 / dz_;
  }
  else if(zSize_ == 1)
  {
    iZm = 0;
    iZ0 = 0;
    iZp = 0;
    WZm = 0.0;
    WZ0 = 1.0;
    WZp = 0.0;
    dWZm = 0.0;
    dWZ0 = 0.0;
    dWZp = 0.0;
  }
}


/** synchronizeMPI */
void Grid1D::synchronizeMPI(pyORBIT_MPI_Comm* pyComm)
{
  // ====== MPI  start ========

  int size_MPI = zSize_;
  int buff_index0 = 0;
  int buff_index1 = 0;
  double* inArr =
    BufferStore::getBufferStore()->getFreeDoubleArr(buff_index0, size_MPI);
  double* outArr =
    BufferStore::getBufferStore()->getFreeDoubleArr(buff_index1, size_MPI);

  for(int i = 0; i < zSize_; i++)
  {
    inArr[i] = arr_[i];
  }

  if(pyComm == NULL)
  {
    ORBIT_MPI_Allreduce(inArr, outArr, size_MPI,
                        MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  }
  else
  {
    ORBIT_MPI_Allreduce(inArr, outArr, size_MPI,
                        MPI_DOUBLE, MPI_SUM, pyComm->comm);
  }

  for(int i = 0; i < zSize_; i++)
  {
    arr_[i] = outArr[i];
  }

  OrbitUtils::BufferStore::getBufferStore()->setUnusedDoubleArr(buff_index0);
  OrbitUtils::BufferStore::getBufferStore()->setUnusedDoubleArr(buff_index1);	
  // ===== MPI end =====
}
