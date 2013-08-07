#include "FieldTracker.hh"

#include <iostream>
#include <fstream>
#include <cmath>
#include <cfloat>
#include <cstdlib>
#include "OrbitConst.hh"
#include "SyncPart.hh"

//////////////////////////////////////////////////////////////////////////
//Constructor for FieldTracker
//   bx     (double) --> _betaX   : The horizontal beta value at the
//                                beginning of the Node [m]
//   by     (double) --> _betaY   : The vertical beta value at the
//                                beginning of the Node [m]
//   ax     (double) --> _alphaX  : The horizontal alpha value at the
//                                beginning of the Node
//   ay     (double) --> _alphaY  : The horizontal alpha  value at the
//                                beginning of the Node
//   ex     (double) --> _etaX    : The horizontal dispersion [m]
//   epx    (double) --> _etaPX   : The horizontal dispersion prime
//   l      (double) --> _length  : The length of the node
//   zi     (double) --> _zi      : Initial tracking position [m]
//   zf     (double) --> _zf      : Final tracking position [m]
//   ds     (double) --> _ds      : Integration step size
//   niters  (int) --> _niters  : Number predictor-corrector iterations
//   resid  (double) --> _resid   : Predictor-corrector residual
//   xrefi  (double) --> _xrefi   : Initial reference particle x [mm]
//   yrefi  (double) --> _yrefi   : Initial reference particle y [mm]
//   eulerai(double) --> _eulerai : Initial reference Euler angle alpha [mr]
//   eulerbi(double) --> _eulerbi : Initial reference Euler angle beta [mr]
//   eulergi(double) --> _eulergi : Initial reference Euler angle gamma [mr]
//   b (bunch)                    : Bunch of particles
//	 filename(string)			  : string name of file to be parsed
//////////////////////////////////////////////////////////////////////////

//Constructor

FieldTracker::FieldTracker(const double &bx, const double &by,
	       const double &ax, const double &ay,
	       const double &ex, const double &epx,
	       const double &l,
	       const double &zi, const double &zf,
	       const double &ds, const int &niters,
	       const double &resid,
	       const double &xrefi, const double &yrefi,
	       const double &eulerai, const double &eulerbi,
	       const double &eulergi, Bunch* b, string &filename) {

	cerr << "Instantiating the 3D field track class \n" ;
	double ZPARSEMIN =  100.0 * zi - 1.0;
	double ZPARSEMAX =  100.0 * zf + 1.0;

	ParseGrid3D(filename, -21.0, 21.0, -13.0, 13.0, ZPARSEMIN, ZPARSEMAX, 1, 1, 1);
initVars();

_length = l;
_zi = zi;
_zf = zf;
_ds = ds;
_niters = niters;
_resid = resid;
_xrefi = xrefi;
_yrefi = yrefi;
_eulerai = eulerai;
_eulerbi = eulerbi;
_eulergi = eulergi;

_betaX = bx;
_betaY = by;
_alphaX = ax;
_alphaY = ay;
_etaX = ex;
_etaPX = epx;

nodeCalculator(b);

}

void FieldTracker::trackBunch(Bunch* b) {
	int i, j;
	double** part_coord_arr = b->coordArr();

	  if(_length < 0.0) return;
	  if(abs(_length) <= OrbitConst::tiny) return;

	  double ds, resid;
	  double dx,  dy,  dz,  dxo,  dyo,  dzo;
	  double dpx, dpy, dpz, dpxo, dpyo, dpzo;

	  double xrefi = 1.0e-03 * _xrefi;
	  double yrefi = 1.0e-03 * _yrefi;
	  double xreff = 1.0e-03 * _xreff;
	  double yreff = 1.0e-03 * _yreff;

	  double ca, sa, cb, sb, cg, sg;

	  double eulera = 1.0e-03 * _eulerai;
	  double eulerb = 1.0e-03 * _eulerbi;
	  double eulerg = 1.0e-03 * _eulergi;

	  ca = cos(eulera);
	  sa = sin(eulera);
	  cb = cos(eulerb);
	  sb = sin(eulerb);
	  cg = cos(eulerg);
	  sg = sin(eulerg);

	  double cxx   =  cb * ca * cg - sa * sg;
	  double cxy   = -cb * ca * sg - sa * cg;
	  double cxz   =  sb * ca;
	  double cyx   =  cb * sa * cg + ca * sg;
	  double cyy   = -cb * sa * sg + ca * cg;
	  double cyz   =  sb * sa;
	  double czx   = -sb * cg;
	  double czy   =  sb * sg;
	  double czz   =  cb;

	  eulera = 1.0e-03 * _euleraf;
	  eulerb = 1.0e-03 * _eulerbf;
	  eulerg = 1.0e-03 * _eulergf;

	  ca = cos(eulera);
	  sa = sin(eulera);
	  cb = cos(eulerb);
	  sb = sin(eulerb);
	  cg = cos(eulerg);
	  sg = sin(eulerg);


	  double dxx =  cb * ca * cg - sa * sg;
	  double dxy =  cb * sa * cg + ca * sg;
	  double dxz = -sb * cg;
	  double dyx = -cb * ca * sg - sa * cg;
	  double dyy = -cb * sa * sg + ca * cg;
	  double dyz =  sb * sg;
	  double dzx =  sb * ca;
	  double dzy =  sb * sa;
	  double dzz =  cb;

	  SyncPart* syncPart = b->getSyncPart();


	  double eTotal, pMomentum, coeff, pmag;
	  double Factor = 2.0 * OrbitConst::PI * 1.0 / _length;
	  double gamma2i = 1.0 / (syncPart->getGamma() * syncPart->getGamma());


	  double mpx, mpy, mpz, mpxp, mpyp, mpzp;
	  double xj, yj, zj, pxj, pyj, pzj, sj;
	  int lost = 0;

	    ofstream fio("Path", std::ios::out);

	  if(getPath == 1)
	  {
	    fio << "#s  x   y   z  px   py  pz" << "\n";
	  }
	  int printonce = 0;
	  for(j = 0; j < b->getSize(); j++)
	  {
	    if(getPath == 1)
	    {
	    	if(printonce<1){
	    	std::cerr << "Printing Path\n";
	    	printonce++;
	    	}
	    	fio<<"\n\n\n\n";
	    }
	    double beta = syncPart->getBeta();

	    double eTotal = syncPart->getEnergy() + syncPart->getMass();
	    double dpp = 1/(beta*beta)*part_coord_arr[j][5]/eTotal;


        pMomentum = sqrt(eTotal*eTotal - (syncPart->getMass()*syncPart->getEnergy())*(syncPart->getMass()*syncPart->getEnergy()));
        double e = 1.602e-19;
	    double amu = 1.66053886e-27;
	    double fac = e * 1.0e9 * pMomentum / OrbitConst::c;

	    coeff = 1.e-09 * OrbitConst::c * OrbitConst::charge_electron / pMomentum;



	    mpx  = 1.0e-03 * b->x(j);
	    mpy  = 1.0e-03 * b->y(j);
	    mpz  = -(czx * mpx + czy * mpy) / czz;
	    mpxp = 1.0e-03 * b->xp(j);
	    mpyp = 1.0e-03 * b->yp(j);
	    mpzp = 1.0 + dpp;
	    mpx  += mpxp * mpz / mpzp;
	    mpy  += mpyp * mpz / mpzp;


	    xj = xrefi + cxx * mpx + cxy * mpy + cxz * mpz;
	    yj = yrefi + cyx * mpx + cyy * mpy + cyz * mpz;
	    zj = _zi;
	    if(FieldTracker::foilAngle != 0)
	    {
	      zj = _zi + (xj - xrefi) / tan(FieldTracker::foilAngle * OrbitConst::PI*2 / 360.0);
	    }
	    pxj = cxx * mpxp + cxy * mpyp + cxz * mpzp;
	    pyj = cyx * mpxp + cyy * mpyp + cyz * mpzp;
	    pzj = czx * mpxp + czy * mpyp + czz * mpzp;
	    pmag = sqrt(pxj * pxj + pyj * pyj + pzj * pzj);
	    pxj /= pmag;
	    pyj /= pmag;
	    pzj /= pmag;

	    sj = 0.0;
	    lost = 0;

	    int iquit = 0;
	    int nreflections = 0;
	    double ppar0 = 0.0; double pperp0 = 0.0; double bmag0 = 0.0;
	    double firststep = 0;

	    double  rj = pow(xj * xj + yj * yj, 0.5);
	    while(iquit == 0)
	    {

	      ds = _ds;
	      dz = ds * pzj;
	      if(dz > (_zf - zj))
	      {
	        dz = _zf - zj;
	        ds = dz / pzj;
	        iquit = 1;
	      }
	      dx = ds * pxj;
	      dy = ds * pyj;

	      xField3D = xj + dx / 2.0;
	      yField3D = yj + dy / 2.0;
	      zField3D = zj + dz / 2.0;

	      BGrid3D();

	      BxField3D *= FieldTracker::BScale;
	      ByField3D *= FieldTracker::BScale;
	      BzField3D *= FieldTracker::BScale;

	      dpx = coeff * (dy * BzField3D - dz * ByField3D);
	      dpy = coeff * (dz * BxField3D - dx * BzField3D);
	      dpz = coeff * (dx * ByField3D - dy * BxField3D);



	      for(i = 1; i < _niters; i++)
	      {
	        dxo = dx;
	        dyo = dy;
	        dzo = dz;
	        dpxo = dpx;
	        dpyo = dpy;
	        dpzo = dpz;
	        ds = _ds;
	        dz = ds * (pzj + dpz / 2.0);


	        if(dz > (_zf - zj))
	        {
	          dz = _zf - zj;
	          ds = dz / (pzj + dpz / 2.0);
	          iquit = 1;
	        }
	        dx = ds * (pxj + dpx / 2.0);
	        dy = ds * (pyj + dpy / 2.0);

	        xField3D = xj + dx / 2.0;
	        yField3D = yj + dy / 2.0;
	        zField3D = zj + dz / 2.0;

	        BGrid3D();

	        BxField3D *= FieldTracker::BScale;
	        ByField3D *= FieldTracker::BScale;
	        BzField3D *= FieldTracker::BScale;

	        dpx = coeff * (dy * BzField3D - dz * ByField3D);
	        dpy = coeff * (dz * BxField3D - dx * BzField3D);
	        dpz = coeff * (dx * ByField3D - dy * BxField3D);
	        resid = sqrt((dx  -  dxo) * (dx - dxo)   +
	                     (dpx - dpxo) * (dpx - dpxo) +
	                     (dy  -  dyo) * (dy - dyo)   +
	                     (dpy - dpyo) * (dpy - dpyo) +
	                     (dz  -  dzo) * (dz - dzo)   +
	                     (dpz - dpzo) * (dpz - dpzo) );

	        if(resid < _resid) {
	        	break;
	        }
	      }

	      xj  +=  dx;
	      yj  +=  dy;
	      zj  +=  dz;
	      pxj += dpx;
	      pyj += dpy;
	      pzj += dpz;
	      sj +=  ds;

	      if(sj > FieldTracker::sMax * _length) iquit = 1;
	      rj = pow(xj * xj + yj * yj, 0.5);
	      if(rj > FieldTracker::rMax) iquit = 1;

	      if(getPath==1)
	      {
	        //Print out the particle coords
	        fio <<  sj << "  " <<  xj << "  " <<  yj << "  " <<  zj << "  "
	                           << pxj << "  " << pyj << "  " << pzj << "\n";
	      }

	    }

	      // If not lost, calc. final coords in ORBITs ref frame upon leaving 3D field.
	      dpx  = dxx * pxj + dxy * pyj + dxz * pzj;
	      dpy  = dyx * pxj + dyy * pyj + dyz * pzj;
	      dpz  = dzx * pxj + dzy * pyj + dzz * pzj;
	      dxo = xj - xreff;
	      dyo = yj - yreff;
	      dzo = dzx * dxo + dzy * dyo;

	      dx = dxx * dxo + dxy * dyo - dpx * dzo / dpz;
	      dy = dyx * dxo + dyy * dyo - dpy * dzo / dpz;
	      b->x(j)  = 1.0e+03 * dx;
	      b->xp(j) = 1.0e+03 * dpx * pmag;
	      b->y(j)  = 1.0e+03 * dy;
	      b->yp(j) = 1.0e+03 * dpy * pmag;
	      part_coord_arr[j][4] += Factor * (sj * (1.0 - gamma2i * dpp) -_sref);
	      if(part_coord_arr[j][4] >  OrbitConst::PI) part_coord_arr[j][4] -= OrbitConst::PI*2;
	      if(part_coord_arr[j][4] < -OrbitConst::PI) part_coord_arr[j][4] += OrbitConst::PI*2;

	      fio << b->x(j) << " " << b->y(j) << " " << b->xp(j) << " " << b->yp(j) << "\n" ;

	  }
	  fio.close();

	  double newtime = syncPart->getTime() + _length/( syncPart->getBeta()*OrbitConst::c );
	  syncPart->setTime(newtime);
}

////////////////////////////////////////////////////////////////////////////////
//
//   Routine to parse 3D grid field file.
//
// PARAMETERS
//   fileName -> File containing input field information.
//   xmin     -> Minimum x for data.
//   xmax     -> Maximum x for data.
//   ymin     -> Minimum y for data.
//   ymax     -> Maximum y for data.
//   zmin     -> Minimum z for data.
//   zmax     -> Maximum z for data.
//   skipX    -> Determines number of x gridpoints used.
//   skipY    -> Determines number of y gridpoints used.
//   skipZ    -> Determines number of z gridpoints used.
//
////////////////////////////////////////////////////////////////////////////////

void FieldTracker::ParseGrid3D(const string &fileName, const double &xmin,
		const double &xmax, const double &ymin, const double &ymax,
		const double &zmin, const double &zmax, const int &skipX,
		const int &skipY, const int &skipZ) {
	int nXTab = 0;
	int nYTab = 0;
	int nZTab = 0;
	int iDummy, i, j, k;
	int xindex, yindex, zindex;
	string Dummy;
	double x, y, z, Bx, By, Bz;


	std::cerr << "Filename: " << fileName << "\n";

	ifstream fio(fileName.c_str(), ios::in);
	if (!fio) {
		std::cerr << "Filename " << fileName << " not found\n";
	}

	fio >> nXTab >> nYTab >> nZTab >> iDummy;
	fio >> iDummy >> Dummy >> Dummy;
	fio >> iDummy >> Dummy >> Dummy;
	fio >> iDummy >> Dummy >> Dummy;
	fio >> iDummy >> Dummy >> Dummy;
	fio >> iDummy >> Dummy >> Dummy;
	fio >> iDummy >> Dummy >> Dummy;
	fio >> iDummy;

	nXGrid = 0;
	for (i = 0; i < nXTab; i++) {
		nYGrid = 0;
		for (j = 0; j < nYTab; j++) {
			nZGrid = 0;
			for (k = 0; k < nZTab; k++) {
				fio >> x >> y >> z >> Bx >> By >> Bz;

				if ((z >= zmin) && (z <= zmax)) {
					if ((k == (k / skipZ) * skipZ) || (k == nZTab - 1)) {
						nZGrid++;
					}
				}
			}
			if ((y >= ymin) && (y <= ymax)) {
				if ((j == (j / skipY) * skipY) || (j == nYTab - 1)) {
					nYGrid++;
				}
			}
		}
		if ((x >= xmin) && (x <= xmax)) {
			if ((i == (i / skipX) * skipX) || (i == nXTab - 1)) {
				nXGrid++;
			}
		}
	}
	fio.close();

	XGrid = new double[nXTab];
	YGrid = new double[nYTab];
	ZGrid = new double[nZTab];
	BXGrid = new Grid3D(nXTab, nYTab, nZTab);
	BYGrid = new Grid3D(nXTab, nYTab, nZTab);
	BZGrid = new Grid3D(nXTab, nYTab, nZTab);
	BMagGrid = new Grid3D(nXTab, nYTab, nZTab);

	ifstream fio2(fileName.c_str(), ios::in);
	if (!fio2) {
		std::cerr << "Filename " << fileName << " not found\n";
	}

	fio2 >> nXTab >> nYTab >> nZTab >> iDummy;
	fio2 >> iDummy >> Dummy >> Dummy;
	fio2 >> iDummy >> Dummy >> Dummy;
	fio2 >> iDummy >> Dummy >> Dummy;
	fio2 >> iDummy >> Dummy >> Dummy;
	fio2 >> iDummy >> Dummy >> Dummy;
	fio2 >> iDummy >> Dummy >> Dummy;
	fio2 >> iDummy;


	xindex = 0;
	for (i = 0; i < nXTab; i++) {
		yindex = 0;
		for (j = 0; j < nYTab; j++) {
			zindex = 0;
			for (k = 0; k < nZTab; k++) {
				fio2 >> x >> y >> z >> Bx >> By >> Bz;

				if ((z >= zmin) && (z <= zmax)) {
					if ((k == (k / skipZ) * skipZ) || (k == nZTab - 1)) {
						if ((y >= ymin) && (y <= ymax)) {
							if ((j == (j / skipY) * skipY)
									|| (j == nYTab - 1)) {
								if ((x >= xmin) && (x <= xmax)) {
									if ((i == (i / skipX) * skipX)
											|| (i == nXTab - 1)) {
										XGrid[xindex] = x / 100.0;
										YGrid[yindex] = y / 100.0;
										ZGrid[zindex] = z / 100.0;
										BXGrid->setValue(Bx / 10000.0, xindex,
												yindex, zindex);
										BYGrid->setValue(By / 10000.0, xindex,
												yindex, zindex);
										BZGrid->setValue(Bz / 10000.0, xindex,
												yindex, zindex);
									}
								}
							}
						}
						zindex++;
					}
				}
			}
			if ((y >= ymin) && (y <= ymax)) {
				if ((j == (j / skipY) * skipY) || (j == nYTab - 1)) {
					yindex++;
				}
			}
		}
		if ((x >= xmin) && (x <= xmax)) {
			if ((i == (i / skipX) * skipX) || (i == nXTab - 1)) {
				xindex++;
			}
		}
	}


	fio2.close();

}

///////////////////////////////////////////////////////////////////////////
//
//
// NAME
//
//   FieldTracker::BGrid3D
//
// DESCRIPTION
//   Routine to interpolate B from 3D grid for 3D Tracker.
//
// PARAMETERS
//   None
//
// RETURNS
//   Nothing
//
///////////////////////////////////////////////////////////////////////////

void FieldTracker::BGrid3D()
{
	  double xF, yF, zF;

	  xF = xField3D;
	  yF = yField3D;
	  zF = zField3D;

	  if((zsymmetry != 0) && (zF < ZGrid[0]))
	  {
	    zF = ZGrid[0] - zF + ZGrid[0];
	  }

	  int ig, jg, kg, igp, jgp, kgp;
	  double dx, dy, dz, xfrac, yfrac, zfrac, xfracc, yfracc, zfracc;

	  dx = (XGrid[nXGrid-1] - XGrid[0]) / double(nXGrid - 1);
	  dy = (YGrid[nYGrid-1] - YGrid[0]) / double(nYGrid - 1);
	  dz = (ZGrid[nZGrid-1] - ZGrid[0]) / double(nZGrid - 1);

	  ig = int((xF - XGrid[0]) / dx) + 1;
	  jg = int((yF - YGrid[0]) / dy) + 1;
	  kg = int((zF - ZGrid[0]) / dz) + 1;

	  if(ig < 0) ig = 1;
	  if(jg < 0) jg = 1;
	  if(kg < 0) kg = 1;
	  //could need to be changed to -2 do to array differences in C and SC
	  if(ig > nXGrid - 2) ig = nXGrid - 2;
	  if(jg > nYGrid - 2) jg = nYGrid - 2;
	  if(kg > nZGrid - 2) kg = nZGrid - 2;


	  ig = ig-1;
	  jg = jg-1;
	  kg = kg-1;

	  igp = ig + 1;
	  jgp = jg + 1;
	  kgp = kg + 1;

	  xfrac = (xF - XGrid[ig]) / dx;
	  yfrac = (yF - YGrid[jg]) / dy;
	  zfrac = (zF - ZGrid[kg]) / dz;

	  xfracc = 1.0 - xfrac;
	  yfracc = 1.0 - yfrac;
	  zfracc = 1.0 - zfrac;

	  BxField3D =   BXGrid->getValueOnGrid(ig , jg , kg ) * xfracc * yfracc * zfracc +
	              + BXGrid->getValueOnGrid(ig , jg , kgp) * xfracc * yfracc * zfrac  +
	              + BXGrid->getValueOnGrid(ig , jgp, kg ) * xfracc * yfrac  * zfracc +
	              + BXGrid->getValueOnGrid(ig , jgp, kgp) * xfracc * yfrac  * zfrac  +
	              + BXGrid->getValueOnGrid(igp, jg , kg ) * xfrac  * yfracc * zfracc +
	              + BXGrid->getValueOnGrid(igp, jg , kgp) * xfrac  * yfracc * zfrac  +
	              + BXGrid->getValueOnGrid(igp, jgp, kg ) * xfrac  * yfrac  * zfracc +
	              + BXGrid->getValueOnGrid(igp, jgp, kgp) * xfrac  * yfrac  * zfrac;

	  ByField3D =   BYGrid->getValueOnGrid(ig , jg , kg ) * xfracc * yfracc * zfracc +
	              + BYGrid->getValueOnGrid(ig , jg , kgp) * xfracc * yfracc * zfrac  +
	              + BYGrid->getValueOnGrid(ig , jgp, kg ) * xfracc * yfrac  * zfracc +
	              + BYGrid->getValueOnGrid(ig , jgp, kgp) * xfracc * yfrac  * zfrac  +
	              + BYGrid->getValueOnGrid(igp, jg , kg ) * xfrac  * yfracc * zfracc +
	              + BYGrid->getValueOnGrid(igp, jg , kgp) * xfrac  * yfracc * zfrac  +
	              + BYGrid->getValueOnGrid(igp, jgp, kg ) * xfrac  * yfrac  * zfracc +
	              + BYGrid->getValueOnGrid(igp, jgp, kgp) * xfrac  * yfrac  * zfrac;

	  BzField3D =   BZGrid->getValueOnGrid(ig , jg , kg ) * xfracc * yfracc * zfracc +
	              + BZGrid->getValueOnGrid(ig , jg , kgp) * xfracc * yfracc * zfrac  +
	              + BZGrid->getValueOnGrid(ig , jgp, kg ) * xfracc * yfrac  * zfracc +
	              + BZGrid->getValueOnGrid(ig , jgp, kgp) * xfracc * yfrac  * zfrac  +
	              + BZGrid->getValueOnGrid(igp, jg , kg ) * xfrac  * yfracc * zfracc +
	              + BZGrid->getValueOnGrid(igp, jg , kgp) * xfrac  * yfracc * zfrac  +
	              + BZGrid->getValueOnGrid(igp, jgp, kg ) * xfrac  * yfrac  * zfracc +
	              + BZGrid->getValueOnGrid(igp, jgp, kgp) * xfrac  * yfrac  * zfrac;

	  if((zsymmetry != 0) && (zF < ZGrid[1]))
	  {
	    BzField3D = -BzField3D;
	  }
}

//////////////////////////////////////////////////////////////////////////////
//Sets initial variables
//////////////////////////////////////////////////////////////////////////////
void FieldTracker::initVars()
{
// set some initial values
  zsymmetry = 0;
  getPath = 0;
  BScale = 1.0;
  sMax = 10.0;
  rMax = 1.0e+10;
  xFoilMin = 0.0;
  xFoilMax = 0.0;
  yFoilMin = 0.0;
  yFoilMax = 0.0;
  zFoilMin = 0.0;
  zFoilMax = 0.0;
  _euleraf = 0;
  _eulerbf = 0;
  _eulergf = 0;
  BxField3D = 0;
  ByField3D = 0;
  BzField3D = 0;
  foilAngle = 0.0;
  rMax = 1.0e+10;
  doRefPath = 0;
  getPath = 0;

}

////////////////////////////////////////////////////////////////////////////////
//
//Tracks the reference particle and computes the final Euler Angles.
//Arguments: Bunch* B - the bunch to be tracked
//
////////////////////////////////////////////////////////////////////////////////
void FieldTracker::nodeCalculator(Bunch* b)
{

 //  If tune calculation is on, check storage sizes:
  int oldSize, i;

  double betaX  = _betaX;
  double betaY  = _betaY;
  double alphaX = _alphaX;
  double alphaY = _alphaY;
  double etaX   = _etaX;
  double etaPX  = _etaPX;
  double gammaX = (1. + (alphaX*alphaX))/betaX;
  double gammaY = (1. + (alphaY*alphaY))/betaY;

	double** part_coord_arr = b->coordArr();
	SyncPart* syncPart = b->getSyncPart();
  if(_length < 0.0) return;
  if(abs(_length) <= OrbitConst::tiny) return;

  double eTotal = syncPart->getEnergy() + syncPart->getMass();


  double pMomentum = sqrt(eTotal*eTotal - (syncPart->getMass()*syncPart->getEnergy())*(syncPart->getMass()*syncPart->getEnergy()));

  double coeff = 1.e-09 * OrbitConst::c * OrbitConst::charge_electron / pMomentum;

  //MIGHT NEED TO BE FIXED ^^^^
  double ds, resid;
  double  dx,  dy,  dz,  dxo,  dyo,  dzo;
  double dpx, dpy, dpz, dpxo, dpyo, dpzo;

  double eulera = 1.0e-03 * _eulerai;
  double eulerb = 1.0e-03 * _eulerbi;
  double eulerg = 1.0e-03 * _eulergi;



  double  xref = 1.0e-03 * _xrefi;
  double  yref = 1.0e-03 * _yrefi;
  double  zref = _zi;
  double pxref = sin(eulerb) * cos(eulera);
  double pyref = sin(eulerb) * sin(eulera);
  double pzref = cos(eulerb);

  _sref = 0.0;
  int iquit = 0;
  double xwidth = FieldTracker::xFoilMax - FieldTracker::xFoilMin;

 ofstream fio("RefPath", ios::out);

  xField3D = xref;
  yField3D = yref;
  zField3D = zref;

  BGrid3D();
  BxField3D *= FieldTracker::BScale;
  ByField3D *= FieldTracker::BScale;
  BzField3D *= FieldTracker::BScale;

  fio << _sref << "  " <<  xref << "  " <<  yref << "  " <<  zref << "  "
                       << pxref << "  " << pyref << "  " << pzref << "  "
                       << BxField3D << "  "
                       << ByField3D << "  "
                       << BzField3D << "\n";

  if(FieldTracker::doRefPath != 0) iquit = 1;


  double  rref = pow(xref * xref + yref * yref, 0.5);

  while(iquit == 0)
  {
    ds = _ds;
    dz = ds * pzref;
    if(dz > (_zf - zref))
    {
      dz = _zf - zref;
      ds = dz / pzref;
      iquit = 1;
    }
    dx = ds * pxref;
    dy = ds * pyref;

    xField3D = xref + dx / 2.0;
    yField3D = yref + dy / 2.0;
    zField3D = zref + dz / 2.0;
    BGrid3D();


    BxField3D *= FieldTracker::BScale;
    ByField3D *= FieldTracker::BScale;
    BzField3D *= FieldTracker::BScale;
    dpx = coeff * (dy * BzField3D - dz * ByField3D);
    dpy = coeff * (dz * BxField3D - dx * BzField3D);
    dpz = coeff * (dx * ByField3D - dy * BxField3D);


    for(i = 1; i < _niters; i++)
    {
      dxo = dx;
      dyo = dy;
      dzo = dz;
      dpxo = dpx;
      dpyo = dpy;
      dpzo = dpz;

      ds = _ds;
      dz = ds * (pzref + dpz / 2.0);
      if(dz > (_zf - zref))
      {
        dz = _zf - zref;
        ds = dz / (pzref + dpz / 2.0);
        iquit = 1;
      }
      dx = ds * (pxref + dpx / 2.0);
      dy = ds * (pyref + dpy / 2.0);

      xField3D = xref + dx / 2.0;
      yField3D = yref + dy / 2.0;
      zField3D = zref + dz / 2.0;
      BGrid3D();
      BxField3D *= FieldTracker::BScale;
      ByField3D *= FieldTracker::BScale;
      BzField3D *= FieldTracker::BScale;
      dpx = coeff * (dy * BzField3D - dz * ByField3D);
      dpy = coeff * (dz * BxField3D - dx * BzField3D);
      dpz = coeff * (dx * ByField3D - dy * BxField3D);
      resid = sqrt((dx  -  dxo) * (dx - dxo)   +
                   (dpx - dpxo) * (dpx - dpxo) +
                   (dy  -  dyo) * (dy - dyo)   +
                   (dpy - dpyo) * (dpy - dpyo) +
                   (dz  -  dzo) * (dz - dzo)   +
                   (dpz - dpzo) * (dpz - dpzo) );

      if(resid < _resid){
    	  break;
      }

    }

    xref  +=  dx;
    pxref += dpx;
    yref  +=  dy;
    pyref += dpy;
    zref  +=  dz;
    pzref += dpz;
    _sref +=  ds;

    if(_sref > FieldTracker::sMax * _length) {
    	iquit = 1;
    }
    rref = pow(xref * xref + yref * yref, 0.5);
    if(rref > FieldTracker::rMax) {
    	iquit = 1;
    }

    fio << _sref << "  " <<  xref << "  " <<  yref << "  " <<  zref << "  "
                         << pxref << "  " << pyref << "  " << pzref << "  "
                         << BxField3D << "  "
                         << ByField3D << "  "
                         << BzField3D << "\n";
  }

  fio.close();

  _xreff = 1.e+03 * xref;
  _yreff = 1.e+03 * yref;

  eulera = atan2(pyref, pxref);
  eulerb = atan2(sqrt(pyref * pyref + pxref * pxref), pzref);
  _euleraf = 1.0e+03 * eulera;
  _eulerbf = 1.0e+03 * eulerb;
  _eulergf = -_euleraf;
}

void FieldTracker::setPathVariable(int i){
	if (i == 1){
	getPath = 1;
	doRefPath = 1;
	}
	else {
		getPath = 0;
		doRefPath = 0;
	}
}
