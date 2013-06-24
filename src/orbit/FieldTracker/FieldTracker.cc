#include "FieldTracker.hh"

#include <iostream>
#include <fstream>
#include <cmath>
#include <cfloat>
#include <cstdlib>


//////////////////////////////////////////////////////////////////////////
//
// Comments and Stuff
//
//
//
//
//////////////////////////////////////////////////////////////////////////

//Constructor

FieldTracker::FieldTracker(double a) {

	cerr << "Instantiating the 3D field track class \n";
	ParseGrid3D("testFile.data", -200.0, 200.0, -200.0, 200.0, 0.0, 2900.0, 1, 1, 1);
//	BGrid3D(150.4732-92.55,46.0-23.0,0.307193,1);


}

void FieldTracker::trackBunch(Bunch* b) {

	cerr << "Ready to track in the 3D field!\n";
//	int i, j;
//
//	  if(_length < 0.0) return;
//	  if(abs(_length) <= Consts::tiny) return;
//
//	  double ds, resid;
//	  double dx,  dy,  dz,  dxo,  dyo,  dzo;
//	  double dpx, dpy, dpz, dpxo, dpyo, dpzo;
//
//	  double xrefi = 1.0e-03 * _xrefi;
//	  double yrefi = 1.0e-03 * _yrefi;
//	  double xreff = 1.0e-03 * _xreff;
//	  double yreff = 1.0e-03 * _yreff;
//
//	  double ca, sa, cb, sb, cg, sg;
//
//	  double eulera = 1.0e-03 * _eulerai;
//	  double eulerb = 1.0e-03 * _eulerbi;
//	  double eulerg = 1.0e-03 * _eulergi;
//
//	  ca = cos(eulera);
//	  sa = sin(eulera);
//	  cb = cos(eulerb);
//	  sb = sin(eulerb);
//	  cg = cos(eulerg);
//	  sg = sin(eulerg);
//
//	  double cxx   =  cb * ca * cg - sa * sg;
//	  double cxy   = -cb * ca * sg - sa * cg;
//	  double cxz   =  sb * ca;
//	  double cyx   =  cb * sa * cg + ca * sg;
//	  double cyy   = -cb * sa * sg + ca * cg;
//	  double cyz   =  sb * sa;
//	  double czx   = -sb * cg;
//	  double czy   =  sb * sg;
//	  double czz   =  cb;
//
//	  eulera = 1.0e-03 * _euleraf;
//	  eulerb = 1.0e-03 * _eulerbf;
//	  eulerg = 1.0e-03 * _eulergf;
//
//	  ca = cos(eulera);
//	  sa = sin(eulera);
//	  cb = cos(eulerb);
//	  sb = sin(eulerb);
//	  cg = cos(eulerg);
//	  sg = sin(eulerg);
//
//	  double dxx =  cb * ca * cg - sa * sg;
//	  double dxy =  cb * sa * cg + ca * sg;
//	  double dxz = -sb * cg;
//	  double dyx = -cb * ca * sg - sa * cg;
//	  double dyy = -cb * sa * sg + ca * cg;
//	  double dyz =  sb * sg;
//	  double dzx =  sb * ca;
//	  double dzy =  sb * sa;
//	  double dzz =  cb;
//
//	  double eTotal, pMomentum, coeff, pmag;
//	  double Factor = 2.0 * Consts::pi * Ring::harmonicNumber / Ring::lRing;
//	  double gamma2i = 1.0 / (mp._syncPart._gammaSync * mp._syncPart._gammaSync);
//
//
//	  double mpx, mpy, mpz, mpxp, mpyp, mpzp;
//	  double xj, yj, zj, pxj, pyj, pzj, sj;
//	  int lost = 0;
//
//	    OFstream fio("Path", std::ios::out);
//
//	  if(FieldTracker::getPath == 1)
//	  {
//	    fio << "#s  x   y   z  px   py  pz" << "\n";
//	  }
//
//	  for(j = 1; j <= mp._nMacros; j++)
//	  {
//	    //std::cerr << "Next particle" << "\n";
//	    if(FieldTracker::getPath == 1)
//	    {
//	      fio<<"\n\n\n\n";
//	    }
//
//	    eTotal = mp._syncPart._eTotal + mp._deltaE(j);
//	    pMomentum = Sqrt(Sqr(eTotal) - Sqr(mp._syncPart._e0));
//	    double e = 1.602e-19;
//	    double amu = 1.66053886e-27;
//	    double fac = e * 1.0e9 * pMomentum / Consts::vLight;
//	    coeff = 1.e-09 * Consts::vLight * mp._syncPart._charge / pMomentum;
//	    mpx  = 1.0e-03 * mp._x(j);
//	    mpy  = 1.0e-03 * mp._y(j);
//	    mpz  = -(czx * mpx + czy * mpy) / czz;
//	    mpxp = 1.0e-03 * mp._xp(j);
//	    mpyp = 1.0e-03 * mp._yp(j);
//	    mpzp = 1.0 + mp._dp_p(j);
//	    mpx  += mpxp * mpz / mpzp;
//	    mpy  += mpyp * mpz / mpzp;
//
//	    xj = xrefi + cxx * mpx + cxy * mpy + cxz * mpz;
//	    yj = yrefi + cyx * mpx + cyy * mpy + cyz * mpz;
//	    zj = _zi;
//	    if(FieldTracker::foilAngle != 0)
//	    {
//	      zj = _zi + (xj - xrefi) / tan(FieldTracker::foilAngle * Consts::twoPi / 360.0);
//	    }
//	    pxj = cxx * mpxp + cxy * mpyp + cxz * mpzp;
//	    pyj = cyx * mpxp + cyy * mpyp + cyz * mpzp;
//	    pzj = czx * mpxp + czy * mpyp + czz * mpzp;
//	    pmag = Sqrt(pxj * pxj + pyj * pyj + pzj * pzj);
//	    pxj /= pmag;
//	    pyj /= pmag;
//	    pzj /= pmag;
//
//	    sj = 0.0;
//	    lost = 0;
//
//	    int iquit = 0;
//	    int nreflections = 0;
//	    double ppar0 = 0.0; double pperp0 = 0.0; double bmag0 = 0.0;
//	    double firststep = 0;
//
//	    double  rj = pow(xj * xj + yj * yj, 0.5);
//
//	    while(iquit == 0)
//	    {
//
//	      ds = _ds;
//	      dz = ds * pzj;
//	      if(dz > (_zf - zj))
//	      {
//	        dz = _zf - zj;
//	        ds = dz / pzj;
//	        iquit = 1;
//	      }
//	      dx = ds * pxj;
//	      dy = ds * pyj;
//
//	      xField3D = xj + dx / 2.0;
//	      yField3D = yj + dy / 2.0;
//	      zField3D = zj + dz / 2.0;
//	      _sub();
//	      BxField3D *= FieldTracker::BScale;
//	      ByField3D *= FieldTracker::BScale;
//	      BzField3D *= FieldTracker::BScale;
//
//	      dpx = coeff * (dy * BzField3D - dz * ByField3D);
//	      dpy = coeff * (dz * BxField3D - dx * BzField3D);
//	      dpz = coeff * (dx * ByField3D - dy * BxField3D);
//
//	      for(i = 1; i < _niters; i++)
//	      {
//	        dxo = dx;
//	        dyo = dy;
//	        dzo = dz;
//	        dpxo = dpx;
//	        dpyo = dpy;
//	        dpzo = dpz;
//	        ds = _ds;
//	        dz = ds * (pzj + dpz / 2.0);
//	        if(dz > (_zf - zj))
//	        {
//	          dz = _zf - zj;
//	          ds = dz / (pzj + dpz / 2.0);
//	          iquit = 1;
//	        }
//	        dx = ds * (pxj + dpx / 2.0);
//	        dy = ds * (pyj + dpy / 2.0);
//
//	        xField3D = xj + dx / 2.0;
//	        yField3D = yj + dy / 2.0;
//	        zField3D = zj + dz / 2.0;
//
//	        _sub();
//	        BxField3D *= FieldTracker::BScale;
//	        ByField3D *= FieldTracker::BScale;
//	        BzField3D *= FieldTracker::BScale;
//	        dpx = coeff * (dy * BzField3D - dz * ByField3D);
//	        dpy = coeff * (dz * BxField3D - dx * BzField3D);
//	        dpz = coeff * (dx * ByField3D - dy * BxField3D);
//	        resid = Sqrt((dx  -  dxo) * (dx - dxo)   +
//	                     (dpx - dpxo) * (dpx - dpxo) +
//	                     (dy  -  dyo) * (dy - dyo)   +
//	                     (dpy - dpyo) * (dpy - dpyo) +
//	                     (dz  -  dzo) * (dz - dzo)   +
//	                     (dpz - dpzo) * (dpz - dpzo) );
//	        if(resid < _resid) break;
//	      }
//
//	      xj  +=  dx;
//	      yj  +=  dy;
//	      zj  +=  dz;
//	      pxj += dpx;
//	      pyj += dpy;
//	      pzj += dpz;
//	      sj +=  ds;
//
//	      if(sj > FieldTracker::sMax * _length) iquit = 1;
//	      rj = pow(xj * xj + yj * yj, 0.5);
//	      if(rj > FieldTracker::rMax) iquit = 1;
//
//	      if(FieldTracker::getPath==1)
//	      {
//	        //Print out the particle coords
//	        fio <<  sj << "  " <<  xj << "  " <<  yj << "  " <<  zj << "  "
//	                           << pxj << "  " << pyj << "  " << pzj << "\n";
//	      }
//
//	    }
//
//
//	      // If not lost, calc. final coords in ORBITs ref frame upon leaving 3D field.
//	      dpx  = dxx * pxj + dxy * pyj + dxz * pzj;
//	      dpy  = dyx * pxj + dyy * pyj + dyz * pzj;
//	      dpz  = dzx * pxj + dzy * pyj + dzz * pzj;
//	      dxo = xj - xreff;
//	      dyo = yj - yreff;
//	      dzo = dzx * dxo + dzy * dyo;
//
//	      dx = dxx * dxo + dxy * dyo - dpx * dzo / dpz;
//	      dy = dyx * dxo + dyy * dyo - dpy * dzo / dpz;
//	      mp._x(j)  = 1.0e+03 * dx;
//	      mp._xp(j) = 1.0e+03 * dpx * pmag;
//	      mp._y(j)  = 1.0e+03 * dy;
//	      mp._yp(j) = 1.0e+03 * dpy * pmag;
//	      mp._phi(j) += Factor * (sj * (1.0 - gamma2i * mp._dp_p(j)) -_sref);
//	      if(mp._phi(j) >  Consts::pi) mp._phi(j) -= Consts::twoPi;
//	      if(mp._phi(j) < -Consts::pi) mp._phi(j) += Consts::twoPi;
//
//	  }
//
//	  // particle tune calculation:
//
//	  if(!TransMap::tuneCalcOn) return; // return if tune calculation is off
//
//	  FTdoTunes(*this, mp);
//	  fio.close();
//




}





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

	std::cerr << nXTab << "  " << nYTab << " " << nZTab << "\n";
	std::cerr << "Filename: " << fileName;

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
	fio >> Dummy;

	std::cerr << "\nnow\n" << Dummy << "\n";
	std::cerr << nXTab << "  " << nYTab << " " << nZTab;

	nXGrid = 0;
	for (i = 0; i < nXTab; i++) {
		nYGrid = 0;
		for (j = 0; j < nYTab; j++) {
			nZGrid = 0;
			for (k = 0; k < nZTab; k++) {
				fio >> x >> y >> z >> Bx >> By >> Bz;
//				std::cerr << x <<", " << y << ", " << z << ", " <<Bx << ", " << By << ", " << Bz << "\n";

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
	fio2 >> Dummy;


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

void FieldTracker::BGrid3D(double xField3D,double yField3D,double zField3D,
		int zsymmetry)
{
	  double xF, yF, zF;

	  xF = xField3D;
	  yF = yField3D;
	  zF = zField3D;

	  std::cerr << "\nBeginning BGrid Field Values: " << xField3D <<", " << yField3D << ", " << zField3D << "\n";

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


	  if(ig < 0) ig = 0;
	  if(jg < 0) jg = 0;
	  if(kg < 0) kg = 0;
	  //could need to be changed to -2 do to array differences in C and SC
	  if(ig > nXGrid - 2) ig = nXGrid - 1;
	  if(jg > nYGrid - 2) jg = nYGrid - 1;
	  if(kg > nZGrid - 2) kg = nZGrid - 1;


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
	  std::cerr << "\nBxField3D: " << BxField3D << " ByField3D: " << ByField3D << " BzField3D: " << BzField3D << "\n";

}

