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

}

void FieldTracker::trackBunch(Bunch* b) {

	cerr << "Ready to track in the 3D field!\n";

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
//				std::cerr << x <<", " << y << ", " << z << ", " <<Bx << ", " << By << ", " << Bz << "\n";


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
										BZGrid->setValue(Bx / 10000.0, xindex,
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

	std::cerr << "VALUE AT START is" << BXGrid->getValueOnGrid(0,0,0);


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

double * FieldTracker::BGrid3D(double xField3D,double yField3D,double zField3D,
		int nXGrid,int nYGrid,int nZGrid,
		double BxField3D,double ByField3D,double BzField3D,
		int zsymmetry)
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

	  if(ig < 0) ig = 0;
	  if(jg < 0) jg = 0;
	  if(kg < 0) kg = 0;
	  //could need to be changed to -2 do to array differences in C and SC
	  if(ig > nXGrid - 1) ig = nXGrid - 1;
	  if(jg > nYGrid - 1) jg = nYGrid - 1;
	  if(kg > nZGrid - 1) kg = nZGrid - 1;

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
	  double FieldArr[3];
	  FieldArr[0] = BxField3D;
	  FieldArr[1] = ByField3D;
	  FieldArr[2] = BzField3D;
	  return FieldArr;


}

