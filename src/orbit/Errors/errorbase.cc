/////////////////////////////////////////////////////////////////////////////
//
// FILE NAME
//   errorbase.cc
//
// AUTHOR
//   Jeff Holmes, ORNL
//   Steven Bunch, UTK
//
// CREATED IN ORBIT: 6/17/02
//
// MODIFIED FOR py_ORBIT: 07/2015
//
//  DESCRIPTION
//    Define functions for magnet alignment errors
//
/////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////
//
// include files
//
/////////////////////////////////////////////////////////////////////////////

#include "errorbase.hh"
#include "OrbitConst.hh"
#include "SyncPart.hh"
#include "Random.hh"

#include <complex>
#include <iostream>
#include <cmath>
#include <cfloat>
#include <cstdlib>

namespace error_base
{

///////////////////////////////////////////////////////////////////////////
//
// NAME
//   CoordDisplacement
//
// DESCRIPTION
//   Performs generalized coordinate displacement
//
// PARAMETERS
//   dx, dxp, dy, dyp, dphi, dE
//
///////////////////////////////////////////////////////////////////////////

void CoordDisplacement(Bunch* bunch,
                       double dx, double dxp,
                       double dy, double dyp,
                       double dz, double dE)
{
  //coordinate array [part. index][x,xp,y,yp,z,dE]
  double** arr = bunch->coordArr();

  for(int i = 0; i < bunch->getSize(); i++)
  {
    arr[i][0] += dx;
    arr[i][1] += dxp;
    arr[i][2] += dy;
    arr[i][3] += dyp;
    arr[i][4] += dz;
    arr[i][5] += dE;
  }
}

///////////////////////////////////////////////////////////////////////////
//
// NAME
//   LongDisplacement
//
// DESCRIPTION
//   Performs longitudinal displacement
//
// PARAMETERS
//   ds
//
///////////////////////////////////////////////////////////////////////////

void LongDisplacement(Bunch* bunch, double ds)
{
  //coordinate array [part. index][x,xp,y,yp,z,dE]
  double** arr = bunch->coordArr();

  for(int i = 0; i < bunch->getSize(); i++)
  {
    drifti(bunch, i, ds);
  }
}

///////////////////////////////////////////////////////////////////////////
//
// NAME
//   StraightRotationXY
//
// DESCRIPTION
//   Performs X-Y coordinate rotation
//
// PARAMETERS
//   anglexy
//
///////////////////////////////////////////////////////////////////////////

void StraightRotationXY(Bunch* bunch, double anglexy)
{

  double xtemp, xptemp, ytemp, yptemp;
  double cs = cos(anglexy);
  double sn = sin(anglexy);

  //coordinate array [part. index][x,xp,y,yp,z,dE]
  double** arr = bunch->coordArr();

  for(int i = 0; i < bunch->getSize(); i++)
  {
    xtemp  = arr[i][0];
    xptemp = arr[i][1];
    ytemp  = arr[i][2];
    yptemp = arr[i][3];

    arr[i][0] =  cs * xtemp  + sn * ytemp;
    arr[i][1] =  cs * xptemp + sn * yptemp;
    arr[i][2] = -sn * xtemp  + cs * ytemp;
    arr[i][3] = -sn * xptemp + cs * yptemp;
  }
}

///////////////////////////////////////////////////////////////////////////
//
// NAME
//   StraightRotationXSI
//
// DESCRIPTION
//   Performs X-S coordinate initial rotation of element
//
// PARAMETERS
//   anglexsi
//   length
//
///////////////////////////////////////////////////////////////////////////

void StraightRotationXSI(Bunch* bunch, double anglexsi, double length)
{
  double xtemp, xptemp, ytemp, yptemp, ztemp, zptemp, z, zp, s;

  double cs = cos(anglexsi);
  double sn = sin(anglexsi);
  double lengtho2 = length / 2.0;

  SyncPart* syncPart = bunch->getSyncPart();

  double dp_p_coeff = 1.0 / (syncPart->getMomentum() * syncPart->getBeta());

  //coordinate array [part. index][x,xp,y,yp,z,dE]
  double** arr = bunch->coordArr();

  for(int i = 0; i < bunch->getSize(); i++)
  {
    xtemp  = arr[i][0];
    xptemp = arr[i][1];
    ytemp  = arr[i][2];
    yptemp = arr[i][3];
    ztemp  = arr[i][4];
    zptemp = 1.0 + arr[i][5] * dp_p_coeff;

    arr[i][0] = cs * xtemp  + sn * lengtho2;
    arr[i][1] = cs * xptemp - sn * zptemp;
    arr[i][4] = sn * xtemp  + cs * ztemp;
    zp        = sn * xptemp + cs * zptemp;
    arr[i][5] = (zp - 1.0) / dp_p_coeff;

    s = -(1.0 - cs) * lengtho2 - sn * xtemp;

    drifti(bunch, i, s);
  }
}

///////////////////////////////////////////////////////////////////////////
//
// NAME
//   StraightRotationXSF
//
// DESCRIPTION
//   Performs X-S coordinate final rotation of element
//
// PARAMETERS
//   anglexsf
//   length
//
///////////////////////////////////////////////////////////////////////////

void StraightRotationXSF(Bunch* bunch, double anglexsf, double length)
{
  double xtemp, xptemp, ytemp, yptemp, ztemp, zptemp, z, zp, s;

  double cs = cos(anglexsf);
  double sn = sin(anglexsf);
  double lengtho2 = length / 2.0;

  SyncPart* syncPart = bunch->getSyncPart();

  double dp_p_coeff = 1.0 / (syncPart->getMomentum() * syncPart->getBeta());

  //coordinate array [part. index][x,xp,y,yp,z,dE]
  double** arr = bunch->coordArr();

  for(int i = 0; i < bunch->getSize(); i++)
  {
    xtemp  = arr[i][0];
    xptemp = arr[i][1];
    ytemp  = arr[i][2];
    yptemp = arr[i][3];
    ztemp  = arr[i][4];
    zptemp = 1.0 + arr[i][5] * dp_p_coeff;

    arr[i][0] =  cs * xtemp  + sn * lengtho2;
    arr[i][1] =  cs * xptemp + sn * zptemp;
    arr[i][4] = -sn * xtemp  + cs * ztemp;
    zp        = -sn * xptemp + cs * zptemp;
    arr[i][5] = (zp - 1.0) / dp_p_coeff;

    s = (1.0 - cs) * lengtho2 + sn * xtemp;

    drifti(bunch, i, s);
  }
}

///////////////////////////////////////////////////////////////////////////
//
// NAME
//   StraightRotationYSI
//
// DESCRIPTION
//   Performs Y-S coordinate initial rotation of element
//
// PARAMETERS
//   angleysi
//   length
//
///////////////////////////////////////////////////////////////////////////

void StraightRotationYSI(Bunch* bunch, double angleysi, double length)
{
  double xtemp, xptemp, ytemp, yptemp, ztemp, zptemp, z, zp, s;

  double cs = cos(angleysi);
  double sn = sin(angleysi);
  double lengtho2 = length / 2.0;

  SyncPart* syncPart = bunch->getSyncPart();

  double dp_p_coeff = 1.0 / (syncPart->getMomentum() * syncPart->getBeta());

  //coordinate array [part. index][x,xp,y,yp,z,dE]
  double** arr = bunch->coordArr();

  for(int i = 0; i < bunch->getSize(); i++)
  {
    xtemp  = arr[i][0];
    xptemp = arr[i][1];
    ytemp  = arr[i][2];
    yptemp = arr[i][3];
    ztemp  = arr[i][4];
    zptemp = 1.0 + arr[i][5] * dp_p_coeff;

    arr[i][2] =  cs * ytemp  - sn * lengtho2;
    arr[i][3] =  cs * yptemp + sn * zptemp;
    arr[i][4] = -sn * ytemp  + cs * ztemp;
    zp        = -sn * yptemp + cs * zptemp;
    arr[i][5] = (zp - 1.0) / dp_p_coeff;

    s = -(1.0 - cs) * lengtho2 + sn * ytemp;

    drifti(bunch, i, s);
  }
}

///////////////////////////////////////////////////////////////////////////
//
// NAME
//   StraightRotationYSF
//
// DESCRIPTION
//   Performs X-S coordinate final rotation of element
//
// PARAMETERS
//   angleysf
//   length
//
///////////////////////////////////////////////////////////////////////////

void StraightRotationYSF(Bunch* bunch, double angleysf, double length)
{
  double xtemp, xptemp, ytemp, yptemp, ztemp, zptemp, z, zp, s;

  double cs = cos(angleysf);
  double sn = sin(angleysf);
  double lengtho2 = length / 2.0;

  SyncPart* syncPart = bunch->getSyncPart();

  double dp_p_coeff = 1.0 / (syncPart->getMomentum() * syncPart->getBeta());

  //coordinate array [part. index][x,xp,y,yp,z,dE]
  double** arr = bunch->coordArr();

  for(int i = 0; i < bunch->getSize(); i++)
  {
    xtemp  = arr[i][0];
    xptemp = arr[i][1];
    ytemp  = arr[i][2];
    yptemp = arr[i][3];
    ztemp  = arr[i][4];
    zptemp = 1.0 + arr[i][5] * dp_p_coeff;

    arr[i][2] = cs * ytemp  - sn * lengtho2;
    arr[i][3] = cs * yptemp - sn * zptemp;
    arr[i][4] = sn * ytemp  + cs * ztemp;
    zp        = sn * yptemp + cs * zptemp;
    arr[i][5] = (zp - 1.0) / dp_p_coeff;

    s = (1.0 - cs) * lengtho2 - sn * ytemp;

    drifti(bunch, i, s);
  }
}

///////////////////////////////////////////////////////////////////////////
//
// NAME
//   BendFieldI
//
// DESCRIPTION
//   Performs bend initial field error
//
// PARAMETERS
//   drho
//
///////////////////////////////////////////////////////////////////////////

void BendFieldI(Bunch* bunch, double drho)
{
  //coordinate array [part. index][x,xp,y,yp,z,dE]
  double** arr = bunch->coordArr();

  for(int i = 0; i < bunch->getSize(); i++)
  {
    arr[i][0] -= drho;
  }
}

///////////////////////////////////////////////////////////////////////////
//
// NAME
//   BendFieldF
//
// DESCRIPTION
//   Performs bend final field error
//
// PARAMETERS
//   drho
//
///////////////////////////////////////////////////////////////////////////

void BendFieldF(Bunch* bunch, double drho)
{
  //coordinate array [part. index][x,xp,y,yp,z,dE]
  double** arr = bunch->coordArr();

  for(int i = 0; i < bunch->getSize(); i++)
  {
    arr[i][0] += drho;
  }
}

///////////////////////////////////////////////////////////////////////////
//
// NAME
//   BendDisplacementXI
//
// DESCRIPTION
//   Performs X initial dipole displacement
//
// PARAMETERS
//   anglexi
//   disp
//
///////////////////////////////////////////////////////////////////////////

void BendDisplacementXI(Bunch* bunch, double anglexi, double disp)
{
  double dx = -disp * cos(anglexi / 2.0);
  double ds =  disp * sin(anglexi / 2.0);

  //coordinate array [part. index][x,xp,y,yp,z,dE]
  double** arr = bunch->coordArr();

  for(int i = 0; i < bunch->getSize(); i++)
  {
    arr[i][0] += dx;
    drifti(bunch, i, ds);
  }
}

///////////////////////////////////////////////////////////////////////////
//
// NAME
//   BendDisplacementXF
//
// DESCRIPTION
//   Performs X final dipole displacement
//
// PARAMETERS
//   anglexf
//   disp
//
///////////////////////////////////////////////////////////////////////////

void BendDisplacementXF(Bunch* bunch, double anglexf, double disp)
{
  double dx = disp * cos(anglexf / 2.0);
  double ds = disp * sin(anglexf / 2.0);

  //coordinate array [part. index][x,xp,y,yp,z,dE]
  double** arr = bunch->coordArr();

  for(int i = 0; i < bunch->getSize(); i++)
  {
    arr[i][0] += dx;
    drifti(bunch, i, ds);
  }
}

///////////////////////////////////////////////////////////////////////////
//
// NAME
//   BendDisplacementYI
//
// DESCRIPTION
//   Performs Y initial dipole displacement
//
// PARAMETERS
//   disp
//
///////////////////////////////////////////////////////////////////////////

void BendDisplacementYI(Bunch* bunch, double disp)
{
  //coordinate array [part. index][x,xp,y,yp,z,dE]
  double** arr = bunch->coordArr();

  for(int i = 0; i < bunch->getSize(); i++)
  {
    arr[i][2] -= disp;
  }
}

///////////////////////////////////////////////////////////////////////////
//
// NAME
//   BendDisplacementYF
//
// DESCRIPTION
//   Performs Y final dipole displacement
//
// PARAMETERS
//   disp
//
///////////////////////////////////////////////////////////////////////////

void BendDisplacementYF(Bunch* bunch, double disp)
{
  //coordinate array [part. index][x,xp,y,yp,z,dE]
  double** arr = bunch->coordArr();

  for(int i = 0; i < bunch->getSize(); i++)
  {
    arr[i][2] += disp;
  }
}

///////////////////////////////////////////////////////////////////////////
//
// NAME
//   BendDisplacementLI
//
// DESCRIPTION
//   Performs longitudinal initial dipole displacement
//
// PARAMETERS
//   angleli
//   disp
//
///////////////////////////////////////////////////////////////////////////

void BendDisplacementLI(Bunch* bunch, double angleli, double disp)
{
  double dx = disp * sin(angleli / 2.0);
  double ds = disp * cos(angleli / 2.0);

  //coordinate array [part. index][x,xp,y,yp,z,dE]
  double** arr = bunch->coordArr();

  for(int i = 0; i < bunch->getSize(); i++)
  {
    arr[i][0] += dx;
    drifti(bunch, i, ds);
  }
}

///////////////////////////////////////////////////////////////////////////
//
// NAME
//   BendDisplacementLF
//
// DESCRIPTION
//   Performs longitudinal final dipole displacement
//
// PARAMETERS
//   anglelf
//   disp
//
///////////////////////////////////////////////////////////////////////////

void BendDisplacementLF(Bunch* bunch, double anglelf, double disp)
{
  double dx =  disp * sin(anglelf / 2.0);
  double ds = -disp * cos(anglelf / 2.0);

  //coordinate array [part. index][x,xp,y,yp,z,dE]
  double** arr = bunch->coordArr();

  for(int i = 0; i < bunch->getSize(); i++)
  {
    arr[i][0] += dx;
    drifti(bunch, i, ds);
  }
}

///////////////////////////////////////////////////////////////////////////
//
// NAME
//   RotationI
//
// DESCRIPTION
//   Performs initial rotation
//
// PARAMETERS
//   anglei
//   rhoi
//   theta
//   length
//   et
//   type
//
///////////////////////////////////////////////////////////////////////////

void RotationI(Bunch* bunch, double anglei, double rhoi, double theta,
               double length, std::string et, std::string type)
{
  double lengtho2, cb, sb, ce, se, s;
  double ROT11, ROT12, ROT13,
         ROT21, ROT22, ROT23,
         ROT31, ROT32, ROT33;
  double v1, v2, v3, vr1, vr2, vr3;

  if((et == "SBEN") || (et == "sbend") || (et == "rbend"))
  {
    lengtho2 = tan(theta / 2.0) / rhoi;
    cb = cos(theta / 2.0);
    sb = sin(theta / 2.0);
  }
  else
  {
    lengtho2 = length / 2.0;
    cb = 1.0;
    sb = 0.0;
  }

  ce = cos(anglei);
  se = sin(anglei);

  if(type == "XS")
  {
    ROT11 =  ce;
    ROT12 =  0.0;
    ROT13 = -se;
    ROT21 =  0.0;
    ROT22 =  1.0;
    ROT23 =  0.0;
    ROT31 =  se;
    ROT32 =  0.0;
    ROT33 =  ce;
  }

  if(type == "YS")
  {
    ROT11 = (cb * cb) + (ce * sb * sb);
    ROT12 =  sb * se;
    ROT13 = (cb * sb) - (cb * ce * sb);
    ROT21 = -sb * se;
    ROT22 =  ce;
    ROT23 =  cb * se;
    ROT31 = (cb * sb) - (cb * ce * sb);
    ROT32 = -cb * se;
    ROT33 = (sb * sb) + (cb * cb * ce);
  }

  if(type == "XY")
  {
    ROT11 =  (sb * sb) + (cb * cb * ce);
    ROT12 =   cb * se;
    ROT13 = -(cb * sb) + (cb * ce * sb);
    ROT21 =  -cb * se;
    ROT22 =   ce;
    ROT23 =  -sb * se;
    ROT31 = -(cb * sb) + (cb * ce * sb);
    ROT32 = sb * se;
    ROT33 =  (cb * cb) + (ce * sb * sb);
  }

  SyncPart* syncPart = bunch->getSyncPart();

  double dp_p_coeff = 1.0 / (syncPart->getMomentum() * syncPart->getBeta());

  //coordinate array [part. index][x,xp,y,yp,z,dE]
  double** arr = bunch->coordArr();

  for(int i = 0; i < bunch->getSize(); i++)
  {
    v1 = arr[i][1];
    v2 = arr[i][3];
    v3 = 1.0 + arr[i][5] * dp_p_coeff;

    vr1 = ROT11 * v1 + ROT12 * v2 + ROT13 * v3;
    vr2 = ROT21 * v1 + ROT22 * v2 + ROT23 * v3;
    vr3 = ROT31 * v1 + ROT32 * v2 + ROT33 * v3;

    arr[i][1] = vr1;
    arr[i][3] = vr2;
    arr[i][5] = (vr3 - 1.0) / dp_p_coeff;

    v1 = arr[i][0];
    v2 = arr[i][2];
    v3 = -lengtho2;

    vr1 = ROT11 * v1 + ROT12 * v2 + ROT13 * v3;
    vr2 = ROT21 * v1 + ROT22 * v2 + ROT23 * v3;
    vr3 = ROT31 * v1 + ROT32 * v2 + ROT33 * v3;

    arr[i][0] = vr1;
    arr[i][2] = vr2;
    s = -lengtho2 - vr3;
    drifti(bunch, i, s);
  }
}

///////////////////////////////////////////////////////////////////////////
//
// NAME
//   RotationF
//
// DESCRIPTION
//   Performs final rotation
//
// PARAMETERS
//   anglef
//   rhoi
//   theta
//   length
//   et
//   type
//
///////////////////////////////////////////////////////////////////////////

void RotationF(Bunch* bunch, double anglef, double rhoi, double theta,
               double length, std::string et, std::string type)
{
  double lengtho2, cb, sb, ce, se, s;
  double ROT11, ROT12, ROT13,
         ROT21, ROT22, ROT23,
         ROT31, ROT32, ROT33;
  double v1, v2, v3, vr1, vr2, vr3;

  if((et == "SBEN") || (et == "sbend") || (et == "rbend"))
  {
    lengtho2 = tan(theta / 2.0) / rhoi;
    cb = cos(theta / 2.0);
    sb = sin(theta / 2.0);
  }
  else
  {
    lengtho2 = length / 2.0;
    cb = 1.0;
    sb = 0.0;
  }

  ce = cos(anglef);
  se = sin(anglef);

  if(type == "XS")
  {
    ROT11 =  ce;
    ROT12 =  0.0;
    ROT13 =  se;
    ROT21 =  0.0;
    ROT22 =  1.0;
    ROT23 =  0.0;
    ROT31 = -se;
    ROT32 =  0.0;
    ROT33 =  ce;
  }

  if(type == "YS")
  {
    ROT11 =  (cb * cb) + (ce * sb * sb);
    ROT12 =   sb * se;
    ROT13 = -(cb * sb) + (cb * ce * sb);
    ROT21 =  -sb * se;
    ROT22 =   ce;
    ROT23 =  -cb * se;
    ROT31 = -(cb * sb) + (cb * ce * sb);
    ROT32 =   cb * se;
    ROT33 =  (sb * sb) + (cb * cb * ce);
  }

  if(type == "XY")
  {
    ROT11 = (sb * sb) + (cb * cb * ce);
    ROT12 = -cb * se;
    ROT13 = (cb * sb) - (cb * ce * sb);
    ROT21 =  cb * se;
    ROT22 =  ce;
    ROT23 = -sb * se;
    ROT31 = (cb * sb) - (cb * ce * sb);
    ROT32 =  sb * se;
    ROT33 = (cb * cb) + (ce * sb * sb);
  }

  SyncPart* syncPart = bunch->getSyncPart();

  double dp_p_coeff = 1.0 / (syncPart->getMomentum() * syncPart->getBeta());

  //coordinate array [part. index][x,xp,y,yp,z,dE]
  double** arr = bunch->coordArr();

  for(int i = 0; i < bunch->getSize(); i++)
  {
    v1 = arr[i][1];
    v2 = arr[i][3];
    v3 = 1.0 + arr[i][5] * dp_p_coeff;

    vr1 = ROT11 * v1 + ROT12 * v2 + ROT13 * v3;
    vr2 = ROT21 * v1 + ROT22 * v2 + ROT23 * v3;
    vr3 = ROT31 * v1 + ROT32 * v2 + ROT33 * v3;

    arr[i][1] = vr1;
    arr[i][3] = vr2;
    arr[i][5] = (vr3 - 1.0) / dp_p_coeff;

    v1 = arr[i][0];
    v2 = arr[i][2];
    v3 = lengtho2;

    vr1 = ROT11 * v1 + ROT12 * v2 + ROT13 * v3;
    vr2 = ROT21 * v1 + ROT22 * v2 + ROT23 * v3;
    vr3 = ROT31 * v1 + ROT32 * v2 + ROT33 * v3;

    arr[i][0] = vr1;
    arr[i][2] = vr2;
    s = lengtho2 - vr3;
    drifti(bunch, i, s);
  }
}

///////////////////////////////////////////////////////////////////////////
//
// NAME
//   DipoleKickerOsc
//
// DESCRIPTION
//   Performs oscillating dipole kick
//
// PARAMETERS
//   k, phaselength, phase
//
///////////////////////////////////////////////////////////////////////////

void DipoleKickerOsc(Bunch* bunch, double k,
                     double phaselength, double phase)
{
  double ztophi = 2.0 * OrbitConst::PI / phaselength;
  double kick;

  //coordinate array [part. index][x,xp,y,yp,z,dE]
  double** arr = bunch->coordArr();

  for(int i = 0; i < bunch->getSize(); i++)
  {
    kick = k * sin(ztophi * arr[i][4] + phase);
    arr[i][1] += kick;
    arr[i][3] -= kick;
  }
}

///////////////////////////////////////////////////////////////////////////
//
// NAME
//   QuadKicker
//
// DESCRIPTION
//   Performs quadrupole kick
//
// PARAMETERS
//   k
//
///////////////////////////////////////////////////////////////////////////

void QuadKicker(Bunch* bunch, double k)
{
  //coordinate array [part. index][x,xp,y,yp,z,dE]
  double** arr = bunch->coordArr();

  for(int i = 0; i < bunch->getSize(); i++)
  {
    arr[i][1] += k * arr[i][0];
    arr[i][3] -= k * arr[i][2];
  }
}

///////////////////////////////////////////////////////////////////////////
//
// NAME
//   QuadKickerOsc
//
// DESCRIPTION
//   Performs oscillating quadrupole kick
//
// PARAMETERS
//   k, phaselength, phase
//
///////////////////////////////////////////////////////////////////////////

void QuadKickerOsc(Bunch* bunch, double k,
                   double phaselength, double phase)
{
  double ztophi = 2.0 * OrbitConst::PI / phaselength;
  double kick;

  //coordinate array [part. index][x,xp,y,yp,z,dE]
  double** arr = bunch->coordArr();

  for(int i = 0; i < bunch->getSize(); i++)
  {
    kick = k * sin(ztophi * arr[i][4] + phase);
    arr[i][1] += kick * arr[i][0];
    arr[i][3] -= kick * arr[i][2];
  }
}

///////////////////////////////////////////////////////////////////////////
// NAME
//   drifti
//
// DESCRIPTION
//   Drifts a single particle. Length < 0 is allowed.
//
// PARAMETERS
//   bunch = reference to the macro-particle bunch
//   i = particle index
//   length = length of the drift
//
///////////////////////////////////////////////////////////////////////////

void drifti(Bunch* bunch, int i, double length)
{
  double KNL, phifac, dp_p;

  SyncPart* syncPart = bunch->getSyncPart();

  double gamma2i = 1.0 / (syncPart->getGamma() * syncPart->getGamma());
  double dp_p_coeff = 1.0 / (syncPart->getMomentum() * syncPart->getBeta());

  //coordinate array [part. index][x,xp,y,yp,z,dE]
  double** arr = bunch->coordArr();

  dp_p = arr[i][5] * dp_p_coeff;
  KNL  = 1.0 / (1.0 + dp_p);
  arr[i][0] += KNL * length * arr[i][1];
  arr[i][2] += KNL * length * arr[i][3];
  phifac = (arr[i][1] * arr[i][1] + arr[i][3] * arr[i][3] +
            dp_p * dp_p * gamma2i) / 2.0;
  phifac = (phifac * KNL - dp_p * gamma2i) * KNL;
  arr[i][4] -= length * phifac;
}

}  //end of namespace error_base
