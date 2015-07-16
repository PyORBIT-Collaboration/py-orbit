///////////////////////////////////////////////////////////////////////////
//
// FILE NAME
//   errorbase.hh
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
//   Define functions for magnet alignment errors
//
///////////////////////////////////////////////////////////////////////////
#ifndef ERROR_BASE_H
#define ERROR_BASE_H

#include "Bunch.hh"

namespace error_base
{
  void drifti(Bunch* bunch, int i, double length);
  void CoordDisplacement(Bunch* bunch,
                         double dx, double dxp,
                         double dy, double dyp,
                         double dz, double dE);
  void QuadKicker(Bunch* bunch, double k);
  void QuadKickerOsc(Bunch* bunch, double k,
                     double phaselength, double phase);
  void DipoleKickerOsc(Bunch* bunch, double k,
                       double phaselength, double phase);
  void LongDisplacement(Bunch* bunch, double ds);
  void StraightRotationXY(Bunch* bunch, double anglexy);
  void StraightRotationXSI(Bunch* bunch, double anglexsi, double length);
  void StraightRotationXSF(Bunch* bunch, double anglexsf, double length);
  void StraightRotationYSI(Bunch* bunch, double angleysi, double length);
  void StraightRotationYSF(Bunch* bunch, double angleysf, double length);
  void BendFieldI(Bunch* bunch, double drho);
  void BendFieldF(Bunch* bunch, double drho);
  void BendDisplacementXI(Bunch* bunch, double anglexi, double disp);
  void BendDisplacementXF(Bunch* bunch, double anglexf, double disp);
  void BendDisplacementYI(Bunch* bunch, double disp);
  void BendDisplacementYF(Bunch* bunch, double disp);
  void BendDisplacementLI(Bunch* bunch, double angleli, double disp);
  void BendDisplacementLF(Bunch* bunch, double anglelf, double disp);
  void RotationI(Bunch* bunch, double anglei, double rhoi,
  	             double theta, double length,
                 std::string et, std::string type);
  void RotationF(Bunch* bunch, double anglef, double rhoi,
  	             double theta, double length,
                 std::string et, std::string type);
  double drand(double r);
  double derf(double x);
  double root_normal(double errtest, double ymin,
                     double ymax, double tol);
  double getGauss(double mean, double sigma, double cutoff);
}

#endif  //ERROR_BASE_H
