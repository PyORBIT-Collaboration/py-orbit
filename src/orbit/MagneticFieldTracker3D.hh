//******************************************************************************
// Program: MagneticFieldTracker3D.h
// Creator: Matt Perkett
// Date: 10/09/2006
// Purpose: header file for MagneticFieldTracker3D class that will track
//    particles for pyORBIT
// Algorithm:
// Limitations/Bugs:
// Modified:
//******************************************************************************

#ifndef _MAGNETID_FIELD_TRACKER_3D_H
#define _MAGNETID_FIELD_TRACKER_3D_H

#include "BaseFieldSource.hh"
#include "Bunch.hh"

class MagneticFieldTracker3D
{
   public:
   MagneticFieldTracker3D();
   MagneticFieldTracker3D(double stepSize);
   ~MagneticFieldTracker3D();
   
   void track(double zmin, double zmax, double x_0, double y_0, double &alpha,
              double &beta, double &gamma, Bunch* bunch, BaseFieldSource* fieldSource,
				  int method_flag);
				  // alpha, beta, and gamma will store the Euler angles to transfer
				  //    back into the tracker coordinate system
   void setStepSize(double stepSize);
   double getStepSize() const;
   
   private:
   static bool isUnitVector(double n1, double n2, double n3);
   
   double ds;  // step size
};

#endif


