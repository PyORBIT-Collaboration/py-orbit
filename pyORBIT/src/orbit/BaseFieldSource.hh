//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    BaseFieldSource.hh
//
// AUTHOR
//    Y. Sato, A. Shishlo
//
// CREATED
//    03/31/2003
//
// DESCRIPTION
//    The base class for field source. Can be sub-classed
//
///////////////////////////////////////////////////////////////////////////

#ifndef BASE_FIELD_SOURCE_H
#define BASE_FIELD_SOURCE_H

class  BaseFieldSource
{
public:

  BaseFieldSource();
  virtual ~BaseFieldSource();

  virtual void getElectricField(double t, double x, double y, double z, double& f_x, double& f_y, double& f_z);
  virtual void getMagneticField(double t, double x, double y, double z, double& f_x, double& f_y, double& f_z);

};

///////////////////////////////////////////////////////////////////////////
//
// END OF FILE
//
///////////////////////////////////////////////////////////////////////////

#endif
