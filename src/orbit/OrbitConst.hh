//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    OrbitConst.hh
//
// AUTHOR
//    A. Shishlo
//
// CREATED
//    06/28/2005
//
// DESCRIPTION
//    Keeps all constants for ORBIT
//
///////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////
//
// include files
//
/////////////////////////////////////////////////////////////////////////////

#ifndef ORBIT_CONSTANTS_H
#define ORBIT_CONSTANTS_H

class  OrbitConst
{
public:

  OrbitConst();
  ~OrbitConst();

  //-------------------------------------------------------------
  //constants
  //-------------------------------------------------------------

  //PI = 3.141592653589793238462643383279502
  const static double PI;

  //speed of light in [m/sec]
  const static double c;

  //elementary charge (CGS system)
  const static double elementary_charge_CGS;

  //coeffitient for magnetic field strength
  //to shift from Tesla to H(CGS)/e(SGS) in 1/m^2
  const static double coeff_Tesla_to_inner;

  //to shift from Volts/m to E in 1/m^2
  const static double coeff_VoltsPerM_to_inner;

  //coeffitient to get volts from our potentials
  const static double coeff_Phi_to_Volts;

  //------------------------------------------------------------
  //electron parameters
  //------------------------------------------------------------

  //electron mass in GeV
  const static double mass_electron;

  //electron classical radius in m
  const static double classicalRadius_electron;

  //charge of the electron
  const static double charge_electron;

  //------------------------------------------------------------
  //proton parameters
  //------------------------------------------------------------

  //proton mass in GeV
   const static double mass_proton;

  //proton classical radius in m
  const static double classicalRadius_proton;

  //charge of the proton
  const static double charge_proton;

};

///////////////////////////////////////////////////////////////////////////
//
// END OF FILE
//
///////////////////////////////////////////////////////////////////////////

#endif
