//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    OrbitConst.cc
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


#include "OrbitConst.hh"

const double OrbitConst::PI = 3.14159265358979323846264;
const double OrbitConst::c = 2.99792458e+8;
const double OrbitConst::elementary_charge_CGS = 4.8032068e-10;
const double OrbitConst::coeff_Tesla_to_inner = (1.0e+4/4.8032068e-10)*0.01*(1.0e+6);
const double OrbitConst::coeff_VoltsPerM_to_inner = 1.0/((1.6021773e-19)*(8.98755179e+9));
const double OrbitConst::coeff_Phi_to_Volts = (1.6021773e-19)*(8.98755179e+9);
const double OrbitConst::mass_electron = 0.51099906e-3;
const double OrbitConst::classicalRadius_electron = 2.8179409e-15;
const double OrbitConst::charge_electron = -1.0;
const double OrbitConst::mass_proton = 0.9382723;
const double OrbitConst::classicalRadius_proton = 1.534698e-18;
const double OrbitConst::charge_proton = +1.0;
const double OrbitConst::permeability = 4*3.14159265358979323846264*1e-007;

OrbitConst::OrbitConst()
{
}

OrbitConst::~OrbitConst()
{
}

///////////////////////////////////////////////////////////////////////////
//
// END OF FILE
//
///////////////////////////////////////////////////////////////////////////
