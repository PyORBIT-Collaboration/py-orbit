//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    cross_sections.hh
//
// CREATED
//    09/13/2011
//
// DESCRIPTION
//   Energy dependent elastic scattering cross sections for ~GeV protons on 
//   various materials. 
///////////////////////////////////////////////////////////////////////////
#ifndef CROSS_SECTIONS_H
#define CROSS_SECTIONS_H

namespace OrbitUtils{

	double get_z(int ma);
	double get_a(int ma);
	double get_rho(int ma);
	double get_radlength(int ma);
	double get_elastic_crosssection(double energyrequest, int ma_index);
	double get_inelastic_crosssection(double energyrequest, int ma_index);
				
}

#endif
