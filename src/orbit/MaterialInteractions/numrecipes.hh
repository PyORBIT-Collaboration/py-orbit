//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    numrecipes.hh
//
// CREATED
//    10/11/2011
//
// DESCRIPTION
//   Some numerical root finding recipes
///////////////////////////////////////////////////////////////////////////
#ifndef NUMRECIPES_H
#define NUMRECIPES_H

namespace OrbitUtils{

	float fstep(float s, float r_o, float pr_o, float theta);
	float rfunc(float x, float p, float fac1);
	int zbrak(float (*fx)(float, float, float, float), float x1, float x2, int n, float xb1[], float xb2[], int &nb, float param1, 
			  float param2, float param3);
	float rtbis(float (*func)(float, float, float, float), float x1, float x2, float xacc, float param1, float param2, float param3);	
	float bessj0(float x);
	float bessj1(float x);
	float qsimp(float (*func)(float, float, float), float a, float b, float p, float fac1);
	float trapzd(float (*rfunc)(float, float, float), float a, float b, int n, float p, float fac1);
}

#endif
