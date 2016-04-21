//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//   TwissFilteringFunctions.cc
//
// AUTHOR
//    A. Shishlo
//
// CREATED
//    04/20/2016
//
// DESCRIPTION
//    A set of functions for filtering macro-particles from the bunch
//    according to their positions relative to the emittance phase space center 
//
///////////////////////////////////////////////////////////////////////////

#include "BunchTwissAnalysis.hh"

#include "orbit_mpi.hh"

namespace OrbitUtils{
	
		/** 
    Bunch filtering according to the Twiss parameters.
    bunch_in is the input bunch, and bunch_bad collects
    the macro-particles that are exluded from the initial bunch.
    The coefficients are from 0 to infinity and define the limit for 
    the ratio
    (gamma*x^2+2*alpha*x*xp+beta*xp^2)/(2*emittance)
    The function returns the total number of removed particles. 
    */
    int bunch_twiss_filtering(Bunch* bunch_in,Bunch* bunch_bad, double coeff_x, double coeff_y, double coeff_z){
      
      BunchTwissAnalysis* twissAnalysis = new BunchTwissAnalysis();
      twissAnalysis->analyzeBunch(bunch_in);
      
      bunch_in->copyBunchTo(bunch_bad);
      
      double coeff_arr[3] = {coeff_x, coeff_y, coeff_z};
      int n_parts =  bunch_in->getSize();
      double emitt_arr[3] = {};
      double gamma_arr[3] = {};
      double alpha_arr[3] = {};
      double beta_arr[3] = {};
      
      for(int ic = 0; ic < 3; ic++){ 
        emitt_arr[ic] =  twissAnalysis->getEffectiveEmittance(ic);
        gamma_arr[ic] =  twissAnalysis->getEffectiveGamma(ic);
        alpha_arr[ic] =  twissAnalysis->getEffectiveAlpha(ic);
        beta_arr[ic]  =  twissAnalysis->getEffectiveBeta(ic);
      }
      
      double coeff = 0;
      double x = 0.;
      double xp = 0.;
      int info = +1;
      int ic2 = 0;
      int ic21 = 0; 
      int wt = 0.;
      for(int ind = 0; ind < n_parts; ind++){
        info = +1;
        for(int ic = 0; ic < 3; ic++){
          if(coeff_arr[ic] <= 0.) continue;
          ic2 = ic*2;
          ic21 = ic*2+1;        
          x = bunch_in->coordArr()[ind][ic2];
          xp = bunch_in->coordArr()[ind][ic21];
          wt = (gamma_arr[ic]*x*x+2*alpha_arr[ic]*x*xp+beta_arr[ic]*xp*xp)/(2*emitt_arr[ic]);
          if(wt >= coeff_arr[ic]) info = -1;
        }
        if(info > 0){
          bunch_bad->deleteParticleFast(ind);
        } else {
          bunch_in->deleteParticleFast(ind);
        }
      }
      bunch_in->compress();
      bunch_bad->compress();
      
      int bad_count = bunch_bad->getSizeGlobal(); 
      return bad_count;
    }
}
