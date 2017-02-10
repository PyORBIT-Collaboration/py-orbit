/////////////////////////////////////////////////////////////////////////////
//
// FILE NAME
//   linac_tracking.cc
//
// AUTHOR
//   Andrei Shishlo, ORNL, shishlo@ornl.gov
//   2017.01.24
//
// DESCRIPTION
//   Define elementary functions for tracking through
//   specific linac accelerator elements without using
//   TEAPOT algorithms. It is necessary for the cases 
//   with a huge energy spread in the bunch.
//
/////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
//
// Include files
//
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
//
// Local Functions:
//
///////////////////////////////////////////////////////////////////////////

#include "OrbitConst.hh"
#include "Bunch.hh"
#include "SyncPart.hh"

#include <complex>

namespace linac_tracking
{
	
	///////////////////////////////////////////////////////////////////////////
	// NAME
	//   linac_drift
	//
	// DESCRIPTION
	//   Drifts a particle bunch. Length < 0 is not allowed.
	//
	// PARAMETERS
	//   bunch = reference to the macro-particle bunch
	//   length = length of the drift
	//
	// RETURNS
	//   Nothing
	//
	///////////////////////////////////////////////////////////////////////////
	
	void linac_drift(Bunch* bunch, double length)
	{
		SyncPart* syncPart = bunch->getSyncPart();
		
		double beta_s = syncPart->getBeta();
		double v_s = OrbitConst::c * beta_s;
		double beta1_s = 1.0/beta_s;
		if(length < 0.) return;
		
		syncPart->setTime(syncPart->getTime() + length / v_s);
		
		double mass = syncPart->getMass();
		double Ekin_s = syncPart->getEnergy();
		double p_s = syncPart->getMomentum();
		double p2_s = p_s*p_s;
		double p2 = 0.;
		double dE = 0.;
		double xp = 0.;
		double yp = 0.;
		double coeff = 0.;
		double Ekin = 0.;
		//double p_z2 = 0.;
		//double beta_z = 0.;
		
		//coordinate array [part. index][x,xp,y,yp,z,dE]
		double** arr = bunch->coordArr();
		for(int i = 0; i < bunch->getSize(); i++)
		{
			dE = arr[i][5];
			xp = arr[i][1];
			yp = arr[i][3];
			Ekin = Ekin_s+dE;
			p2 = Ekin*(Ekin+2.0*mass);
			//p_z2 = p2 - (xp*xp + yp*yp)*p2_s;
			//beta_z = sqrt(p_z2)/(Ekin+mass);
			coeff = 1.0/sqrt(p2/p2_s - xp*xp - yp*yp);
			arr[i][0] += length * xp * coeff;
			arr[i][2] += length * yp * coeff;
			arr[i][4] += length * beta_s * (beta1_s - (beta1_s + dE/p_s)*coeff);
		}
	}
	
	////////////////////////////
	// NAME
	//   linac_quad1
	//
	// DESCRIPTION
	//   Quadrupole element one: linear transport matrix for linac
	//
	// PARAMETERS
	//   bunch  = reference to the macro-particle bunch
	//   length = length of transport
	//   kq = quadrupole field strength [m^(-2)]
	//
	// RETURNS
	//   Nothing
	//
	///////////////////////////////////////////////////////////////////////////
	
	void linac_quad1(Bunch* bunch, double length, double kq, int useCharge)
	{
    double charge = +1.0;
    if(useCharge == 1) charge = bunch->getCharge();
    
    double kqc = kq * charge;
    if(kqc == 0.)
    {
    	linac_drift(bunch,length);
    	return;
    }
    
    double x_init, xp_init, y_init, yp_init;
    double sqrt_kq, kqlength;
    double cx, sx, cy, sy;
    double m11 = 0., m12 = 0., m21 = 0., m22 = 0.;
    double m33 = 0., m34 = 0., m43 = 0., m44 = 0.;
    
    SyncPart* syncPart = bunch->getSyncPart();
    double beta_s = syncPart->getBeta();
		double beta1_s = 1.0/beta_s;    
		double v_s = OrbitConst::c * beta_s;    
    
    if(length > 0.)
    {
    	syncPart->setTime(syncPart->getTime() + length / v_s);
    }
    
 		// ==== B*rho = 3.335640952*momentum [T*m] if momentum in GeV/c ===
 		double dB_dr = kqc*3.335640952*syncPart->getMomentum();
 		
 		double mass = syncPart->getMass();
		double Ekin_s = syncPart->getEnergy();
		double p_s = syncPart->getMomentum();
		double p2_s = p_s*p_s;
		double p2 = 0.;
		double dE = 0.;
		double p = 0.;
		double Ekin = 0.;    
		double coeff = 0.; 
   
    //coordinate array [part. index][x,xp,y,yp,z,dE]
    double** arr = bunch->coordArr();
    int nParts = bunch->getSize();
    
    for(int i = 0; i < nParts; i++)
    {
    	
 			dE = arr[i][5];
			Ekin = Ekin_s+dE;
			p2 = Ekin*(Ekin+2.0*mass); 
			p = sqrt(p2);
    	kqc = dB_dr/(3.335640952*p);
			
    	sqrt_kq  = sqrt(fabs(kqc));
    	kqlength = sqrt_kq * length; 
    	
    	if(kqc > 0.)
    	{
    		cx = cos(kqlength);
    		sx = sin(kqlength);
    		cy = cosh(kqlength);
    		sy = sinh(kqlength);
    		m21 = -sx * sqrt_kq;
    		m43 = sy * sqrt_kq;
    	}
    	else
    	{
    		cx = cosh(kqlength);
    		sx = sinh(kqlength);
    		cy = cos(kqlength);
    		sy = sin(kqlength);
    		m21 = sx * sqrt_kq;
    		m43 = -sy * sqrt_kq;
    	}
    	
    	m11 = cx;
    	m12 = sx / sqrt_kq;
    	m22 = cx;
    	m33 = cy;
    	m34 = sy / sqrt_kq;
    	m44 = cy;    
    	
     	xp_init = arr[i][1];
    	yp_init = arr[i][3];
     	
    	// coeff = p_s/p_z
  		coeff = 1.0/sqrt(p2/p2_s - xp_init*xp_init - yp_init*yp_init);
  		
    	x_init  = arr[i][0];
    	xp_init = xp_init*coeff;
    	y_init  = arr[i][2];
    	yp_init = yp_init*coeff;
    	
    	arr[i][0]  = x_init * m11 + xp_init * m12;
    	arr[i][1]  = (x_init * m21 + xp_init * m22)/coeff;
    	arr[i][2]  = y_init * m33 + yp_init * m34;
    	arr[i][3]  = (y_init * m43 + yp_init * m44)/coeff;
 			arr[i][4] += length * beta_s * (beta1_s - (beta1_s + dE/p_s)*coeff);
    }
  }
  
  ///////////////////////////////////////////////////////////////////////////
  // NAME
  //  linac_ quad2
  //
  // DESCRIPTION
  //   Quadrupole element two: nonlinear piece for linac
  //
  // PARAMETERS
  //   bunch  = reference to the macro-particle bunch
  //   length = length of the element
  //
  // RETURNS
  //   Nothing
  //
  ///////////////////////////////////////////////////////////////////////////
  
  void linac_quad2(Bunch* bunch, double length)
  {
  	//there are no non-linear components
  }
}  //end of namespace linac_tracking
