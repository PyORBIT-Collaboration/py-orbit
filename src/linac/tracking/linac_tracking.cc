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
		if(length <= 0.) return;
		
		double delta_t = length / v_s;
		double delta_tc = delta_t*OrbitConst::c;
		
		syncPart->setTime(syncPart->getTime() + delta_t);
		
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
		double Etotal = 0.;
		double p_z2 = 0.;
		double beta_z = 0.;
		double beta_x = 0.;
		double beta_y = 0.;
		
		//coordinate array [part. index][x,xp,y,yp,z,dE]
		//sqrt(fabs(p_z2)) - fabs function is a protection for case when xp and yp are too big for dE
		//It means the particles are nonphysical and simulations do not make sense 
		double** arr = bunch->coordArr();
		for(int i = 0; i < bunch->getSize(); i++)
		{
			dE = arr[i][5];
			xp = arr[i][1];
			yp = arr[i][3];
			Ekin = Ekin_s+dE;
			Etotal = Ekin + mass;
			p2 = Ekin*(Etotal + mass);
			beta_x = (xp*p_s)/Etotal;
			beta_y = (yp*p_s)/Etotal;
			p_z2 = p2 - (xp*xp + yp*yp)*p2_s;
			beta_z = sqrt(fabs(p_z2))/Etotal;
			arr[i][0] += beta_x * delta_tc;
			arr[i][2] += beta_y * delta_tc;
			arr[i][4] += (beta_z - beta_s) * delta_tc;
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
	//   kq = G/(B*rho) B*rho = 3.33564*P0/charge [T*m] where P0 in GeV/c
	//   G = dB/dr [T/m] quad gradient
	//
	// RETURNS
	//   Nothing
	//
	///////////////////////////////////////////////////////////////////////////
	
	void linac_quad1(Bunch* bunch, double length, double kq, int useCharge)
	{
    if(kq == 0. || bunch->getCharge() == 0.)
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
		
 		double delta_t = length / v_s;
		double delta_tc = delta_t*OrbitConst::c;   
		
    if(length > 0.)
    {
    	syncPart->setTime(syncPart->getTime() + delta_t);
    }
    
 		// ==== B*rho = 3.335640952*momentum/charge [T*m] if momentum in GeV/c ===
 		double dB_dr = kq*bunch->getB_Rho();
 		
 		double mass = syncPart->getMass();
		double Ekin_s = syncPart->getEnergy();
		double p_s = syncPart->getMomentum();
		double p2_s = p_s*p_s;
		double p2 = 0.;
		double dE = 0.;
		double p = 0.;
		double Ekin = 0.;  
		double Etotal = 0.;
		double coeff = 0.; 
		
		double p_z2 = 0.;
		double beta_z = 0.;
   
    //coordinate array [part. index][x,xp,y,yp,z,dE]
    double** arr = bunch->coordArr();
    int nParts = bunch->getSize();
    
    double kqc = 0.;
    
    for(int i = 0; i < nParts; i++)
    {
    	
 			dE = arr[i][5];
			Ekin = Ekin_s+dE;
			Etotal = Ekin + mass;
			p2 = Ekin*(Ekin+2.0*mass); 
			p = sqrt(p2);
    	kqc = dB_dr/(3.335640952*p/bunch->getCharge());
			beta_z = p/Etotal;
			
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
     	
    	// coeff = p_s/p
  		coeff = p_s/p;
  		
    	x_init  = arr[i][0];
    	xp_init = xp_init*coeff;
    	y_init  = arr[i][2];
    	yp_init = yp_init*coeff;
    	
    	arr[i][0]  = x_init * m11 + xp_init * m12;
    	arr[i][1]  = (x_init * m21 + xp_init * m22)/coeff;
    	arr[i][2]  = y_init * m33 + yp_init * m34;
    	arr[i][3]  = (y_init * m43 + yp_init * m44)/coeff;

    	arr[i][4] += (beta_z - beta_s) * delta_tc;
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
  

	////////////////////////////
	// NAME
	//   linac_quad3
	//
	// DESCRIPTION
	//   Quadrupole element 3: longitudinal field component of the quad field
	//   Bz(x,y,z) = x*y*dG(z)/dz
	//   This function performs the transverse kicks correction, so the length
	//   is just a parameter.
	//
	// PARAMETERS
	//   bunch  = reference to the macro-particle bunch
	//   length = length of transport
	//   G = dB/dr [T/m] quad gradient
	//   0.299792458 = (1/v_light) * 10^9
	//
	// RETURNS
	//   Nothing
	//
	///////////////////////////////////////////////////////////////////////////
	
	void linac_quad3(Bunch* bunch, double length, double dB_dz)
	{
    double charge = bunch->getCharge();

    if(dB_dz == 0. || charge == 0.)
    {
    	return;
    }
    
    SyncPart* syncPart = bunch->getSyncPart();
    //momentum in GeV/c
    double momentum = syncPart->getMomentum();
    double kick_coeff = 0.299792458*charge*dB_dz*length/momentum;
     
    //coordinate array [part. index][x,xp,y,yp,z,dE]
    double** arr = bunch->coordArr();
    int nParts = bunch->getSize();
    
    double coef_xy = 0.;
    
    for(int i = 0; i < nParts; i++){
    	coef_xy = kick_coeff*arr[i][0]*arr[i][2];
    	arr[i][1]  +=  arr[i][3]*coef_xy;
    	arr[i][3]  += -arr[i][1]*coef_xy;
    }
  }
   
///////////////////////////////////////////////////////////////////////////
// NAME
//   kick
//
// DESCRIPTION
//   Kicks a particle bunch
//
// PARAMETERS
//   bunch = reference to the macro-particle bunch
//   kx = strength of the horizontal kick in rad
//   ky = strength of the vertical kick in rad
//   kE = strength of the energy kick in GeV
//
// RETURNS
//   Nothing
//
///////////////////////////////////////////////////////////////////////////

	void kick(Bunch* bunch, double kx, double ky, double kE, int useCharge)
	{
		SyncPart* syncPart = bunch->getSyncPart();
		
		double mass = syncPart->getMass();
		double Ekin_s = syncPart->getEnergy();  	
		double p_s = syncPart->getMomentum();
		double p2_s = p_s*p_s;
		double dE = 0.;
		double p2 = 0.;
		double p_z = 0.;
		double Ekin = 0.;    
		double coeff = 0.;
		
		//coordinate array [part. index][x,xp,y,yp,z,dE]
		double** arr = bunch->coordArr();
		if(kx != 0.)
		{
			for(int i = 0; i < bunch->getSize(); i++)
			{
				dE = arr[i][5];
				Ekin = Ekin_s+dE;
				p2 = Ekin*(Ekin+2.0*mass);
				p_z = sqrt(p2 - (arr[i][1]*arr[i][1] + arr[i][3]*arr[i][3])*p2_s);
				coeff = p_s/p_z;
				arr[i][1] += kx*coeff;
			}
		}
		if(ky != 0.)
		{
			for(int i = 0; i < bunch->getSize(); i++)
			{
				dE = arr[i][5];
				Ekin = Ekin_s+dE;
				p2 = Ekin*(Ekin+2.0*mass);
				p_z = sqrt(p2 - (arr[i][1]*arr[i][1] + arr[i][3]*arr[i][3])*p2_s);
				coeff = p_s/p_z;        	
				arr[i][3] += ky*coeff;
			}
		}
		if(kE != 0.)
		{
			for(int i = 0; i < bunch->getSize(); i++)
			{
				arr[i][5] += kE;
			}
		}
	}
	
}  //end of namespace linac_tracking
