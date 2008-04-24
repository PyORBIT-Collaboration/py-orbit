/////////////////////////////////////////////////////////////////////////////
//
// FILE NAME
//   teapotbase.cc
//
// AUTHOR
//   Joshua Abrams, Knox College, jabrams@knox.edu
//   Steven Bunch, University of Tennessee, sbunch2@utk.edu
//
// modified by Andrei Shishlo
//   12/30/05
//
// DESCRIPTION
//   define elementary function for different elements
//
/////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
//
// include files
//
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
//
// Local Functions:
//
///////////////////////////////////////////////////////////////////////////
#include "teapotbase.hh"

#include "OrbitConst.hh"
#include "Bunch.hh"
#include "SyncPart.hh"

#include <complex>

namespace teapot_base{

	static double* factorial = NULL;

	void init_factorial(){
		if(factorial == NULL){
			int n = 50;
			factorial = new double[n];
			factorial[0] = 1.0;
			for(int i = 1; i < n; i++){
				factorial[i] = i*factorial[i-1];
			}
		}
	}

	void delete_factorial(){
		delete [] factorial;
	}

///////////////////////////////////////////////////////////////////////////
// NAME
//
//  rotatexy
//
// DESCRIPTION
//  rotates particles coordinates
//
// PARAMETERS
//  bunch = reference to the macro-particle bunch
//  anglexy = rotation angle
//
// RETURNS
//    Nothing
//
///////////////////////////////////////////////////////////////////////////

void rotatexy(Bunch* bunch, double anglexy){
	double xtemp, pxtemp, ytemp, pytemp;
	double cs = cos(anglexy);
	double sn = sin(anglexy);

	//coordinate array [part. index][x,xp,y,yp,z,dE]
	double** arr = bunch->coordArr();

	for(int i = 0, n = bunch->getSize(); i < n; i++)
	{
		xtemp  = arr[i][0];
		pxtemp = arr[i][1];
		ytemp  = arr[i][2];
		pytemp = arr[i][3];

		arr[i][0] =  cs * xtemp  - sn * ytemp;
		arr[i][1] =  cs * pxtemp - sn * pytemp;
		arr[i][2] =  sn * xtemp  + cs * ytemp;
		arr[i][3] =  sn * pxtemp + cs * pytemp;
	}
}

///////////////////////////////////////////////////////////////////////////
// NAME
//  drifti
//
// DESCRIPTION
//  drifts a single particle
//
// PARAMETERS
//  bunch = reference to the macro-particle bunch
//  i = particle index
//  length = length of the drift
//
// RETURNS
//    Nothing
//
///////////////////////////////////////////////////////////////////////////

void drifti(Bunch* bunch, int i, double length){

	if(length <= 0.) return;

	//coordinate array [part. index][x,xp,y,yp,z,dE]
	double** arr = bunch->coordArr();

	double KNL, px, py, phifac;

	SyncPart* syncPart = bunch->getSyncPart();

	double gamma = syncPart->getGamma();
	double gamma2i = 1.0 / (gamma*gamma);
	double dp_p = (arr[i][5] / syncPart->getMomentum()) * gamma;

	KNL = 1.0 / (1.0 + dp_p);
	px = arr[i][1];
	py = arr[i][3];

	arr[i][0] += KNL * length * px;
	arr[i][2] += KNL * length * py;
	phifac = (px * px +  py * py + dp_p * dp_p * gamma2i ) / 2.0;
	phifac = (phifac * KNL - dp_p * gamma2i) * KNL;
	arr[i][4] += length * phifac;
}

///////////////////////////////////////////////////////////////////////////
// NAME
//
//  drift
//
// DESCRIPTION
//  drifts a particle bunch
//
// PARAMETERS
//  bunch = reference to the macro-particle bunch
//  length = length of the drift
//
// RETURNS
//    Nothing
//
///////////////////////////////////////////////////////////////////////////

void drift(Bunch* bunch, double length){

	if(length <= 0.) return;

	double KNL, phifac, dp_p;

	SyncPart* syncPart = bunch->getSyncPart();
	double v = OrbitConst::c * syncPart->getBeta();
	if(length > 0.){
	   syncPart->setTime( syncPart->getTime() + length/v);
	}

	double gamma2i = 1.0 / (syncPart->getGamma() * syncPart->getGamma());
	double dp_p_coeff = 1./(syncPart->getMomentum()*syncPart->getBeta());

	//coordinate array [part. index][x,xp,y,yp,z,dE]
	double** arr = bunch->coordArr();

	for(int i = 0, n_part =  bunch->getSize(); i < n_part; i++){
		dp_p = arr[i][5] * dp_p_coeff;
		KNL = 1.0 / (1.0 + dp_p);
		arr[i][0] += KNL * length * arr[i][1];
		arr[i][2] += KNL * length * arr[i][3];
		phifac = (arr[i][1]*arr[i][1] + arr[i][3]*arr[i][3] + dp_p * dp_p * gamma2i) / 2.0;
		phifac = (phifac * KNL - dp_p * gamma2i) * KNL;
		arr[i][4] += length * phifac;
	}
}

///////////////////////////////////////////////////////////////////////////
// NAME
//
//  kick
//
// DESCRIPTION
//  kicks a particle bunch
//
// PARAMETERS
//  bunch = reference to the macro-particle bunch
//  kx = strength of the horizontal kick in rad
//  ky = strength of the vertical kick in rad
//  kE = strength of the energy kick in rad
//
// RETURNS
//    Nothing
//
///////////////////////////////////////////////////////////////////////////

void kick(Bunch* bunch, double kx, double ky, double kE)
{
	//coordinate array [part. index][x,xp,y,yp,z,dE]
	double** arr = bunch->coordArr();

	if(kx != 0.) {
		for(int i = 0; i < bunch->getSize(); i++){
      arr[i][1] += kx;
		}
  }
	if(ky != 0.) {
		for(int i = 0; i < bunch->getSize(); i++){
      arr[i][3] += ky;
		}
  }
	if(kE != 0.) {
		for(int i = 0; i < bunch->getSize(); i++){
      arr[i][5] += kE;
		}
  }
}

///////////////////////////////////////////////////////////////////////////
// NAME
//
//  multpi
//
// DESCRIPTION
//  gives particle a multipole momentum kick
//
// PARAMETERS
//  bunch =  reference to the macro-particle bunch
//  i = particle index
//  pole = multipole number
//  kl = integrated strength of the kick [m^(-pole)]
//  skew = 0 - normal, 1 - skew
//
// RETURNS
//  Nothing
//
///////////////////////////////////////////////////////////////////////////

void multpi(Bunch* bunch, int i, int pole, double kl, int skew){

  double z_re, z_im, zn_re, zn_im;
	double zn_re0, zn_im0;
  double kl1;

	//coordinate array [part. index][x,xp,y,yp,z,dE]
	double** arr = bunch->coordArr();

  kl1 = kl / factorial[pole];
  z_re = arr[i][0];
	z_im = arr[i][2];

  // take power of z to the n
	zn_re = 1.0;
	zn_im = 0.;
  for (int k = 0; k < pole; k++)
  {
    zn_re0 = zn_re * z_re - zn_im*z_im;
		zn_im0 = zn_im * z_re + zn_re*z_im;
		zn_re = zn_re0;
		zn_im = zn_im0;
  }

  // MAD Conventions on signs of multipole terms

  if(skew)
  {
    arr[i][1] += kl1 * zn_im;
    arr[i][3] += kl1 * zn_re;
  }
  else
  {
    arr[i][1] -= kl1 * zn_re;
    arr[i][3] += kl1 * zn_im;
  }
}

///////////////////////////////////////////////////////////////////////////
// NAME
//
//  multp
//
// DESCRIPTION
//  gives particle a multipole momentum kick
//
// PARAMETERS
//  bunch =  reference to the macro-particle bunch
//  pole = multipole number
//  kl = integrated strength of the kick [m^(-pole)]
//  skew = 0 - normal, 1 - skew
//
// RETURNS
//  Nothing
//
///////////////////////////////////////////////////////////////////////////

void multp(Bunch* bunch, int pole, double kl, int skew){

  double z_re, z_im, zn_re, zn_im;
	double zn_re0, zn_im0;
  double kl1;

	//coordinate array [part. index][x,xp,y,yp,z,dE]
	double** arr = bunch->coordArr();

  kl1 = kl / factorial[pole];

  // MAD Conventions on signs of multipole terms
	for(int i = 0; i < bunch->getSize(); i++){
		z_re = arr[i][0];
		z_im = arr[i][2];

		// take power of z to the n
		zn_re = 1.0;
		zn_im = 0.;
		for (int k = 0; k < pole; k++)
		{
			zn_re0 = zn_re * z_re - zn_im*z_im;
			zn_im0 = zn_im * z_re + zn_re*z_im;
			zn_re = zn_re0;
			zn_im = zn_im0;
		}

		if(skew)
		{
			arr[i][1] += kl1 * zn_im;
			arr[i][3] += kl1 * zn_re;
		}
		else
		{
			arr[i][1] -= kl1 * zn_re;
			arr[i][3] += kl1 * zn_im;
		}
	}
}


///////////////////////////////////////////////////////////////////////////
// NAME
//
//  multpfringeIN
//
// DESCRIPTION
//  Hard edge fringe field for a multipole
//
// PARAMETERS
//  bunch  = reference to the macro-particle bunch
//  pole = multipole number
//  kl = multipole strength
//  skew = multipole skew
//
// RETURNS
//  Nothing
//
///////////////////////////////////////////////////////////////////////////

void multpfringeIN(Bunch* bunch, int pole, double kl, int skew)
{
  int nm2, nm1, n, np1, np2;
  double x, px, y, py, cs, sn;
  double rfac, ifac, xyp2, xp2y, xy2;
	std::complex<double> z, pz, zn1, zn2, skewfacm, skewfacp;
	std::complex<double> xfac, pxfac, yfac, pyfac, phifac;
	std::complex<double> cdum1, cdum2;
  double kl1;
	double dp_p;

	SyncPart* syncPart = bunch->getSyncPart();

	double dp_p_coeff = 1./(syncPart->getMomentum()*syncPart->getBeta());

  nm2 = pole - 1;
  nm1 = pole;
  n = pole + 1;
  np1 = pole + 2;
  np2 = pole + 3;
  kl1 = kl / (4.0 * factorial[np1]);

	//coordinate array [part. index][x,xp,y,yp,z,dE]
	double** arr = bunch->coordArr();

  // MAD Conventions on signs of multipole terms

  if(skew)
  {
    cs = cos(OrbitConst::PI / (2. * n));
    sn = sin(OrbitConst::PI / (2. * n));
		skewfacm = std::complex<double>(cs, -sn);
    skewfacp = std::complex<double>(cs,  sn);
  }

  for(int i = 0, n_part =  bunch->getSize(); i < n_part; i++){

		dp_p = arr[i][5] * dp_p_coeff;

    x = arr[i][0];
    y = arr[i][2];
		z = std::complex<double>(x ,y);
    px = arr[i][1];
    py = arr[i][3];
		pz = std::complex<double>(px,py);

    if(skew)
    {
			z = z*skewfacp;
      pz = pz * skewfacp;
      x = std::real(z);
      y = std::imag(z);
      px = std::real(pz);
      py = std::imag(pz);
    }

    // take power of z to the n-1, n-2

    zn2 = std::complex<double>(1.,0.);
    for (int k = 0; k < nm2; k++)
    {
			zn2 = zn2 * z;
    }
		zn1 = zn2 * z;

    xyp2 = n * x * x + np2 * y * y;
    xp2y = np2 * x * x + n * y * y;
    xy2 = 2.0 * x * y;

    rfac = xyp2;
    ifac = -xy2;
    xfac = std::complex<double>(rfac, ifac);
    xfac = zn1 * xfac;

    rfac = nm1 * (px * xyp2 - py * xy2);
    ifac = nm1 * (-px * xy2 + py * xp2y);
    cdum1 = std::complex<double>(rfac, ifac);
    rfac = 2.0 * (n * px * x - py * y);
    ifac = 2.0 * (-px * y + np2 * py * x);
    cdum2 = std::complex<double>(rfac, ifac);
    pxfac = cdum1 + cdum2 * z;
    pxfac = zn2 * pxfac;

    rfac = -xy2;
    ifac = xp2y;
    yfac = std::complex<double>(rfac, ifac);
    yfac = zn1 * yfac;

    rfac = nm1 * (px * xy2 - py * xp2y);
    ifac = nm1 * (px * xyp2 - py * xy2);
    cdum1 = std::complex<double>(rfac, ifac);
    rfac = 2.0 * (np2 * px * y - py * x);
    ifac = 2.0 * (-px * x + n * py * y);
    cdum2 = std::complex<double>(rfac, ifac);
    pyfac = cdum1 + cdum2 * z;
    pyfac = zn2 * pyfac;

    rfac = px * xyp2 - py * xy2;
    ifac = -px * xy2 + py * xp2y;
    phifac = std::complex<double>(rfac, ifac);
    phifac = zn1 * phifac;

    if(skew)
    {
      xfac = xfac * skewfacm;
      pxfac = pxfac * skewfacm;
      yfac = yfac * skewfacm;
      pyfac = pyfac * skewfacm;
      phifac = phifac * skewfacm;
    }

    arr[i][0] += kl1 * std::real(xfac) / (1.0 + dp_p);
    arr[i][1] -= kl1 * std::real(pxfac) / (1.0 + dp_p);
    arr[i][2] += kl1 * std::real(yfac) / (1.0 + dp_p);
    arr[i][3] -= kl1 * std::real(pyfac) / (1.0 + dp_p);
    arr[i][4] += kl1 * std::real(yfac) /((1.0 + dp_p) * (1.0 + dp_p));
  }
}

///////////////////////////////////////////////////////////////////////////
// NAME
//
//  multpfringeOUT
//
// DESCRIPTION
//  Hard edge fringe field for a multipole
//
// PARAMETERS
//  bunch  = reference to the macro-particle bunch
//  pole = multipole number
//  kl = multipole strength
//  skew = multipole skew
//
// RETURNS
//  Nothing
//
///////////////////////////////////////////////////////////////////////////

void multpfringeOUT(Bunch* bunch, int pole, double kl,int skew)
{
  int nm2, nm1, n, np1, np2;
  double x, px, y, py, cs, sn;
  double rfac, ifac, xyp2, xp2y, xy2;
  std::complex<double> z, pz, zn1, zn2, skewfacm, skewfacp;
  std::complex<double> xfac, pxfac, yfac, pyfac, phifac;
  std::complex<double> cdum1, cdum2;
  double kl1;
	double dp_p;

	SyncPart* syncPart = bunch->getSyncPart();

	double dp_p_coeff = 1./(syncPart->getMomentum()*syncPart->getBeta());

  nm2 = pole - 1;
  nm1 = pole;
  n = pole + 1;
  np1 = pole + 2;
  np2 = pole + 3;
  kl1 = kl / (4.0 * factorial[np1]);

	//coordinate array [part. index][x,xp,y,yp,z,dE]
	double** arr = bunch->coordArr();

  // MAD Conventions on signs of multipole terms

  if(skew)
  {
    cs = cos(OrbitConst::PI / (2. * n));
    sn = sin(OrbitConst::PI / (2. * n));
    skewfacm = std::complex<double>(cs, -sn);
    skewfacp = std::complex<double>(cs,  sn);
  }

  for(int i = 1, n_part =  bunch->getSize(); i < n_part; i++){

		dp_p = arr[i][5] * dp_p_coeff;

    x = arr[i][0];
    y = arr[i][2];
    z = std::complex<double>(x, y);
    px = arr[i][1];
    py = arr[i][3];
    pz = std::complex<double>(px, py);

    if(skew)
    {
      z = z * skewfacp;
      pz = pz * skewfacp;
      x = std::real(z);
      y = std::imag(z);
      px = std::real(pz);
      py = std::imag(pz);
    }

    // take power of z to the n-1, n-2

    zn2 = std::complex<double>(1., 0.);
    for (int k = 0; k < nm2; k++)
    {
      zn2 = zn2 * z;
    }
    zn1 = zn2 * z;

    xyp2 = n * x * x + np2 * y * y;
    xp2y = np2 * x * x + n * y * y;
    xy2 = 2.0 * x * y;

    rfac = xyp2;
    ifac = -xy2;
    xfac = std::complex<double>(rfac, ifac);
    xfac = zn1 * xfac;

    rfac = nm1 * (px * xyp2 - py * xy2);
    ifac = nm1 * (-px * xy2 + py * xp2y);
    cdum1 = std::complex<double>(rfac, ifac);
    rfac = 2.0 * (n * px * x - py * y);
    ifac = 2.0 * (-px * y + np2 * py * x);
    cdum2 = std::complex<double>(rfac, ifac);
    pxfac = cdum1 + cdum2 * z;
    pxfac = zn2 * pxfac;

    rfac = -xy2;
    ifac = xp2y;
    yfac = std::complex<double>(rfac, ifac);
    yfac = zn1 * yfac;

    rfac = nm1 * (px * xy2 - py * xp2y);
    ifac = nm1 * (px * xyp2 - py * xy2);
    cdum1 = std::complex<double>(rfac, ifac);
    rfac = 2.0 * (np2 * px * y - py * x);
    ifac = 2.0 * (-px * x + n * py * y);
    cdum2 = std::complex<double>(rfac, ifac);
    pyfac = cdum1 + cdum2 * z;
    pyfac = zn2 * pyfac;

    rfac = px * xyp2 - py * xy2;
    ifac = -px * xy2 + py * xp2y;
    phifac = std::complex<double>(rfac, ifac);
    phifac = zn1 * phifac;

    if(skew)
    {
      xfac = xfac * skewfacm;
      pxfac = pxfac * skewfacm;
      yfac = yfac * skewfacm;
      pyfac = pyfac * skewfacm;
      phifac = phifac * skewfacm;
    }

    arr[i][0] -= kl1 * std::real(xfac) / (1.0 + dp_p);
    arr[i][1] += kl1 * std::real(pxfac) / (1.0 + dp_p);
    arr[i][2] -= kl1 * std::real(yfac) / (1.0 + dp_p);
    arr[i][3] += kl1 * std::real(pyfac) / (1.0 + dp_p);
    arr[i][4] -= kl1 * std::real(yfac)/((1.0 + dp_p) * (1.0 + dp_p));
  }
}


///////////////////////////////////////////////////////////////////////////
// NAME
//
//  quad1
//
// DESCRIPTION
//  quadrupole element one: linear transport matrix
//
// PARAMETERS
//  bunch  = reference to the macro-particle bunch
//  length = length of transport
//  kq = quadrupole field strength [m^(-2)]
//
// RETURNS
//  Nothing
//
///////////////////////////////////////////////////////////////////////////

void quad1(Bunch* bunch, double length, double kq)
{
  double x_init, xp_init, y_init, yp_init, z_init;
  double sqrt_kq, kqlength;
  double cx, sx, cy, sy, m11 = 0., m12 = 0., m21 = 0., m22 = 0.;
	double m33 = 0., m34 = 0., m43 = 0., m44 = 0.;

	SyncPart* syncPart = bunch->getSyncPart();
	double v = OrbitConst::c * syncPart->getBeta();
	if(length > 0.){
	   syncPart->setTime( syncPart->getTime() + length/v);
	}

	double gamma2i = 1.0 / (syncPart->getGamma() * syncPart->getGamma());
	double dp_p_coeff = 1./(syncPart->getMomentum()*syncPart->getBeta());

  if(kq >= 0.)
  {
    sqrt_kq = pow(kq, 0.5);
    kqlength = sqrt_kq * length;
    cx = cos(kqlength);
    sx = sin(kqlength);
    cy = cosh(kqlength);
    sy = sinh(kqlength);
    m11 = cx;
    m12 = sx / sqrt_kq;
    m21 = -sx * sqrt_kq;
    m22 = cx;
    m33 = cy;
    m34 = sy / sqrt_kq;
    m43 = sy * sqrt_kq;
    m44 = cy;
  }
  else if(kq < 0.)
  {
    sqrt_kq = pow(-kq, 0.5);
    kqlength = sqrt_kq * length;
    cx = cosh(kqlength);
    sx = sinh(kqlength);
    cy = cos(kqlength);
    sy = sin(kqlength);
    m11 = cx;
    m12 = sx / sqrt_kq;
    m21 = sx * sqrt_kq;
    m22 = cx;
    m33 = cy;
    m34 = sy / sqrt_kq;
    m43 = -sy * sqrt_kq;
    m44 = cy;
  }

	double dp_p;

	//coordinate array [part. index][x,xp,y,yp,z,dE]
	double** arr = bunch->coordArr();

	for(int i = 0, n_part =  bunch->getSize(); i < n_part; i++){
		dp_p = arr[i][5] * dp_p_coeff;
    x_init = arr[i][0];
    xp_init = arr[i][1];
    y_init = arr[i][2];
    yp_init = arr[i][3];
    z_init = arr[i][4];

    arr[i][0] = x_init * m11 + xp_init * m12;
    arr[i][1] = x_init * m21 + xp_init * m22;
    arr[i][2] = y_init * m33 + yp_init * m34;
    arr[i][3] = y_init * m43 + yp_init * m44;
    arr[i][4] = z_init - dp_p * gamma2i * length;
  }
}


///////////////////////////////////////////////////////////////////////////
// NAME
//
//  quad2
//
// DESCRIPTION
//  quadrupole transport through drift
//
// PARAMETERS
//  bunch   = reference to the macro-particle bunch
//  length = length of the element
//
// RETURNS
//    Nothing
//
///////////////////////////////////////////////////////////////////////////

void quad2(Bunch* bunch, double length)
{
  double x_init, y_init, z_init;
  double KNL, px, py, phifac;

	SyncPart* syncPart = bunch->getSyncPart();

	double gamma2i = 1.0 / (syncPart->getGamma() * syncPart->getGamma());
	double dp_p_coeff = 1./(syncPart->getMomentum()*syncPart->getBeta());

	double dp_p;

	//coordinate array [part. index][x,xp,y,yp,z,dE]
	double** arr = bunch->coordArr();

	for(int i = 0, n_part =  bunch->getSize(); i < n_part; i++){
		dp_p = arr[i][5] * dp_p_coeff;
    x_init = arr[i][0];
    y_init = arr[i][2];
    z_init = arr[i][4];
    KNL = 1.0 / (1.0 + dp_p);
    px = arr[i][1];
    py = arr[i][3];

    arr[i][0] = x_init - KNL * length * dp_p * arr[i][1];
    arr[i][2] = y_init - KNL * length * dp_p * arr[i][3];
    phifac = (px * px +
              py * py + dp_p * dp_p * gamma2i) / 2.0;
    phifac = (phifac * KNL + dp_p * dp_p * gamma2i) * KNL;
    arr[i][4] = z_init + length * phifac;
  }
}

///////////////////////////////////////////////////////////////////////////
// NAME
//
//  quadfringeIN
//
// DESCRIPTION
//  Hard edge fringe field for a quad
//
// PARAMETERS
//  bunch =  reference to the macro-particle bunch
//  kq  = strength of quad
//
// RETURNS
//  Nothing
//
///////////////////////////////////////////////////////////////////////////

void quadfringeIN(Bunch* bunch, double kq)
{
  double x_init, xp_init, y_init, yp_init, z_init, dp_init;

	SyncPart* syncPart = bunch->getSyncPart();

	double dp_p_coeff = 1./(syncPart->getMomentum()*syncPart->getBeta());

	double dp_p;

	//coordinate array [part. index][x,xp,y,yp,z,dE]
	double** arr = bunch->coordArr();

  for(int i = 0, n_part =  bunch->getSize(); i < n_part; i++)
  {
		dp_p = arr[i][5] * dp_p_coeff;
    x_init = arr[i][0];
    xp_init = arr[i][1];
    y_init = arr[i][2];
    yp_init = arr[i][3];
    z_init = arr[i][4];
    dp_init = dp_p;

    arr[i][0] = (x_init +
                        kq / (12. * (1. + dp_init)) * x_init
                        * (x_init * x_init + 3. * y_init * y_init));

    arr[i][1] = (xp_init -
                         kq / (4. * (1. + dp_init))
                         * (xp_init * (x_init * x_init + y_init * y_init)
	                    - 2. * yp_init * x_init * y_init));

    arr[i][2] = (y_init -
                        kq / (12. * (1. + dp_init)) * y_init
                        * (y_init * y_init + 3. * x_init * x_init));

    arr[i][3] = (yp_init -
                         kq / (4. * (1. + dp_init))
                         * (-yp_init * (x_init * x_init + y_init * y_init)
	                    + 2. * xp_init * x_init * y_init));

    arr[i][4] = z_init +
                        kq / (12. * (1. + dp_init) * (1. + dp_init))
                        * (xp_init * x_init *
                           (x_init * x_init + 3. * y_init * y_init) -
                           yp_init * y_init *
                           (y_init * y_init + 3. * x_init * x_init));
  }
}

///////////////////////////////////////////////////////////////////////////
// NAME
//
//  quadfringeOUT
//
// DESCRIPTION
//  Hard edge fringe field for a quad
//
// PARAMETERS
//  bunch  = reference to the macro-particle bunch
//  kq  = strength of quad
//
// RETURNS
//  Nothing
//
///////////////////////////////////////////////////////////////////////////

void quadfringeOUT(Bunch* bunch, double kq)
{

  double x_init, xp_init, y_init, yp_init, z_init, dp_init;

	SyncPart* syncPart = bunch->getSyncPart();

	double dp_p_coeff = 1./(syncPart->getMomentum()*syncPart->getBeta());

	double dp_p;

	//coordinate array [part. index][x,xp,y,yp,z,dE]
	double** arr = bunch->coordArr();

  for(int i = 0, n_part =  bunch->getSize(); i < n_part; i++)
  {
		dp_p = arr[i][5] * dp_p_coeff;
    x_init = arr[i][0];
    xp_init = arr[i][1];
    y_init = arr[i][2];
    yp_init = arr[i][3];
    z_init = arr[i][4];
    dp_init = dp_p;

    arr[i][0] = (x_init -
                        kq / (12. * (1. + dp_init)) * x_init
                        * (x_init * x_init + 3. * y_init * y_init));

    arr[i][1] = (xp_init +
                         kq / (4. * (1. + dp_init))
                         * (xp_init * (x_init * x_init + y_init * y_init)
	                    - 2. * yp_init * x_init * y_init));

    arr[i][2] = (y_init +
                        kq / (12. * (1. + dp_init)) * y_init
                        * (y_init * y_init + 3. * x_init * x_init));

    arr[i][3] = (yp_init +
                         kq / (4. * (1. + dp_init))
                         * (-yp_init * (x_init * x_init + y_init * y_init)
	                    + 2. * xp_init * x_init * y_init));

    arr[i][4] = z_init -
                        kq / (12. * (1. + dp_init) * (1. + dp_init))
                        * (xp_init * x_init *
                           (x_init * x_init + 3. * y_init * y_init) -
                           yp_init * y_init *
                           (y_init * y_init + 3. * x_init * x_init));
  }
}

///////////////////////////////////////////////////////////////////////////
// NAME
//
//  wedgerotate
//
// DESCRIPTION
//  Rotates coordinates by e for fringe fields at non-SBEND
//
// PARAMETERS
//  bunch  = reference to the macro-particle bunch
//  e = rotation angle
//  frinout = 0 before fringe, 1 after fringe
//
// RETURNS
//  Nothing
//
///////////////////////////////////////////////////////////////////////////

void wedgerotate(Bunch* bunch, double e, int frinout)
{
  double cs, sn;
  double x_init, xp_init,  z_init, p0_init, p0;


	SyncPart* syncPart = bunch->getSyncPart();

	double dp_p_coeff = 1./(syncPart->getMomentum()*syncPart->getBeta());

	double dp_p;

  cs = cos(e);
  sn = sin(e);

	//coordinate array [part. index][x,xp,y,yp,z,dE]
	double** arr = bunch->coordArr();

	int info = 1;

  for(int i = 0, n_part =  bunch->getSize(); i < n_part; i++)
  {
		info = 1;
		if(arr[i][5] == 0.) info = 0;
		dp_p = arr[i][5] * dp_p_coeff;
    x_init = arr[i][0];
    xp_init = arr[i][1];
    z_init = arr[i][4];
    p0_init = (1.0 + dp_p);

    if(frinout == 0)
    {
      arr[i][0] = x_init / cs;
      arr[i][1] = xp_init * cs + p0_init * sn;
      arr[i][4] = cs * z_init + sn * arr[i][0];
      p0 = -xp_init * sn + p0_init * cs;
      dp_p = p0 - 1.0;
    }
    else
    {
      arr[i][0] = x_init * cs;
      arr[i][1] = xp_init * cs - p0_init * sn;
      arr[i][4] = (z_init - sn * arr[i][0]) / cs;
      p0 = xp_init * sn + p0_init * cs;
      dp_p = p0 - 1.0;
    }
		arr[i][5] = dp_p/dp_p_coeff;
		if(info == 0) arr[i][5] = 0.;
  }
}

///////////////////////////////////////////////////////////////////////////
// NAME
//
//  wedgedrift
//
// DESCRIPTION
//  Drifts particles through wedge for non-SBEND
//
// PARAMETERS
//  bunch =  reference to the macro-particle bunch
//  e = wedge angle
//  inout = 0 for in, 1 for out
//
// RETURNS
//  Nothing
//
///////////////////////////////////////////////////////////////////////////

void wedgedrift(Bunch* bunch, double e, int inout)
{
  double cs, sn, ct, tn;
  double x_init, p0_init, s;

	SyncPart* syncPart = bunch->getSyncPart();
	double dp_p_coeff = 1./(syncPart->getMomentum()*syncPart->getBeta());
	double dp_p;

  cs = cos(e);
  sn = sin(e);
  ct = cs / sn;

	//coordinate array [part. index][x,xp,y,yp,z,dE]
	double** arr = bunch->coordArr();

  for(int i = 0, n_part =  bunch->getSize(); i < n_part; i++)
  {
		dp_p = arr[i][5] * dp_p_coeff;
    x_init = arr[i][0];

    if(inout == 0)
    {
      p0_init = (1.0 + dp_p);
      tn = arr[i][1] / p0_init;
      s = x_init / (ct - tn);
    }
    else
    {
      s = x_init / ct;
    }

    drifti(bunch, i, s);
  }
}

///////////////////////////////////////////////////////////////////////////
// NAME
//
//  wedgebend
//
// DESCRIPTION
//  Straight bends particles through wedge for non-SBEND
//
// PARAMETERS
//  bunch  = reference to the macro-particle bunch
//  e = wedge angle
//  inout = 0 for in, 1 for out
//  rho = radius of curvature
//  nsteps = number of integraton steps
//
// RETURNS
//  Nothing
//
///////////////////////////////////////////////////////////////////////////

void wedgebend(Bunch* bunch, double e, int inout, double rho, int nsteps)
{
  double kappa, cs, sn, ct, tn;
  double x_init, p0_init, s, sm, sm2;
  int nst;

	SyncPart* syncPart = bunch->getSyncPart();
	double dp_p_coeff = 1./(syncPart->getMomentum()*syncPart->getBeta());
	double dp_p;

  nst = nsteps / 2;
  if(nst < 1) nst = 1;
  kappa = 1.0 / rho;
  cs = cos(e);
  sn = sin(e);
  ct = cs / sn;

	//coordinate array [part. index][x,xp,y,yp,z,dE]
	double** arr = bunch->coordArr();

  for(int i = 0, n_part =  bunch->getSize(); i < n_part; i++)
  {
		dp_p = arr[i][5] * dp_p_coeff;
    x_init = arr[i][0];

    if(inout == 0)
    {
      s = -x_init / ct;
    }
    else
    {
      p0_init = (1.0 + dp_p);
      tn = arr[i][1] / p0_init;
      s = -x_init / (ct + tn);
    }

    sm = s / nst;
    sm2 = sm / 2.0;

    drifti(bunch, i, sm2);
    arr[i][1] = arr[i][1] - sm * kappa;
    for(int j = 1; j < nst; j++)
    {
      drifti(bunch, i, sm);
      arr[i][1] = arr[i][1] - sm * kappa;
    }
    drifti(bunch, i, sm2);
  }
}

///////////////////////////////////////////////////////////////////////////
// NAME
//
//  bend1
//
// DESCRIPTION
//  Linear bend transport
//
// PARAMETERS
//  bunch  = reference to the macro-particle bunch
//  length = length of transport
//  th = bending angle
//
// RETURNS
//  Nothing
//
///////////////////////////////////////////////////////////////////////////

void bend1(Bunch* bunch, double length, double th)
{
  double x_init, xp_init, y_init, yp_init, z_init;
  double cx, sx, rho;
  double m11, m12, m16, m21, m22, m26,
       m33, m34, m43, m44,
       m51, m52, m56;

	SyncPart* syncPart = bunch->getSyncPart();
	double v = OrbitConst::c * syncPart->getBeta();
	if(length > 0.){
	   syncPart->setTime( syncPart->getTime() + length/v);
	}

	double betasq = syncPart->getBeta() * syncPart->getBeta();
	double dp_p_coeff = 1./(syncPart->getMomentum()*syncPart->getBeta());

	double dp_p;

  rho = length / th;
  cx = cos(th);
  sx = sin(th);
  m11 = cx;
  m12 = rho * sx;
  m16 = rho * (1.0 - cx);
  m21 = -sx / rho;
  m22 = cx;
  m26 = sx;
  m33 = 1.0;
  m34 = length;
  m43 = 0.0;
  m44 = 1.0;
  m51 = sx;
  m52 = rho * (1.0 - cx);
  m56 = (betasq * length - rho * sx);

	//coordinate array [part. index][x,xp,y,yp,z,dE]
	double** arr = bunch->coordArr();

  for(int i = 0, n_part =  bunch->getSize(); i < n_part; i++)
  {
		dp_p = arr[i][5] * dp_p_coeff;
    x_init = arr[i][0];
    xp_init = arr[i][1];
    y_init = arr[i][2];
    yp_init = arr[i][3];
    z_init = arr[i][4];

    arr[i][0] = x_init * m11 + xp_init * m12 + dp_p * m16;
    arr[i][1] =x_init * m21 + xp_init * m22 + dp_p * m26;
    arr[i][2] = y_init + length * arr[i][3];

    arr[i][4] = z_init + x_init * m51 + xp_init * m52 + dp_p * m56;
  }
}

///////////////////////////////////////////////////////////////////////////
// NAME
//
//  bend2
//
// DESCRIPTION
//  Kinetic bend transport (same as nonlinear quad transport - quad2)
//
// PARAMETERS
//  bunch =  reference to the macro-particle bunch
//  length = length of element (either full of half step)
//
// RETURNS
//  Nothing
//
///////////////////////////////////////////////////////////////////////////

void bend2(Bunch* bunch, double length)
{
  double x_init, y_init, z_init;
  double KNL, px, py, phifac;

	SyncPart* syncPart = bunch->getSyncPart();

	double gamma2i = 1.0 / (syncPart->getGamma() * syncPart->getGamma());
	double dp_p_coeff = 1./(syncPart->getMomentum()*syncPart->getBeta());

	double dp_p;

	//coordinate array [part. index][x,xp,y,yp,z,dE]
	double** arr = bunch->coordArr();

  for(int i = 0, n_part =  bunch->getSize(); i < n_part; i++)
  {
		dp_p = arr[i][5] * dp_p_coeff;
    x_init = arr[i][0];
    y_init = arr[i][2];
    z_init = arr[i][4];
    KNL = 1.0 / (1.0 + dp_p);
    px = arr[i][1];
    py = arr[i][3];

    arr[i][0] = x_init - KNL * length * dp_p * arr[i][1];
    arr[i][2] = y_init - KNL * length * dp_p * arr[i][3];
    phifac = (px * px + py * py + dp_p * dp_p * gamma2i) / 2.0;
    phifac = (phifac * KNL + dp_p * dp_p * gamma2i) * KNL;
    arr[i][4] = z_init + length * phifac;
  }
}

///////////////////////////////////////////////////////////////////////////
// NAME
//
//  bend3
//
// DESCRIPTION
//  Nonlinear curvature bend transport
//  depending on py and dE in Hamiltonian
//
// PARAMETERS
//  bunch =  reference to the macro-particle bunch
//  th = bending angle
//
// RETURNS
//  Nothing
//
///////////////////////////////////////////////////////////////////////////

void bend3(Bunch* bunch, double th)
{
  double xp_init, y_init, z_init;
  double KNL, py, phifac;

	SyncPart* syncPart = bunch->getSyncPart();

	double gamma2i = 1.0 / (syncPart->getGamma() * syncPart->getGamma());
	double dp_p_coeff = 1./(syncPart->getMomentum()*syncPart->getBeta());

	double dp_p;

	//coordinate array [part. index][x,xp,y,yp,z,dE]
	double** arr = bunch->coordArr();

  for(int i = 0, n_part =  bunch->getSize(); i < n_part; i++)
  {
		dp_p = arr[i][5] * dp_p_coeff;
    xp_init = arr[i][1];
    y_init = arr[i][2];
    z_init = arr[i][4];
    KNL = 1.0 / (1.0 + dp_p);
    py = arr[i][3];

    phifac = (py * py +
              dp_p * dp_p * gamma2i
             ) / 2.0;
    arr[i][1] = xp_init - phifac * KNL * th;
    arr[i][2] = y_init + KNL * py * arr[i][0] * th;
    phifac = (phifac * KNL - dp_p * gamma2i) * KNL;
    arr[i][4] = z_init + th * phifac * arr[i][0];
  }
}

///////////////////////////////////////////////////////////////////////////
// NAME
//
//  bend4
//
// DESCRIPTION
//  Nonlinear curvature bend transport
//  depending on px in Hamiltonian
//
// PARAMETERS
//  bunch =  reference to the macro-particle bunch
//  th = bending angle
//
// RETURNS
//  Nothing
//
///////////////////////////////////////////////////////////////////////////

void bend4(Bunch* bunch, double th)
{
  double x_init, xp_init, z_init;
  double KNL, px, phifac, xfac;

	SyncPart* syncPart = bunch->getSyncPart();

	double dp_p_coeff = 1./(syncPart->getMomentum()*syncPart->getBeta());

	double dp_p;

	//coordinate array [part. index][x,xp,y,yp,z,dE]
	double** arr = bunch->coordArr();

  for(int i = 0, n_part =  bunch->getSize(); i < n_part; i++)
  {
		dp_p = arr[i][5] * dp_p_coeff;
    x_init = arr[i][0];
    xp_init = arr[i][1];
    z_init = arr[i][4];
    KNL = 1.0 / (1.0 + dp_p);
    px = arr[i][1];

    xfac = 1.0 + KNL * px * th / 2.0;
    phifac = (KNL * KNL * px * px) / 2.0;
    arr[i][0] = x_init * xfac * xfac;
    arr[i][1] = xp_init / xfac;
    arr[i][4] = z_init + th * phifac * arr[i][0];
  }
}

///////////////////////////////////////////////////////////////////////////
// NAME
//
//  bendfringeIN
//
// DESCRIPTION
//  Hard edge fringe field for a bend
//
// PARAMETERS
//  bunch =  reference to the macro-particle bunch
//  rho = radius of curvature for bending
//
// RETURNS
//  Nothing
//
///////////////////////////////////////////////////////////////////////////

void bendfringeIN(Bunch* bunch, double rho)
{

  double x_init, xp_init, y_init, yp_init, z_init, dp_init, kappa;

	SyncPart* syncPart = bunch->getSyncPart();

	double dp_p_coeff = 1./(syncPart->getMomentum()*syncPart->getBeta());

	double dp_p;

  kappa = 1.0 / rho;

	//coordinate array [part. index][x,xp,y,yp,z,dE]
	double** arr = bunch->coordArr();

  for(int i = 0, n_part =  bunch->getSize(); i < n_part; i++)
  {
		dp_p = arr[i][5] * dp_p_coeff;
    x_init = arr[i][0];
    xp_init = arr[i][1];
    y_init = arr[i][2];
    yp_init = arr[i][3];
    z_init = arr[i][4];
    dp_init = dp_p;

    arr[i][0] = (x_init
               + (kappa * y_init * y_init) / (2. * (1. + dp_init)));
    arr[i][3] = (yp_init
                - (kappa * xp_init * y_init) / (1. + dp_init));
    arr[i][4] = z_init + kappa * xp_init * y_init * y_init
                 / (2. * (1. + dp_init) * (1. + dp_init));
  }
}

///////////////////////////////////////////////////////////////////////////
// NAME
//
//  bendfringeOUT
//
// DESCRIPTION
//  Hard edge fringe field for a bend
//
// PARAMETERS
//  bunch =  reference to the macro-particle bunch
//  rho = radius of curvature for bending
//
// RETURNS
//  Nothing
//
///////////////////////////////////////////////////////////////////////////

void bendfringeOUT(Bunch* bunch, double rho)
{

  double x_init, xp_init, y_init, yp_init, z_init, dp_init, kappa;

	SyncPart* syncPart = bunch->getSyncPart();

	double dp_p_coeff = 1./(syncPart->getMomentum()*syncPart->getBeta());

	double dp_p;

  kappa = 1.0 / rho;

	//coordinate array [part. index][x,xp,y,yp,z,dE]
	double** arr = bunch->coordArr();

  for(int i = 0, n_part =  bunch->getSize(); i < n_part; i++)
  {
		dp_p = arr[i][5] * dp_p_coeff;
    x_init = arr[i][0];
    xp_init = arr[i][1];
    y_init = arr[i][2];
    yp_init = arr[i][3];
    z_init = arr[i][4];
    dp_init = dp_p;

    arr[i][0] = (x_init
               - (kappa * y_init * y_init) / (2. * (1. + dp_init)));
    arr[i][3] = (yp_init
                + (kappa * xp_init * y_init) / (1. + dp_init));
    arr[i][4] = z_init - kappa * xp_init * y_init * y_init
                 / (2. * (1. + dp_init) * (1. + dp_init));
  }
}

///////////////////////////////////////////////////////////////////////////
// NAME
//
//  soln
//
// DESCRIPTION
//  Integration through a solenoid
//
// PARAMETERS
//  bunch =  reference to the macro-particle bunch
//  length = integration length
//  B   = magnetic field (1/m)
//
// RETURNS
//  Nothing
//
///////////////////////////////////////////////////////////////////////////

void soln(Bunch* bunch, double length, double B)
{

  double x_init, px_init, y_init, py_init, z_init;
  double KNL, phase, cs, sn;
  double cu, cpu, u_init, pu_init, u, pu, phifac;

	SyncPart* syncPart = bunch->getSyncPart();
	double v = OrbitConst::c * syncPart->getBeta();
	if(length > 0.){
	   syncPart->setTime( syncPart->getTime() + length/v);
	}

	double gamma2i = 1.0 / (syncPart->getGamma() * syncPart->getGamma());
	double dp_p_coeff = 1./(syncPart->getMomentum()*syncPart->getBeta());

	double dp_p;

	//coordinate array [part. index][x,xp,y,yp,z,dE]
	double** arr = bunch->coordArr();

  for(int i = 0, n_part =  bunch->getSize(); i < n_part; i++)
  {
		dp_p = arr[i][5] * dp_p_coeff;
    x_init = arr[i][0];
    px_init = arr[i][1];
    y_init = arr[i][2];
    py_init = arr[i][3];
    z_init = arr[i][4];
    KNL = 1.0 / (1.0 + dp_p);

    cu = y_init / 2. - px_init / B;
    cpu = x_init * B / 2. + py_init;
    u_init = y_init / 2. + px_init / B;
    pu_init = -x_init * B / 2. + py_init;
    phase = KNL * B * length;
    cs = cos(phase);
    sn = sin(phase);

    u = u_init * cs + pu_init * sn / B;
    pu = -u_init * B * sn + pu_init * cs;

    arr[i][0] =  (-pu + cpu) / B;
    arr[i][1] = 0.5 * (u - cu) * B;
    arr[i][2] =  (u + cu);
    arr[i][3] = 0.5 * (pu + cpu);

    phifac = (pu_init * pu_init +
              B * B * u_init * u_init +
              dp_p * dp_p * gamma2i
             ) / 2.0;
    phifac = (phifac * KNL - dp_p * gamma2i) * KNL;
    arr[i][4] = z_init + length * phifac;
  }
}

///////////////////////////////////////////////////////////////////////////
// NAME
//
//  solnfringeIN
//
// DESCRIPTION
//  Hard edge fringe field for a solenoid
//
// PARAMETERS
//  bunch =  reference to the macro-particle bunch
//  B   = magnetic field (1/m)
//
// RETURNS
//  Nothing
//
///////////////////////////////////////////////////////////////////////////

void solnfringeIN(Bunch* bunch, double B)
{
}

///////////////////////////////////////////////////////////////////////////
// NAME
//
//  solnfringeOUT
//
// DESCRIPTION
//  Hard edge fringe field for a solenoid
//
// PARAMETERS
//  bunch =  reference to the macro-particle bunch
//  B   = magnetic field (1/m)
//
// RETURNS
//  Nothing
//
///////////////////////////////////////////////////////////////////////////

void solnfringeOUT(Bunch* bunch, double B)
{
}

///////////////////////////////////////////////////////////////////////////
// NAME
//
//  wedgebendCF
//
// DESCRIPTION
//  Straight bends particles through wedge for Combined Function non-SBEND
//
// PARAMETERS
//  bunch =  reference to the macro-particle bunch
//  e = wedge angle
//  inout = 0 for in, 1 for out
//  rho = radius of curvature
//  vecnum = number of multipole terms
//  pole = multipolarities of multipole terms
//  kl = integrated strengths of multipole terms
//  skew = skewness  of multipole terms
//  nsteps = number of integraton steps
//
// RETURNS
//  Nothing
//
///////////////////////////////////////////////////////////////////////////

void wedgebendCF(Bunch* bunch, double e, int inout,
                 double rho,
                 int vecnum,
                 std::vector<int>& pole,
                 std::vector<double>& kl,
                 std::vector<int>& skew,
                 int nsteps)
{
  double kappa, cs, sn, ct, tn;
  double x_init, p0_init, s, sm, sm2, klint;
  int nst;

	SyncPart* syncPart = bunch->getSyncPart();
	double dp_p_coeff = 1./(syncPart->getMomentum()*syncPart->getBeta());
	double dp_p;

  nst = nsteps / 2;
  if(nst < 1) nst = 1;
  kappa = 1.0 / rho;
  cs = cos(e);
  sn = sin(e);
  ct = cs / sn;

	//coordinate array [part. index][x,xp,y,yp,z,dE]
	double** arr = bunch->coordArr();

  for(int i = 0, n_part =  bunch->getSize(); i < n_part; i++)
  {
		dp_p = arr[i][5] * dp_p_coeff;
    x_init = arr[i][0];

    if(inout == 0)
    {
      s = -x_init / ct;
    }
    else
    {
      p0_init = (1.0 + dp_p);
      tn = arr[i][1] / p0_init;
      s = -x_init / (ct + tn);
    }

    sm = s / nst;
    sm2 = sm / 2.0;

    drifti(bunch, i, sm2);
    arr[i][1] = arr[i][1] - sm * kappa;
    for (int l = 0; l < vecnum; l++)
    {
      klint = kl[l] * sm;
      multpi(bunch, i, pole[l], klint, skew[l]);
    }
    for(int j = 1; j < nst; j++)
    {
      drifti(bunch, i, sm);
      arr[i][1] = arr[i][1] - sm * kappa;
      for (int l = 0; l < vecnum; l++)
      {
        klint = kl[l] * sm;
        multpi(bunch, i, pole[l], klint, skew[l]);
      }
    }
    drifti(bunch, i, sm2);
  }
}


///////////////////////////////////////////////////////////////////////////
// NAME
//
//  RingRF
//
// DESCRIPTION
//  The ring type of RF cavity. Transition time factor T(k) = const = T(k0).
//  There is no need for a simplectic phase correction.
//
// PARAMETERS
//  bunch =  reference to the macro-particle bunch
//  harmonic_numb = harmonics number
//  voltage = voltage in Giga Volts
//  phase_s = synchronous phase in Rad
//
// RETURNS
//  Nothing
//
///////////////////////////////////////////////////////////////////////////

void RingRF(Bunch* bunch, double ring_length, int harmonic_numb, double voltage, double phase_s)
{
	double deltaV = 0.;
	double charge = bunch->getCharge();
	double coeff =  charge;

	double Factor = 2.0*OrbitConst::PI/ring_length;

	SyncPart* syncPart = bunch->getSyncPart();
	if(phase_s != 0.){
		double kin_e = syncPart->getEnergy();
		kin_e = kin_e + coeff * voltage * sin(phase_s);
		syncPart->setMomentum(syncPart->energyToMomentum(kin_e));
	}

	//coordinate array [part. index][x,xp,y,yp,z,dE]
	double** arr = bunch->coordArr();

  for(int i = 0, n_part =  bunch->getSize(); i < n_part; i++)
  {
		deltaV = voltage * ( sin(harmonic_numb*Factor*arr[i][4] + phase_s));
		arr[i][5] += coeff * deltaV;
	}
}


}  //end of namespace teapot_base
