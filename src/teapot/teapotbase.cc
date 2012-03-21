/////////////////////////////////////////////////////////////////////////////
//
// FILE NAME
//   teapotbase.cc
//
// AUTHOR
//   Jeff Holmes, ORNL, jzh@ornl.gov
//   Joshua Abrams, Knox College, jabrams@knox.edu
//   Steven Bunch, University of Tennessee, sbunch2@utk.edu
//
// Modified by Andrei Shishlo
//   12/30/05
//
// Checked by Jeff Holmes
//   02/2012
//
// DESCRIPTION
//   Define elementary functions for different elements
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

#include "teapotbase.hh"
#include "OrbitConst.hh"
#include "Bunch.hh"
#include "SyncPart.hh"

#include <complex>

namespace teapot_base
{
    static double* factorial = NULL;

    void init_factorial()
    {
        if(factorial == NULL)
        {
            int n = 50;
            factorial = new double[n];
            factorial[0] = 1.0;
            for(int i = 1; i < n; i++)
            {
                factorial[i] = i * factorial[i - 1];
            }
        }
    }

    void delete_factorial()
    {
        delete [] factorial;
    }

///////////////////////////////////////////////////////////////////////////
// NAME
//   rotatexy
//
// DESCRIPTION
//   Rotates particle coordinates
//
// PARAMETERS
//   bunch = reference to the macro-particle bunch
//   anglexy = rotation angle
//
// RETURNS
//   Nothing
//
///////////////////////////////////////////////////////////////////////////

void rotatexy(Bunch* bunch, double anglexy)
{
    double xtemp, pxtemp, ytemp, pytemp;
    double cs = cos(anglexy);
    double sn = sin(anglexy);

    //coordinate array [part. index][x,xp,y,yp,z,dE]
    double** arr = bunch->coordArr();

    for(int i = 0; i < bunch->getSize(); i++)
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
//   drifti
//
// DESCRIPTION
//   Drifts a single particle. Length < 0 is allowed.
//
// PARAMETERS
//   bunch = reference to the macro-particle bunch
//   i = particle index
//   length = length of the drift
//
// RETURNS
//   Nothing
//
///////////////////////////////////////////////////////////////////////////

void drifti(Bunch* bunch, int i, double length)
{
    double KNL, phifac, dp_p;

    SyncPart* syncPart = bunch->getSyncPart();

    double gamma2i = 1.0 / (syncPart->getGamma() * syncPart->getGamma());
    double dp_p_coeff = 1.0 / (syncPart->getMomentum() * syncPart->getBeta());

    //coordinate array [part. index][x,xp,y,yp,z,dE]
    double** arr = bunch->coordArr();

    dp_p = arr[i][5] * dp_p_coeff;
    KNL  = 1.0 / (1.0 + dp_p);
    arr[i][0] += KNL * length * arr[i][1];
    arr[i][2] += KNL * length * arr[i][3];
    phifac = (arr[i][1] * arr[i][1] + arr[i][3] * arr[i][3] +
              dp_p * dp_p * gamma2i) / 2.0;
    phifac = (phifac * KNL - dp_p * gamma2i) * KNL;
    arr[i][4] -= length * phifac;
}

///////////////////////////////////////////////////////////////////////////
// NAME
//   drift
//
// DESCRIPTION
//   Drifts a particle bunch. Length < 0 is allowed.
//
// PARAMETERS
//   bunch = reference to the macro-particle bunch
//   length = length of the drift
//
// RETURNS
//   Nothing
//
///////////////////////////////////////////////////////////////////////////

void drift(Bunch* bunch, double length)
{
    double KNL, phifac, dp_p;

    SyncPart* syncPart = bunch->getSyncPart();
    
    double v = OrbitConst::c * syncPart->getBeta();
    if(length > 0.)
    {
	   syncPart->setTime(syncPart->getTime() + length / v);
    }

    double gamma2i = 1.0 / (syncPart->getGamma() * syncPart->getGamma());
    double dp_p_coeff = 1.0 / (syncPart->getMomentum() * syncPart->getBeta());

    //coordinate array [part. index][x,xp,y,yp,z,dE]
    double** arr = bunch->coordArr();

    for(int i = 0; i < bunch->getSize(); i++)
    {
        dp_p = arr[i][5] * dp_p_coeff;
        KNL  = 1.0 / (1.0 + dp_p);
        arr[i][0] += KNL * length * arr[i][1];
        arr[i][2] += KNL * length * arr[i][3];
        phifac = (arr[i][1] * arr[i][1] + arr[i][3] * arr[i][3] +
                  dp_p * dp_p * gamma2i) / 2.0;
        phifac = (phifac * KNL - dp_p * gamma2i) * KNL;
        arr[i][4] -= length * phifac;
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

void kick(Bunch* bunch, double kx, double ky, double kE)
{
    //coordinate array [part. index][x,xp,y,yp,z,dE]
    double** arr = bunch->coordArr();

    if(kx != 0.)
    {
        for(int i = 0; i < bunch->getSize(); i++)
        {
            arr[i][1] += kx;
        }
    }
    if(ky != 0.)
    {
        for(int i = 0; i < bunch->getSize(); i++)
        {
            arr[i][3] += ky;
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

///////////////////////////////////////////////////////////////////////////
// NAME
//   multpi
//
// DESCRIPTION
//   Gives particle a multipole momentum kick
//
// PARAMETERS
//   bunch =  reference to the macro-particle bunch
//   i = particle index
//   pole = multipole number
//   kl = integrated strength of the kick [m^(-pole)]
//   skew = 0 - normal, 1 - skew
//
// RETURNS
//   Nothing
//
///////////////////////////////////////////////////////////////////////////

void multpi(Bunch* bunch, int i, int pole, double kl, int skew)
{
    std::complex<double> z, zn;
    double kl1;

    //coordinate array [part. index][x,xp,y,yp,z,dE]
    double** arr = bunch->coordArr();

    kl1 = kl / factorial[pole];

    z = (arr[i][0], arr[i][2]);

    // take power of z to the n
    zn = (1.0, 0.0);
    for (int k = 0; k < pole; k++)
    {
        zn *= z;
    }

    // MAD Conventions on signs of multipole terms
    if(skew)
    {
        arr[i][1] += kl1 * std::imag(zn);
        arr[i][3] += kl1 * std::real(zn);
    }
    else
    {
        arr[i][1] -= kl1 * std::real(zn);
        arr[i][3] += kl1 * std::imag(zn);
    }
}

///////////////////////////////////////////////////////////////////////////
// NAME
//   multp
//
// DESCRIPTION
//   Gives particles multipole momentum kicks
//
// PARAMETERS
//   bunch =  reference to the macro-particle bunch
//   pole = multipole number
//   kl = integrated strength of the kick [m^(-pole)]
//   skew = 0 - normal, 1 - skew
//
// RETURNS
//   Nothing
//
///////////////////////////////////////////////////////////////////////////

void multp(Bunch* bunch, int pole, double kl, int skew)
{
    std::complex<double> z, zn;
    double kl1;

    //coordinate array [part. index][x,xp,y,yp,z,dE]
    double** arr = bunch->coordArr();

    kl1 = kl / factorial[pole];
    
    for(int i = 0; i < bunch->getSize(); i++)
    {
        z = (arr[i][0], arr[i][2]);

        // take power of z to the n
        zn = (1.0, 0.0);
        for (int k = 0; k < pole; k++)
        {
            zn *= z;
        }

        // MAD Conventions on signs of multipole terms
        if(skew)
        {
            arr[i][1] += kl1 * std::imag(zn);
            arr[i][3] += kl1 * std::real(zn);
        }
        else
        {
            arr[i][1] -= kl1 * std::real(zn);
            arr[i][3] += kl1 * std::imag(zn);
        }
    }
}

///////////////////////////////////////////////////////////////////////////
// NAME
//   multpfringeIN
//
// DESCRIPTION
//   Hard edge fringe field for a multipole
//
// PARAMETERS
//   bunch  = reference to the macro-particle bunch
//   pole = multipole number
//   kl = multipole strength
//   skew = multipole skew
//
// RETURNS
//   Nothing
//
///////////////////////////////////////////////////////////////////////////

void multpfringeIN(Bunch* bunch, int pole, double kl, int skew)
{
    std::complex<double> rootm1 = std::complex<double>(0.0, 1.0);

    SyncPart* syncPart = bunch->getSyncPart();

    double dp_p_coeff = 1.0 / (syncPart->getMomentum() * syncPart->getBeta());

    int lm1 = pole;
    int l   = pole + 1;
    int lp1 = pole + 2;
    int lp2 = pole + 3;
    std::complex<double> cxlm1 = std::complex<double>(lm1, 0.0);
    std::complex<double> cxl   = std::complex<double>(l  , 0.0);
    std::complex<double> cxlp1 = std::complex<double>(lp1, 0.0);
    std::complex<double> cxlp2 = std::complex<double>(lp2, 0.0);

    double klfactlp1 = kl / (4.0 * factorial[lp1]);

    // MAD Conventions on signs of multipole terms

    std::complex<double> kterm;
    if(skew)
    {
        kterm = std::complex<double>(0.0, klfactlp1);
    }
    else
    {
        kterm = std::complex<double>(klfactlp1, 0.0);
    }

    double dp_p, KNL;

    //coordinate array [part. index][x,xp,y,yp,z,dE]
    double** arr = bunch->coordArr();

    for(int i = 0; i < bunch->getSize(); i++)
    {
        double x = arr[i][0];
        double y = arr[i][2];
        std::complex<double> z = std::complex<double>(x, y);
        double px = arr[i][1];
        double py = arr[i][3];
        dp_p = arr[i][5] * dp_p_coeff;
        KNL  = 1.0 / (1.0 + dp_p);

        // take power of z to the lm1, l

        std::complex<double> zlm1 = std::complex<double>(1., 0.);
        for (int k = 0; k < lm1; k++)
        {
            zlm1 *= z;
        }
        std::complex<double> zl = zlm1 * z;

        std::complex<double> fxterm   = std::complex<double>(l * x, -lp2 * y);
        std::complex<double> dxfxterm = cxl * (fxterm + z);
        std::complex<double> dyfxterm = rootm1 * (cxl * fxterm - cxlp2 * z);
        std::complex<double> fyterm   = std::complex<double>(l * y, lp2 * x);
        std::complex<double> dxfyterm = cxl * fyterm + rootm1 * cxlp2 * z;
        std::complex<double> dyfyterm = cxl * (rootm1 * fyterm + z);

        std::complex<double> fxcx   = -kterm * zl   * fxterm;
        std::complex<double> dxfxcx = -kterm * zlm1 * dxfxterm;
        std::complex<double> dyfxcx = -kterm * zlm1 * dyfxterm;
        std::complex<double> fycx   = -kterm * zl   * fyterm;
        std::complex<double> dxfycx = -kterm * zlm1 * dxfyterm;
        std::complex<double> dyfycx = -kterm * zlm1 * dyfyterm;

        arr[i][0] -= std::real(fxcx) * KNL;
        arr[i][2] -= std::real(fycx) * KNL;

        double M11 = 1.0 - std::real(dxfxcx) * KNL;
        double M12 =     - std::real(dxfycx) * KNL;
        double M21 =     - std::real(dyfxcx) * KNL;
        double M22 = 1.0 - std::real(dyfycx) * KNL;
        double detM = M11 * M22 - M12 * M21;

        double pxnew = ( M22 * px - M12 * py) / detM;
        double pynew = (-M21 * px + M11 * py) / detM;

        arr[i][1] = pxnew;
        arr[i][3] = pynew;

        arr[i][4] -= (pxnew * std::real(fxcx) + pynew * std::real(fycx)) *
                     KNL * KNL;
    }
}

///////////////////////////////////////////////////////////////////////////
// NAME
//   multpfringeOUT
//
// DESCRIPTION
//   Hard edge fringe field for a multipole
//
// PARAMETERS
//   bunch  = reference to the macro-particle bunch
//   pole = multipole number
//   kl = multipole strength
//   skew = multipole skew
//
// RETURNS
//   Nothing
//
///////////////////////////////////////////////////////////////////////////

void multpfringeOUT(Bunch* bunch, int pole, double kl, int skew)
{
    std::complex<double> rootm1 = std::complex<double>(0.0, 1.0);

    SyncPart* syncPart = bunch->getSyncPart();

    double dp_p_coeff = 1.0 / (syncPart->getMomentum() * syncPart->getBeta());

    int lm1 = pole;
    int l   = pole + 1;
    int lp1 = pole + 2;
    int lp2 = pole + 3;
    std::complex<double> cxlm1 = std::complex<double>(lm1, 0.0);
    std::complex<double> cxl   = std::complex<double>(l  , 0.0);
    std::complex<double> cxlp1 = std::complex<double>(lp1, 0.0);
    std::complex<double> cxlp2 = std::complex<double>(lp2, 0.0);

    double klfactlp1 = kl / (4.0 * factorial[lp1]);

    // MAD Conventions on signs of multipole terms

    std::complex<double> kterm;
    if(skew)
    {
        kterm = std::complex<double>(0.0, klfactlp1);
    }
    else
    {
        kterm = std::complex<double>(klfactlp1, 0.0);
    }

    double dp_p, KNL;

    //coordinate array [part. index][x,xp,y,yp,z,dE]
    double** arr = bunch->coordArr();

    for(int i = 0; i < bunch->getSize(); i++)
    {
        double x = arr[i][0];
        double y = arr[i][2];
        std::complex<double> z = std::complex<double>(x, y);
        double px = arr[i][1];
        double py = arr[i][3];
        dp_p = arr[i][5] * dp_p_coeff;
        KNL  = 1.0 / (1.0 + dp_p);

        // take power of z to the lm1, l

        std::complex<double> zlm1 = std::complex<double>(1., 0.);
        for (int k = 0; k < lm1; k++)
        {
            zlm1 *= z;
        }
        std::complex<double> zl = zlm1 * z;

        std::complex<double> fxterm   = std::complex<double>(l * x, -lp2 * y);
        std::complex<double> dxfxterm = cxl * (fxterm + z);
        std::complex<double> dyfxterm = rootm1 * (cxl * fxterm - cxlp2 * z);
        std::complex<double> fyterm   = std::complex<double>(l * y, lp2 * x);
        std::complex<double> dxfyterm = cxl * fyterm + rootm1 * cxlp2 * z;
        std::complex<double> dyfyterm = cxl * (rootm1 * fyterm + z);

        std::complex<double> fxcx   = kterm * zl   * fxterm;
        std::complex<double> dxfxcx = kterm * zlm1 * dxfxterm;
        std::complex<double> dyfxcx = kterm * zlm1 * dyfxterm;
        std::complex<double> fycx   = kterm * zl   * fyterm;
        std::complex<double> dxfycx = kterm * zlm1 * dxfyterm;
        std::complex<double> dyfycx = kterm * zlm1 * dyfyterm;

        arr[i][0] -= std::real(fxcx) * KNL;
        arr[i][2] -= std::real(fycx) * KNL;

        double M11 = 1.0 - std::real(dxfxcx) * KNL;
        double M12 =     - std::real(dxfycx) * KNL;
        double M21 =     - std::real(dyfxcx) * KNL;
        double M22 = 1.0 - std::real(dyfycx) * KNL;
        double detM = M11 * M22 - M12 * M21;

        double pxnew = ( M22 * px - M12 * py) / detM;
        double pynew = (-M21 * px + M11 * py) / detM;

        arr[i][1] = pxnew;
        arr[i][3] = pynew;

        arr[i][4] -= (pxnew * std::real(fxcx) + pynew * std::real(fycx)) *
                     KNL * KNL;
    }
}

////////////////////////////
// NAME
//   quad1
//
// DESCRIPTION
//   Quadrupole element one: linear transport matrix
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

void quad1(Bunch* bunch, double length, double kq)
{
    if(kq == 0.)
    {
        drift(bunch,length);
        return;
    }

    double x_init, xp_init, y_init, yp_init;
    double sqrt_kq, kqlength;
    double cx, sx, cy, sy;
    double m11 = 0., m12 = 0., m21 = 0., m22 = 0.;
    double m33 = 0., m34 = 0., m43 = 0., m44 = 0.;

    SyncPart* syncPart = bunch->getSyncPart();

    double v = OrbitConst::c * syncPart->getBeta();
    if(length > 0.)
    {
        syncPart->setTime(syncPart->getTime() + length / v);
    }

    double gamma2i = 1.0 / (syncPart->getGamma() * syncPart->getGamma());
    double dp_p_coeff = 1.0 /(syncPart->getMomentum() * syncPart->getBeta());

    if(kq > 0.)
    {
        sqrt_kq  = pow(kq, 0.5);
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
        sqrt_kq  = pow(-kq, 0.5);
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

    for(int i = 0; i < bunch->getSize(); i++)
    {
        dp_p    = arr[i][5] * dp_p_coeff;
        x_init  = arr[i][0];
        xp_init = arr[i][1];
        y_init  = arr[i][2];
        yp_init = arr[i][3];

        arr[i][0]  = x_init * m11 + xp_init * m12;
        arr[i][1]  = x_init * m21 + xp_init * m22;
        arr[i][2]  = y_init * m33 + yp_init * m34;
        arr[i][3]  = y_init * m43 + yp_init * m44;
        arr[i][4] += dp_p * gamma2i * length;
    }
}

///////////////////////////////////////////////////////////////////////////
// NAME
//   quad2
//
// DESCRIPTION
//   Quadrupole element two: nonlinear piece
//
// PARAMETERS
//   bunch  = reference to the macro-particle bunch
//   length = length of the element
//
// RETURNS
//   Nothing
//
///////////////////////////////////////////////////////////////////////////

void quad2(Bunch* bunch, double length)
{
    double KNL, phifac;

    SyncPart* syncPart = bunch->getSyncPart();

    double gamma2i = 1.0 / (syncPart->getGamma() * syncPart->getGamma());
    double dp_p_coeff = 1.0 / (syncPart->getMomentum() * syncPart->getBeta());
    double dp_p;

    //coordinate array [part. index][x,xp,y,yp,z,dE]
    double** arr = bunch->coordArr();

    for(int i = 0; i < bunch->getSize(); i++)
    {
        dp_p = arr[i][5] * dp_p_coeff;
        KNL = 1.0 / (1.0 + dp_p);

        arr[i][0] -= KNL * length * dp_p * arr[i][1];
        arr[i][2] -= KNL * length * dp_p * arr[i][3];
        phifac = (arr[i][1] * arr[i][1] + arr[i][3] * arr[i][3] +
                  dp_p * dp_p * gamma2i) / 2.0;
        phifac = (phifac * KNL + dp_p * dp_p * gamma2i) * KNL;
        arr[i][4] -= length * phifac;
    }
}

///////////////////////////////////////////////////////////////////////////
// NAME
//   quadfringeIN
//
// DESCRIPTION
//   Hard edge fringe field for a quad
//
// PARAMETERS
//   bunch =  reference to the macro-particle bunch
//   kq  = strength of quad
//
// RETURNS
//   Nothing
//
///////////////////////////////////////////////////////////////////////////

void quadfringeIN(Bunch* bunch, double kq)
{
    double KNL, x_init, xp_init, y_init, yp_init;

    SyncPart* syncPart = bunch->getSyncPart();

    double dp_p_coeff = 1.0 / (syncPart->getMomentum() * syncPart->getBeta());
    double dp_p;

    //coordinate array [part. index][x,xp,y,yp,z,dE]
    double** arr = bunch->coordArr();

    for(int i = 0; i < bunch->getSize(); i++)
    {
        dp_p    = arr[i][5] * dp_p_coeff;
        KNL     = 1.0 / (1.0 + dp_p);
        x_init  = arr[i][0];
        xp_init = arr[i][1];
        y_init  = arr[i][2];
        yp_init = arr[i][3];

        arr[i][0] += (kq * KNL / 12.) * x_init *
                     (x_init * x_init + 3. * y_init * y_init);

        arr[i][1] -= (kq * KNL / 4.) *
                     (xp_init * (x_init * x_init + y_init * y_init) -
                      2. * yp_init * x_init * y_init);

        arr[i][2] -= (kq * KNL / 12.) * y_init *
                     (y_init * y_init + 3. * x_init * x_init);

        arr[i][3] -= (kq * KNL / 4.) *
                     (-yp_init * (x_init * x_init + y_init * y_init) +
                      2. * xp_init * x_init * y_init);

        arr[i][4] += (kq * KNL * KNL / 12.) *
                     (xp_init * x_init *
                      (x_init * x_init + 3. * y_init * y_init) -
                      yp_init * y_init *
                      (y_init * y_init + 3. * x_init * x_init));
    }
}

///////////////////////////////////////////////////////////////////////////
// NAME
//   quadfringeOUT
//
// DESCRIPTION
//   Hard edge fringe field for a quad
//
// PARAMETERS
//   bunch  = reference to the macro-particle bunch
//   kq  = strength of quad
//
// RETURNS
//   Nothing
//
///////////////////////////////////////////////////////////////////////////

void quadfringeOUT(Bunch* bunch, double kq)
{
    double KNL, x_init, xp_init, y_init, yp_init;

    SyncPart* syncPart = bunch->getSyncPart();

    double dp_p_coeff = 1.0 / (syncPart->getMomentum() * syncPart->getBeta());
    double dp_p;

    //coordinate array [part. index][x,xp,y,yp,z,dE]
    double** arr = bunch->coordArr();

    for(int i = 0; i < bunch->getSize(); i++)
    {
        dp_p    = arr[i][5] * dp_p_coeff;
        KNL     = 1.0 / (1.0 + dp_p);
        x_init  = arr[i][0];
        xp_init = arr[i][1];
        y_init  = arr[i][2];
        yp_init = arr[i][3];

        arr[i][0] -= (kq * KNL / 12.) * x_init *
                     (x_init * x_init + 3. * y_init * y_init);

        arr[i][1] += (kq * KNL / 4.) *
                     (xp_init * (x_init * x_init + y_init * y_init) -
                      2. * yp_init * x_init * y_init);

        arr[i][2] += (kq * KNL / 12.) * y_init *
                     (y_init * y_init + 3. * x_init * x_init);

        arr[i][3] += (kq * KNL / 4.) *
                     (-yp_init * (x_init * x_init + y_init * y_init) +
                      2. * xp_init * x_init * y_init);

        arr[i][4] -= (kq * KNL * KNL / 12.) *
                     (xp_init * x_init *
                      (x_init * x_init + 3. * y_init * y_init) -
                      yp_init * y_init *
                      (y_init * y_init + 3. * x_init * x_init));
    }
}

///////////////////////////////////////////////////////////////////////////
// NAME
//   wedgerotate
//
// DESCRIPTION
//   Rotates coordinates by e for fringe fields at non-SBEND
//
// PARAMETERS
//   bunch  = reference to the macro-particle bunch
//   e = rotation angle
//   frinout = 0 before fringe, 1 after fringe
//
// RETURNS
//   Nothing
//
///////////////////////////////////////////////////////////////////////////

void wedgerotate(Bunch* bunch, double e, int frinout)
{
    double cs, sn;
    double xp_temp, p0_temp, p0;

    SyncPart* syncPart = bunch->getSyncPart();

    double dp_p_coeff = 1.0 / (syncPart->getMomentum() * syncPart->getBeta());
    double dp_p;

    cs = cos(e);
    sn = sin(e);

    //coordinate array [part. index][x,xp,y,yp,z,dE]
    double** arr = bunch->coordArr();

    for(int i = 0; i < bunch->getSize(); i++)
    {
        if(frinout == 0)
        {
            dp_p    = arr[i][5] * dp_p_coeff;
            xp_temp = arr[i][1];
            p0_temp = 1.0 + dp_p;

            arr[i][0] /=  cs;
            arr[i][1]  =  xp_temp * cs + p0_temp * sn;
            p0         = -xp_temp * sn + p0_temp * cs;
            dp_p       =  p0 - 1.0;
            arr[i][4]  =  (-arr[i][0] * sn + arr[i][4]) / cs;
        }
        else
        {
            dp_p = arr[i][5] * dp_p_coeff;
            p0   = 1.0 + dp_p;

            arr[i][4]  = arr[i][0] * sn + arr[i][4] * cs;
            arr[i][0] *= cs;
            xp_temp    = arr[i][1] * cs - p0 * sn;
            p0_temp    = arr[i][1] * sn + p0 * cs;
            arr[i][1]  = xp_temp;
            dp_p       = p0_temp - 1.0;
        }
        arr[i][5] = dp_p / dp_p_coeff;
    }
}

///////////////////////////////////////////////////////////////////////////
// NAME
//   wedgedrift
//
// DESCRIPTION
//   Drifts particles through wedge for non-SBEND
//
// PARAMETERS
//   bunch =  reference to the macro-particle bunch
//   e = wedge angle
//   inout = 0 for in, 1 for out
//
// RETURNS
//   Nothing
//
///////////////////////////////////////////////////////////////////////////

void wedgedrift(Bunch* bunch, double e, int inout)
{
    double ct, tn;
    double s;

    SyncPart* syncPart = bunch->getSyncPart();

    double dp_p_coeff = 1.0 / (syncPart->getMomentum() * syncPart->getBeta());
    double dp_p;

    ct = cos(e) / sin(e);

    //coordinate array [part. index][x,xp,y,yp,z,dE]
    double** arr = bunch->coordArr();

    for(int i = 0; i < bunch->getSize(); i++)
    {
        if(inout == 0)
        {
            dp_p = arr[i][5] * dp_p_coeff;
            tn   = arr[i][1] / (1.0 + dp_p);
            s    = arr[i][0] / (ct - tn);
        }
        else
        {
            s    = arr[i][0] / ct;
        }

        drifti(bunch, i, s);
    }
}

///////////////////////////////////////////////////////////////////////////
// NAME
//   wedgebend
//
// DESCRIPTION
//   Straight bends particles through wedge for non-SBEND
//
// PARAMETERS
//   bunch  = reference to the macro-particle bunch
//   e = wedge angle
//   inout = 0 for in, 1 for out
//   rho = radius of curvature
//   nsteps = number of integraton steps
//
// RETURNS
//   Nothing
//
///////////////////////////////////////////////////////////////////////////

void wedgebend(Bunch* bunch, double e, int inout, double rho, int nsteps)
{
    double ct, tn;
    double s, sm, sm2;
    int nst;

    SyncPart* syncPart = bunch->getSyncPart();

    double dp_p_coeff = 1.0 / (syncPart->getMomentum() * syncPart->getBeta());
    double dp_p;

    nst = nsteps / 2;
    if(nst < 1) nst = 1;
    ct = cos(e) / sin(e);

    //coordinate array [part. index][x,xp,y,yp,z,dE]
    double** arr = bunch->coordArr();

    for(int i = 0; i < bunch->getSize(); i++)
    {
        if(inout == 0)
        {
            s    = -arr[i][0] / ct;
        }
        else
        {
            dp_p =  arr[i][5] * dp_p_coeff;
            tn   =  arr[i][1] / (1.0 + dp_p);
            s    = -arr[i][0] / (ct + tn);
        }

        sm  = s / nst;
        sm2 = sm / 2.0;

        drifti(bunch, i, sm2);
        arr[i][1] -= sm / rho;
        for(int j  = 1; j < nst; j++)
        {
            drifti(bunch, i, sm);
            arr[i][1] -= sm / rho;
        }
        drifti(bunch, i, sm2);
    }
}

///////////////////////////////////////////////////////////////////////////
// NAME
//   bend1
//
// DESCRIPTION
//   Linear bend transport
//
// PARAMETERS
//   bunch  = reference to the macro-particle bunch
//   length = length of transport
//   th = bending angle
//
// RETURNS
//   Nothing
//
///////////////////////////////////////////////////////////////////////////

void bend1(Bunch* bunch, double length, double th)
{
    double x_init, xp_init;
    double cx, sx, rho;
    double m11, m12, m16;
    double m21, m22, m26;
    double m51, m52, m56;

    SyncPart* syncPart = bunch->getSyncPart();

    double v = OrbitConst::c * syncPart->getBeta();
    if(length > 0.)
    {
	   syncPart->setTime( syncPart->getTime() + length/v);
    }

    double betasq = syncPart->getBeta() * syncPart->getBeta();
    double dp_p_coeff = 1.0 / (syncPart->getMomentum() * syncPart->getBeta());
    double dp_p;

    rho = length / th;
    cx  = cos(th);
    sx  = sin(th);
    m11 = cx;
    m12 = rho * sx;
    m16 = rho * (1.0 - cx);
    m21 = -sx / rho;
    m22 = cx;
    m26 = sx;
    m51 = -sx;
    m52 = -rho * (1.0 - cx);
    m56 = -betasq * length + rho * sx;

    //coordinate array [part. index][x,xp,y,yp,z,dE]
    double** arr = bunch->coordArr();

    for(int i = 0; i < bunch->getSize(); i++)
    {
        dp_p   = arr[i][5] * dp_p_coeff;
        x_init = arr[i][0];
        xp_init = arr[i][1];

        arr[i][0]  = x_init * m11 + xp_init * m12 + dp_p * m16;
        arr[i][1]  = x_init * m21 + xp_init * m22 + dp_p * m26;
        arr[i][2] += length * arr[i][3];
        arr[i][4] += x_init * m51 + xp_init * m52 + dp_p * m56;
    }
}

///////////////////////////////////////////////////////////////////////////
// NAME
//   bend2
//
// DESCRIPTION
//   Kinetic bend transport (same as nonlinear quad transport - quad2)
//
// PARAMETERS
//   bunch  = reference to the macro-particle bunch
//   length = length of element (either full of half step)
//
// RETURNS
//   Nothing
//
///////////////////////////////////////////////////////////////////////////

void bend2(Bunch* bunch, double length)
{
    double KNL, phifac;

    SyncPart* syncPart = bunch->getSyncPart();

    double gamma2i = 1.0 / (syncPart->getGamma() * syncPart->getGamma());
    double dp_p_coeff = 1.0 / (syncPart->getMomentum() * syncPart->getBeta());
    double dp_p;

    //coordinate array [part. index][x,xp,y,yp,z,dE]
    double** arr = bunch->coordArr();

    for(int i = 0; i < bunch->getSize(); i++)
    {
        dp_p = arr[i][5] * dp_p_coeff;
        KNL = 1.0 / (1.0 + dp_p);

        arr[i][0] -= KNL * length * dp_p * arr[i][1];
        arr[i][2] -= KNL * length * dp_p * arr[i][3];
        phifac = (arr[i][1] * arr[i][1] + arr[i][3] * arr[i][3] +
                  dp_p * dp_p * gamma2i) / 2.0;
        phifac = (phifac * KNL + dp_p * dp_p * gamma2i) * KNL;
        arr[i][4] -= length * phifac;
    }
}

///////////////////////////////////////////////////////////////////////////
// NAME
//   bend3
//
// DESCRIPTION
//   Nonlinear curvature bend transport
//   depending on py and dE in Hamiltonian
//
// PARAMETERS
//   bunch =  reference to the macro-particle bunch
//   th = bending angle
//
// RETURNS
//   Nothing
//
///////////////////////////////////////////////////////////////////////////

void bend3(Bunch* bunch, double th)
{
    double KNL, phifac;

    SyncPart* syncPart = bunch->getSyncPart();

    double gamma2i = 1.0 / (syncPart->getGamma() * syncPart->getGamma());
    double dp_p_coeff = 1.0 / (syncPart->getMomentum() * syncPart->getBeta());
    double dp_p;

    //coordinate array [part. index][x,xp,y,yp,z,dE]
    double** arr = bunch->coordArr();

    for(int i = 0; i < bunch->getSize(); i++)
    {
        dp_p = arr[i][5] * dp_p_coeff;
        KNL  = 1.0 / (1.0 + dp_p);

        phifac     = (arr[i][3] * arr[i][3] + dp_p * dp_p * gamma2i) / 2.0;
        arr[i][1] -= phifac * KNL * th;
        arr[i][2] += KNL * arr[i][3] * arr[i][0] * th;
        phifac     = (phifac * KNL - dp_p * gamma2i) * KNL;
        arr[i][4] -= th * phifac * arr[i][0];
    }
}

///////////////////////////////////////////////////////////////////////////
// NAME
//   bend4
//
// DESCRIPTION
//   Nonlinear curvature bend transport
//   depending on px in Hamiltonian
//
// PARAMETERS
//   bunch =  reference to the macro-particle bunch
//   th = bending angle
//
// RETURNS
//   Nothing
//
///////////////////////////////////////////////////////////////////////////

void bend4(Bunch* bunch, double th)
{
    double KNL, phifac, xfac;

    SyncPart* syncPart = bunch->getSyncPart();

    double dp_p_coeff = 1.0 / (syncPart->getMomentum() * syncPart->getBeta());
    double dp_p;

    //coordinate array [part. index][x,xp,y,yp,z,dE]
    double** arr = bunch->coordArr();

    for(int i = 0; i < bunch->getSize(); i++)
    {
        dp_p = arr[i][5] * dp_p_coeff;
        KNL  = 1.0 / (1.0 + dp_p);

        xfac   = 1.0 + KNL * arr[i][1] * th / 2.0;
        phifac = KNL * KNL * arr[i][1] * arr[i][1] / 2.0;
        arr[i][0] *= xfac * xfac;
        arr[i][1] /= xfac;
        arr[i][4] -= th * phifac * arr[i][0];
    }
}

///////////////////////////////////////////////////////////////////////////
// NAME
//   bendfringeIN
//
// DESCRIPTION
//   Hard edge fringe field for a bend
//
// PARAMETERS
//   bunch =  reference to the macro-particle bunch
//   rho = radius of curvature for bending
//
// RETURNS
//   Nothing
//
///////////////////////////////////////////////////////////////////////////

void bendfringeIN(Bunch* bunch, double rho)
{
    SyncPart* syncPart = bunch->getSyncPart();

    double dp_p_coeff = 1.0 / (syncPart->getMomentum() * syncPart->getBeta());
    double dp_p, KNL;

    //coordinate array [part. index][x,xp,y,yp,z,dE]
    double** arr = bunch->coordArr();

    for(int i = 0; i < bunch->getSize(); i++)
    {
        dp_p    = arr[i][5] * dp_p_coeff;
        KNL  = 1.0 / (1.0 + dp_p);

        arr[i][0] += KNL * arr[i][2] * arr[i][2] / (2. * rho);
        arr[i][3] -= KNL * arr[i][1] * arr[i][2] / rho;
        arr[i][4] -= KNL * KNL * arr[i][1] * arr[i][2] * arr[i][2] / (2. * rho);
    }
}

///////////////////////////////////////////////////////////////////////////
// NAME
//   bendfringeOUT
//
// DESCRIPTION
//   Hard edge fringe field for a bend
//
// PARAMETERS
//   bunch =  reference to the macro-particle bunch
//   rho = radius of curvature for bending
//
// RETURNS
//   Nothing
//
///////////////////////////////////////////////////////////////////////////

void bendfringeOUT(Bunch* bunch, double rho)
{
    SyncPart* syncPart = bunch->getSyncPart();

    double dp_p_coeff = 1.0 / (syncPart->getMomentum() * syncPart->getBeta());
    double dp_p, KNL;

    //coordinate array [part. index][x,xp,y,yp,z,dE]
    double** arr = bunch->coordArr();

    for(int i = 0; i < bunch->getSize(); i++)
    {
        dp_p    = arr[i][5] * dp_p_coeff;
        KNL  = 1.0 / (1.0 + dp_p);

        arr[i][0] -= KNL * arr[i][2] * arr[i][2] / (2. * rho);
        arr[i][3] += KNL * arr[i][1] * arr[i][2] / rho;
        arr[i][4] += KNL * KNL * arr[i][1] * arr[i][2] * arr[i][2] / (2. * rho);
    }
}

///////////////////////////////////////////////////////////////////////////
// NAME
//   soln
//
// DESCRIPTION
//   Integration through a solenoid
//
// PARAMETERS
//   bunch  =  reference to the macro-particle bunch
//   length = integration length
//   B      = magnetic field (1/m)
//
// RETURNS
//   Nothing
//
///////////////////////////////////////////////////////////////////////////

void soln(Bunch* bunch, double length, double B)
{
    double KNL, phase, cs, sn;
    double cu, cpu, u_init, pu_init, u, pu, phifac;

    SyncPart* syncPart = bunch->getSyncPart();

    double v = OrbitConst::c * syncPart->getBeta();
    if(length > 0.)
    {
	   syncPart->setTime( syncPart->getTime() + length/v);
    }

    double gamma2i = 1.0 / (syncPart->getGamma() * syncPart->getGamma());
    double dp_p_coeff = 1.0 / (syncPart->getMomentum() * syncPart->getBeta());
    double dp_p;

    //coordinate array [part. index][x,xp,y,yp,z,dE]
    double** arr = bunch->coordArr();

    for(int i = 0; i < bunch->getSize(); i++)
    {
        dp_p = arr[i][5] * dp_p_coeff;
        KNL  = 1.0 / (1.0 + dp_p);

        cu      =  arr[i][2] / 2.     - arr[i][1] / B;
        cpu     =  arr[i][0] * B / 2. + arr[i][3];
        u_init  =  arr[i][2] / 2.     + arr[i][1] / B;
        pu_init = -arr[i][0] * B / 2. + arr[i][3];
        phase = KNL * B * length;
        cs = cos(phase);
        sn = sin(phase);

        u =   u_init * cs     + pu_init * sn / B;
        pu = -u_init * B * sn + pu_init * cs;

        arr[i][0] = (-pu + cpu) / B;
        arr[i][1] = 0.5 * (u - cu) * B;
        arr[i][2] = u + cu;
        arr[i][3] = 0.5 * (pu + cpu);

        phifac = (pu_init * pu_init +
                  B * B * u_init * u_init +
                  dp_p * dp_p * gamma2i
                 ) / 2.0;
        phifac = (phifac * KNL - dp_p * gamma2i) * KNL;
        arr[i][4] -= length * phifac;
    }
}

///////////////////////////////////////////////////////////////////////////
// NAME
//   wedgebendCF
//
// DESCRIPTION
//   Straight bends particles through wedge for Combined Function non-SBEND
//
// PARAMETERS
//   bunch =  reference to the macro-particle bunch
//   e = wedge angle
//   inout = 0 for in, 1 for out
//   rho = radius of curvature
//   vecnum = number of multipole terms
//   pole = multipolarities of multipole terms
//   kl = integrated strengths of multipole terms
//   skew = skewness  of multipole terms
//   nsteps = number of integraton steps
//
// RETURNS
//   Nothing
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
    double ct, tn;
    double s, sm, sm2, klint;
    int nst;

    SyncPart* syncPart = bunch->getSyncPart();

    double dp_p_coeff = 1.0 / (syncPart->getMomentum() * syncPart->getBeta());
    double dp_p;

    nst = nsteps / 2;
    if(nst < 1) nst = 1;
    ct = cos(e) / sin(e);

    //coordinate array [part. index][x,xp,y,yp,z,dE]
    double** arr = bunch->coordArr();

    for(int i = 0; i < bunch->getSize(); i++)
    {
        if(inout == 0)
        {
            s = -arr[i][0] / ct;
        }
        else
        {
            dp_p =  arr[i][5] * dp_p_coeff;
            tn   =  arr[i][1] / (1.0 + dp_p);
            s    = -arr[i][0] / (ct + tn);
        }

        sm = s / nst;
        sm2 = sm / 2.0;

        drifti(bunch, i, sm2);
        arr[i][1] -= sm / rho;
        for (int l = 0; l < vecnum; l++)
        {
            klint = kl[l] * sm;
            multpi(bunch, i, pole[l], klint, skew[l]);
        }
        for(int j = 1; j < nst; j++)
        {
            drifti(bunch, i, sm);
            arr[i][1] -= sm / rho;
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
//   RingRF
//
// DESCRIPTION
//   Ring type RF cavity. Transition time factor T(k) = const = T(k0).
//   No need for symplectic phase correction.
//
// PARAMETERS
//   bunch =  reference to the macro-particle bunch
//   harmonic_numb = harmonic number
//   voltage = voltage in Giga Volts
//   phase_s = synchronous phase in Rad
//
// RETURNS
//   Nothing
//
///////////////////////////////////////////////////////////////////////////

void RingRF(Bunch* bunch, double ring_length, int harmonic_numb,
            double voltage, double phase_s)
{
    double deltaV = 0.;
    double charge = bunch->getCharge();
    double coeff  = charge;

    double Factor = 0.;
    if(ring_length > 0.)
    {
        Factor = 2.0 * OrbitConst::PI/ring_length;
    }

    SyncPart* syncPart = bunch->getSyncPart();

    if(phase_s != 0.)
    {
        double kin_e = syncPart->getEnergy();
        kin_e += coeff * voltage * sin(phase_s);
        syncPart->setMomentum(syncPart->energyToMomentum(kin_e));
    }

    //coordinate array [part. index][x,xp,y,yp,z,dE]
    double** arr = bunch->coordArr();

    for(int i = 0; i < bunch->getSize(); i++)
    {
        deltaV = voltage * ( sin(harmonic_numb*Factor*arr[i][4] + phase_s));
        arr[i][5] += coeff * deltaV;
    }
}

}  //end of namespace teapot_base
