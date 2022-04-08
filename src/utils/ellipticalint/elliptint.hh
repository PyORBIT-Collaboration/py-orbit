//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//   elliptint.hh
//
// AUTHOR
//    A. Shishlo
//
// CREATED
//    04/05/2022
//
// DESCRIPTION
//    Special Elliptical Integrals.
//    Numerical Recipes.
//    6.12 Elliptic Integrals and Jacobian Elliptic Functions
///////////////////////////////////////////////////////////////////////////

#ifndef ELLIPTICAL_INTEGRALS_H
#define ELLIPTICAL_INTEGRALS_H

/*
Chapter 6. Special Functions elliptint.hh
**/

#include <limits>
#include <cstdlib>

namespace EllipticalIntegrals{
	
	/*
	Computes Carlson’s degenerate elliptic integral, RC .x; y/. x 
	must be nonnegative and y must be nonzero. If y < 0, 
	the Cauchy principal value is returned.
	**/
	double rc(const double x, const double y) {
		static const double ERRTOL=0.0012, THIRD=1.0/3.0, C1=0.3, C2=1.0/7.0,
		C3=0.375, C4=9.0/22.0;
		static const double TINY=5.0*std::numeric_limits<double>::min(),
		BIG=0.2*std::numeric_limits<double>::max(), COMP1=2.236/sqrt(TINY),
		COMP2=(TINY*BIG)*(TINY*BIG)/25.0;
		double alamb,ave,s,w,xt,yt;
		if (x < 0.0 || y == 0.0 || (x+fabs(y)) < TINY || (x+fabs(y)) > BIG ||
			(y<-COMP1 && x > 0.0 && x < COMP2)) throw("invalid arguments in rc");
			if (y > 0.0) {
				xt=x;
				yt=y;
				w=1.0;
			} else {
				xt=x-y;
				yt= -y;
				w=sqrt(x)/sqrt(xt);
			}
			do {
				alamb=2.0*sqrt(xt)*sqrt(yt)+yt;
				xt=0.25*(xt+alamb);
				yt=0.25*(yt+alamb);
				ave=THIRD*(xt+yt+yt);
				s=(yt-ave)/ave;
			} while (fabs(s) > ERRTOL);
			return w*(1.0+s*s*(C1+s*(C2+s*(C3+s*C4))))/sqrt(ave);
	}
	
	/*
	Computes Carlson’s elliptic integral of the first kind, RF .x; y; z/. x, y, and z 
	must be non-negative, and at most one can be zero.
	**/
	double rf(const double x, const double y, const double z) {
		static const double ERRTOL=0.0025, THIRD=1.0/3.0,C1=1.0/24.0, C2=0.1, C3=3.0/44.0, C4=1.0/14.0;
		static const double TINY=5.0*std::numeric_limits<double>::min(),
		BIG=0.2*std::numeric_limits<double>::max();
		double alamb,ave,delx,dely,delz,e2,e3,sqrtx,sqrty,sqrtz,xt,yt,zt;
		if (std::min(std::min(x,y),z) < 0.0 || std::min(std::min(x+y,x+z),y+z) < TINY ||
			std::max(std::max(x,y),z) > BIG) throw("invalid arguments in rf");
			xt=x;
			yt=y;
			zt=z;
			do {
				sqrtx=sqrt(xt);
				sqrty=sqrt(yt);
				sqrtz=sqrt(zt);
				alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz;
				xt=0.25*(xt+alamb);
				yt=0.25*(yt+alamb);
				zt=0.25*(zt+alamb);
				ave=THIRD*(xt+yt+zt);
				delx=(ave-xt)/ave;
				dely=(ave-yt)/ave;
				delz=(ave-zt)/ave;
			} while (std::max(std::max(fabs(delx),fabs(dely)),fabs(delz)) > ERRTOL);
			e2=delx*dely-delz*delz;
			e3=delx*dely*delz;
			return (1.0+(C1*e2-C2-C3*e3)*e2+C4*e3)/sqrt(ave);
	}
	
	/*
	Computes Carlson’s elliptic integral of the second kind, RD .x; y; z/. x and y 
	must be nonnegative, and at most one can be zero. z must be positive.	
	**/
	
	double rd(const double x, const double y, const double z) {
		static const double ERRTOL=0.0015, C1=3.0/14.0, C2=1.0/6.0, C3=9.0/22.0,
		C4=3.0/26.0, C5=0.25*C3, C6=1.5*C4;
		static const double TINY=2.0*pow(std::numeric_limits<double>::max(),-2./3.),
		BIG=0.1*ERRTOL*pow(std::numeric_limits<double>::min(),-2./3.);
		double alamb,ave,delx,dely,delz,ea,eb,ec,ed,ee,fac,sqrtx,sqrty,
		sqrtz,sum,xt,yt,zt;
		if (std::min(x,y) < 0.0 || std::min(x+y,z) < TINY || std::max(std::max(x,y),z) > BIG)
			throw("invalid arguments in rd");
		xt=x;
		yt=y;
		zt=z;
		sum=0.0;
		fac=1.0;
		do {
			sqrtx=sqrt(xt);
			sqrty=sqrt(yt);
			sqrtz=sqrt(zt);
			alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz;
			sum += fac/(sqrtz*(zt+alamb));
			fac=0.25*fac;
			xt=0.25*(xt+alamb);
			yt=0.25*(yt+alamb);
			zt=0.25*(zt+alamb);
			ave=0.2*(xt+yt+3.0*zt);
			delx=(ave-xt)/ave;
			dely=(ave-yt)/ave;
			delz=(ave-zt)/ave;
		} while (std::max(std::max(abs(delx),abs(dely)),abs(delz)) > ERRTOL);
		ea=delx*dely;
		eb=delz*delz;
		ec=ea-eb;
		ed=ea-6.0*eb;
		ee=ed+ec+ec;
		return 3.0*sum+fac*(1.0+ed*(-C1+C5*ed-C6*delz*ee)
			+delz*(C2*ee+delz*(-C3*ec+delz*C4*ea)))/(ave*sqrt(ave));
	}
	
	/*
	Computes Carlson’s elliptic integral of the third kind, RJ .x; y; z; p/. x, y, and z 
	must be nonnegative, and at most one can be zero. p must be nonzero. If p < 0, 
	the Cauchy principal value is returned.
	**/
	double rj(const double x, const double y, const double z, const double p) {
		static const double ERRTOL=0.0015, C1=3.0/14.0, C2=1.0/3.0, C3=3.0/22.0,
		C4=3.0/26.0, C5=0.75*C3, C6=1.5*C4, C7=0.5*C2, C8=C3+C3;
		static const double TINY=pow(5.0*std::numeric_limits<double>::min(),1./3.),
		BIG=0.3*pow(0.2*std::numeric_limits<double>::max(),1./3.);
		double a,alamb,alpha,ans,ave,b,beta,delp,delx,dely,delz,ea,eb,ec,ed,ee,
		fac,pt,rcx,rho,sqrtx,sqrty,sqrtz,sum,tau,xt,yt,zt;
		if (std::min(std::min(x,y),z) < 0.0 || std::min(std::min(x+y,x+z),std::min(y+z,fabs(p))) < TINY
			|| std::max(std::max(x,y),std::max(z,fabs(p))) > BIG) throw("invalid arguments in rj");
			sum=0.0;
			fac=1.0;
			if (p > 0.0) {
				xt=x;
				yt=y;
				zt=z;
				pt=p;
			} else {
				xt=std::min(std::min(x,y),z);
				zt=std::max(std::max(x,y),z);
				yt=x+y+z-xt-zt;
				a=1.0/(yt-p);
				b=a*(zt-yt)*(yt-xt);
				pt=yt+b;
				rho=xt*zt/yt;
				tau=p*pt/yt;
				rcx=rc(rho,tau);
			}
			do {
				sqrtx=sqrt(xt);
				sqrty=sqrt(yt);
				sqrtz=sqrt(zt);
				alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz;
				alpha=pow((pt*(sqrtx+sqrty+sqrtz)+sqrtx*sqrty*sqrtz),2);
				beta=pt*pow((pt+alamb),2);
				sum += fac*rc(alpha,beta);
				fac=0.25*fac;
				xt=0.25*(xt+alamb);
				yt=0.25*(yt+alamb);
				zt=0.25*(zt+alamb);
				pt=0.25*(pt+alamb);
				ave=0.2*(xt+yt+zt+pt+pt);
				delx=(ave-xt)/ave;
				dely=(ave-yt)/ave;
				delz=(ave-zt)/ave;
				delp=(ave-pt)/ave;
			} while (std::max(std::max(fabs(delx),fabs(dely)),
				std::max(fabs(delz),fabs(delp))) > ERRTOL);
			ea=delx*(dely+delz)+dely*delz;
			eb=delx*dely*delz;
			ec=delp*delp;
			ed=ea-3.0*ec;
			ee=eb+2.0*delp*(ea-ec);
			ans=3.0*sum+fac*(1.0+ed*(-C1+C5*ed-C6*ee)+eb*(C7+delp*(-C8+delp*C4))
				+delp*ea*(C2-delp*C3)-C2*delp*ec)/(ave*sqrt(ave));
				if (p <= 0.0) ans=a*(b*ans+3.0*(rcx-rf(xt,yt,zt)));
				return ans;
	}
	
	/*
	Legendre elliptic integral of the first kind F(phi,ak),
	evaluated using Carlson’s function RF . 
	The argument ranges are 0<= phi <= pi/2, 0 <= k*sin(phi) <= 1.
	**/
	double ellf(const double phi, const double ak) {
		double s=sin(phi);
		return s*rf(pow(cos(phi),2),(1.0-s*ak)*(1.0+s*ak),1.0);
	}

	/*
	Legendre elliptic integral of the second kind E(phi,ak),
	evaluated using Carlson’s function RF and RD. 
	The argument ranges are 0<= phi <= pi/2, 0 <= k*sin(phi) <= 1.	
	**/
	double elle(const double phi, const double ak) {
		double cc,q,s;
		s=sin(phi);
		cc=pow(cos(phi),2);
		q=(1.0-s*ak)*(1.0+s*ak);
		return s*(rf(cc,q,1.0)-(pow(s*ak,2))*rd(cc,q,1.0)/3.0);
	}
		
	/*
	Legendre elliptic integral of the third kind P(phi,ak),
	evaluated using Carlson’s function RJ and RD.
	(Note that the sign convention on n is opposite that of Abramowitz and Stegun.)
	The argument ranges are 0<= phi <= pi/2, 0 <= k*sin(phi) <= 1.
	**/
	double ellpi(const double phi, const double en, const double ak) {
		double cc,enss,q,s;
		s=sin(phi);
		enss=en*s*s;
		cc=pow(cos(phi),2);
		q=(1.0-s*ak)*(1.0+s*ak);
		return s*(rf(cc,q,1.0)-enss*rj(cc,q,1.0,1.0+enss)/3.0);
	}
};

#endif

