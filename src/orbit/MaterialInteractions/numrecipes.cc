#include <cstdlib>
#include <iostream>
#include <cmath>
#include <cfloat>
#include "numrecipes.hh"

namespace OrbitUtils{
	
	
	float fstep(float s, float r_o, float pr_o, float theta)
	{
		float a, b, c, d;
		float f_x;
		float n_o = 3.;
		
		a=pow((n_o*theta)/3., 2.0); 
		b=pr_o * pr_o; 
		c=2*r_o*pr_o; 
		d=r_o * r_o;
		
		f_x=a*s*s*s - b*s*s - c*s - d;
		//std::cerr<<"fx, a, s, b, c, d"<<f_x<<" "<<a<<" "<<s<<" "<<b<<" "<<c<<" "<<d<<"\n";
		return f_x;
	}
	
	
	float rfunc(float x, float p, float fac1)
	{
		return sin(x)*exp(-fac1*pow(2*p*sin(x/2), 2))/pow(sin(x/2),4);
	}
	
	
	int zbrak(float (*fx)(float, float, float, float), float x1, float x2, 
			  int n, float xb1[], float xb2[], int &nb, float param1, 
			  float param2, float param3)
	{
		int nbb, i;
		float x, fp, fc, dx;
		
		nbb=0;
		dx=(x2-x1)/n;
		fp=(*fx)(x=x1, param1, param2, param3);
		for (i=1; i<=n; i++) {
			fc=(*fx)(x += dx, param1, param2, param3);
			if(fc*fp <= 0.0) {
				xb1[++nbb]=x-dx;
				xb2[nbb]=x;
				if(nb == nbb) return 0;
			}
			fp=fc;
		}
		nb=nbb;
	}
	
	#undef JMAX
	#define JMAX 40
	
	float rtbis(float (*func)(float, float, float, float), float x1, float x2, 
				float xacc, float param1, float param2, float param3)
	{
		int j;
		float dx, f, fmid, xmid, rtb;
		
		f=(*func)(x1, param1, param2, param3);
		fmid=(*func)(x2, param1, param2, param3);
		if(f*fmid >= 0.0) std::cerr<<"Root must be bracketed for bisection in rtbis\n";
		rtb = f < 0.0 ? (dx=x2-x1, x1) : (dx = x1-x2, x2);
		for (j=1; j <= JMAX; j++) {
			fmid = (*func)(xmid=rtb+(dx *= 0.5), param1, param2, param3);
			if (fmid <= 0.0) rtb = xmid; 
			if (fabs(dx) < xacc || fmid == 0.0) return rtb;
		}
		std::cerr<<"Too many bisections in rtbis \n";
		return 0.0;
	}
	
	
	
	float bessj0(float x)
	{
		float ax,z;
		double xx,y,ans,ans1,ans2;
		
		if ((ax=fabs(x)) < 8.0) {
			y=x*x;
			ans1=57568490574.0+y*(-13362590354.0+y*(651619640.7
													+y*(-11214424.18+y*(77392.33017+y*(-184.9052456)))));
			ans2=57568490411.0+y*(1029532985.0+y*(9494680.718
												  +y*(59272.64853+y*(267.8532712+y*1.0))));
			ans=ans1/ans2;
		} else {
			z=8.0/ax;
			y=z*z;
			xx=ax-0.785398164;
			ans1=1.0+y*(-0.1098628627e-2+y*(0.2734510407e-4
											+y*(-0.2073370639e-5+y*0.2093887211e-6)));
			ans2=-0.1562499995e-1+y*(0.1430488765e-3
									   +y*(-0.6911147651e-5+y*(0.7621095161e-6
															   -y*0.934935152e-7)));
			ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
		}
		return ans;
	}
	
	
	
	float bessj1(float x)
	{
		float ax,z;
		double xx,y,ans,ans1,ans2;
		
		if ((ax=fabs(x)) < 8.0) {
			y=x*x;
			ans1=x*(72362614232.0+y*(-7895059235.0+y*(242396853.1
													  +y*(-2972611.439+y*(15704.48260+y*(-30.16036606))))));
			ans2=144725228442.0+y*(2300535178.0+y*(18583304.74
												   +y*(99447.43394+y*(376.9991397+y*1.0))));
			ans=ans1/ans2;
		} else {
			z=8.0/ax;
			y=z*z;
			xx=ax-2.356194491;
			ans1=1.0+y*(0.183105e-2+y*(-0.3516396496e-4
									   +y*(0.2457520174e-5+y*(-0.240337019e-6))));
			ans2=0.04687499995+y*(-0.2002690873e-3
								  +y*(0.8449199096e-5+y*(-0.88228987e-6
														 +y*0.105787412e-6)));
			ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
			if (x < 0.0) ans = -ans;
		}
		return ans;
	}
	
#undef EPS
#undef JMAX
#define EPS 1.0e-6
#define JMAX 20
	
	float qsimp(float (*func)(float, float, float), float a, float b, 
				float p, float fac1)
	{
		float trapzd(float (*func)(float, float, float), float a, 
					 float b, int n, float p, float fac1);
		void nrerror(char error_text[]);
		int j;
		float s,st,ost,os;
		
		ost = os = -1.0e30;
		for (j=1;j<=JMAX;j++) {
			st=trapzd(func,a,b,j,p,fac1);
			s=(4.0*st-ost)/3.0;
			if (j > 5)
				if (fabs(s-os) < EPS*fabs(os) ||
					(s == 0.0 && os == 0.0)) return s;
			os=s;
			ost=st;
		}
		//nrerror("Too many steps in routine qsimp");
		return 0.0;
	}
#undef EPS
#undef JMAX
	
	
	/* note #undef's at end of file */
#define FUNC(x, p, fac1) ((*rfunc)(x, p, fac1))
	
	float trapzd(float (*rfunc)(float, float, float), float a, 
				 float b, int n, float p, float fac1)
	{
		float x,tnm,sum,del;
		static float s;
		int it,j;
		
		if (n == 1) {
			return (s=0.5*(b-a)*(FUNC(a, p, fac1)+FUNC(b, p, fac1)));
		} else {
			for (it=1,j=1;j<n-1;j++) it <<= 1;
			tnm=it;
			del=(b-a)/tnm;
			x=a+0.5*del;
			for (sum=0.0,j=1;j<=it;j++,x+=del) sum += FUNC(x, p, fac1);
			s=0.5*(s+(b-a)*sum/tnm);
			return s;
		}
	}
#undef FUNC
			

	
		
}

	
	
	
