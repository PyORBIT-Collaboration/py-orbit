//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    bessj0.hh
//
//    This subroutine calculates the First and Second Kind Bessel Function of
//    order 0,1,n for any real number x. The polynomial approximation by
//    series of Chebyshev polynomials is used for 0<x<8 and 0<8/x<1.
//    REFERENCES:
//    M.ABRAMOWITZ,I.A.STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS, 1965.
//    C.W.CLENSHAW, NATIONAL PHYSICAL LABORATORY MATHEMATICAL TABLES,
//    VOL.5, 1962.
//
///////////////////////////////////////////////////////////////////////////

#include "bessel.hh"

namespace OrbitUtils{
	
	double bessj0 (double x) {
		
		const double
		P1=1.0, P2=-0.1098628627E-2, P3=0.2734510407E-4, P4=-0.2073370639E-5, P5= 0.2093887211E-6,
		Q1=-0.1562499995E-1, Q2= 0.1430488765E-3, Q3=-0.6911147651E-5, Q4= 0.7621095161E-6, Q5=-0.9349451520E-7,
		R1= 57568490574.0, R2=-13362590354.0, R3=651619640.7, R4=-11214424.18, R5= 77392.33017, R6=-184.9052456,
		S1= 57568490411.0, S2=1029532985.0, S3=9494680.718, S4= 59272.64853, S5=267.8532712, S6=1.0;
		
		double AX,FR,FS,Z,FP,FQ,XX,Y, TMP;
		
		if (x==0.0) return 1.0;
		
		AX = fabs(x);
		if (AX < 8.0) {
			Y = x*x;
			FR = R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6))));
			FS = S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6))));
			TMP = FR/FS;
		}
		else {
			Z = 8./AX;
			Y = Z*Z;
			XX = AX-0.785398164;
			FP = P1+Y*(P2+Y*(P3+Y*(P4+Y*P5)));
			FQ = Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)));
			TMP = sqrt(0.636619772/AX)*(FP*cos(XX)-Z*FQ*sin(XX));
		}
		return TMP;
	}
	
	double BSign(double X, double Y) {
		if (Y<0.0) return (-fabs(X));
		else return (fabs(X));
	}
	
	double bessj1 (double X) {
		
    const double  
		P1=1.0, P2=0.183105E-2, P3=-0.3516396496E-4, P4=0.2457520174E-5,
		P5=-0.240337019E-6,  P6=0.636619772,
		Q1= 0.04687499995, Q2=-0.2002690873E-3, Q3=0.8449199096E-5,
		Q4=-0.88228987E-6, Q5= 0.105787412E-6,
		R1= 72362614232.0, R2=-7895059235.0, R3=242396853.1,
		R4=-2972611.439,   R5=15704.48260,  R6=-30.16036606,
		S1=144725228442.0, S2=2300535178.0, S3=18583304.74,
		S4=99447.43394,    S5=376.9991397,  S6=1.0;
		
	  double AX,FR,FS,Y,Z,FP,FQ,XX, TMP;
		
		AX = fabs(X);
		if (AX < 8.0) {
			Y = X*X;
			FR = R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6))));
			FS = S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6))));
			TMP = X*(FR/FS);
		}
		else {
			Z = 8.0/AX;
			Y = Z*Z;
			XX = AX-2.35619491;
			FP = P1+Y*(P2+Y*(P3+Y*(P4+Y*P5)));
			FQ = Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)));
			TMP = sqrt(P6/AX)*(cos(XX)*FP-Z*sin(XX)*FQ)*BSign(S6,X);
		}
	  return TMP;
	}
	
	double bessj (int N, double X) {
		
		const int IACC = 40; 
		const double BIGNO = 1e10,  BIGNI = 1e-10;
    
	  double TOX,BJM,BJ,BJP,SUM,TMP;
		int J, JSUM, M;
		
		if (N == 0) return bessj0(X);
		if (N == 1) return bessj1(X);
		if (X == 0.0) return 0.0;
		
		TOX = 2.0/X;
		if (X > 1.0*N) {
			BJM = bessj0(X);
			BJ  = bessj1(X);
			for (J=1; J<N; J++) {
				BJP = J*TOX*BJ-BJM;
				BJM = BJ;
				BJ  = BJP;
			}
			return BJ;
		}
		else {
			M = (int) (2*((N+floor(sqrt(1.0*(IACC*N))))/2));
			TMP = 0.0;
			JSUM = 0;
			SUM = 0.0;
			BJP = 0.0;
			BJ  = 1.0;
			for (J=M; J>0; J--) {
				BJM = J*TOX*BJ-BJP;
				BJP = BJ;
				BJ  = BJM;
				if (fabs(BJ) > BIGNO) {
					BJ  = BJ*BIGNI;
					BJP = BJP*BIGNI;
					TMP = TMP*BIGNI;
					SUM = SUM*BIGNI;
				}
				if (JSUM != 0)  SUM += BJ;
				JSUM = 1-JSUM;
				if (J == N)  TMP = BJP;
			}
			SUM = 2.0*SUM-BJ;
			return (TMP/SUM);
		}
	}
	
  double bessi0(double X) {
      double Y,AX,BX;
      const double P1=1.0, P2=3.5156229, P3=3.0899424, P4=1.2067429, P5=0.2659732, 
			P6=0.360768e-1, P7=0.45813e-2, Q1=0.39894228, Q2=0.1328592e-1, Q3=0.225319e-2,
      Q4=-0.157565e-2, Q5=0.916281e-2, Q6=-0.2057706e-1, Q7=0.2635537e-1, 
			Q8=-0.1647633e-1, Q9=0.392377e-2;
			
      if (fabs(X) < 3.75) {
        Y=(X/3.75)*(X/3.75);
        return (P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))));
      }
      else {
        AX=fabs(X);
        Y=3.75/AX;
        BX=exp(AX)/sqrt(AX);
        AX=Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))));
        return (AX*BX);
      }
  }

  double bessi1(double X) {
      double Y,AX,BX;
			
      const double P1=0.5, P2=0.87890594, P3=0.51498869, P4=0.15084934,
      P5=0.2658733e-1, P6=0.301532e-2, P7=0.32411e-3,
      Q1=0.39894228, Q2=-0.3988024e-1, Q3=-0.362018e-2,
      Q4=0.163801e-2, Q5=-0.1031555e-1, Q6=0.2282967e-1,
      Q7=-0.2895312e-1, Q8=0.1787654e-1, Q9=-0.420059e-2;
			
      if (fabs(X) < 3.75) {
        Y=(X/3.75)*(X/3.75);
        return(X*(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7)))))));
      }
      else {
        AX=fabs(X);
        Y=3.75/AX;
        BX=exp(AX)/sqrt(AX);
        AX=Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))));
        return (AX*BX);
      }
  }

  double bessi(int N, double X) {
		int IACC = 40; 
	  double BIGNO = 1e10, BIGNI = 1e-10;
		double TOX, BIM, BI, BIP, BSI;
		int J, M;
		
		if (N==0)  return (bessi0(X));
		if (N==1)  return (bessi1(X));
		if (X==0.0) return 0.0;
		
		TOX = 2.0/X;
		BIP = 0.0;
		BI  = 1.0;
		BSI = 0.0;
		M = (int) (2*((N+floor(sqrt(IACC*N)))));
		for (J = M; J>0; J--) {
			BIM = BIP+J*TOX*BI;
			BIP = BI;
			BI  = BIM;
			if (fabs(BI) > BIGNO) {
				BI  = BI*BIGNI;
				BIP = BIP*BIGNI;
				BSI = BSI*BIGNI;
			}
			if (J==N)  BSI = BIP;
		}
		return (BSI*bessi0(X)/BI);
  }
}
