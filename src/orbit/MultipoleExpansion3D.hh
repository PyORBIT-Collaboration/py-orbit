//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    MultipoleExpansion3D.hh
//
// AUTHOR
//    M. Perkett
//
// CREATED
//    09/19/2006
//
// MODIFIED
// 
//
// DESCRIPTION
//    Subclass of BaseFieldSource for 3D multipole expansion
//
// INPUT FILE FORMAT
//    It may be useful to look at a 3D Multiple Expansion explicit form
//        (from JG Wang's paper in this case)
//     Basically, the file must be ordered with all the coefficients
//     delimited by some sort of spacing (shouldn't matter what type) The
//     order will be z value then all the m=0 coeffiecients (from l=0 to l=# of l) 
//     next all the coefficients for m=1 (from l=0 to l=# of l) and each l will
//     have 2 coefficients (the first for sine, the second for cosine)
//     It is important that everything in the table carry a value that 
//     corresponds to the necessary coefficient.  Failure to comply will
//     result in an incorrect calculation.  Step size should NOT matter.
//     It is also important to discuss the m=0 series... since m=0, there
//     will not be a sine component (sin(0)=0) and it should therefore be
//     omitted from the remainder of the list.
//     It is also assumed that although an even step size is not needed,
//     the numbers go from least to greatest
//
//     Note: # of m assumes that you are starting from 0 and going up to
//           value m  (same with # of l)
//
//    line 1:  <# of m>    <# of l>    <# coefficients(file length-1)>
//    line 2:  z    .....
//      .
//      .
//      .
//    line n:
//
///////////////////////////////////////////////////////////////////////////

#ifndef _MULTIPOLE_EXPANSION_3D
#define _MULTIPOLE_EXPANSION_3D

#include<string>
#include "BaseFieldSource.hh"
using namespace std;


class MultipoleExpansion3D : public BaseFieldSource
{ 
   public:
   MultipoleExpansion3D();
   MultipoleExpansion3D(string inputFile);
   ~MultipoleExpansion3D();
   void getElectricField(double t, double x, double y, double z,
                         double &f_x, double &f_y, double &f_z);
   void getMagneticField(double t, double x, double y, double z,
                         double &f_x, double &f_y, double &f_z);
   double getCoefficient(const double &z, const int &m, const int &l,
        int sORc) const;
   int factorial(const int &number) const;   // used by callBesselFunction
        
   private:
   // helper functions


   
   // coefficients
   double ****coefficients;
   int m,kappa,length;   // kappa is distinct from l right now -> need better
                        // notation in future so that it is more readable
   double *zValues;      // array to store z values in
  
   // TEMPORARY CODE USED FOR TESTING PURPOSES
//   double tempArray[2401][3];
};

#endif
///////////////////////////////////////////////////////////////////////////
//
// END OF FILE
//
///////////////////////////////////////////////////////////////////////////
