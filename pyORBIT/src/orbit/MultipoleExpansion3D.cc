//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    MultipoleExpansion3D.cc
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
//
///////////////////////////////////////////////////////////////////////////


#include<iostream>
#include<math.h>
#include<fstream>
#include<iomanip>
#include "MultipoleExpansion3D.hh"
using namespace std;

//*********************************************************************************
// Function: MultipoleExpansion3D()
// Purpose: default constructor
// Precondition: none
// Postcondition: everything is setup so that any calculations that occur should
//    return zero in the getMagneticField()
//*********************************************************************************
MultipoleExpansion3D:: MultipoleExpansion3D()
{
   m = 3;
   length = 3;
   kappa = 3;
   
   zValues = new double[3];
   zValues[0] = -1000000.0;
   zValues[1] = 0.0;
   zValues[2] = 1000000.0;
   
   // since there are so few zValues... it is most convenient at this time to
   //    set up a larger number of values in the coefficients array
   coefficients = new double***[length];
   for (int i = 0; i < length; i++)
   {
      coefficients[i] = new double**[m + 1];
      for (int j = 0; j <= m; j++)
      {
         coefficients[i][j] = new double*[kappa + 1];
         for (int k = 0; k <= kappa; k++)
         {
            // this conditional ensures that there is no memory leak when the
            // ~MultipoleExpansion3D is called
            if ( !(j == 0 && k == 0))
            {
               coefficients[i][j][k] = new double[2];
               coefficients[i][j][k][0] = 0.0;
               coefficients[i][j][k][1] = 0.0;
            }
         }
      }
   }
}


//*********************************************************************************
// Function: MultipoleExpansion3D(string inputFileName)
// Purpose: constructor
// Precondition: the inputFileName1 and inputFileName2 must exist
// Postcondition: coeffients are all read into proper arrays
//*********************************************************************************
MultipoleExpansion3D:: MultipoleExpansion3D(string inputFileName)
{
   ifstream input(inputFileName.c_str());
   
   // if the precondition is violated, exit
   if (!input.is_open())
   {
       cout << "MultipoleExpansion3D:: readInCoefficients : " << endl;
       cout << "   " << inputFileName << "  could not be opened.  Error!  Exiting..."
            << endl;
       exit(1);
   }
        
   //!!! later add a test to ensure correct values were read in and that file doesn't terminate prematurely

   // read in stuff at top of file
   input >> m >> kappa >> length;
// DEBUGGING   cout << " m = " << m << "  kappa = " << kappa << "  length = " << length << endl;
   
   
   // dynamically allocated memory for array coefficients
   //    1st index: line number (starting at zero)
   //    2nd index: m
   //    3rd index: l
   //    4th index: 0 for sin and 1 for cosine
   coefficients = new double***[length];
   zValues = new double[length];
   
   // read everything else in assuming correct file format outlined in header
   for(int i = 0; i < length; i++)
   {
         input >> zValues[i];
         
// DEBUGGING         cout << zValues[i];
         
         coefficients[i] = new double**[m + 1];

      for(int j = 0; j <= m; j++)
      {
         coefficients[i][j] = new double*[kappa + 1];

         for(int k = 0; k <= kappa; k++)
         {
             // this ensures that there is no wasted memory for both sine and
             //   cosine coefficients for the m = 0 (only needs the cosine)
             // Note: this means that the cosine values will be stored at index
             //   0 instead of 1 like they usually are!
             if(j == 0)
             {
                if (k == 0)
                {
                   // don't read in m=0 and l=0
                }
                else
                {
                   coefficients[i][j][k] = new double[1];
                   input >> coefficients[i][j][k][0];
// DEBUGGING                   cout << setw(15) << coefficients[i][j][k][0];
                }
             }
             else
             {
                coefficients[i][j][k] = new double[2];
                input >> coefficients[i][j][k][0] >> coefficients[i][j][k][1];
// DEBUGGING                cout << setw(15) << coefficients[i][j][k][0] << setw(15)
// DEBUGGING                     << coefficients[i][j][k][1];
             }

         }
      }
// DEBUGGING         cout << endl;
      
   }
   
   // close files opened at the beginning
   input.close();
   
   //readInCoefficients(inputFileName1, inputFileName2);

   // TEMPORARY CODE FOR TESTING PURPOSES
//   ifstream temp("Chicane_other_4c.txt");
//   for (int i = 0; i < 2401; i++)
//   {
//       temp >> tempArray[i][0] >> tempArray[i][1] >> tempArray[i][2];
//   }
   

}


//*********************************************************************************
// Function: ~MultipoleExpansion3D()
// Purpose: destructor
// Precondition: none
// Postcondition: deallocates dynamically allocated memory (if any)
//*********************************************************************************
MultipoleExpansion3D:: ~MultipoleExpansion3D()
{
           delete [] zValues;
   // iterate through deleting all the dynamically allocated memory from the
   //   constructor
   for(int i = 0; i < length; i++)
   {
      for(int j = 0; j <= m; j++)
      {
         for(int k = 0; k <= kappa; k++)
         {
            if (!(j == 0 && k == 0))
            {
               delete [] coefficients[i][j][k];
            }

         }
         delete [] coefficients[i][j];
      }
      delete [] coefficients[i];
      
   }
  // cout << "here" << endl;   
   delete [] coefficients;
   
   
}


//*********************************************************************************
// Function: void getElectricField(double t, double x, double y, double z,
//                                 double &f_x, double &f_y, double &f_z)
// Purpose: allows user to get electric field componenets in the x,y, and z dir.
// Precondition: a valid time and position is specified
// Postcondition: x,y,z electric field components are stored in f_x,f_y,f_z
//*********************************************************************************
void MultipoleExpansion3D:: getElectricField(double t, double x, double y, double z,
   double &f_x, double &f_y, double &f_z)
{
   f_x = f_y = f_z = 0.0;
}


//*********************************************************************************
// Function: void getMagneticField(double t, double x, double y, double z,
//                            double &f_x, double &f_y, double &f_z)
// Purpose: allows user to return the magnetic field components in the x, y, and z
//          directions (f_x, f_y, f_z)
// Precondition: a valid time and position is specified
// Postcondition: x, y, and z magnetic field components are stored in f_x, f_y, f_z
//*********************************************************************************
void MultipoleExpansion3D:: getMagneticField(double t, double x, double y, double z,
   double &f_x, double &f_y, double &B_z)
{
 /*       
   // if out of bounds, exit
   if (z < -200.0 || z > 400.0)
   {
      cerr << "MultipoleExpansion3D::getMagneticField() : z < -200.0 or z > 400.0"
           << "  Error!  Exiting..." << endl;
      exit(1);
   }
*/
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ OLD METHOD FOR DEBUGGING @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


/*
            double r = sqrt( pow(x,2) + pow(y,2) );
            double phi = atan2( y,x );

            double B_r = -0.50*getCoefficient(z,0,2,0)*r //+ (1.0/16.0)*(~~~~~)*pow(r,3) 
               //- (1.0/384.0)*C_0_6[2*k]*pow(r,5)
               + (getCoefficient(z,1,0,0) - (3.0/8.0)*getCoefficient(z,1,2,0)*pow(r,2))*sin(phi)
               + (getCoefficient(z,1,0,1) - (3.0/8.0)*getCoefficient(z,1,2,1)*pow(r,2))*cos(phi)
               + (2.0*getCoefficient(z,2,0,0)*r - (1.0/3.0)*getCoefficient(z,2,2,0)*pow(r,3))*sin(2.0*phi)
               + (2.0*getCoefficient(z,2,0,1)*r - (1.0/3.0)*getCoefficient(z,2,2,1)*pow(r,3))*cos(2.0*phi)
               + (3.0*getCoefficient(z,3,0,0)*pow(r,2) - 5.0/16.0*getCoefficient(z,3,2,0)*pow(r,4))*sin(3.0*phi)
               + (3.0*getCoefficient(z,3,0,1)*pow(r,2) - 5.0/16.0*getCoefficient(z,3,2,1)*pow(r,4))*cos(3.0*phi);
 
            double B_phi = (getCoefficient(z,1,0,0) - (1.0/8.0)*getCoefficient(z,1,2,0)*pow(r,2))*cos(phi)
                 - (getCoefficient(z,1,0,1) - (1.0/8.0)*getCoefficient(z,1,2,1)*pow(r,2))*sin(phi)
                 + (2.0*getCoefficient(z,2,0,0)*r - (1.0/6.0)*getCoefficient(z,2,2,0)*pow(r,3))*cos(2.0*phi)
                 - (2.0*getCoefficient(z,2,0,1)*r - (1.0/6.0)*getCoefficient(z,2,2,1)*pow(r,3))*sin(2.0*phi)
                 + (3.0*getCoefficient(z,3,0,0)*pow(r,2) - (3.0/16.0)*getCoefficient(z,3,2,0)*pow(r,4))*cos(3.0*phi)
                 - (3.0*getCoefficient(z,3,0,1)*pow(r,2) - (3.0/16.0)*getCoefficient(z,3,2,1)*pow(r,4))*sin(3.0*phi);

            B_z = getCoefficient(z,0,1,0) - 0.25*getCoefficient(z,0,3,0)*pow(r,2) 
               //+ (1.0/64.0)*C_0_5[2*k]*pow(r,4)

               + (getCoefficient(z,1,1,0)*r - (1.0/8.0)*getCoefficient(z,1,3,0)*pow(r,3))*sin(phi)
               + (getCoefficient(z,1,1,1)*r - (1.0/8.0)*getCoefficient(z,1,3,1)*pow(r,3))*cos(phi)
               + (getCoefficient(z,2,1,0)*pow(r,2) - (1.0/12.0)*getCoefficient(z,2,3,0)*pow(r,4))*sin(2.0*phi)
               + (getCoefficient(z,2,1,1)*pow(r,2) - (1.0/12.0)*getCoefficient(z,2,3,1)*pow(r,4))*cos(2.0*phi)
               + (getCoefficient(z,3,1,0)*pow(r,3) - (1.0/16.0)*getCoefficient(z,3,3,0)*pow(r,5))*sin(3.0*phi)
               + (getCoefficient(z,3,1,1)*pow(r,3) - (1.0/16.0)*getCoefficient(z,3,3,1)*pow(r,5))*cos(3.0*phi);

            f_x = B_r*cos(phi) - B_phi*sin(phi);
            f_y = B_r*sin(phi) + B_phi*cos(phi);
*/            
            
// This should be more efficient code that is the same as above
            double r = sqrt( pow(x,2) + pow(y,2) );
            double phi = atan2( y,x );
            double sinphi = sin(phi);
            double cosphi = cos(phi);
            double sin2phi = sin(2.0*phi);
            double cos2phi = cos(2.0*phi);
            double sin3phi = sin(3.0*phi);
            double cos3phi = cos(3.0*phi);
            
            double B_r = -0.50*getCoefficient(z,0,2,0)*r //+ (1.0/16.0)*(~~~~~)*pow(r,3) 
               //- (1.0/384.0)*C_0_6[2*k]*pow(r,5)
               + (getCoefficient(z,1,0,0) - (3.0/8.0)*getCoefficient(z,1,2,0)*r*r)*sinphi
               + (getCoefficient(z,1,0,1) - (3.0/8.0)*getCoefficient(z,1,2,1)*r*r)*cosphi
               + (2.0*getCoefficient(z,2,0,0)*r - (1.0/3.0)*getCoefficient(z,2,2,0)*r*r*r)*sin2phi
               + (2.0*getCoefficient(z,2,0,1)*r - (1.0/3.0)*getCoefficient(z,2,2,1)*r*r*r)*cos2phi
               + (3.0*getCoefficient(z,3,0,0)*r*r - 5.0/16.0*getCoefficient(z,3,2,0)*r*r*r*r)*sin3phi
               + (3.0*getCoefficient(z,3,0,1)*r*r - 5.0/16.0*getCoefficient(z,3,2,1)*r*r*r*r)*cos3phi;
 
            double B_phi = (getCoefficient(z,1,0,0) - (1.0/8.0)*getCoefficient(z,1,2,0)*r*r)*cosphi
                 - (getCoefficient(z,1,0,1) - (1.0/8.0)*getCoefficient(z,1,2,1)*r*r)*sinphi
                 + (2.0*getCoefficient(z,2,0,0)*r - (1.0/6.0)*getCoefficient(z,2,2,0)*r*r*r)*cos2phi
                 - (2.0*getCoefficient(z,2,0,1)*r - (1.0/6.0)*getCoefficient(z,2,2,1)*r*r*r)*sin2phi
                 + (3.0*getCoefficient(z,3,0,0)*r*r - (3.0/16.0)*getCoefficient(z,3,2,0)*r*r*r*r)*cos3phi
                 - (3.0*getCoefficient(z,3,0,1)*r*r - (3.0/16.0)*getCoefficient(z,3,2,1)*r*r*r*r)*sin3phi;

            B_z = getCoefficient(z,0,1,0) - 0.25*getCoefficient(z,0,3,0)*r*r 
               //+ (1.0/64.0)*C_0_5[2*k]*pow(r,4)

               + (getCoefficient(z,1,1,0)*r - (1.0/8.0)*getCoefficient(z,1,3,0)*r*r*r)*sinphi
               + (getCoefficient(z,1,1,1)*r - (1.0/8.0)*getCoefficient(z,1,3,1)*r*r*r)*cosphi
               + (getCoefficient(z,2,1,0)*r*r - (1.0/12.0)*getCoefficient(z,2,3,0)*r*r*r*r)*sin2phi
               + (getCoefficient(z,2,1,1)*r*r - (1.0/12.0)*getCoefficient(z,2,3,1)*r*r*r*r)*cos2phi
               + (getCoefficient(z,3,1,0)*r*r*r - (1.0/16.0)*getCoefficient(z,3,3,0)*r*r*r*r*r)*sin3phi
               + (getCoefficient(z,3,1,1)*r*r*r - (1.0/16.0)*getCoefficient(z,3,3,1)*r*r*r*r*r)*cos3phi;

            f_x = (B_r*cosphi - B_phi*sinphi)/10000.0;
            f_y = (B_r*sinphi + B_phi*cosphi)/10000.0;
            B_z = B_z/10000.0;

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ OLD METHOD FOR DEBUGGING @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@



//f_x = 0.0;
//f_y = 1.0e-5;
//B_z = 0.0;





/*
//                    ACTUAL VERSION
   double B_r, B_phi;
   const double r = sqrt( x*x + y*y);
   const double phi = atan2(y,x);
   
   //
   int l = 1;
   
      
   B_r = 0.0;
   B_phi = 0.0;
   B_z = 0.0;
//cout << "::::x = " << x << "  y = " << y << "  z = " << z << endl;
//cout << "::::r = " << r << "  phi = " << phi << endl;
//cout << "::::B_x = " << f_x << "  B_y = " << f_y << "  B_z = " << B_z << endl;
   
cout << endl;

// something wrong below this point
   for (int i = 0; i <= m; i++)
   {
      for (int j = 0; j <= l; j++)
      {
//cout << " i = " << i << "  j = " << j << endl;
         if ( !(i == 0 && j == 0) )
         {
//cout << "hit!" << endl;
            // these are to cut down on the largest computation time wasters
            int mFact = factorial(i);
            int lFact = factorial(j);
            int lmFact = factorial(i + j);
            double sin_m_phi = sin(i*phi);
            double cos_m_phi = cos(i*phi);
            double number = pow(-1.0,j)*(double)mFact/(pow(2.0,2*j)*(double)lFact*(double)lmFact)
                         *pow(r,(2*j + i - 1));
//cout << "mFact = " << mFact << "  lFact = " << lFact << "  lmFact = " << lmFact << endl;
//cout << "sin_m_phi = " << sin_m_phi << "  cos_m_phi = " << cos_m_phi << "  number = " << number << endl;
//cout << "r^(2*j+i-1) = " << 2*j+i-1 << endl;
            if ( i == 0 )
            {
cout << " m = " << i << "  l = " << j << "  sORc = " << 0 << endl;
//cout << "again!" << endl;
                    double coefficient1 = getCoefficient(z,i,2*j,0);
               double coefficient1z = getCoefficient(z,i,2*j+1,0);

//cout << "%%% z = " << z << "  m = " << i << "  l = " << j << endl;
//cout << "getCoefficient(" << z << "," << i << "," << 2*j << ",0) = " << coefficient1 << endl;
//cout << "getCoefficient(" << z << "," << i << "," << 2*j+1 << ",0) = " << coefficient1z << endl;
               
               B_r = B_r + number*(2*j+i)*coefficient1;
               B_phi = B_phi + number*i*coefficient1;
               B_z = B_z + number*coefficient1z*r;
            }
            else
            {
//cout << "again!" << endl;
cout << " m = " << i << "  l = " << j << "  sORc = " << 0 << endl;
cout << " m = " << i << "  l = " << j << "  sORc = " << 1 << endl;
               double coefficient1 = getCoefficient(z,i,2*j,0);
               double coefficient2 = getCoefficient(z,i,2*j,1);
               double coefficient1z = getCoefficient(z,i,2*j+1,0);
               double coefficient2z = getCoefficient(z,i,2*j+1,1);
//cout << "%%% z = " << z << "  m = " << i << "  l = " << j << endl;
//cout << "getCoefficient(" << z << "," << i << "," << 2*j << ",0) = " << coefficient1 << endl;
//cout << "getCoefficient(" << z << "," << i << "," << 2*j << ",1) = " << coefficient2 << endl;
//cout << "getCoefficient(" << z << "," << i << "," << 2*j+1 << ",0) = " << coefficient1z << endl;
//cout << "getCoefficient(" << z << "," << i << "," << 2*j+1 << ",1) = " << coefficient2z << endl;



               B_r = B_r + number*(2*j+i)*coefficient1*sin_m_phi + 
                           number*(2*j+i)*coefficient2*cos_m_phi;
               B_phi = B_phi - number*i*coefficient1*sin_m_phi +
                               number*i*coefficient2*cos_m_phi;
               B_z = B_z + number*coefficient1z*sin_m_phi*r +
                           number*coefficient2z*cos_m_phi*r;
            }
         }
         else
         {
cout << " m = " << i << "  l = " << j << "  sORc = " << 0 << endl;
            B_z = B_z + getCoefficient(z,i,2*j+1,0);
         }
         
      }
   }
//   something wrong above this point ^^^

cout << endl;

// TEMPORARY CODE FOR THE PURPOSE OF FIXING BUG
//   int tempIndex = (int)((z + 200.0)/0.25);
//  if (tempIndex < 2400 && tempIndex > -1)
//   {
//      B_r = B_r + 1.0/16.0*tempArray[tempIndex][0]*pow(r,3) -1.0/384.0*tempArray[tempIndex][2]*pow(r,5);
//      B_z = B_z + 1.0/64.0*tempArray[tempIndex][1]*pow(r,4);
//   }
// TEMPORARY CODE FOR THE PURPOSE OF FIXING BUG ^^^


   f_x = B_r*cos(phi) - B_phi*sin(phi);
   f_y = B_r*sin(phi) + B_phi*cos(phi);

   
//cout << "::::B_x = " << f_x << "  B_y = " << f_y << "  B_z = " << B_z << endl;
*/
}


//*********************************************************************************
// Function: int factorial( const int &number ) const
// Purpose: returns the factorial of a number to the user
// Precondition: number >= 0
// Postcondition: factorial of number is returned
//*********************************************************************************
int MultipoleExpansion3D:: factorial (const int &number ) const
{
   // if precondition is violated exit the program
   if ( number < 0 )
   {
      cerr << "factorial: number < 0  Error!  Exiting..." << endl;
      exit(1);
   }

   int result = 1;

   for ( int i = number; i > 0; i-- )
   {
      result = result * i;
   }

   // this is to ensure that if the result is too large to fit in an int, it will
   //    theoretically roll back over into the negative region giving a negative #
   if ( result < 0 )
   {
      cerr << "factorial: result < 0 Error!  Exiting..." << endl;
      exit(1);
   }

   return result;
}



//*********************************************************************************
// Function:
// Purpose: 
// Precondition: no one element array
// Postcondition: 
//*********************************************************************************
double MultipoleExpansion3D:: getCoefficient(const double &z, const int &m,
   const int &l, int sORc) const
{
   double x1 = 0.0, x2 = 0.0, y1 = 0.0, y2 = 0.0, xn = z;
   int minIndex = 0;
   int maxIndex = length - 1;
   int guessIndex = length/2;
   
   // this loop should definitely not go through length-1 iterations...if so, 
   //    there is something wrong
   for (int i = 0; i < length; i++)
   {

//cerr << "length = " << length << endl;
//cerr << "loop # " << i << endl;

      // the top level conditionals test for when the correct indeces have
      //    been found
      if (maxIndex < minIndex)
      {
//cerr << "2) minIndex = " << minIndex << "  maxIndex = " << maxIndex << endl;
         if (maxIndex < 0 || minIndex >= length)
         {
                 /*
            // check to see if equal to the boundary point
            if (z = zValues[0])
            {
               return coefficients[0][m][l][sORc];
            }
            else if (z == zValues[length - 1])
            {
               return coefficients[length - 1][m][l][sORc];
            }
            */
            cerr << "MultipoleExpansion3D:: getCoefficient: maxIndex = 0"
                 << " or minIndex >= length  Error!  Exiting..." << endl
                 << "minIndex = " << minIndex << "  maxIndex = " << maxIndex
                 << "  guessIndex = " << guessIndex  << "  z = " << z << endl;
            exit(1);
         }
         x1 = zValues[maxIndex];
         x2 = zValues[minIndex];
         y1 = coefficients[maxIndex][m][l][sORc];
         y2 = coefficients[minIndex][m][l][sORc];
         xn = z;
         break;
      }
      
      if (z == zValues[guessIndex])
      {
         return coefficients[guessIndex][m][l][sORc];
      }
      else if (z > zValues[guessIndex])
      {
//cerr << " z > coefficients" << endl;
         minIndex = guessIndex + 1;
      }
      else if (z < zValues[guessIndex])
      {
//cerr << " z < coefficients" << endl;
         maxIndex = guessIndex - 1;
      }
      
      guessIndex = (minIndex + maxIndex)/2;

      if (i==length-1)
      {
         cerr << "MultipoleExpansion3D:: getCoefficient:  i == length-1  Error!"
              << "  Exiting..." << endl
              << "minIndex = " << minIndex << "  maxIndex = " << maxIndex
              << "  guessIndex = " << guessIndex << endl;
         exit(1);
      }
   }
//cerr << " z = " << z << "   zValues[minIndex] = " << zValues[minIndex] << "   zValues[maxIndex] = " << zValues[maxIndex] << endl;

   if (x2 - x1 == 0)
   {
      cerr << "MultipoleExpansion3D:: getCoefficient:  Cannot divide by zero!"
           << " Error!  Exiting..." << endl
           << "minIndex = " << minIndex << "  maxIndex = " << maxIndex
           << "  guessIndex = " << guessIndex << endl;
      exit(1);
   }
   
   // same point cannot be repeated and (x2-x1)!=0
   return ((y2 - y1)/(x2 - x1)*(xn - x2) + y2);     
}





///////////////////////////////////////////////////////////////////////////
//
// END OF FILE
//
///////////////////////////////////////////////////////////////////////////
