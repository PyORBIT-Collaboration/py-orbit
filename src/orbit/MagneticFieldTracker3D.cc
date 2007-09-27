//******************************************************************************
// Program: 3DMagneticFieldTracker.cpp
// Creator: Matt Perkett
// Date: 10/09/2006
// Purpose: tracks particles through a magnetic field
// Input:
// Output:
// Algorithm:
// Limitations/Bugs:
// Modified: 10/31/2006 -> I just made a change that could possibly do something
//   that I didn't think of, I changed the orbit_xi,y,z initialization values to
//   include the reference particle x and y values so that the method_flag = 1
//   will work properly
//******************************************************************************

#include "MagneticFieldTracker3D.hh"
#include "ParticleMassCharge.hh"
#include "BaseFieldSource.hh"
#include "OrbitConst.hh"
#include<iomanip>
#include<fstream>
#include<string>
#include <math.h>

//******************************************************************************
// Function: MagneticFieldTracker3D()
// Purpose: default constructor
// Precondition:
// Postcondition:
//******************************************************************************
MagneticFieldTracker3D:: MagneticFieldTracker3D()
{
   ds = 0.0001;
}


//******************************************************************************
// Function: MagneticFieldTracker3D(double stepSize)
// Purpose: alternate constructor
// Precondition:
// Postcondition:
//******************************************************************************
MagneticFieldTracker3D:: MagneticFieldTracker3D(double stepSize)
{
   ds = stepSize;
}


//******************************************************************************
// Function: ~MagneticFieldTracker3D()
// Purpose: destructor
// Precondition:
// Postcondition:
//******************************************************************************
MagneticFieldTracker3D:: ~MagneticFieldTracker3D()
{
   
}


//******************************************************************************
// Function: track(double zmin, double zmax)
// Purpose: tracks a particle from zmin to zmax
//    - method_flag = 1 then do NOT transfer back to orbit coordinates
//    - method_flag = 0 then do transfer back to orbit coordinates
//    - make sure to call anything except the initial call of track() with
//        alpha = beta = gamma = x_0 = y_0 = 0 or it will not work properly
//    - ensure that the reference particle has x and y set to 0 initially or it
//        will do something that you don't want
//
// Modifications to be made:  
//      -> use orbit constants file wherever possible
//      -> what is the actual charge of the reference particle?
//      -> I'm still not sure on the transfer back to ORBIT coordinates
//
// Precondition: zmin < zmax, "masscharge" is an attribute of bunch
// Postcondition:
//******************************************************************************
void MagneticFieldTracker3D::track(double zmin, double zmax, double x_0, double y_0,
double &a, double &b, double &g, Bunch* bunch, BaseFieldSource* fieldSource, int method_flag)
{
   // make sure that region tracking over makes sense
   if (zmin >= zmax)
   {
      cerr << "MagneticFieldTracker3D::track : zmin >= zmax  Error!  Exiting..." << endl;
      exit(1);
   }
   else if ( !(method_flag == 0 || method_flag == 1) )
   {
      cerr << "MagneticFieldTracker3D::track: method_flag does not equal either 0"
           << " or 1  Error!  Exiting..." << endl;
      exit(1);
   }
   
   // later add a conditional that checks to ensure that "masscharge" is an attribute

   const double c = OrbitConst::c;  // speed of light
   
   // get array from bunch containing all the particles, the masscharge particle attribute
   //    and syncronous particle
   double** values = bunch->coordArr();
   ParticleMassCharge* attributes = (ParticleMassCharge*) bunch->getParticleAttributes("masscharge");
   SyncPart* refParticle = bunch->getSyncPart();

   // variables used by every particle in loop
   double refPX = 0.0, refPY = 0.0, refPZ = 0.0;
   double refX = 0.0, refY = 0.0;
   double theta1 = -9.99E-99, theta2 = -9.99E-99;

   ///////////////////////////////////////////////////////
   // for file output
/*
   ofstream output;
   output.open("output.dat",ios::app);
*/
   ofstream before;
   before.open("before.dat");
   ofstream after;
   after.open("after.dat");
   
   // before.dat
   for (int i = -1; i < bunch->getSize(); i++)
   {
      if (i == -1)
      {
         before << refParticle->getX() << '\t' << refParticle->getY() << endl;
      }
      else
      {
         before << values[i][0] << '\t' << values[i][2] << endl;
      }
   }
   ///////////////////////////////////////////////////////
   
   // track reference particle first, then track every other particle in bunch
   for (int i = -1; i < bunch->getSize(); i++)
   {
      ///////////////////////////////////////////////////////
      // for file output
      char str[5];
      sprintf(str, "%d", i);
      string fileName = "particle_";
      fileName = fileName + str + ".dat";
      //ofstream output;
      //output.open(fileName.c_str(),ios::app);
      ///////////////////////////////////////////////////////
      
      // variable initializations
      double mass = 0.0, q = 0.0, x_i = 0.0, y_i = 0.0, z_i = 0.0, 
             speedx = 0.0, speedy = 0.0, speedz = 0.0,
             orbit_speedx = 0.0, orbit_speedy = 0.0, orbit_speedz = 0.0,
             orbit_x_i = 0.0, orbit_y_i = 0.0, orbit_z_i = 0.0;
      double p = 0.0, W = 0.0, beta = 0.0, gamma = 0.0, magnitude = 0.0;
      double dx = -9.99E-99, dy = -9.99E-99, dz = -9.99E-99;
      double dv_x = -9.99E-99, dv_y = -9.99E-99, dv_z = -9.99E-99;
      double B_x = -9.99E-99, B_y = -9.99E-99, B_z = -9.99E-99;
      double sDist = 0.0;
      int counter = 0;
             
      // FOR REFERENCE PARTICLE
      if (i == -1)
      {
         // mass, p, and W are in MeV
         mass = refParticle->getMass();
         refPX = refParticle->getPX();
         refPY = refParticle->getPY();
         refPZ = refParticle->getPZ();
         refX = refParticle->getX();
         refY = refParticle->getY();
         p = sqrt(refPX*refPX + refPY*refPY + refPZ*refPZ);
         W = sqrt(mass*mass + p*p);
         beta = p/W;
         gamma = W/mass;
         
         // speed in m/s
         orbit_speedx = beta*c*refParticle->getPX()/p;
         orbit_speedy = beta*c*refParticle->getPY()/p;
         orbit_speedz = beta*c*refParticle->getPZ()/p;
         
         // x_i, y_i, and z_i in m
         orbit_x_i = (x_0 + refParticle->getX())*0.001;
         orbit_y_i = (y_0 + refParticle->getY())*0.001;
         orbit_z_i = zmin;
         
         // convert mass and charge into mks units
         mass = mass*1.78266173e-30;      // mass (kg)
         
//!!! what is the actual charge of the reference particle?

         q = 1.0;                         // charge (0,+/-1,+/-2...)
         q = q*1.60217733e-19;            // charge (Coulombs)

      }
      // FOR ALL OTHER PARTICLES
      else
      {
         // mass, p, and W are in MeV
         mass = attributes->getMass(i);
         p = sqrt( pow( (refPX + values[i][1]),2 ) + 
                   pow( (refPY + values[i][3]),2 ) + 
                   pow( (refPZ + values[i][5]),2 ) );
         W = sqrt(mass*mass + p*p);
         beta = p/W;
         gamma = W/mass;

         // speeds in m/s
         orbit_speedx = beta*c*(refPX + values[i][1])/p;
         orbit_speedy = beta*c*(refPY + values[i][3])/p;
         orbit_speedz = beta*c*(refPZ + values[i][5])/p;

         // x_i, y_i, and z_i in m
         orbit_x_i = (x_0 + refX + values[i][0])*0.001;
         orbit_y_i = (y_0 + refY + values[i][2])*0.001;
         orbit_z_i = zmin;

         // convert mass and charge to mks units
         q = attributes->getCharge(i);      // charge (0,+/-1,+/-2...)
         q = q*1.60217733e-19;              // charge (Coulombs)              
         mass = mass*1.78266173e-30;        // mass (kg)
      }
      
      // EVERYTHING IS IN MKS UNITS
      
      // transfer from ORBIT coordinate system to cartesian coordinates
      //   Euler angles a (alpha) - rotation about the z-axis, b (beta) - rotation
      //   about the y-axis, and g (gamma) - rotation about the z-axis
      speedx = ( ( orbit_speedx*cos(a) + orbit_speedy*sin(a) )*cos(b)
               - orbit_speedz*sin(b) )*cos(g) + ( -1.0*orbit_speedx*sin(a)
               + orbit_speedy*cos(a) )*sin(g);
      speedy = -1.0*( (orbit_speedx*cos(a) + orbit_speedy*sin(a) )*cos(b)
               - orbit_speedz*sin(b) )*sin(g) + (-1.0*orbit_speedx*sin(a) 
               + orbit_speedy*cos(a) )*cos(g);
      speedz = ( orbit_speedx*cos(a) + orbit_speedy*sin(a) )*sin(b)
               + (orbit_speedz)*cos(b);

      x_i = ( ( orbit_x_i*cos(a) + orbit_y_i*sin(a) )*cos(b)
               - orbit_z_i*sin(b) )*cos(g) + ( -1.0*orbit_x_i*sin(a)
               + orbit_y_i*cos(a) )*sin(g);
      y_i = -1.0*( (orbit_x_i*cos(a) + orbit_y_i*sin(a) )*cos(b)
               - orbit_z_i*sin(b) )*sin(g) + (-1.0*orbit_x_i*sin(a) 
               + orbit_y_i*cos(a) )*cos(g);
      z_i = ( orbit_x_i*cos(a) + orbit_y_i*sin(a) )*sin(b)
               + (orbit_z_i)*cos(b);

      // speed of particle
      magnitude = beta*c;

      // unit vector of speed
      double v_ix = speedx/magnitude, 
             v_iy = speedy/magnitude, 
             v_iz = speedz/magnitude;
             
      // a constant used in tracking portion
      const double kappa = q/(mass*gamma*magnitude);

      // 'guess' at position
      double x_g = x_i,
             y_g = y_i,
             z_g = z_i;
             
      // 'guess' at speed
      double v_gx = v_ix, 
             v_gy = v_iy, 
             v_gz = v_iz;    

/*
      // print relevant information to screen
      cout << endl;
      cout << "~~~~~~~~ Particle " << i << " Information ~~~~~~~~" << endl;
      cout << "mass (kg)      -> " << mass << endl;
      cout << "charge (C)     -> " << q << endl;
      cout << "magnitude (m/s)-> " << magnitude << endl;
      cout << "p (MeV)        -> " << p << endl;
//      cout << "W (MeV)        -> " << W << endl;
//      cout << "beta           -> " << beta << endl;
//      cout << "gamma          -> " << gamma << endl;
//      cout << "kappa          -> " << kappa << endl;
      cout << "x_i = " << x_i << "  y_i = " << y_i << "  z_i = " << z_i << endl;
//      cout << "x_g = " << x_g << "  y_g = " << y_g << "  z_g = " << z_g << endl;
      cout << "speedx = " << speedx << "  speedy = " << speedy << "  speedz = " 
           << speedz << endl;
      cout << "v_ix = " << v_ix << "  v_iy = " << v_iy << "  v_iz = " << v_iz << endl;
//      cout << "v_gx = " << v_gx << "  v_gy = " << v_gy << "  v_gz = " << v_gz << endl;
      cout << "-" << endl;
*/

      // TRACKING
      //    continue until particle reaches zmax
      while (z_i < zmax)
      {
         counter++;

/* 
         // FOR DEBUGGING PURPOSES         
         if (counter%100000 == 0)
         {
            cerr << "iteration " << counter << endl;
         }
*/
/*
         ///////////////////////////////////////////////////////
         // for file output
         if (counter%500 == 0)
         {
            output << x_i << '\t' << y_i << '\t' << z_i << endl;
         }
         ///////////////////////////////////////////////////////
*/
         
         // Prevents an infinite loop
         if ( counter > 5000000 )
         {
            cerr << "counter > 5000000   break!" << endl;
            break;
         }
              
         // initializations
         x_g = x_i;
         y_g = y_i;
         z_g = z_i;
         v_gx = v_ix;
         v_gy = v_iy;
         v_gz = v_iz;
         
         for ( int j = 1; j <= 2; j++ )
         {
            dx = v_gx*ds;
            dy = v_gy*ds;
            dz = v_gz*ds;

            x_g = x_i + dx/2.0;
            y_g = y_i + dy/2.0;
            z_g = z_i + dz/2.0;

            // getMagneticfield requires arguments in cm
            fieldSource->getMagneticField(0.0,x_g*100.0,y_g*100.0,z_g*100.0,B_x,B_y,B_z);

            dv_x = (v_gy*B_z - v_gz*B_y)*kappa*ds;
            dv_y = (v_gz*B_x - v_gx*B_z)*kappa*ds;
            dv_z = (v_gx*B_y - v_gy*B_x)*kappa*ds;

            v_gx = v_ix + dv_x/2.0;
            v_gy = v_iy + dv_y/2.0;
            v_gz = v_iz + dv_z/2.0;
         }

         x_i = x_i + dx;
         y_i = y_i + dy;
         z_i = z_i + dz;
         v_ix = v_ix + dv_x;
         v_iy = v_iy + dv_y;
         v_iz = v_iz + dv_z; 

         // if velocity no longer makes a unit vector then exit program
         if (!isUnitVector (v_ix,v_iy,v_iz))
         {
            cerr << "MagneticFieldTracker3D::track:  Unit Vector Condition Failed!  "
                 << "Exiting..." << endl;
            exit(1);
         }

         // increment total distance traveled along curve
         sDist += ds;
      }
      
      // multiply original magnitude of momentum by v_ix,v_iy,&v_iz to get
      //     momentum in respective directions
      double px = p*v_ix,
             py = p*v_iy,
             pz = p*v_iz;

      // FOR REFERENCE PARTICLE
      if (i == -1)
      {
         // THE REFERENCE PARTICLE WILL RETAIN THE RELEVANT INFORMATION NEEDED
         //   BY OTHER PARTICLES TO TRANSFER INTO ORBIT COORDINATES UNTIL THE
         //   ENTIRE TRACKING PROCESS IS COMPLETE
         
         // angles to tranfer back to ORBIT
         theta1 = atan2(v_ix,v_iz);
         theta2 = atan2(v_iy,v_iz);
           
         // method_flag == 1 means do NOT transfer back into ORBIT coordinates
         if (method_flag == 1)
         {
            theta1 = 0.0;
            theta2 = 0.0;
         }
         
         refParticle->setPX( px*cos(theta1) - pz*sin(theta1) );
         refParticle->setPY( (px*sin(theta1) + pz*cos(theta1))*sin(theta2)
                            + py*cos(theta2) );
         refParticle->setPZ( (px*sin(theta1) + pz*cos(theta1))*cos(theta2)
                            - py*sin(theta2) );

         // multiply by 1000.0 for the x and y components to convert back to mm
         refParticle->setX( (x_i*cos(theta1) - z_i*sin(theta1))*1000.0 );
         refParticle->setY( ((x_i*sin(theta1) + z_i*cos(theta1))*sin(theta2)
                            + y_i*cos(theta2))*1000.0 );
         refParticle->setZ( sDist );
      }
      // FOR ALL OTHER PARTICLES
      else
      {
         // momentum
         values[i][1] = px*cos(theta1) - pz*sin(theta1) - refParticle->getPX();
         values[i][3] = (px*sin(theta1) + pz*cos(theta1))*sin(theta2)
                        + py*cos(theta2) - refParticle->getPY();
         values[i][5] = (px*sin(theta1) + pz*cos(theta1))*cos(theta2) 
                        - py*sin(theta2) - refParticle->getPZ();
   
         // x,y,z
         // multiply by 1000.0 for the x and y components to convert back to mm
         values[i][0] = (x_i*cos(theta1) - z_i*sin(theta1))*1000.0
                        - refParticle->getX();
         values[i][2] = ((x_i*sin(theta1) + z_i*cos(theta1))*sin(theta2) 
                        + y_i*cos(theta2))*1000.0 - refParticle->getY();
         values[i][4] = sDist - refParticle->getZ();
      }

/*
      // finish printing relevant information
      cout << "-" << endl;
      cout << "x_f = " << x_i << "  y_f = " << y_i << "  z_f = " << z_i << endl;
      cout << "v_fx = " << v_ix << "  v_fy = " << v_iy << "  v_fz = " << v_iz << endl;
      cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
      cout << endl;
*/

   }
   
/*
      // THIS SECTION DOES NOT WORK PROPERLY YET
      
      // set a, b, and g to Euler angles needed to transfer from the resulting
      //    ORBIT coordinate system back into the tracking coordinate system
      double mag = sqrt(refParticle->getPX()*refParticle->getPX() + refParticle->getPY()*
                        refParticle->getPY() + refParticle->getPZ()*refParticle->getPZ());
                        
      // A, B, and C are the i,j,k components of vector pointing in the final direction
      //   of reference particle (in the tracker's coordinate system) - it's normalized
      double A = refParticle->getPX()/mag;
      double B = refParticle->getPY()/mag;
      double C = refParticle->getPZ()/mag;
      
      // see notebook for justification
      a = acos( (-1.0*C*B)/(sqrt(B*B + A*A)*sqrt(C*C + A*A)) );
      b = acos(C);
      g = acos(A/sqrt(B*B + A*A));
*/

      // set reference particle attributes correctly if method_flag = 0
      //   (indicating a transfer back to ORBIT coordinates)
   if (method_flag == 0)
   {
      refParticle->setXYZ(0.0,0.0,0.0);
   }

   // after.dat
   for (int i = -1; i < bunch->getSize(); i++)
   {
      if (i == -1)
      {
         after << refParticle->getX() << '\t' << refParticle->getY() << endl;
      }
      else
      {
         after << values[i][0] << '\t' << values[i][2] << endl;
      }
   }
   
   // close files
   before.close();
   after.close();
//   output.close();
   
}


//******************************************************************************
// Function: void setStepSize(double stepSize)
// Purpose: sets the step size
// Precondition: step size > 0.0
// Postcondition: ds = stepSize
//******************************************************************************
void MagneticFieldTracker3D:: setStepSize(double stepSize)
{
   if (ds <= 0.0)
   {
      cerr << "MagneticFieldTracker3D:: setStepSize: stepSize <= 0.0  Error!"
           << "  Exiting..." << endl;
      exit(1);
   }
   ds = stepSize;
}


//******************************************************************************
// Function: double getStepSize() const
// Purpose: returns the step size
// Precondition: none
// Postcondition: returns the step size
//******************************************************************************
double MagneticFieldTracker3D:: getStepSize() const
{
   return ds;
}


//*********************************************************************************
// Function: bool isUnitVector( double n1, double n2, double n3 )
// Purpose: tests to see if input is a unit vector (i.e. n1^2 + n2^2 + n3^2 = 1)
// Precondition: none
// Postcondition: returns True if it is a unit vector and False if it isn't
//*********************************************************************************
bool MagneticFieldTracker3D::isUnitVector (double n1, double n2, double n3)
{
//   cerr << fabs(n1*n1 + n2*n2 + n3*n3 - 1.0) << endl;
   return fabs(n1*n1 + n2*n2 + n3*n3 - 1.0) < 1.0E-3;
}


//******************************************************************************
// Function: 
// Purpose:
// Precondition:
// Postcondition:
//******************************************************************************



//EOF
