#ifndef CPPEXTERNALEFFECTS_HH_
#define CPPEXTERNALEFFECTS_HH_



#include "Python.h"

#include "ExternalEffects.hh"
#include <complex>

typedef std::complex<double>	tcomplex;

namespace Tracker3DField{
	
	class  LasStripExternalEffects: public ExternalEffects
	{
		public:
		
			/** Constructor. */
			LasStripExternalEffects(double a,double b,double c, char *addressEG,int states);
			
			/** Destructor. */
			~LasStripExternalEffects();
		
		/** It initializes effects. */
		void setupEffects(Bunch* bunch);
		
		/** It finalizes effects. */
		void finalizeEffects(Bunch* bunch);

		/** It applies the external effects to a particle with certain index. */
		void applyEffects(Bunch* bunch, int index, 
	                            double* y_in_vct, double* y_out_vct, 
														  double t, double t_step, 
														  OrbitUtils::BaseFieldSource* fieldSource,
															RungeKuttaTracker* tracker);	
		
		/** It defines parameters of the laser beam**/
		void	setLaserHalfAngle(double a);
		double	getLaserHalfAngle();
		void	setLaserPower(double a);
		double	getLaserPower();
		void	setLaser_lambda(double a);
		double	getLaser_lambda();
		
		void	read_transitions(char* adressEG,int states);
		
		/** It gets laser field in the laser frame**/
		void GetOwnLaserField(double x, double y, double z, double t, double& E_x, double& E_y, double& E_z,double& H_x, double& H_y, double& H_z);
		/** It gets laser field in laboratory frame**/

		

		

		  private:

			  
			  //parameters of dipole transition, energy, lifetime, spontaneous relaxation, 
			  
			  double*** dipole_transition_x;
			  double*** dipole_transition_y;
			  double*** dipole_transition_z;
			  double*** gamma_spontaneous_relax;
			  double** energy;
			  double** gamma_autoionization;
			  
			  //this array is used on each step of solution of density matrix equation at definite field  
			  tcomplex*** exp_mu_El;
			  tcomplex*** k_RungeKutt;
			  double* Gamma_i;
			  double* E_i;
			  double** gamma_ij;
			  bool** cond;

			  
			  double delta_F;
			  int n_data;
			  int levels;

			  
			  //Laser parameters
			  double Laser_lambda;
			  double LaserPower;
			  double Laser_half_angle;
			  
		  
			  //Parameters of orientation of laser in space
			  double param1;
			  double param2;
			  double param3;
			  double param4;
			  
			  //Normized direction of the Poiting (or wave) vector
			  double nx_lab;
			  double ny_lab;
			  double nz_lab;
			 
			  //time and frequensy of laser in frame of particle
			  double omega_part;
			  double part_t_step;
			  double t_part;
			  double phasa_part;
			  					  
			  //Fields in the particle frame  
			  double Ex_las;
			  double Ey_las;
			  double Ez_las;
			  
			  double Ex_stat;
			  double Ey_stat;
			  double Ez_stat;
			  			 

			  
			  			  
				/** It converts laser field from laser frame to the laboratory frame**/
				void GetLabLaserField(double x, double y, double z, double t, double& E_x, double& E_y, double& E_z,double& B_x, double& B_y, double& B_z);
					
				/**Solver for Amplitudes**/
				void AmplSolver4step(int i,Bunch* bunch);
				
				/*this method gives koefficient of autoionization as a function of field  */
				void GetEnergyAutoionization(int k,double& E, double& Gamma);
				
				/*this method gives frequendy transitions as a function of field */
				void GetDipoleTransition(int k,int ks,double& relax,double& mu_x,double& mu_y,double& mu_z);				
				
				/*this method gives laser and static fields transformed by rotation relatively z axes  */
				void GetFrameParticleParameters(int i, double t_step, Bunch* bunch);
				

				
	};
};



#endif /*CPPEXTERNALEFFECTS_HH_*/


