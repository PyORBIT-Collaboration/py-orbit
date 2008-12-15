//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    CppExternalEffects.cc
//
// CREATED
//    06/27/2008
//
// DESCRIPTION
//    The base class for C++ implementation of a external effects 
//    during the transport of particles through the external field. 
//    It should be sub-classed on Python level and implement
//    setupEffects(Bunch* bunch)
//    finalizeEffects(Bunch* bunch)
//    applyEffects(Bunch* bunch, int index, 
//	                            double* y_in_vct, double* y_out_vct, 
//														  double t, double t_step, 
//														  OrbitUtils::BaseFieldSource* fieldSource)
//    methods.
//    The results of these methods will be available from the c++ level.
//    This is an example of embedding Python in C++ Orbit level.
//
//     
//         
//       
//     
//
///////////////////////////////////////////////////////////////////////////
#include "orbit_mpi.hh"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>

#include "LasStripExternalEffects.hh"
#include "RungeKuttaTracker.hh"
#include "OrbitConst.hh"
#include "LorentzTransformationEM.hh"


#define Re(i,n,m) bunch->arrAttr[i][(n-1)*levels+m]		//i-part index, n,m-attr index
#define Im(i,n,m) bunch->arrAttr[i][(n-1)*levels+m+levels*levels]
#define phasa(i)  bunch->arrAttr[i][0]
#define time(i)  bunch->arrAttr[i][2*levels*levels+1]
#define x0(i)  bunch->arrAttr[i][2*levels*levels+2]
#define y0(i)  bunch->arrAttr[i][2*levels*levels+3]
#define z0(i)  bunch->arrAttr[i][2*levels*levels+4]
#define px0(i)  bunch->arrAttr[i][2*levels*levels+5]
#define py0(i)  bunch->arrAttr[i][2*levels*levels+6]
#define pz0(i)  bunch->arrAttr[i][2*levels*levels+7]
#define dm(i,n,m) tcomplex(bunch->arrAttr[i][(n-1)*levels+m],bunch->arrAttr[i][(n-1)*levels+m+levels*levels])
#define J tcomplex(0.,1.)
#define k_rk(j,n,m) k_RungeKutt[j][(n-1)*levels+m]

typedef std::complex<double>	tcomplex;

using namespace Tracker3DField;
using namespace OrbitUtils;

LasStripExternalEffects::LasStripExternalEffects(double a,double b,double c,char* addressEG,int states)
{
	setName("unnamed");

	Laser_half_angle=a;
	LaserPower=b;
	Laser_lambda=c;

	
	read_transitions(addressEG,states);	

//allocating memory for koefficients of 4-th order Runge-Kutta method and other koeeficients of the master equation

k_RungeKutt=new tcomplex**[levels+1];	for (int i=0;i<levels+1;i++)	k_RungeKutt[i]=new tcomplex*[levels+1]; for (int i=0;i<levels+1;i++) for (int j=0;j<levels+1;j++)	k_RungeKutt[i][j]=new tcomplex[5];
exp_mu_El=new tcomplex**[levels+1];	for (int i=0;i<levels+1;i++)	exp_mu_El[i]=new tcomplex*[levels+1]; for (int i=0;i<levels+1;i++) for (int j=0;j<levels+1;j++)	exp_mu_El[i][j]=new tcomplex[5];
gamma_ij=new double*[levels+1];	for (int i=0;i<levels+1;i++)	gamma_ij[i]=new double[levels+1];
cond=new bool*[levels+1];	for (int i=0;i<levels+1;i++)	cond[i]=new bool[levels+1];
Gamma_i=new double[levels+1];
E_i=new double[levels+1];


}


LasStripExternalEffects::~LasStripExternalEffects()
{

	for (int i=0;i<levels+1;i++) for (int j=0;j<levels+1;j++)	delete [] k_RungeKutt[i][j]; for (int i=0;i<levels+1;i++)	delete [] k_RungeKutt[i];	delete	[]	k_RungeKutt;
	for (int i=0;i<levels+1;i++) for (int j=0;j<levels+1;j++)	delete [] exp_mu_El[i][j]; for (int i=0;i<levels+1;i++)	delete [] exp_mu_El[i];	delete	[]	exp_mu_El;		
	delete [] E_i;
	delete [] Gamma_i;
	for (int i=0;i<levels+1;i++)	delete	[]	gamma_ij[i];	delete	[]	gamma_ij;
	for (int i=0;i<levels+1;i++)	delete	[]	cond[i];		delete	[]	cond;
	for (int i=0;i<levels+1;i++)	delete	[]	energy[i];		delete	[]	energy;	
	for (int i=0;i<levels+1;i++)	delete	[]	gamma_autoionization[i];		delete	[]	gamma_autoionization;
	for (int i=0;i<levels+1;i++) for (int j=0;j<levels+1;j++)	delete [] gamma_spontaneous_relax[i][j]; for (int i=0;i<levels+1;i++)	delete [] gamma_spontaneous_relax[i];	delete	[]	gamma_spontaneous_relax;
	for (int i=0;i<levels+1;i++) for (int j=0;j<levels+1;j++)	delete [] dipole_transition_x[i][j]; for (int i=0;i<levels+1;i++)	delete [] dipole_transition_x[i];	delete	[]	dipole_transition_x;
	for (int i=0;i<levels+1;i++) for (int j=0;j<levels+1;j++)	delete [] dipole_transition_y[i][j]; for (int i=0;i<levels+1;i++)	delete [] dipole_transition_y[i];	delete	[]	dipole_transition_y;
	for (int i=0;i<levels+1;i++) for (int j=0;j<levels+1;j++)	delete [] dipole_transition_z[i][j]; for (int i=0;i<levels+1;i++)	delete [] dipole_transition_z[i];	delete	[]	dipole_transition_z;

}




void LasStripExternalEffects::setupEffects(Bunch* bunch){
}
		

	
void LasStripExternalEffects::finalizeEffects(Bunch* bunch){
}



void LasStripExternalEffects::applyEffects(Bunch* bunch, int index, 
	                            double* y_in_vct, double* y_out_vct, 
														  double t, double t_step, 
														  BaseFieldSource* fieldSource,
															RungeKuttaTracker* tracker)
{
	double** xyz = bunch->coordArr();
	double B_x,B_y,B_z;
	double mass = bunch->getMass();
	
		for (int i=0; i<bunch->getSizeGlobal();i++)	{

		GetLabLaserField(xyz[i][0],xyz[i][2],xyz[i][4],t,Ex_las,Ey_las,Ez_las,B_x,B_y,B_z);
		LorentzTransformationEM::transform(mass,
			                                 xyz[i][1],xyz[i][3],xyz[i][5],
																			 Ex_las,Ey_las,Ez_las,
																			 B_x,B_y,B_z);	
		
			
		fieldSource->getElectricField(xyz[i][0],xyz[i][2],xyz[i][4],t,Ex_stat,Ey_stat,Ez_stat);	
		fieldSource->getMagneticField(xyz[i][0],xyz[i][2],xyz[i][4],t,B_x,B_y,B_z);		
		LorentzTransformationEM::transform(mass,
			                                 xyz[i][1],xyz[i][3],xyz[i][5],
																			 Ex_stat,Ey_stat,
																			 Ez_stat,B_x,B_y,B_z);	
			
			GetFrameParticleParameters(i,t_step,bunch);	//	This function gives parameters Ez_stat	Ex_las	Ey_las	Ez_las	t_part	omega_part	part_t_step phasa_part

			ofstream file("data_ampl.txt",ios::app);
			file<<t<<"\t";
			for(int n=1;n<levels+1;n++)	file<<Re(i,n,n)<<"\t";
			double sum=0;	for(int n=1;n<levels+1;n++)	sum+=Re(i,n,n);
			file<<sum<<"\n";
			file.close();			
		
//				cout<<scientific<<setprecision(10)<<bunch->x(0)<<"\t"<<bunch->y(0)<<"\t"<<bunch->z(0)<<"\n";			
				
			
		AmplSolver4step(i,bunch);	
		
	
		
		}
	


	
	
	
//	if (getName()=="second_effect")	{
//	cout<<scientific<<setprecision(20)<<bunch->x(0)<<"\t"<<bunch->y(0)<<"\t"<<bunch->z(0)<<"\n";
//	}
	

	
}

void	LasStripExternalEffects::setLaserHalfAngle(double a)	{	Laser_half_angle=a;}

double	LasStripExternalEffects::getLaserHalfAngle()	{	return Laser_half_angle;}

void	LasStripExternalEffects::setLaserPower(double a)	{	LaserPower=a;}

double	LasStripExternalEffects::getLaserPower()	{	return LaserPower;}

void	LasStripExternalEffects::setLaser_lambda(double a)	{	Laser_lambda=a;}

double	LasStripExternalEffects::getLaser_lambda()	{	return Laser_lambda;}








void	LasStripExternalEffects::read_transitions(char* addressEG,int states)	{

	ifstream file;
	double F,alpha=7.297352570e-3;
	int k,ks,fi;
	char nameEG[1024];

	levels=states*(1+states)*(1+2*states)/6;

	
//this function measures parameter of input files (length, delta_F) using groung level file 000.txt
sprintf(nameEG,"%s000.txt",addressEG);	
file.open(nameEG);file>>F>>F>>F>>delta_F; file.close();
file.open(nameEG);fi=0;	while(!file.eof())	{file>>F>>F>>F;fi++;} file.close();n_data=fi-1;


//allocating memory for dynamic massive of data that will be read fron files
energy=new double*[levels+1];	for (int i=0;i<levels+1;i++)	energy[i]=new double[n_data+10];
gamma_autoionization=new double*[levels+1];	for (int i=0;i<levels+1;i++)	gamma_autoionization[i]=new double[n_data+10];
gamma_spontaneous_relax=new double**[levels+1];	for (int i=0;i<levels+1;i++)	gamma_spontaneous_relax[i]=new double*[levels+1]; for (int i=0;i<levels+1;i++) for (int j=0;j<levels+1;j++)	gamma_spontaneous_relax[i][j]=new double[n_data+10];
dipole_transition_x=new double**[levels+1];	for (int i=0;i<levels+1;i++)	dipole_transition_x[i]=new double*[levels+1]; for (int i=0;i<levels+1;i++) for (int j=0;j<levels+1;j++)	dipole_transition_x[i][j]=new double[n_data+10];
dipole_transition_y=new double**[levels+1];	for (int i=0;i<levels+1;i++)	dipole_transition_y[i]=new double*[levels+1]; for (int i=0;i<levels+1;i++) for (int j=0;j<levels+1;j++)	dipole_transition_y[i][j]=new double[n_data+10];
dipole_transition_z=new double**[levels+1];	for (int i=0;i<levels+1;i++)	dipole_transition_z[i]=new double*[levels+1]; for (int i=0;i<levels+1;i++) for (int j=0;j<levels+1;j++)	dipole_transition_z[i][j]=new double[n_data+10];	





	//this loop reads energies and autoionization coefficients
	for(int n=1;n<states+1;n++)
	for(int m=-(n-1);m<(n-1)+1;m++)
	for(int n1=0;n1<n-abs(m)-1+1;n1++)	{
		
		k=1+m*n+n1+(n*n*n-n)/3-m*abs(m-1)/2;
		
sprintf(nameEG,"%s%i%i%i.txt",addressEG,n1,n-n1-abs(m)-1,abs(m));		
file.open(nameEG);for (fi=0; fi<n_data; fi++)	{file>>F>>energy[k][fi]>>gamma_autoionization[k][fi];}	file.close();

	}
	

	
	
	//this double loop reads dipole transitions anf fills spontaneous relaxation coefficients
	for(int n=1;n<states+1;n++)
	for(int m=-(n-1);m<(n-1)+1;m++)
	for(int n1=0;n1<n-abs(m)-1+1;n1++)

	for(int ns=1;ns<states+1;ns++)
	for(int ms=-(ns-1);ms<(ns-1)+1;ms++)
	for(int n1s=0;n1s<ns-abs(ms)-1+1;n1s++)		{
			
			k=1+m*n+n1+(n*n*n-n)/3-m*abs(m-1)/2;
			ks=1+ms*ns+n1s+(ns*ns*ns-ns)/3-ms*abs(ms-1)/2;
			
	sprintf(nameEG,"%s%i%i%i---%i%i%i.txt",addressEG,n1,n-n1-abs(m)-1,m,n1s,ns-n1s-abs(ms)-1,ms);
	file.open(nameEG);
	for (fi=0; fi<n_data; fi++)	{file>>F>>dipole_transition_x[k][ks][fi]>>dipole_transition_y[k][ks][fi]>>dipole_transition_z[k][ks][fi];
	
	// this condition assumes that probability of spontaneous (see next line) and indused tansition  between levels with the same principal quantum number n is sero
	if(ns==n)	{dipole_transition_x[k][ks][fi]=0;dipole_transition_y[k][ks][fi]=0;dipole_transition_z[k][ks][fi]=0;}
	
	//this loop fills spontaneous relaxation (transition) of atom in relative atomic units
	gamma_spontaneous_relax[k][ks][fi]=fabs((4*alpha*alpha*alpha/3)*pow(energy[k][fi]-energy[ks][fi],3)*(pow(dipole_transition_x[k][ks][fi],2)+pow(dipole_transition_y[k][ks][fi],2)+pow(dipole_transition_z[k][ks][fi],2)));
	
	}
	file.close();
	

//		cout<<setprecision(20)<<n1<<n-n1-abs(m)-1<<m<<"---"<<n1s<<ns-n1s-abs(ms)-1<<ms<<"  dx="<<dipole_transition_x[k][ks][1]<<"  dy="<<dipole_transition_y[k][ks][1]<<"  dz="<<dipole_transition_z[k][ks][1]<<"\n";
		}


}






void LasStripExternalEffects::GetDipoleTransition(int k,int ks,double& relax,double& mu_x,double& mu_y,double& mu_z){
	
	int i=(int)(Ez_stat/delta_F);
	double c=(Ez_stat-i*delta_F)/delta_F;
	
mu_x=dipole_transition_x[k][ks][i]+c*(dipole_transition_x[k][ks][i+1]-dipole_transition_x[k][ks][i]);
mu_y=dipole_transition_y[k][ks][i]+c*(dipole_transition_y[k][ks][i+1]-dipole_transition_y[k][ks][i]);
mu_z=dipole_transition_z[k][ks][i]+c*(dipole_transition_z[k][ks][i+1]-dipole_transition_z[k][ks][i]);

relax=fabs(gamma_spontaneous_relax[k][ks][i]+c*(gamma_spontaneous_relax[k][ks][i+1]-gamma_spontaneous_relax[k][ks][i]));

}




void LasStripExternalEffects::GetEnergyAutoionization(int k,double& E, double& Gamma){
	
	int i=(int)(Ez_stat/delta_F);
	double c=(Ez_stat-i*delta_F)/delta_F;

	
	if ((gamma_autoionization[k][i]>1e-20)&&(gamma_autoionization[k][i+1]>1e-20))
	Gamma=gamma_autoionization[k][i]*pow(gamma_autoionization[k][i+1]/gamma_autoionization[k][i],c);
	else	Gamma=0;
	
	E=energy[k][i]+c*(energy[k][i+1]-energy[k][i]);
	

}







void 	LasStripExternalEffects::GetOwnLaserField(double x, double y, double z,double t, double& E_x, double& E_y, double& E_z,double& B_x, double& B_y, double& B_z)
	
{
	//Some constant
	double eps_0=8.854187817e-012;
	double mu_0=1.2566370614e-006;
	
	double z0=Laser_lambda/(OrbitConst::PI*Laser_half_angle*Laser_half_angle);
	double k=2*OrbitConst::PI/Laser_lambda;
	
	//non-normized Direction of Poiting vector
	double	nx=x*z;
	double	ny=y*z;
	double	nz=z*z+z0*z0-(x*x+y*y)*(0.5+z*z/(z*z+z0*z0))-(z0/k);
	
	double normn=1./sqrt(nx*nx+ny*ny+nz*nz);
	nx_lab=nx*normn;
	ny_lab=ny*normn;
	nz_lab=nz*normn;
	
	//non-normized direction of H-field 
	double nBx=0;
	double nBy=nz;
	double nBz=-ny;
	
	//non-normized direction of E-field 
	double nEx=ny*ny+nz*nz;
	double nEy=-nx*ny;
	double nEz=-nx*nz;
	
	//Absolute value of Poiting vector
	double absS=(k*LaserPower/(OrbitConst::PI*z0))*(z0*z0/(z0*z0+z*z))*exp(-k*z0*(x*x+y*y)/(z0*z0+z*z));
	//Absolute value of E vector
	double absE=sqrt(2*absS*sqrt(mu_0/eps_0));
	//Absolute value of H vector
	double absB=absE*sqrt(eps_0*mu_0);
	
	//components of H-field
	B_x=absB*nBx/sqrt(nEx);
	B_y=absB*nBy/sqrt(nEx);
	B_z=absB*nBz/sqrt(nEx);
	
	//components of E-field
	double a=1./sqrt(nEx*nEx+nEy*nEy+nEz*nEz);
	
	E_x=absE*nEx*a;
	E_y=absE*nEy*a;
	E_z=absE*nEz*a;
	

			}
	

void 	LasStripExternalEffects::GetLabLaserField(double x, double y, double z, double t, double& E_x, double& E_y, double& E_z,double& B_x, double& B_y, double& B_z){
	
	/**Here must be the code of conversion of fields and coordinates from laser frame to lab frame**/
	
	GetOwnLaserField(x,y,z,t,E_x,E_y,E_z,B_x,B_y,B_z);
	
	/**Here must be the code of conversion of fields and coordinates from laser frame to lab frame**/
	
}








void	LasStripExternalEffects::GetFrameParticleParameters(int i, double t_step, Bunch* bunch)	{
	
	
	double** xyz = bunch->coordArr();		//definition of coordinates an impulses of particles
	
	double ta=2.418884326505e-17;			//atomic unit of time
	double Ea=5.14220642e11;				//Atomic unit of electric field
	double mp=bunch->getMass();

	
//next four lines calculate relyativistic factor-Gamma
double p2=xyz[i][1]*xyz[i][1]+xyz[i][3]*xyz[i][3]+xyz[i][5]*xyz[i][5];
double mp2=mp*mp;
double E=sqrt(mp2+p2);
double gamma=E/mp;
double beta_cos_pn=(xyz[i][1]*nx_lab+xyz[i][3]*ny_lab+xyz[i][5]*nz_lab)/E;


//next nine lines convert vector of laser and static field in basis where static fiels is parallel to z axes	(in natural (non-atomic) units)
double Estat=sqrt(Ex_stat*Ex_stat+Ey_stat*Ey_stat+Ez_stat*Ez_stat);	
double Elas2=Ex_las*Ex_las+Ey_las*Ey_las+Ez_las*Ez_las;
double ElasEsat=Ex_las*Ex_stat+Ey_las*Ey_stat+Ez_las*Ez_stat;
Ex_las=ElasEsat/Estat;
Ey_las=0;
Ez_las=sqrt(Elas2-Ex_las*Ex_las);		Ex_las/=Ea;	Ey_las/=Ea;	Ez_las/=Ea;
Ez_stat=sqrt(Ex_stat*Ex_stat+Ey_stat*Ey_stat+Ez_stat*Ez_stat)/Ea;
Ex_stat=0;	// In fact parameters Ex_stat and Ey_stat are not used and they simply can be deleted
Ey_stat=0;	// In fact parameters Ex_stat and Ey_stat are not used and they simply can be deleted

part_t_step=t_step/gamma/ta;	//time step in frame of particle (in atomic units)
t_part=time(i);					//time  in frame of particle (in atomic units). Attention.!!! It is supposed that each particle moves with constant velocity
phasa_part=phasa(i);			//phasa of laser field in frame of particle
omega_part=2*OrbitConst::PI*OrbitConst::c*gamma*(1-beta_cos_pn)*ta/Laser_lambda;		// frequensy of laser in particle frame (in atomic units)


//cout<<"Ex_las="<<Ex_las<<"\n";
//cout<<"Ey_las="<<Ey_las<<"\n";
//cout<<"Ez_las="<<Ez_las<<"\n\n";


			

//Ez_stat	Ex_las	Ey_las	Ez_las	t_part	omega_part	part_t_step

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////froissar stora test      froissar stora test      froissar stora test      froissar stora test    ///////////////   
//next eight lines are testing parameters of laser electric field (in atomic units)

double Rabi=1e+12*ta;					//Rabi frequensy (in atomic units)
double Elas=64*Rabi*sqrt(2.)/27;		//Amplitude of laser field (in atonic units)
Ex_las=0.5773502691896258*Elas;			//components of laser field (in atomic units)
Ey_las=0.5773502691896258*Elas;	
Ez_las=0.5773502691896258*Elas;

Ex_las=Elas;
Ey_las=0;
Ez_las=0;


Ez_stat=0.00005;


double gamma_sweep=-OrbitConst::PI*Rabi*Rabi/2/log(1-0.9);

omega_part=0*0.44421741823240407099+4/9.+gamma_sweep*(t_part-10000*(2*OrbitConst::PI/Rabi)/2);				// frequensy of laser in particle frame (in atomic units) 

/////////////froissar stora test      froissar stora test      froissar stora test      froissar stora test ///////////          
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	


//cout<<setprecision(20)<<"gamma="<<gamma_sweep<<"\n";


	
}














void LasStripExternalEffects::AmplSolver4step(int i, Bunch* bunch)	{
	
	
	
	
	

	double** xyz = bunch->coordArr();
	tcomplex z,mu_Elas;
	double dt,mu_x,mu_y,mu_z;
		

		for(int n=1; n<levels+1;n++)	GetEnergyAutoionization(n,E_i[n], Gamma_i[n]);
		
//		 cout<<setprecision(20)<<fabs(E_i[1]-E_i[7])<<"\n";
			
		for(int n=2; n<levels+1;n++)
		for(int m=1; m<n;m++)	{
			
			GetDipoleTransition(n,m,gamma_ij[n][m],mu_x,mu_y,mu_z);
			mu_Elas=mu_x*Ex_las+J*mu_y*Ey_las+mu_z*Ez_las;							
			cond[n][m]=fabs(fabs(E_i[n]-E_i[m])-omega_part)<1000*abs(mu_Elas);
			
			if (cond[n][m])	{
				
			exp_mu_El[n][m][1]=exp(J*(t_part*fabs(E_i[n]-E_i[m])-phasa_part))*mu_Elas;
			exp_mu_El[n][m][2]=exp(J*((t_part+part_t_step/2)*fabs(E_i[n]-E_i[m])-(phasa_part+omega_part*part_t_step/2)))*mu_Elas;
			exp_mu_El[n][m][3]=exp_mu_El[n][m][2];
			exp_mu_El[n][m][4]=exp(J*((t_part+part_t_step)*fabs(E_i[n]-E_i[m])-(phasa_part+omega_part*part_t_step)))*mu_Elas;	
			}
			
		}

					
		

		for(int j=1; j<5; j++)
		for(int n=1; n<levels+1;n++)
		for(int m=1; m<n+1;m++)	{	
			
			if (j==4)	dt=part_t_step;	else dt=part_t_step/2.;

			k_RungeKutt[n][m][j]*=0.;
			for(int k=1;k<n;k++)			if (cond[n][k]) 	k_RungeKutt[n][m][j]+=exp_mu_El[n][k][j]*(dm(i,k,m)+k_RungeKutt[k][m][j-1]*dt);	
			for(int k=n+1;k<levels+1;k++)	if (cond[k][n]) 	k_RungeKutt[n][m][j]+=conj(exp_mu_El[k][n][j])*(dm(i,k,m)+k_RungeKutt[k][m][j-1]*dt);
			for(int k=1;k<m;k++)			if (cond[m][k]) 	k_RungeKutt[n][m][j]-=conj(exp_mu_El[m][k][j])*(dm(i,n,k)+k_RungeKutt[n][k][j-1]*dt);
			for(int k=m+1;k<levels+1;k++)	if (cond[k][m]) 	k_RungeKutt[n][m][j]-=exp_mu_El[k][m][j]*(dm(i,n,k)+k_RungeKutt[n][k][j-1]*dt);
			k_RungeKutt[n][m][j]*=J/2.;
			



			if(n==m){
			for(int k=m+1;k<levels+1;k++)	k_RungeKutt[n][m][j]+=gamma_ij[k][m]*(dm(i,k,k)+k_RungeKutt[k][k][j-1]*dt);
			for(int k=1;k<m;k++)			k_RungeKutt[n][m][j]-=gamma_ij[m][k]*(dm(i,n,m)+k_RungeKutt[n][m][j-1]*dt);
			}
			
			if(n!=m){
			for(int k=1;k<m;k++)			k_RungeKutt[n][m][j]-=gamma_ij[m][k]*(dm(i,n,m)+k_RungeKutt[n][m][j-1]*dt)/2.;
			for(int k=1;k<n;k++)			k_RungeKutt[n][m][j]-=gamma_ij[n][k]*(dm(i,n,m)+k_RungeKutt[n][m][j-1]*dt)/2.;
			}
		
			
			k_RungeKutt[n][m][j]-=(Gamma_i[n]+Gamma_i[m])*(dm(i,n,m)+k_RungeKutt[n][m][j-1]*dt)/2.;


			k_RungeKutt[m][n][j]=conj(k_RungeKutt[n][m][j]);
			
		}

		
		
		
		
		
		
	for(int n=1;n<levels+1;n++)
	for(int m=1;m<n+1;m++)	{
		
		
		z=(k_RungeKutt[n][m][1]+2.*k_RungeKutt[n][m][2]+2.*k_RungeKutt[n][m][3]+k_RungeKutt[n][m][4])/6.;
		Re(i,n,m)+=z.real()*part_t_step;
		Im(i,n,m)+=z.imag()*part_t_step;
		
		Re(i,m,n)=Re(i,n,m);
		Im(i,m,n)=-Im(i,n,m);

		
	}

	time(i)+=part_t_step;
	phasa(i)+=omega_part*part_t_step;
	
	x0(i)=xyz[i][0];
	y0(i)=xyz[i][2];
	z0(i)=xyz[i][4];
	
	px0(i)=xyz[i][1];
	py0(i)=xyz[i][3];
	pz0(i)=xyz[i][5];
		
	
}




