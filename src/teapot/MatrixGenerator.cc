/////////////////////////////////////////////////////////////////////////////
//
// FILE NAME
//   atrixGenerator.cc
//
// AUTHOR
//   Andrei Shishlo
//
// TIME
//  03/18/2008
//
// DESCRIPTION
//   A generator of the linear transport matrix from the bunch tracking.
//   It solves the system of equations with respect to a,b,c:
//   y0 = a*x0^2 + b*x0 + c
//   y1 = a*x1^2 + b*x1 + c
//   y2 = a*x2^2 + b*x2 + c
//   for x0 = 0., we need only b coefficient.
//   Solution:
//   b = ((y1-y0)*x2^2 - (y2-y0)*x1^2)/(x1*x2*(x2 -x1))
//
/////////////////////////////////////////////////////////////////////////////
#include "orbit_mpi.hh"
#include "MatrixGenerator.hh"

using namespace OrbitUtils;

namespace teapot_base{
	
	MatrixGenerator::MatrixGenerator(){
		step_arr = new double[6];
		step_arr_init = new double[6];
		step_arr_init[0] = 0.000001;
		step_arr_init[1] = 0.000001;
		step_arr_init[2] = 0.000001;
		step_arr_init[3] = 0.000001;
		step_arr_init[4] = 0.000001;
		step_arr_init[5] = 0.000001;
		for(int i = 0; i < 6; i++){
			step_arr[i] = step_arr_init[i];
		}		
		step_reduce = 20.;
	}
	
	MatrixGenerator::~MatrixGenerator(){
		delete [] step_arr;
		delete [] step_arr_init;
	}	
	
	double& MatrixGenerator::step(int index){
		return step_arr_init[index];
	}
	
	void MatrixGenerator::initBunch(Bunch* bunch){
		for(int i = 0; i < 6; i++){
			step_arr[i] = step_arr_init[i];
		}
		double kin_energy = bunch->getSyncPart()->getEnergy();
		step_arr[5] = kin_energy*step_arr[5];
		bunch->deleteAllParticles();
		bunch->addParticle(0.,0.,0.,0.,0.,0.);
		bunch->addParticle(step_arr[0]/step_reduce,0.,0.,0.,0.,0.);
		bunch->addParticle(0.,step_arr[1]/step_reduce,0.,0.,0.,0.);
		bunch->addParticle(0.,0.,step_arr[2]/step_reduce,0.,0.,0.);
		bunch->addParticle(0.,0.,0.,step_arr[3]/step_reduce,0.,0.);
		bunch->addParticle(0.,0.,0.,0.,step_arr[4]/step_reduce,0.);
		bunch->addParticle(0.,0.,0.,0.,0.,step_arr[5]/step_reduce);		
		bunch->addParticle(step_arr[0],0.,0.,0.,0.,0.);
		bunch->addParticle(0.,step_arr[1],0.,0.,0.,0.);
		bunch->addParticle(0.,0.,step_arr[2],0.,0.,0.);
		bunch->addParticle(0.,0.,0.,step_arr[3],0.,0.);
		bunch->addParticle(0.,0.,0.,0.,step_arr[4],0.);
		bunch->addParticle(0.,0.,0.,0.,0.,step_arr[5]);
	}
	
	void MatrixGenerator::calculateMatrix(Bunch* bunch,Matrix* mtrx){
		if(mtrx->rows() != mtrx->columns() || mtrx->rows() < 6){
			ORBIT_MPI_Finalize("MatrixGenerator::calculateMatrix: Matrix has a wrong size.");
		}
		if(bunch->getSize() != 13){
			ORBIT_MPI_Finalize("MatrixGenerator::calculateMatrix: Bunch should have 13 macro-particles.");
		}
		double x1,x2,y1,y2,y0;
		double** coord_arr = bunch->coordArr();
		double** arr = mtrx->getArray();
		for(int i = 0; i < 6; i++){
			for(int j = 0; j < 6; j++){
				x1 = step_arr[j]/step_reduce;
				x2 = step_arr[j];
				y0 = coord_arr[0][j];
				y1 = coord_arr[i+1][j];
				y2 = coord_arr[i+1+6][j];
				arr[j][i] = ((y1-y0)*x2*x2 - (y2-y0)*x1*x1)/(x1*x2*(x2 -x1));
			}
		}
		if(mtrx->rows() == 7){
			for(int i = 0; i < 6; i++){
				arr[i][6] = coord_arr[0][i];
				arr[6][i] = 0.;
			}
			arr[6][6] = 1.0;
		}
	}
	
	
	void MatrixGenerator::initBunchForChromaticityCoeff(Bunch* bunch){
		for(int i = 0; i < 6; i++){
			step_arr[i] = step_arr_init[i];
		}
		double kin_energy = bunch->getSyncPart()->getEnergy();
		step_arr[5] = kin_energy*step_arr[5];
		bunch->deleteAllParticles();
		//initial point
		bunch->addParticle(0.,0.,0.,0.,0.,0.);
		bunch->addParticle(0.,0.,0.,0.,0.,step_arr[5]/step_reduce);	
		bunch->addParticle(0.,0.,0.,0.,0.,step_arr[5]);
		//for x-dE coefficient
		bunch->addParticle(step_arr[0]/step_reduce,0.,0.,0.,0.,0.);
		bunch->addParticle(step_arr[0]/step_reduce,0.,0.,0.,0.,step_arr[5]/step_reduce);
		bunch->addParticle(step_arr[0]/step_reduce,0.,0.,0.,0.,step_arr[5]);
		bunch->addParticle(step_arr[0],0.,0.,0.,0.,0.);
		bunch->addParticle(step_arr[0],0.,0.,0.,0.,step_arr[5]/step_reduce);
		bunch->addParticle(step_arr[0],0.,0.,0.,0.,step_arr[5]);
		//for xp-dE coefficient
		bunch->addParticle(0.,step_arr[1]/step_reduce,0.,0.,0.,0.);
		bunch->addParticle(0.,step_arr[1]/step_reduce,0.,0.,0.,step_arr[5]/step_reduce);
		bunch->addParticle(0.,step_arr[1]/step_reduce,0.,0.,0.,step_arr[5]);
		bunch->addParticle(0.,step_arr[1],0.,0.,0.,0.);
		bunch->addParticle(0.,step_arr[1],0.,0.,0.,step_arr[5]/step_reduce);
		bunch->addParticle(0.,step_arr[1],0.,0.,0.,step_arr[5]);
		//for y-dE coefficient
		bunch->addParticle(0.,0.,step_arr[2]/step_reduce,0.,0.,0.);
		bunch->addParticle(0.,0.,step_arr[2]/step_reduce,0.,0.,step_arr[5]/step_reduce);
		bunch->addParticle(0.,0.,step_arr[2]/step_reduce,0.,0.,step_arr[5]);
		bunch->addParticle(0.,0.,step_arr[2],0.,0.,0.);
		bunch->addParticle(0.,0.,step_arr[2],0.,0.,step_arr[5]/step_reduce);
		bunch->addParticle(0.,0.,step_arr[2],0.,0.,step_arr[5]);
		//for yp-dE coefficient
		bunch->addParticle(0.,0.,0.,step_arr[3]/step_reduce,0.,0.);
		bunch->addParticle(0.,0.,0.,step_arr[3]/step_reduce,0.,step_arr[5]/step_reduce);
		bunch->addParticle(0.,0.,0.,step_arr[3]/step_reduce,0.,step_arr[5]);
		bunch->addParticle(0.,0.,0.,step_arr[3],0.,0.);
		bunch->addParticle(0.,0.,0.,step_arr[3],0.,step_arr[5]/step_reduce);
		bunch->addParticle(0.,0.,0.,step_arr[3],0.,step_arr[5]);		
	}	
	
	void MatrixGenerator::calculateChromaticityCoeff(Bunch* bunch,
		                            double& coeff_x_dE, double& coeff_xp_dE,
		                            double& coeff_y_dE, double& coeff_yp_dE){
		if(bunch->getSize() != 27){
			ORBIT_MPI_Finalize("MatrixGenerator::calculateChromaticityCoeff: Bunch should have 27 macro-particles.");
		}
		double x1,x2,y1,y2,y0;
		double** coord_arr = bunch->coordArr();
		//calculation of dE coeff. at 0 point for (x,xp,y,yp)
		x1 = step_arr[5]/step_reduce;
		x2 = step_arr[5];
		y0 = coord_arr[0][0];
		y1 = coord_arr[1][0];
		y2 = coord_arr[2][0];		
		double coeff_x_dE_0 = ((y1-y0)*x2*x2 - (y2-y0)*x1*x1)/(x1*x2*(x2 -x1));
		y0 = coord_arr[0][1];
		y1 = coord_arr[1][1];
		y2 = coord_arr[2][1];		
		double coeff_xp_dE_0 = ((y1-y0)*x2*x2 - (y2-y0)*x1*x1)/(x1*x2*(x2 -x1));
		y0 = coord_arr[0][2];
		y1 = coord_arr[1][2];
		y2 = coord_arr[2][2];		
		double coeff_y_dE_0 = ((y1-y0)*x2*x2 - (y2-y0)*x1*x1)/(x1*x2*(x2 -x1));
		y0 = coord_arr[0][3];
		y1 = coord_arr[1][3];
		y2 = coord_arr[2][3];		
		double coeff_yp_dE_0 = ((y1-y0)*x2*x2 - (y2-y0)*x1*x1)/(x1*x2*(x2 -x1));
		//calculation of x-dE coefficients
		y0 = coord_arr[3][0];
		y1 = coord_arr[4][0];
		y2 = coord_arr[5][0];		
		double coeff_x_dE_1 = ((y1-y0)*x2*x2 - (y2-y0)*x1*x1)/(x1*x2*(x2 -x1));
		y0 = coord_arr[6][0];
		y1 = coord_arr[7][0];
		y2 = coord_arr[8][0];		
		double coeff_x_dE_2 = ((y1-y0)*x2*x2 - (y2-y0)*x1*x1)/(x1*x2*(x2 -x1));
		//calculation of xp-dE coefficients
		y0 = coord_arr[9][1];
		y1 = coord_arr[10][1];
		y2 = coord_arr[11][1];		
		double coeff_xp_dE_1 = ((y1-y0)*x2*x2 - (y2-y0)*x1*x1)/(x1*x2*(x2 -x1));
		y0 = coord_arr[12][1];
		y1 = coord_arr[13][1];
		y2 = coord_arr[14][1];		
		double coeff_xp_dE_2 = ((y1-y0)*x2*x2 - (y2-y0)*x1*x1)/(x1*x2*(x2 -x1));
		//calculation of y-dE coefficients
		y0 = coord_arr[15][2];
		y1 = coord_arr[16][2];
		y2 = coord_arr[17][2];		
		double coeff_y_dE_1 = ((y1-y0)*x2*x2 - (y2-y0)*x1*x1)/(x1*x2*(x2 -x1));
		y0 = coord_arr[18][2];
		y1 = coord_arr[19][2];
		y2 = coord_arr[20][2];		
		double coeff_y_dE_2 = ((y1-y0)*x2*x2 - (y2-y0)*x1*x1)/(x1*x2*(x2 -x1));
		//calculation of yp-dE coefficients
		y0 = coord_arr[21][3];
		y1 = coord_arr[22][3];
		y2 = coord_arr[23][3];		
		double coeff_yp_dE_1 = ((y1-y0)*x2*x2 - (y2-y0)*x1*x1)/(x1*x2*(x2 -x1));
		y0 = coord_arr[24][3];
		y1 = coord_arr[25][3];
		y2 = coord_arr[26][3];		
		double coeff_yp_dE_2 = ((y1-y0)*x2*x2 - (y2-y0)*x1*x1)/(x1*x2*(x2 -x1));
		//----------------------------------------
		//combination of 3 coefficients
		x1 = step_arr[0]/step_reduce;
		x2 = step_arr[0];
		y0 = coeff_x_dE_0;
		y1 = coeff_x_dE_1;
		y2 = coeff_x_dE_2;		
		coeff_x_dE = ((y1-y0)*x2*x2 - (y2-y0)*x1*x1)/(x1*x2*(x2 -x1));		
		x1 = step_arr[1]/step_reduce;
		x2 = step_arr[1];
		y0 = coeff_xp_dE_0;
		y1 = coeff_xp_dE_1;
		y2 = coeff_xp_dE_2;		
		coeff_xp_dE = ((y1-y0)*x2*x2 - (y2-y0)*x1*x1)/(x1*x2*(x2 -x1));		
		x1 = step_arr[2]/step_reduce;
		x2 = step_arr[2];
		y0 = coeff_y_dE_0;
		y1 = coeff_y_dE_1;
		y2 = coeff_y_dE_2;		
		coeff_y_dE = ((y1-y0)*x2*x2 - (y2-y0)*x1*x1)/(x1*x2*(x2 -x1));		
		x1 = step_arr[3]/step_reduce;
		x2 = step_arr[3];
		y0 = coeff_yp_dE_0;
		y1 = coeff_yp_dE_1;
		y2 = coeff_yp_dE_2;		
		coeff_yp_dE = ((y1-y0)*x2*x2 - (y2-y0)*x1*x1)/(x1*x2*(x2 -x1));
	}	
	
}  //end of namespace teapot_base



