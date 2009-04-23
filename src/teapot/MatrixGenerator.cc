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
		step_arr[0] = 0.000001;
		step_arr[1] = 0.000001;
		step_arr[2] = 0.000001;
		step_arr[3] = 0.000001;
		step_arr[4] = 0.000001;
		step_arr[5] = 0.000001;
		step_reduce = 20.;
	}
	
	MatrixGenerator::~MatrixGenerator(){
		delete [] step_arr;
	}	
	
	double& MatrixGenerator::step(int index){
		return step_arr[index];
	}
	
	void MatrixGenerator::initBunch(Bunch* bunch){
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
			ORBIT_MPI_Finalize("MatrixGenerator::calculateMatrix: Bunch should have 7 macro-particles.");
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
	
}  //end of namespace teapot_base



