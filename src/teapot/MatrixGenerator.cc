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
//
/////////////////////////////////////////////////////////////////////////////
#include "orbit_mpi.hh"
#include "MatrixGenerator.hh"

using namespace OrbitUtils;

namespace teapot_base{
	
	MatrixGenerator::MatrixGenerator(){
		step_arr = new double[6];
		step_arr[0] = 0.0001;
		step_arr[1] = 0.0001;
		step_arr[2] = 0.0001;
		step_arr[3] = 0.0001;
		step_arr[4] = 0.0001;
		step_arr[5] = 0.0001;
	}
	
	MatrixGenerator::~MatrixGenerator(){
		delete [] step_arr;
	}	
	
	double& MatrixGenerator::step(int index){
		return step_arr[index];
	}
	
	void MatrixGenerator::initBunch(Bunch* bunch){
		if(bunch->getSize() != 7){
			bunch->deleteAllParticles();
			for(int i = 0; i < 7; i++){
				bunch->addParticle(0.,0.,0.,0.,0.,0.);
			}
		}
		double tmp[6];
		for(int i = 0; i < 6; i++){
			tmp[i] = 0.;
		}
		double** coord_arr = bunch->coordArr();
		for(int i = 0; i < 6; i++){
			coord_arr[0][i] = 0.;
			tmp[i] = step_arr[i];
			for(int j = 0; j < 6; j++){
			coord_arr[i+1][j] = tmp[j];
			}
			tmp[i] = 0.;
		}
	}
	
	void MatrixGenerator::calculateMatrix(Bunch* bunch,Matrix* mtrx){
		if(mtrx->rows() != mtrx->columns() || mtrx->rows() < 6){
			ORBIT_MPI_Finalize("MatrixGenerator::calculateMatrix: Matrix has a wrong size.");
		}
		if(bunch->getSize() != 7){
			ORBIT_MPI_Finalize("MatrixGenerator::calculateMatrix: Bunch should have 7 macro-particles.");
		}
		double** coord_arr = bunch->coordArr();
		double** arr = mtrx->getArray();
		for(int i = 0; i < 6; i++){
			for(int j = 0; j < 6; j++){
				arr[i][j] = (coord_arr[j+1][i] - coord_arr[0][i])/step_arr[j];
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



