/////////////////////////////////////////////////////////////////////////////
//
// FILE NAME
//   MatrixGenerator.hh
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
#ifndef TEAPOT_BASE_MATRIX_GENERATOR_H
#define TEAPOT_BASE_MATRIX_GENERATOR_H

#include "Bunch.hh"
#include "Matrix.hh"

using namespace OrbitUtils;

namespace teapot_base{
	class MatrixGenerator
	{
	public:
		MatrixGenerator();
		~MatrixGenerator();
		
		double& step(int index);
		void initBunch(Bunch* bunch);
		void calculateMatrix(Bunch* bunch,Matrix* mtrx);
		
	private:
		double* step_arr;		
	};
	
}  //end of namespace teapot_base

#endif  //TEAPOT_BASE_MATRIX_GENERATOR_H
