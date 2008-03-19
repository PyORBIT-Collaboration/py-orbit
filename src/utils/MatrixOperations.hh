#ifndef __MATRIX_OPERATIONS_H_
#define __MATRIX_OPERATIONS_H_

#include "Matrix.hh"
#include "PhaseVector.hh"
#include "Bunch.hh"

namespace OrbitUtils{
	class MatrixOperations{
	public:
		static int invert(double **a, int n);
		static int invert(Matrix* matrix);
		static int mult(PhaseVector* v, Matrix* mtrx, PhaseVector* v_res);
		static int mult(Matrix* mtrx, PhaseVector* v, PhaseVector* v_res);
		static void track(Bunch* bunch,Matrix* mtrx);		
	};
};

#endif
