#ifndef __MATRIX_OPERATIONS_H_
#define __MATRIX_OPERATIONS_H_

#include "Matrix.hh"
#include "PhaseVector.hh"
#include "Bunch.hh"

namespace OrbitUtils{
	
	/**
	  The collection of the static methods for matrices. 
	*/
		
	class MatrixOperations{
	public:
		
		/** Inverts the 2D array as a matrix. */
		static int invert(double **a, int n);
		
		/** Inversts the matrix in place. */
		static int invert(Matrix* matrix);
		
		/** Multiplies a vector with the transposed matrix. */
		static int mult(PhaseVector* v, Matrix* mtrx, PhaseVector* v_res);
		
		/** Multiplies a vector with the matrix. */
		static int mult(Matrix* mtrx, PhaseVector* v, PhaseVector* v_res);
		
		/** 
		  Tracks the bunch through the transport matrix. 
			The matrix should be 6x6 or 7x7.
		*/
		static void track(Bunch* bunch,Matrix* mtrx);		
	};
};

#endif
