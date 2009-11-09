//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//   Matrix.hh
//
// AUTHOR
//    A. Shishlo
//
// CREATED
//    03/12/2008
//
// DESCRIPTION
//    A class for plain NxM matrices with double values. 
//
///////////////////////////////////////////////////////////////////////////

#ifndef PLAIN_MATRIX_H
#define PLAIN_MATRIX_H

namespace OrbitUtils{
	
	/** A class for plain NxM matrices with double values. */
	
	class Matrix
	{
		public:
		
			/** Constructor. */
			Matrix(int n, int m);
			
			/** Destructor */
			Matrix(Matrix* mtrx);
			~Matrix();
		
			/** Returns the pointer to the two dimensional array. */
			double** getArray();
			
			/** Returns the number of rows. */
			int rows();
			
			/** Returns the number of columns. */
			int columns();
			double& value(int i, int j);
			
			/** Copies the matrix to another. */
			int copyTo(Matrix* m_child);
			
			/** Transposes the matrix in place. */
			void transpose();
			
			/** Makes the unit matrix in place. */
			int unit();
			
			/** Sets to zero matrix elements */
			void zero();
			
			/** Adds the value to each element of the matrix. */
			int add(double val);
			
			/** Adds a matrix. */ 
			int add(Matrix* mtrx);
			
			/** Multiplies each element of the matrix by the value. */
			int mult(double val);
			
			/** Multiplies two matrices. */ 
			int mult(Matrix* mtrx);
		
		private:
		
			int n,m;
			
			double** a;
		
	};
};
///////////////////////////////////////////////////////////////////////////
//
// END OF FILE
//
///////////////////////////////////////////////////////////////////////////

#endif
