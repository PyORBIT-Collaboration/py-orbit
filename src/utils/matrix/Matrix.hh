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
	class Matrix
	{
		public:
		
			Matrix(int n, int m);
			Matrix(Matrix* mtrx);
			~Matrix();
		
			double** getArray();
			int rows();
			int columns();
			double& value(int i, int j);
			
			int copyTo(Matrix* m_child);
			void transpose();
			int unit();
			void zero();
			int add(double val);
			int add(Matrix* mtrx);
			int mult(double val);
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
