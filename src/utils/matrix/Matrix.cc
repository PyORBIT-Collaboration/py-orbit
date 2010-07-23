//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//   Matrix.cc
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

#include "orbit_mpi.hh"
#include "Matrix.hh"
#include "BufferStore.hh"

#include <cstdlib>

using namespace OrbitUtils;

Matrix::Matrix(int n_in, int m_in)
{
	n = n_in;
	m = m_in;
		
	a = (double** ) malloc (sizeof(double*)*n);
	for(int i=0; i < n; i++){
		a[i] = (double* ) malloc (sizeof(double)*m);
	}	
	zero();
}

Matrix::Matrix(Matrix* mtrx)
{
	n = mtrx->rows();
	m = mtrx->columns();
	a = (double** ) malloc (sizeof(double*)*n);
	for(int i=0; i < n; i++){
		a[i] = (double* ) malloc (sizeof(double)*m);
	}	
	mtrx->copyTo(this);
}

Matrix::~Matrix()
{
	for(int i=0; i < n; i++){
    free(a[i]);
  }
  free(a);
}

void Matrix::zero(){
	for(int i = 0; i < n; i++){
		for(int j = 0; j < m; j++){
			a[i][j] = 0.;
		}
	}
}

int Matrix::unit(){
	if(n != m){
		ORBIT_MPI_Finalize("Matrix: unit operation only for a square matrix.");
		return 0;
	}
	for(int i = 0; i < n; i++){
		for(int j = 0; j < m; j++){
		  a[i][j] = 0.;
		}
	}
	for(int i = 0; i < n; i++){
		a[i][i] = 1.0;
	}
	return 1;
}

int Matrix::copyTo(Matrix* m_child){
	if(n != m_child->rows() || m != m_child->columns()){
		ORBIT_MPI_Finalize("Matrix: You try to copy Matrices with wrong size.");
		return 0;
	}
	double** m_child_arr = m_child->getArray();
	for(int i = 0; i < n; i++){
		for(int j = 0; j < m; j++){
			m_child_arr[i][j] = a[i][j];
		}
	}
	return 1;
}

void Matrix::transpose(){
	if(n == m){
		double val = 0.;
		for(int i = 0; i < n; i++){
			for(int j = 0; j < i; j++){
				val = a[i][j];
				a[i][j] = a[j][i];
				a[j][i] = val;
			}
		}
		return;
	}
	//make new arrays
  double** arr = new double*[m];
	for(int i=0; i < m; i++){
		arr[i] = new double[n];
	}		
	//copy
	for(int i = 0; i < n; i++){
		for(int j = 0; j < m; j++){
			arr[j][i] = a[i][j];
		}
	} 
	//remove old arrays
	for(int i=0; i < n; i++){
    delete [] a[i];
  }
  delete [] a;
	int nn = n;
	n = m;
	m = nn;
	a = arr;
}

double& Matrix::value(int i, int j){ return a[i][j];}

double** Matrix::getArray(){ return a;}

int Matrix::rows(){ return n;}

int Matrix::columns(){ return m;}

int Matrix::add(double val){
	for(int i = 0; i < n; i++){
		for(int j = 0; j < m; j++){
		  a[i][j] += val;
		}
	}
	return 1;
}

int Matrix::add(Matrix* mtrx){
	if(n != mtrx->rows() || m != mtrx->columns()){
		ORBIT_MPI_Finalize("Matrix: You try to add a matrix with wrong size.");
		return 0;
	}
	double** arr = mtrx->getArray();
	for(int i = 0; i < n; i++){
		for(int j = 0; j < m; j++){
		  a[i][j] += arr[i][j];
		}
	}
	return 1;
}

int Matrix::mult(double val){
	for(int i = 0; i < n; i++){
		for(int j = 0; j < m; j++){
		  a[i][j] *= val;
		}
	}
	return 1;
}

int Matrix::mult(Matrix* mtrx){
	if(m != mtrx->rows()){
		ORBIT_MPI_Finalize("Matrix: You try to multiply Matrices with wrong size.");
		return 0;
	}
	double** arr = mtrx->getArray();
	int buff_index = -1;
	int m_new = mtrx->columns();
  double* tmp = BufferStore::getBufferStore()->getFreeDoubleArr(buff_index,m_new);
	for(int i = 0; i < n; i++){
		for(int j = 0; j < m_new; j++){
			tmp[j] = 0.;
			for(int jj = 0; jj < m; jj++){
				tmp[j] += a[i][jj]*arr[jj][j];
			}
		}
		a[i] = (double *) realloc(a[i],m_new*sizeof(double));
		for(int j = 0; j < m_new; j++){
		  a[i][j] = tmp[j];
		}
	}
	m = m_new;
	BufferStore::getBufferStore()->setUnusedDoubleArr(buff_index);
	return 1;
}


