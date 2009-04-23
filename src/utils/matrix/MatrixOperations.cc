#include "orbit_mpi.hh"

#include "MatrixOperations.hh"
#include "BufferStore.hh"

#include <cmath>

using namespace OrbitUtils;

int MatrixOperations::invert(double **a, int n){
	
	//   Taken from Nrecipes and slightly modified.
	//   Get matrix A(nxn) and transform it into A^(-1).
	//   Returns -1 if A^(-1) doesn't exist.
	
	int i,icol,irow,j,k,l,ll;
	double big,dum,pivinv,temp;
	
	icol = 0;
	irow = 0;
	
	// indexr and indexc track column permutation
	int buff_index0 = -1;
	int buff_index1 = -1;
	int buff_index2 = -1;
	int *indxc= BufferStore::getBufferStore()->getFreeIntArr(buff_index0,n);
	int *indxr= BufferStore::getBufferStore()->getFreeIntArr(buff_index1,n);
	int *ipiv = BufferStore::getBufferStore()->getFreeIntArr(buff_index2,n);
	
	for (j=0;j<n;j++) ipiv[j]=0;
	for (i=0;i<n;i++) {
		big=0.0;
		// Looking for pivot
		for (j=0;j<n;j++)
			if (ipiv[j] != 1)
			for (k=0;k<n;k++) {
				if (ipiv[k] == 0) {
					if (fabs(a[j][k]) >= big) {
						big=fabs(a[j][k]);
						irow=j;
						icol=k;
					}
 		    } else if (ipiv[k] > 1) {
					BufferStore::getBufferStore()->setUnusedIntArr(buff_index0);
					BufferStore::getBufferStore()->setUnusedIntArr(buff_index1);
					BufferStore::getBufferStore()->setUnusedIntArr(buff_index2);					
					return 0;
				};
			}
			++(ipiv[icol]);
			// Pivot found - interchanging rows
			if (irow != icol) {
				for (l=0;l<n;l++) {
					temp = a[irow][l];
					a[irow][l] = a[icol][l];
					a[icol][l] = temp;
				}
			}
			indxr[i]=irow;
			indxc[i]=icol;
			if (a[icol][icol] == 0.0) return 0;
			pivinv=1.0/a[icol][icol];
			a[icol][icol]=1.0;
			for (l=0;l<n;l++) a[icol][l] *= pivinv;
			for (ll=0;ll<n;ll++)
			if (ll != icol) {
				dum=a[ll][icol];
				a[ll][icol]=0.0;
				for (l=0;l<n;l++) a[ll][l] -= a[icol][l]*dum;
			}
	}
	for (l=n-1;l>=0;l--) {
		if (indxr[l] != indxc[l])
		for (k=0;k<n;k++){
			temp = a[k][indxr[l]];
			a[k][indxr[l]] = a[k][indxc[l]];
			a[k][indxc[l]] = temp;
		}
	}
	BufferStore::getBufferStore()->setUnusedIntArr(buff_index0);
	BufferStore::getBufferStore()->setUnusedIntArr(buff_index1);
	BufferStore::getBufferStore()->setUnusedIntArr(buff_index2);
	return 1;
}

int MatrixOperations::invert(Matrix* matrix){
	if(matrix->rows() != matrix->columns()){
		return 0;
	}
	return MatrixOperations::invert(matrix->getArray(),matrix->rows());
}

int MatrixOperations::mult(PhaseVector* v, Matrix* mtrx, PhaseVector* v_res){
	int n = mtrx->rows();
	int m = mtrx->columns();	
	if(v->size() != n || v_res->size() != m){
		ORBIT_MPI_Finalize("MatrixOperations:You try to multiply PhaseVector by Matrix with wrong size.");
		return 0;
	}
	double* vArr = v->getArray();
	double** arr = mtrx->getArray();
	double* vArr_res = v_res->getArray();	
	for(int i = 0; i < m; i++){
		vArr_res[i] = 0.;
		for(int j = 0; j < n; j++){
			vArr_res[i] += vArr[j]*arr[j][i];
		}
	}
	return 1;;
}

int MatrixOperations::mult(Matrix* mtrx, PhaseVector* v, PhaseVector* v_res){
	int n = mtrx->rows();
	int m = mtrx->columns();
	if(v->size() != m || v_res->size() != n){
		ORBIT_MPI_Finalize("MatrixOperations:You try to multiply Matrix by PhaseVector with wrong size.");
		return 0;
	}
	double* vArr = v->getArray();
	double** arr = mtrx->getArray();
	double* vArr_res = v_res->getArray();	
	for(int i = 0; i < n; i++){
		vArr_res[i] = 0.;
		for(int j = 0; j < m; j++){
			vArr_res[i] += arr[i][j]*vArr[j];
		}
	}
	return 1;
}

void MatrixOperations::track(Bunch* bunch,Matrix* mtrx){
	int n = mtrx->rows();
	int m = mtrx->columns();
	if( n != m || (m != 6 && m != 7)){
		ORBIT_MPI_Finalize("MatrixOperations:track(Bunch,Matrix) - Matrix has a wrong size.");
	}
	double tmp[5];
	double** bunch_arr = bunch->coordArr();
	double** arr = mtrx->getArray();
	for(int ip = 0, nParts = bunch->getSize(); ip < nParts; ip++){
		for(int i = 0; i < 6; i++){
			tmp[i] = 0.;
			for(int j = 0; j < 6; j++){
				tmp[i] += arr[i][j]*bunch_arr[ip][j];
			}
		}
		for(int i = 0; i < 6; i++){
			bunch_arr[ip][i] = tmp[i];
		}
		if(n == 7){
			for(int i = 0; i < 6; i++){
				bunch_arr[ip][i] += arr[i][6];
			}
		}
	}
}

