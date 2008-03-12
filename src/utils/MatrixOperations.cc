#include "MatrixOperations.hh"

#include <cmath>

#include "BufferStore.hh"

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
 		    } else if (ipiv[k] > 1) {return 0;};
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

