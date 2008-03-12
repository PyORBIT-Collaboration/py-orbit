/////////////////////////////////////////////////////////////////////////////
//
// FILE NAME
//   BufferStore.hh
//
// AUTHORS
//    A. Shishlo
//
// CREATED
//    07/01/2005
//
// DESCRIPTION
//    BufferStore class - singleton to keep references 
//    to integer and double arrays that will be used as temporary ones.  
//
///////////////////////////////////////////////////////////////////////////// 

#ifndef BUFFER_STORE_H
#define BUFFER_STORE_H

/**
	 This is a class based on a singleton pattern 
	 that keeps sets of one dimensional arrays.
 */

class BufferStore
{
 public:

	 /** Returns the buffers store. */ 
  static BufferStore* getBufferStore();

	/** We do not need a destructor actually. */ 
  ~BufferStore();


	/** 
	   Returns a reference to a free double array and an index of this array.
		 After usage the array should be freed by setUnusedDoubleArr method.
	 */
  double* getFreeDoubleArr(int &index, int size);	

	/** Returns a size of the array with particular index. */
  int getDoubleArrSize(int index);
	
	/** Returns a number of the arrays. */
  int getDoubleArrSize();
	
	/** Sets that the array is free. */
	void setUnusedDoubleArr(int index);

	/** 
	   Returns a reference to a free int array and an index of this array.
		 After usage the array should be freed by setUnusedDoubleArr method.
	 */
	int* getFreeIntArr(int &index, int size);
	
	/** Returns a size of the array with particular index. */
  int getIntArrSize(int index);
	
	/** Returns a number of the arrays. */
	int getIntArrSize();
	
	/** Sets that the array is free. */
	void setUnusedIntArr(int index);	
	
	/** 
	   Returns a reference to a free char array and an index of this array.
		 After usage the array should be freed by setUnusedDoubleArr method.
	 */	
	char* getFreeCharArr(int &index, int size);
	
	/** Returns a size of the array with particular index. */
  int getCharArrSize(int index);
	
	/** Returns a number of the arrays. */
  int getCharArrSize();	
	
	/** Sets that the array is free. */
	void setUnusedCharArr(int index);	

 private:

   /** Constructor. */
  BufferStore();
	
  /**
	   Returns a reference to the array of doubles with particular 
		 index and mark this array as used. 
  */
  double* getDoubleArr(int index, int size);
	
	/** Returns 1 if the array is used somewhere, and 0 otherwise. */
	int isInUseDoubleArr(int index);
	
	/** Sets that the array is used somewhere. */
	void setInUseDoubleArr(int index);  
	
	
  int* getIntArr(int index, int size);		
	int isInUseIntArr(int index);
	void setInUseIntArr(int index);

  char* getCharArr(int index, int size);
	int isInUseCharArr(int index);
	void setInUseCharArr(int index);	
	
	static  BufferStore* bStore;

  //data arrays
  int** int_data;
  int* int_size;
	int* int_in_use;
  int nInt;

  double** dbl_data;
  int* dbl_size;
  int* dbl_in_use;
  int nDbl;

  char** char_data;
  int* char_size;
	int* char_in_use;
  int nChar;

};


///////////////////////////////////////////////////////////////////////////
//
// END OF FILE
//
///////////////////////////////////////////////////////////////////////////

#endif   // BUFFER_STORE_H


