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
//    to integer and double arrays that will be used in parallel exchange.
//
///////////////////////////////////////////////////////////////////////////// 

#ifndef BUFFER_STORE_H
#define BUFFER_STORE_H

/////////////////////////////////////////////////////////////////////////////
//
// include files
//
/////////////////////////////////////////////////////////////////////////////


class BufferStore
{
 public:

  static BufferStore* getBufferStore();

  ~BufferStore();

  //returns the array of doubles that can be used as 
  //temporary one for MPI operations
  double* getDoubleArr(int index, int size);

  int getDoubleArrSize(int index);
  int getDoubleArrSize();

  //returns the array of integers that can be used as 
  //temporary one for MPI operations
  int* getIntArr(int index, int size);

  int getIntArrSize(int index);
  int getIntArrSize();

  //returns the array of chars that can be used as 
  //temporary one for MPI operations
  char* getCharArr(int index, int size);

  int getCharArrSize(int index);
  int getCharArrSize();

 private:

   //constructor
  BufferStore();

  static  BufferStore* bStore;

  //data arrays
  int** int_data;
  int* int_size;
  int nInt;

  double** dbl_data;
  int* dbl_size;
  int nDbl;

  char** char_data;
  int* char_size;
  int nChar;

};


///////////////////////////////////////////////////////////////////////////
//
// END OF FILE
//
///////////////////////////////////////////////////////////////////////////

#endif   // BUFFER_STORE_H


