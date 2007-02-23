/////////////////////////////////////////////////////////////////////////////
//
// FILE NAME
//   BufferStore.cc
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
/////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////
//
// include files
//
/////////////////////////////////////////////////////////////////////////////
#include <cstdlib>
#include <iostream>

#include "BufferStore.hh"


BufferStore* BufferStore::bStore = NULL;

/**
//This is a test program for the BufferStore class
int main (int argc, char **argv)
{
  BufferStore* bs = BufferStore::getBufferStore();
  double* da_0 = bs->getDoubleArr(0,100);
  int* ia_0 = bs->getIntArr(0,200);

  int count = 0;
  int size = 100;

  while(1){
    if(size > 1000000) size = 1000000;
    for(int j = 0; j < 5; j++){
       da_0 = bs->getDoubleArr(j,size+j*1000);
       ia_0 = bs->getIntArr(j,size+j*1000);
       for(int i = 0; i < bs->getDoubleArrSize(j); i++){
	 da_0[i] = i*1.0;
       }
       for(int i = 0; i < bs->getIntArrSize(j); i++){
	 ia_0[i] = i;
       }
    }
    if(count % 10000 == 0) {
      std::cout<<"count="<<count<<" size="<< bs->getDoubleArrSize(0)
	       <<" size="<< bs->getIntArrSize(4)
	       <<" s="<<size
	       <<std::endl;
    }
    count++;
    size += 1000;
  }
  std::cout<<"Stop"<<std::endl;
  return 0;
}
*/


BufferStore::BufferStore()
{

  int_data = (int** ) malloc (sizeof(int*));
  int_data[0] = (int* ) malloc (sizeof(int));
  int_size = (int* ) malloc (sizeof(int));
  int_size[0] = 1;
  nInt = 1;

  dbl_data = (double** ) malloc (sizeof(double*));
  dbl_data[0] = (double* ) malloc (sizeof(double));
  dbl_size = (int* ) malloc (sizeof(int));
  dbl_size[0] = 1;
  nDbl = 1;


  char_data = (char** ) malloc (sizeof(char*));
  char_data[0] = (char* ) malloc (sizeof(char));
  char_size = (int* ) malloc (sizeof(int));
  char_size[0] = 1;
  nChar = 1;

}

BufferStore::~BufferStore()
{
  for(int i = 0; i < nInt; i++){
    free(int_data[i]);
  }
  free(int_data);
  free(int_size);
  
  for(int i = 0; i < nDbl; i++){
    free(dbl_data[i]);
  }
  free(dbl_data);
  free(dbl_size);  

  for(int i = 0; i < nChar; i++){
    free(char_data[i]);
  }
  free(char_data);
  free(char_size); 
}

BufferStore* BufferStore::getBufferStore()
{

  if(!bStore){
    bStore = new  BufferStore();
    return  bStore;
  }
  return bStore;
}

double* BufferStore::getDoubleArr(int index, int size){

  if(size <= 0) size = 1;

  if(index > (nDbl - 1)){
    int nDbl_new = index + 1;
    dbl_data = (double **) realloc(dbl_data, nDbl_new*sizeof(double*));
    dbl_size = (int*) realloc(dbl_size,nDbl_new*sizeof(int));
    for(int i = nDbl; i <= index; i++){
      dbl_data[i] =  (double *) malloc (size*sizeof(double));
      dbl_size[i] = size;
    } 
    nDbl = nDbl_new;   
  }

  if(size > dbl_size[index]){
    dbl_data[index] = (double *) realloc(dbl_data[index],size*sizeof(double));
    dbl_size[index] = size;
  }

  return dbl_data[index];
}

int BufferStore::getDoubleArrSize(int index){
  return dbl_size[index];
}


int BufferStore::getDoubleArrSize(){
  return nDbl;
}

int* BufferStore::getIntArr(int index, int size){

  if(size <= 0) size = 1;

  if(index > (nInt - 1)){
    int nInt_new = index + 1;
    int_data = (int **) realloc(int_data, nInt_new*sizeof(int*));
    int_size = (int*) realloc(int_size,nInt_new*sizeof(int));
    for(int i = nInt; i <= index; i++){
      int_data[i] =  (int *) malloc (size*sizeof(int));
      int_size[i] = size;
    } 
    nInt = nInt_new;   
  }

  if(size > int_size[index]){
    int_data[index] = (int *) realloc(int_data[index],size*sizeof(int));
    int_size[index] = size;
  }

  return int_data[index];
}

int BufferStore::getIntArrSize(int index){
  return int_size[index];
}

int BufferStore::getIntArrSize(){
  return nInt;
}






char* BufferStore::getCharArr(int index, int size){

  if(size <= 0) size = 1;

  if(index > (nInt - 1)){
    int nChar_new = index + 1;
    char_data = (char **) realloc(char_data, nChar_new*sizeof(char*));
    char_size = (int*) realloc(char_size,nChar_new*sizeof(int));
    for(int i = nChar; i <= index; i++){
      char_data[i] =  (char *) malloc (size*sizeof(char));
      char_size[i] = size;
    } 
    nChar = nChar_new;   
  }

  if(size > int_size[index]){
    char_data[index] = (char *) realloc(char_data[index],size*sizeof(char));
    char_size[index] = size;
  }

  return char_data[index];
}

int BufferStore::getCharArrSize(int index){
  return char_size[index];
}

int BufferStore::getCharArrSize(){
  return nChar;
}
