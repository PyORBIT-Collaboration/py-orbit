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
//    to integer, double, and char arrays that will be used in parallel exchange.
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

BufferStore::BufferStore()
{

  int_data = (int** ) malloc (sizeof(int*));
  int_data[0] = (int* ) malloc (sizeof(int));
  int_size = (int* ) malloc (sizeof(int));
  int_size[0] = 1;
	int_in_use = (int* ) malloc (sizeof(int));
	int_in_use[0] = 0;
  nInt = 1;

  dbl_data = (double** ) malloc (sizeof(double*));
  dbl_data[0] = (double* ) malloc (sizeof(double));
  dbl_size = (int* ) malloc (sizeof(int));
  dbl_size[0] = 1;
	dbl_in_use = (int* ) malloc (sizeof(int));
	dbl_in_use[0] = 0;	
  nDbl = 1;


  char_data = (char** ) malloc (sizeof(char*));
  char_data[0] = (char* ) malloc (sizeof(char));
  char_size = (int* ) malloc (sizeof(int));
  char_size[0] = 1;
	char_in_use = (int* ) malloc (sizeof(int));
	char_in_use[0] = 0;	
  nChar = 1;
}

BufferStore::~BufferStore()
{
  for(int i = 0; i < nInt; i++){
    free(int_data[i]);
  }
  free(int_data);
  free(int_size);
	free(int_in_use);
  
  for(int i = 0; i < nDbl; i++){
    free(dbl_data[i]);
  }
  free(dbl_data);
  free(dbl_size);  
	free(dbl_in_use);
	
  for(int i = 0; i < nChar; i++){
    free(char_data[i]);
  }
  free(char_data);
  free(char_size); 
	free(char_in_use);
}

BufferStore* BufferStore::getBufferStore()
{
  if(!bStore){
    bStore = new  BufferStore();
  }
  return bStore;
}

double* BufferStore::getDoubleArr(int index, int size){

  if(size <= 0) size = 1;

  if(index > (nDbl - 1)){
    int nDbl_new = index + 1;
    dbl_data = (double **) realloc(dbl_data, nDbl_new*sizeof(double*));
    dbl_size = (int*) realloc(dbl_size,nDbl_new*sizeof(int));
		dbl_in_use = (int*) realloc(dbl_in_use,nDbl_new*sizeof(int));
    for(int i = nDbl; i <= index; i++){
      dbl_data[i] =  (double *) malloc (size*sizeof(double));
      dbl_size[i] = size;
			dbl_in_use[i] = 0;
    } 
    nDbl = nDbl_new;   
  }

  if(size > dbl_size[index]){
    dbl_data[index] = (double *) realloc(dbl_data[index],size*sizeof(double));
    dbl_size[index] = size;
  }
	
	dbl_in_use[index] = 1;
  return dbl_data[index];
}

double* BufferStore::getFreeDoubleArr(int &index, int size){
	for(int i = 0; i < getDoubleArrSize(); i++){
		int used = isInUseDoubleArr(i);
		if(used == 0){
			index = i;
			return getDoubleArr(i,size);
		}
	}
	index = getDoubleArrSize();
	return getDoubleArr(index,size);
}

int BufferStore::getDoubleArrSize(int index){
  return dbl_size[index];
}

int BufferStore::getDoubleArrSize(){
  return nDbl;
}

int BufferStore::isInUseDoubleArr(int index){
	return dbl_in_use[index];
}
	
void BufferStore::setInUseDoubleArr(int index){
	dbl_in_use[index] = 1;
}

void BufferStore::setUnusedDoubleArr(int index){
	dbl_in_use[index] = 0;
}


int* BufferStore::getIntArr(int index, int size){

  if(size <= 0) size = 1;

  if(index > (nInt - 1)){
    int nInt_new = index + 1;
    int_data = (int **) realloc(int_data, nInt_new*sizeof(int*));
    int_size = (int*) realloc(int_size,nInt_new*sizeof(int));
		int_in_use = (int*) realloc(int_in_use,nInt_new*sizeof(int));
    for(int i = nInt; i <= index; i++){
      int_data[i] =  (int *) malloc (size*sizeof(int));
      int_size[i] = size;
			int_in_use[i] = 0;
    } 
    nInt = nInt_new;   
  }

  if(size > int_size[index]){
    int_data[index] = (int *) realloc(int_data[index],size*sizeof(int));
    int_size[index] = size;
  }

	int_in_use[index] = 1;
  return int_data[index];
}

int* BufferStore::getFreeIntArr(int &index, int size){
	for(int i = 0; i < getIntArrSize(); i++){
		int used = isInUseIntArr(i);
		if(used == 0){
			index = i;
			return getIntArr(i,size);
		}
	}
	index = getIntArrSize();
	return getIntArr(index,size);
}

int BufferStore::getIntArrSize(int index){
  return int_size[index];
}

int BufferStore::getIntArrSize(){
  return nInt;
}

int BufferStore::isInUseIntArr(int index){
	return int_in_use[index];
}
	
void BufferStore::setInUseIntArr(int index){
	int_in_use[index] = 1;
}

void BufferStore::setUnusedIntArr(int index){
	int_in_use[index] = 0;
}


char* BufferStore::getCharArr(int index, int size){

  if(size <= 0) size = 1;

  if(index > (nInt - 1)){
    int nChar_new = index + 1;
    char_data = (char **) realloc(char_data, nChar_new*sizeof(char*));
    char_size = (int*) realloc(char_size,nChar_new*sizeof(int));
		char_in_use = (int*) realloc(char_in_use,nChar_new*sizeof(int));
    for(int i = nChar; i <= index; i++){
      char_data[i] =  (char *) malloc (size*sizeof(char));
      char_size[i] = size;
			char_in_use[i] = 0;
    } 
    nChar = nChar_new;   
  }

  if(size > int_size[index]){
    char_data[index] = (char *) realloc(char_data[index],size*sizeof(char));
    char_size[index] = size;
  }
	
	char_in_use[index] = 1;
  return char_data[index];
}

char* BufferStore::getFreeCharArr(int &index, int size){
	for(int i = 0; i < getCharArrSize(); i++){
		int used = isInUseCharArr(i);
		if(used == 0){
			index = i;
			return getCharArr(i,size);
		}
	}
	index = getCharArrSize();
	return getCharArr(index,size);
}

int BufferStore::getCharArrSize(int index){
  return char_size[index];
}

int BufferStore::getCharArrSize(){
  return nChar;
}

int BufferStore::isInUseCharArr(int index){
	return char_in_use[index];
}
	
void BufferStore::setInUseCharArr(int index){
	char_in_use[index] = 1;
}

void BufferStore::setUnusedCharArr(int index){
	char_in_use[index] = 0;
}


