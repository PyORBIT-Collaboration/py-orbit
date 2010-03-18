//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    Function.cc
//
// AUTHOR
//    Y. Sato, A. Shishlo
//
// CREATED
//    12/31/2003
//
// DESCRIPTION
//    Specification and inline functions for a class that keeps
//    table y(x) and does some operation with the tables.
//    It is using linear interpolation.
//
///////////////////////////////////////////////////////////////////////////
#include "orbit_mpi.hh"
#include "OU_Function.hh"

using namespace OrbitUtils;

Function::Function(): CppPyWrapper(NULL)
{
  clean();
  size = 0;
  sizeChunk = 50;
  maxSize = sizeChunk;
  x_arr = new double[maxSize];
  y_arr = new double[maxSize];

  //MPI stuffs
  rank_MPI = 0;
  size_MPI = 1;
  iMPIini  = 0;
  ORBIT_MPI_Initialized(&iMPIini);

  if(iMPIini > 0){
    ORBIT_MPI_Comm_size(MPI_COMM_WORLD, &size_MPI);
    ORBIT_MPI_Comm_rank(MPI_COMM_WORLD, &rank_MPI);
  }
}


Function::~Function()
{
  delete [] x_arr;
  delete [] y_arr;
}

void Function::resize()
{
  double* x_tmp = x_arr;
  double* y_tmp = y_arr;

  x_arr = new double[maxSize+sizeChunk];
  y_arr = new double[maxSize+sizeChunk];

  for(int i = 0; i < size; i++){
    x_arr[i] = x_tmp[i];
    y_arr[i] = y_tmp[i];
  }

  maxSize += sizeChunk;

  delete [] x_tmp;
  delete [] y_tmp;
}

void Function::finalize(const char* message)
{
  if(iMPIini > 0){
    ORBIT_MPI_Finalize(message);
  }
}

void Function::add(double x, double y)
{
  inf_const_step = 0;
  if((size+1) ==  maxSize){
    resize();
  }

  x_arr[size] = x;
  y_arr[size] = y;
  size++;

  if(xMin > x) xMin = x;
  if(yMin > y) yMin = y;
  if(xMax < x) xMax = x;
  if(yMax < y) yMax = y;
}

int Function::getSize()
{
  return size;
}

double Function::x(int ind)
{
  if(ind < size){
    return x_arr[ind];
  }
  finalize("ORBIT Utils Function class:The index in x(int ind) more than size");
  return 0.0;

}

double Function::y(int ind)
{
  if(ind < size){
    return y_arr[ind];
  }
  finalize("ORBIT Utils Function class:The index in y(int ind) more than size");
  return 0.0;
}

double* Function::xArr(){
	return x_arr;
}
		
double* Function::yArr(){
	return y_arr;
}


double Function::getMinX()
{
  return xMin;
}

double Function::getMinY()
{
  return yMin;
}

double Function::getMaxX()
{
  return xMax;
}

double Function::getMaxY()
{
  return yMax;
}

void Function::clean()
{
  inf_const_step = 0;
  x_step = 0.;
  size = 0;
  xMin = 1.0e+300;
  xMax = -1.0e+300;
  yMin = 1.0e+300;
  yMax = -1.0e+300;
}

double Function::getY(double x)
{

  if(size < 2){
    finalize("ORBIT Utils Function class:The Function method  getY(double x)  (size<2)");
  }

  if(x <= xMin) return y_arr[0];
  if(x >= xMax) return y_arr[size-1];

  int ind = 0;
  double yy = 0.;

  if(inf_const_step > 0){
    ind = (int)((x-xMin)/x_step);
    if(ind < 0 || ind > (size-2)){
      finalize("ORBIT Utils Function class: The Function method  y(double x)  ind < 0 or ind >= (size-2)");
    }
    yy = y_arr[ind] + (y_arr[ind+1] - y_arr[ind])*
      ((x - x_arr[ind])/x_step);
    return yy;
  }

  while(ind <  (size-1) && x > x_arr[ind]){
    ind++;
  }

  ind--;

  yy = y_arr[ind] + (y_arr[ind+1] - y_arr[ind])*
    ((x - x_arr[ind])/(x_arr[ind+1] - x_arr[ind]));
  return yy;
}

//this method should be used only for monotonic function
// f(x1) < f(x2) if x1 < x2
double Function::getX(double y)
{

  if(size < 2){
    finalize("ORBIT Utils Function class: The Function method  getX(double y)  (size<2)");
  }

  if(y <= yMin) return x_arr[0];
  if(y >= yMax) return x_arr[size-1];

  int ind = 0;
  double xx = 0.;

  int ind_start = 0;
  int ind_stop = size-1;
	int count = 0;
  while((ind_stop - ind_start) > 1){
		count++;
    ind = (ind_stop + ind_start)/2;
    if(y > y_arr[ind]){
      ind_start = ind;
    }
    else{
      if(y == y_arr[ind]){
				ind_stop = (size-1);
        ind_start = ind;
      }
      else{
				ind_stop =  ind;
      }
    }
		if(count > 10000){
			finalize("ORBIT Utils Function class: The Function method  getX(double y) has unlimited loop. Check data.");
		}
  }

	ind = ind_stop;

  ind--;

	if(y_arr[ind+1] != y_arr[ind]){
    xx = x_arr[ind] + (x_arr[ind+1] - x_arr[ind])*
    ((y - y_arr[ind])/(y_arr[ind+1] - y_arr[ind]));
	}
	else{
		xx = x_arr[ind];
	}

  return xx;
}

int Function::setConstStep(int info)
{
  if(info == 0){
    inf_const_step = 0;
    x_step = 0.;
		return 0;
  }
  else{
    inf_const_step = 1;

    if(size < 2){
      return 0;
    }

    //check that step is const
    x_step = x_arr[1] - x_arr[0];
    for(int i = 0; i < (size-1); i++){
      if(abs(x_step - (x_arr[i+1] - x_arr[i]))/x_step > 1.0e-11){
	      return 0;
      }
    }
  }
	return 1;
}

//return 1 if step on x is constant and 0 - otherwise
int Function::isStepConst()
{
  return inf_const_step;
}

// Sets the inverse Function. The x-coordinates
// of f_inv could be defined already.
// It will return 1 if it was a success and 0 otherwise
int Function::setInverse(Function* f_inv)
{
  //check that inverse function can be made
	
	if(size < 2){
		f_inv->clean();
		return 0;
	}
	
	//create x-arr if it is not ready
	if(f_inv->getSize() < 2){
		f_inv->clean();
		for(int i = 0; i < size; i++){
			double xx = getMinY() + i*(getMaxY() - getMinY())/(size - 1);
			f_inv->add(xx,0.);
		}
		f_inv->setConstStep(1);
	}
	
	if(y_arr[1] > y_arr[0]){
		for(int i = 0; i < (size-1); i++){
			if(y_arr[i] >= y_arr[i+1]){
				f_inv->clean();
				return 0;
			}
		}
	}
	else{
		for(int i = 0; i < (size-1); i++){
			if(y_arr[i] <= y_arr[i+1]){
				f_inv->clean();
				return 0;
			}
		}
	}
	
	if(f_inv->getMaxX() > getMaxY()){
		finalize("ORBIT Utils Function class: The Function method  setInverse(Function* f_inv) f_inv->getMaxX() > getMaxY()");
		return 0;
	}
	
	if(f_inv->getMinX() < getMinY()){
		finalize("ORBIT Utils Function class: The Function method  setInverse(Function* f_inv) f_inv->getMinX() < getMinY()");
		return 0;
	}
	
	Function* f_tmp = new Function();
	
	int nP = f_inv->getSize();
	
	double xx = 0.;
	double yy = 0.;
	double coeff = 0.;
	int ind = 0;
	
	for(int i = 0; i < nP; i++){
		xx = f_inv->x(i);
		if(y_arr[1] > y_arr[0]){
			ind = 0;
			while(ind <  (size-1) && xx >= y_arr[ind]){
				ind++;
			}
			ind--;
		}
		else{
    	ind = size-1;
			while(ind > 0 && xx <= y_arr[ind]){
				ind--;
			}
		}
		coeff = (xx - y_arr[ind])/(y_arr[ind+1]-y_arr[ind]);
		yy = x_arr[ind] + (x_arr[ind+1] - x_arr[ind])*coeff;
		f_tmp->add(xx,yy);
	}
	
	int inf_tmp = f_inv->isStepConst();
	
	nP = f_tmp->getSize();
	f_inv->clean();
	for(int i = 0; i < nP; i++){
		f_inv->add(f_tmp->x(i),f_tmp->y(i));
	}
	f_inv->setConstStep(inf_tmp);
	
	//remove tmp Function
	delete f_tmp;
	
	return 1;
}


void Function::print(ostream& Out)
{
  if(rank_MPI == 0){
    Out<<"% size = "<< getSize() <<std::endl
       <<"% minX = "<< getMinX() <<std::endl
       <<"% maxX = "<< getMaxX() <<std::endl
       <<"% minY = "<< getMinY() <<std::endl
       <<"% maxY = "<< getMaxY() <<std::endl
       <<"% x-step const = "<< isStepConst() <<std::endl;

      Out<<"% #i      x     y"<<std::endl;

        for(int i = 0; i < size; i++){
          Out<<" "<< i
	     <<"   \t"<< x_arr[i]
	     <<"   \t"<< y_arr[i]
             <<std::endl;
	  if(i % 1000 == 0) Out.flush();
	}
  }
}

void Function::print(const char* fileName)
{
  ofstream F_dump;
  if(rank_MPI == 0)F_dump.open (fileName, ios::out);
  print(F_dump);
  if(rank_MPI == 0){F_dump.close();}
  return;
}

//auxiliary method to create normalize cumulative function
//for probability distribution with y_min = 0 and y_max = 1.0
//It returns 1 if it was a success and 0 otherwise 
int Function::normalize()
{
  if(size < 2){
    return 0;
  }

  for(int i = 0; i < (size-1); i++){
    if(y_arr[i] > y_arr[i+1]){
      return 0;
    }
  }

  for(int i = 0; i < size; i++){
    y_arr[i] /=y_arr[size-1];
  }

  y_arr[size-1] = 1.0;
  y_arr[0] = 0.;

  yMin = 1.0e+300;
  yMax = -1.0e+300;
  for(int i = 0; i < size; i++){
    if(y_arr[i] < yMin) yMin = y_arr[i];
    if(y_arr[i] > yMax) yMax = y_arr[i];
  }
	return 1;
}

