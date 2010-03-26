//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    SplineCH.cc
//
// AUTHOR
//    A. Shishlo
//
// CREATED
//    03/25/2010
//
// DESCRIPTION
//    A cubic Hermite spline. 
//    http://en.wikipedia.org/wiki/Cubic_Hermite_spline
//    y(t) = h00(t)*y0 + h10(t)*m0 + h01(t)*y1 + h11(t)*m1
//    t = (x-x0)/(x1-x0)
//    h00(t)=2*t^3 - 3*t^2 +1
//    h10(t)=t^3 - 2*t^2 +t
//    h01(t)=-2*t^3 + 3*t^2
//    h11(t)=t^3 - t^2
//    m_k = 0.5*((y[k+1]-y[k])/(x[k+1]-x[k])+(y[k]-y[k-1])/(x[k]-x[k-1]))
//    m0 = (y[1] - y[0])/(x[1]-x[0])
//    m[n-2] = (y[n-1] - y[n-2])/(x[n-1]-x[n-2])
//
///////////////////////////////////////////////////////////////////////////
#include "orbit_mpi.hh"
#include "OU_SplineCH.hh"

#include <iomanip>

namespace OrbitUtils{
	
	
	SplineCH::SplineCH(): CppPyWrapper(NULL)
	{  
		x_arr = NULL;
		y_arr = NULL;
		m_arr = NULL;
		size = 0;
		
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
	
	
	SplineCH::~SplineCH()
	{
		if(x_arr != NULL) delete [] x_arr;
		if(y_arr != NULL) delete [] y_arr;
		if(m_arr != NULL) delete [] m_arr;
	}
	
	int SplineCH::compile(OrbitUtils::Function* f)
	{
		if(x_arr != NULL) delete [] x_arr;
		if(y_arr != NULL) delete [] y_arr;
		if(m_arr != NULL) delete [] m_arr;
				
		size = f->getSize();
		if(size < 3){
			size = 0;
			return 0;
		}
		
		//check the x_arr
		for(int i = 0; i < (size-1); i++){
			if(f->x(i) >= f->x(i+1)){
				size = 0;
				return 0;
			}
		}
		
		x_arr = new double[size];
		y_arr = new double[size];
		m_arr = new double[size];
		
		for(int i = 0; i < size; i++){
			x_arr[i] = f->x(i);
			y_arr[i] = f->y(i);
		}
		
		m_arr[0] = (y_arr[1] - y_arr[0])/(x_arr[1] - x_arr[0]);
		m_arr[size-1] = (y_arr[size-1] - y_arr[size-2])/(x_arr[size-1] - x_arr[size-2]);
		
		for(int i = 1, n = size-1; i < n; i++){
			m_arr[i] = 0.5*((y_arr[i+1] - y_arr[i])/(x_arr[i+1] - x_arr[i]) + 
			           (y_arr[i] - y_arr[i-1])/(x_arr[i] - x_arr[i-1])); 
		}				
	}
	
	void SplineCH::finalize(const char* message)
	{
		if(iMPIini > 0){
			ORBIT_MPI_Finalize(message);
		}
	}
	
	int SplineCH::getSize()
	{
		return size;
	}
	
	double SplineCH::x(int ind)
	{
		if(ind < size){
			return x_arr[ind];
		}
		finalize("ORBIT Utils SplineCH class:The index in x(int ind) more than size");
		return 0.0;
		
	}
	
	double SplineCH::y(int ind)
	{
		if(ind < size){
			return y_arr[ind];
		}
		finalize("ORBIT Utils SplineCH class:The index in y(int ind) more than size");
		return 0.0;
	}
	
	double SplineCH::getY(double x)
	{
		
		if(size < 3){
			finalize("ORBIT Utils SplineCH class: number of points less than 3, no spline.");
			return 0.;
		}		
		
		if(x <= x_arr[0]) return y_arr[0];
		if(x >= x_arr[size-1]) return y_arr[size-1];
		
		int ind = 0;

		int ind_start = 0;
		int ind_stop = size-1;
		int count = 0;
		while((ind_stop - ind_start) > 1){
			count++;
			ind = (ind_stop + ind_start)/2;
			if(x > x_arr[ind]){
				ind_start = ind;
			}
			else{
				ind_stop =  ind;
			}
			if(count > 200){
				finalize("ORBIT Utils SplineCH class: The SplineCH method  getX(double y) has unlimited loop. Check data.");
			}
		}	
		
		ind = ind_start;
		double dx = x_arr[ind+1] - x_arr[ind];	
		double t = (x - x_arr[ind])/dx0;
		double t2 = t*t;
		double t3 = t2*t;
		double yy = y_arr[ind]*(2*t3-3*t2+1.0) + m_arr[ind]*(t3-2*t2+t)*dx + y_arr[ind+1]*(-2*t3+3*t2) + m_arr[ind+1]*(t3-t2)*dx; 

		return yy;
	}
	
	double SplineCH::getYP(double x)
	{
		
		if(size < 3){
			finalize("ORBIT Utils SplineCH class: number of points less than 3, no spline.");
			return 0.;
		}		
		
		if(x <= x_arr[0]) return m_arr[0];
		if(x >= x_arr[size-1]) return m_arr[size-1];
		
		int ind = 0;

		int ind_start = 0;
		int ind_stop = size-1;
		int count = 0;
		while((ind_stop - ind_start) > 1){
			count++;
			ind = (ind_stop + ind_start)/2;
			if(x > x_arr[ind]){
				ind_start = ind;
			}
			else{
				ind_stop =  ind;
			}
			if(count > 200){
				finalize("ORBIT Utils SplineCH class: The SplineCH method  getX(double y) has unlimited loop. Check data.");
			}
		}	
		
		ind = ind_start;
		double dx = x_arr[ind+1] - x_arr[ind];
		double t = (x - x_arr[ind])/dx;
		double t2 = t*t;
		double yyp = 6*y_arr[ind]*(t2-t)/dx + m_arr[ind]*(3*t2-4*t+1) + 6*y_arr[ind+1]*(t-t2)/dx + m_arr[ind+1]*(3*t2-2*t); 

		return yyp;
	}	
	
	void SplineCH::print(ostream& Out)
	{
		if(rank_MPI == 0){
		  Out<<std::setprecision(15)<< std::setiosflags(std::ios::scientific);			
			Out<<"% size = "<< getSize() <<std::endl;
			Out<<"% #i      x     y    m(Cubic Hermite)"<<std::endl;
			
			for(int i = 0; i < size; i++){
				Out<<" "<< i
				<<"   \t"<< x_arr[i]
				<<"   \t"<< y_arr[i]
				<<"   \t"<< m_arr[i]
				<<std::endl;
				if(i % 1000 == 0) Out.flush();
			}
		}
	}
	
	void SplineCH::print(const char* fileName)
	{
		ofstream F_dump;
		if(rank_MPI == 0)F_dump.open (fileName, ios::out);
		print(F_dump);
		if(rank_MPI == 0){F_dump.close();}
		return;
	}
}

