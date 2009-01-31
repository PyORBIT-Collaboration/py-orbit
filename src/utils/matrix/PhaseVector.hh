//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//   PhaseVector.hh
//
// AUTHOR
//    A. Shishlo
//
// CREATED
//    03/12/2008
//
// DESCRIPTION
//     A class for a plain double values vector.
//
///////////////////////////////////////////////////////////////////////////

#ifndef PHASE_VECTOR_H
#define PHASE_VECTOR_H

namespace OrbitUtils{
	class PhaseVector
	{
		public:
		
			PhaseVector(int n);
			PhaseVector(PhaseVector* vIn);
			~PhaseVector();
		
			double* getArray();
			int size();
			void zero();			
			double& value(int i);
			
			int copyTo(PhaseVector* vIn);
			int add(double val);
			int add(PhaseVector* mtrx);
			double norm();
			int mult(double val);
			double dot(PhaseVector* vctr);
		
		private:
		
			int n;
			
			double* v;
		
	};
};

#endif
