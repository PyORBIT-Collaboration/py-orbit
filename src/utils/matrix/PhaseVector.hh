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
	
	/** A double values vector. */
	
	class PhaseVector
	{
		public:
		
			/** Constructor. */
			PhaseVector(int n);
			
			/** A copy-constructor. */
			PhaseVector(PhaseVector* vIn);
			
			/** Destructor. */
			~PhaseVector();
		
			/** Returns the pointer to the inner array. */
			double* getArray();
			
			/** Returns the size of the vector. */
			int size();
			
			/** Sets the vector to zero. */
			void zero();		
			
			/** Returns one component with index i. */
			double& value(int i);
			
			/** Copies the existing vector to another. */
			int copyTo(PhaseVector* vIn);
			
			/** Adds a value to each component. */
			int add(double val);
			
			/** Adds a vector to the existing one. */
			int add(PhaseVector* mtrx);
			
			/** Returns the norm. */
			double norm();
			
			/** Multiplies the vector's components by a number. */
			int mult(double val);
			
			/** Returns scalar multiplication of two vectors. */
			double dot(PhaseVector* vctr);
		
		private:
		
			int n;
			
			double* v;
		
	};
};

#endif
