//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    PyExternalEffects.hh
//
// CREATED
//    05/14/2008
//
// DESCRIPTION
//    The base class for Python implementation of a external effects 
//    during the transport of particles through the external field. 
//    It should be sub-classed on Python level and implement
//    setupEffects(Bunch* bunch)
//    finalizeEffects(Bunch* bunch)
//    applyEffects(Bunch* bunch, int index, 
//	                            double* y_in_vct, double* y_out_vct, 
//														  double t, double t_step, 
//														  OrbitUtils::BaseFieldSource* fieldSource)
//    methods.
//    The results of these methods will be available from the c++ level.
//    This is an example of embedding Python in C++ Orbit level.
//
//    ATTENTION: Using this class in real calculations is not wise! 
//               It is slow, because it delegates the field calculations
//               to the Python level. It is for prototyping and 
//               debugging only.
//
///////////////////////////////////////////////////////////////////////////
#ifndef PY_EXTERNAL_EFFECTS_H
#define PY_EXTERNAL_EFFECTS_H

#include "Python.h"

#include "ExternalEffects.hh"

namespace TrackerRK4{
	
	class  PyExternalEffects: public ExternalEffects
	{
		public:
		
			/** Constructor. */
			PyExternalEffects(PyObject* py_wrapperIn);
			
			/** Destructor. */
			~PyExternalEffects();
		
		/** It initializes effects. */
		void setupEffects(Bunch* bunch);
		
		/** It finalizes effects. */
		void finalizeEffects(Bunch* bunch);

		/** It applies the external effects to a particle with certain index. */
		void applyEffects(Bunch* bunch, int index, 
	                            double* y_in_vct, double* y_out_vct, 
														  double t, double t_step, 
														  OrbitUtils::BaseFieldSource* fieldSource,
															RungeKuttaTracker* tracker);	

	};
};

#endif
