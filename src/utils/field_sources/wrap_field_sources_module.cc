#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include <iostream>

#include "wrap_bunch.hh"
#include "wrap_utils.hh"

#include "wrap_magnetic_field_source_grid3d.hh"
#include "wrap_quad_field_source.hh"
#include "wrap_dipole_field_source.hh"

namespace wrap_field_sources_module{

  void error(const char* msg){ ORBIT_MPI_Finalize(msg); }

#ifdef __cplusplus
extern "C" {
#endif

	
	// defenition of the memebers of the python wrapper module for functions
	// they will be vailable from python level
	static PyMethodDef FieldSourcesFunctionMethods[] = { 		
		/* {"name_of+function_on_python_level",  wrap_function_c    , METH_VARARGS, "help description"}, */
		{NULL, NULL, 0, NULL}        /* Sentinel */
	};

	//--------------------------------------------------
	//Initialization function of the module will be called 
	//from utils wrapper initialization.
	//--------------------------------------------------
	
  void initFieldSourcesModule(PyObject* module, const char* field_source_module_name){
    //create == field source == module
    PyObject* module_local = Py_InitModule(field_source_module_name,FieldSourcesFunctionMethods);
		PyModule_AddObject(module,field_source_module_name,module_local);
		Py_INCREF(module_local);
		wrap_field_source_grid3d::initMagnetFieldSourceGrid3D(module_local);
		wrap_quad_field_source::initQuadFieldSource(module_local);
		wrap_dipole_field_source::initDipoleFieldSource(module_local);
	}

#ifdef __cplusplus
}
#endif

//end of namespace wrap_field_sources_module
}
