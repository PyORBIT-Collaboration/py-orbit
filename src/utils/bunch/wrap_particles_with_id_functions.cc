#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include <iostream>

#include "wrap_bunch.hh"

#include "ParticlesWithIdFunctions.hh"

using namespace OrbitUtils;

namespace wrap_utils_bunch_functions{

  void error(const char* msg){ ORBIT_MPI_Finalize(msg); }

#ifdef __cplusplus
extern "C" {
#endif
	//---------------------------------------------------------
	//Python ParticlesWithIdFunctions functions definition
	//---------------------------------------------------------
	
	static PyObject* part_wit_id_sort(PyObject* self, PyObject* args)
	{

		PyObject *pyIn;
		if(!PyArg_ParseTuple(args,"O:bunchSortId",&pyIn)){
			error("bunchSortId(Bunch) - Bunch is needed.");
		}			
		PyObject* pyBunchType = wrap_orbit_bunch::getBunchType("Bunch");
		if((!PyObject_IsInstance(pyIn,pyBunchType))){
			error("bunchSortId(Bunch) - input parameter is not Bunch");
		}		
		Bunch* bunch = (Bunch*) ((pyORBIT_Object*) pyIn)->cpp_obj;
		bunch_sort_id(bunch);
    Py_INCREF(Py_None);
    return Py_None;	
	}	

	// defenition of the memebers of the python ParticlesWithIdFunctions wrapper module for functions
	// they will be vailable from python level
	static PyMethodDef NumrecipesModuleMethods[] = { 		
		{"bunchSortId", part_wit_id_sort , METH_VARARGS, "bunchSortId(bunch) - Sorting bunch according to the Id particles attributes."},
		{NULL, NULL, 0, NULL}        /* Sentinel */
	};

	//--------------------------------------------------
	//Initialization function of the pyParticlesWithIdFunctions module
	//It will be called from Bunch wrapper initialization
	//--------------------------------------------------
	
  void initParticlesWithIdFunctions(PyObject* module, const char* part_with_id_module_name){
    //create numrecipes module
    PyObject* module_partId = Py_InitModule(part_with_id_module_name,NumrecipesModuleMethods);
		PyModule_AddObject(module,part_with_id_module_name,module_partId);	
	}

#ifdef __cplusplus
}
#endif

//end of namespace wrap_utils_bunch_functions
}
