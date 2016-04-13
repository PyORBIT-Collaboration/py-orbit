#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include <iostream>

#include "wrap_bunch.hh"
#include "wrap_utils.hh"

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
	
	static PyObject* wrap_transport_mtrx(PyObject* self, PyObject* args)
	{

		PyObject *pyIn;
		PyObject *pyOut;
		PyObject *pyM;
		if(!PyArg_ParseTuple(args,"OOO:trasportMtrx",&pyIn,&pyOut,&pyM)){
			error("trasportMtrx(Bunch in,Bunch out,Matrix) - Bunches and Matrix are needed.");
		}			
		PyObject* pyBunchType = wrap_orbit_bunch::getBunchType("Bunch");
		if((!PyObject_IsInstance(pyIn,pyBunchType)) || (!PyObject_IsInstance(pyOut,pyBunchType)) ){
			error("trasportMtrx(Bunch in,Bunch out,Matrix) - 1st or 2nd input parameter is not Bunch");
		}	
		PyObject* pyORBIT_Matrix_Type = wrap_orbit_utils::getOrbitUtilsType("Matrix");
		if(!PyObject_IsInstance(pyM,pyORBIT_Matrix_Type)){
			error("trasportMtrx(Bunch in,Bunch out,Matrix) - function needs a matrix.");
		}		
		Bunch* bunch_in = (Bunch*) ((pyORBIT_Object*) pyIn)->cpp_obj;
		Bunch* bunch_out = (Bunch*) ((pyORBIT_Object*) pyOut)->cpp_obj;
		Matrix* mtrx = (Matrix*) ((pyORBIT_Object*) pyM)->cpp_obj;
		
		int n_stat = transport_mtrx(bunch_in,bunch_out,mtrx);
    return Py_BuildValue("i", n_stat);	
	}		

	// defenition of the memebers of the python ParticlesWithIdFunctions wrapper module for functions
	// they will be vailable from python level
	static PyMethodDef NumrecipesModuleMethods[] = { 		
		{"bunchSortId",  part_wit_id_sort    , METH_VARARGS, "bunchSortId(bunch) - Sorting bunch according to the Id particles attributes."},
		{"trasportMtrx", wrap_transport_mtrx , METH_VARARGS, "trasportMtrx(bunch in, bunch out, matrix) - calculates transport matrix."},
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
		Py_INCREF(module_partId);
	}

#ifdef __cplusplus
}
#endif

//end of namespace wrap_utils_bunch_functions
}
