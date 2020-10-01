#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include <iostream>

#include "wrap_bunch.hh"
#include "wrap_utils.hh"

#include "ParticlesWithIdFunctions.hh"
#include "TwissFilteringFunctions.hh"
#include "InitialCoordsAttrFunctions.hh"

using namespace OrbitUtils;

namespace wrap_utils_bunch_functions{

  void error(const char* msg){ ORBIT_MPI_Finalize(msg); }

#ifdef __cplusplus
extern "C" {
#endif
	//---------------------------------------------------------
	//Python Bunch Utils Functions definition
	//---------------------------------------------------------
	
	static PyObject* wrap_part_wit_id_sort(PyObject* self, PyObject* args)
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
		int info_x = 0, info_y = 0, info_z = 0;
		if(!PyArg_ParseTuple(args,"OOO|iii:transport_with_twiss_Mtrx",&pyIn,&pyOut,&pyM,&info_x,&info_y,&info_z)){
			error("transport_with_twiss_Mtrx(Bunch in,Bunch out,Matrix, info_x, info_y, info_z) - Bunches and Matrix are needed.");
		}			
		PyObject* pyBunchType = wrap_orbit_bunch::getBunchType("Bunch");
		if((!PyObject_IsInstance(pyIn,pyBunchType)) || (!PyObject_IsInstance(pyOut,pyBunchType)) ){
			error("transport_with_twiss_Mtrx(Bunch in,Bunch out,Matrix, info_x, info_y, info_z) - 1st or 2nd input parameter is not Bunch");
		}	
		PyObject* pyORBIT_Matrix_Type = wrap_orbit_utils::getOrbitUtilsType("Matrix");
		if(!PyObject_IsInstance(pyM,pyORBIT_Matrix_Type)){
			error("transport_with_twiss_Mtrx(Bunch in,Bunch out,Matrix, info_x, info_y, info_z) - function needs a matrix.");
		}		
		Bunch* bunch_in = (Bunch*) ((pyORBIT_Object*) pyIn)->cpp_obj;
		Bunch* bunch_out = (Bunch*) ((pyORBIT_Object*) pyOut)->cpp_obj;
		Matrix* mtrx = (Matrix*) ((pyORBIT_Object*) pyM)->cpp_obj;
		
		int n_stat = transport_mtrx(bunch_in,bunch_out,mtrx,info_x,info_y,info_z);
    return Py_BuildValue("i", n_stat);	
	}			
	
	static PyObject* wrap_copyCoordsToInitCoordsAttr(PyObject* self, PyObject* args)
	{

		PyObject *pyIn;
		if(!PyArg_ParseTuple(args,"O:copyCoordsToInitCoordsAttr",&pyIn)){
			error("copyCoordsToInitCoordsAttr(Bunch) - Bunch is needed.");
		}			
		PyObject* pyBunchType = wrap_orbit_bunch::getBunchType("Bunch");
		if((!PyObject_IsInstance(pyIn,pyBunchType))){
			error("copyCoordsToInitCoordsAttr(Bunch) - input parameter is not Bunch");
		}		
		Bunch* bunch = (Bunch*) ((pyORBIT_Object*) pyIn)->cpp_obj;
		copyCoordsToInitCoordsAttr(bunch);
    Py_INCREF(Py_None);
    return Py_None;	
	}		
	
	static PyObject* wrap_swapInitCoordsAttrAndCoords(PyObject* self, PyObject* args)
	{

		PyObject *pyIn;
		if(!PyArg_ParseTuple(args,"O:swapInitCoordsAttrAndCoords",&pyIn)){
			error("swapInitCoordsAttrAndCoords(Bunch) - Bunch is needed.");
		}			
		PyObject* pyBunchType = wrap_orbit_bunch::getBunchType("Bunch");
		if((!PyObject_IsInstance(pyIn,pyBunchType))){
			error("swapInitCoordsAttrAndCoords(Bunch) - input parameter is not Bunch");
		}		
		Bunch* bunch = (Bunch*) ((pyORBIT_Object*) pyIn)->cpp_obj;
		swapInitCoordsAttrAndCoords(bunch);
    Py_INCREF(Py_None);
    return Py_None;	
	}			

	static PyObject* wrap_transport_mtrx_from_init_coords(PyObject* self, PyObject* args)
	{
		PyObject *pyIn;
		PyObject *pyM;
		int info_x = 0, info_y = 0, info_z = 0;
		if(!PyArg_ParseTuple(args,"OO|iii:transport_with_twiss_Mtrx",&pyIn,&pyM,&info_x,&info_y,&info_z)){
			error("transport_with_twiss_Mtrx_FromInitCoords(Bunch in,Matrix, info_x, info_y, info_z) - Bunches and Matrix are needed.");
		}			
		PyObject* pyBunchType = wrap_orbit_bunch::getBunchType("Bunch");
		if((!PyObject_IsInstance(pyIn,pyBunchType))){
			error("transport_with_twiss_Mtrx_FromInitCoords(Bunch in,Matrix, info_x, info_y, info_z) - 1st input parameter is not Bunch");
		}	
		PyObject* pyORBIT_Matrix_Type = wrap_orbit_utils::getOrbitUtilsType("Matrix");
		if(!PyObject_IsInstance(pyM,pyORBIT_Matrix_Type)){
			error("transport_with_twiss_Mtrx_FromInitCoords(Bunch in,Matrix, info_x, info_y, info_z) - function needs a matrix.");
		}		
		Bunch* bunch = (Bunch*) ((pyORBIT_Object*) pyIn)->cpp_obj;
		Matrix* mtrx = (Matrix*) ((pyORBIT_Object*) pyM)->cpp_obj;
		
		int n_stat = transport_mtrx_from_init_coords(bunch,mtrx,info_x,info_y,info_z);
    return Py_BuildValue("i", n_stat);	
	}		
	
	static PyObject* wrap_bunch_twiss_filtering(PyObject* self, PyObject* args)
	{
		PyObject *pyIn;
		PyObject *pyOut;
		double coeff_x = -1.0, coeff_y = -1.0, coeff_z = -1.0;
		if(!PyArg_ParseTuple(args,"OO|ddd:bunch_twiss_filtering",&pyIn,&pyOut,&coeff_x,&coeff_y,&coeff_z)){
			error("bunchTwissFiltering(Bunch in,Bunch out, coeff_x, coeff_y, coeff_z) - Bunches and Matrix are needed.");
		}			
		PyObject* pyBunchType = wrap_orbit_bunch::getBunchType("Bunch");
		if((!PyObject_IsInstance(pyIn,pyBunchType)) || (!PyObject_IsInstance(pyOut,pyBunchType)) ){
			error("bunchTwissFiltering(Bunch in,Bunch out, coeff_x, coeff_y, coeff_z) - 1st or 2nd input parameter is not Bunch");
		}	
		Bunch* bunch_in = (Bunch*) ((pyORBIT_Object*) pyIn)->cpp_obj;
		Bunch* bunch_bad = (Bunch*) ((pyORBIT_Object*) pyOut)->cpp_obj;
		
		int n_removed = bunch_twiss_filtering(bunch_in,bunch_bad, coeff_x, coeff_y, coeff_z);
    return Py_BuildValue("i", n_removed);	
	}			
	
	// defenition of the memebers of the python bunch_utils_functions wrapper module for functions
	// they will be vailable from python level
	static PyMethodDef BunchUtilsFunctionMethods[] = { 		
		{"bunchSortId",  wrap_part_wit_id_sort    , METH_VARARGS, "bunchSortId(bunch) - Sorting bunch according to the Id particles attributes."},
		{"transportMtrx", wrap_transport_mtrx , METH_VARARGS, "transportMtrx(bunch in, bunch out, matrix) - calculates transport matrix."},
		{"copyCoordsToInitCoordsAttr",wrap_copyCoordsToInitCoordsAttr, METH_VARARGS,"copyCoordsToInitCoordsAttr(bunch) - copy coords to Init Coords Attr."},
		{"swapInitCoordsAttrAndCoords",wrap_swapInitCoordsAttrAndCoords, METH_VARARGS,"swapInitCoordsAttrAndCoords(bunch) - swap coords to Init Coords Attr."},
		{"transportMtrxFromInitCoords", wrap_transport_mtrx_from_init_coords , METH_VARARGS, "transportMtrx(bunch in, matrix) - calculates transport matrix."},
		{"bunchTwissFiltering",  wrap_bunch_twiss_filtering, METH_VARARGS, "bunchTwissFiltering(bunch_in,bunch_bad, x_lim, y_lim, z_lim) - bunch Twiss filtering"},		
		{NULL, NULL, 0, NULL}        /* Sentinel */
	};

	//--------------------------------------------------
	//Initialization function of the Bunch utilities module
	//It will be called from utils wrapper initialization
	//--------------------------------------------------
	
  void initBunchUtilsFunctions(PyObject* module, const char* bunch_utils_module_name){
    //create ==operations with bunches== module
    PyObject* module_local = Py_InitModule(bunch_utils_module_name,BunchUtilsFunctionMethods);
		PyModule_AddObject(module,bunch_utils_module_name,module_local);
		Py_INCREF(module_local);
	}

#ifdef __cplusplus
}
#endif

//end of namespace wrap_utils_bunch_functions
}
