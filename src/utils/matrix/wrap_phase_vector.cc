#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_utils.hh"
#include "wrap_phase_vector.hh"

#include <iostream>

#include "PhaseVector.hh"
#include "Matrix.hh"
#include "MatrixOperations.hh"

using namespace OrbitUtils;
using namespace wrap_orbit_utils;

namespace wrap_utils_phase_vector{

  void error(const char* msg){ ORBIT_MPI_Finalize(msg); }

#ifdef __cplusplus
extern "C" {
#endif

	//---------------------------------------------------------
	//Python PhaseVector class definition
	//---------------------------------------------------------

	//constructor for python class wrapping PhaseVector instance
	//It never will be called directly
	static PyObject* PhaseVector_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
	{
		pyORBIT_Object* self;
		self = (pyORBIT_Object *) type->tp_alloc(type, 0);
		self->cpp_obj = NULL;
		return (PyObject *) self;
	}

  //initializator for python  PhaseVector class
  //this is implementation of the __init__ method
  static int PhaseVector_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){
		PyObject* pyIn;
		if(!PyArg_ParseTuple(args,"O:__init__",&pyIn)){
			error("PyPhaseVector - PhaseVector(vector or size) - constructor needs a parameter.");
		}		
		if(PyNumber_Check(pyIn)){
			int size;
			if(!PyArg_ParseTuple(args,"i:__init__",&size)){
				error("PyPhaseVector - __init__(size) - input parameter is needed.");
			}		
			self->cpp_obj = new PhaseVector(size);
			return 0;
		}			
		PyObject* pyORBIT_PhaseVector_Type = getOrbitUtilsType("PhaseVector");
		if(!PyObject_IsInstance(pyIn,pyORBIT_PhaseVector_Type)){
			error("PyPhaseVector - PhaseVector(vector_parent) - constructor needs a parent matrix.");
		}
		PhaseVector* v = (PhaseVector*)(((pyORBIT_Object*) pyIn)->cpp_obj);
		self->cpp_obj = new PhaseVector(v->size());
		v->copyTo((PhaseVector*) self->cpp_obj);
    return 0;
  }

  //-----------------------------------------------------
  //destructor for python PhaseVector class (__del__ method).
  //-----------------------------------------------------
  static void PhaseVector_del(pyORBIT_Object* self){
		//std::cerr<<"The PhaseVector __del__ has been called!"<<std::endl;
		delete ((PhaseVector*)self->cpp_obj);
		self->ob_type->tp_free((PyObject*)self);
  }

	//Sets or returns the value of the phase vector element with particular index. 
  static PyObject* PhaseVector_get_set(PyObject *self, PyObject *args){
    //if nVars == 1 this is get value
    //if nVars == 2 this is set value
    int nVars = PyTuple_Size(args);
    pyORBIT_Object* pyPhaseVector = (pyORBIT_Object*) self;
    double val = 0.;
		int i;
		PhaseVector* cpp_PhaseVector = (PhaseVector*) pyPhaseVector->cpp_obj;
		if(!PyArg_ParseTuple(	args,"i|d:set",&i,&val)){
			error("PyPhaseVector - get/set(i[,value]) - something is missing.");
		}	
		if(i >= cpp_PhaseVector->size() || i < 0){
			error("PyPhaseVector - get/set(i[,value]) - wrong i or j.");
		}		
    if(nVars == 1 ||  nVars == 2){
      if(nVars == 1){			
        val = cpp_PhaseVector->value(i);
      }
      else{
				cpp_PhaseVector->value(i) = val;
      }
			return Py_BuildValue("d",val);
    }
    else{
      error("PyPhaseVector. You should call get(i) or set(i,value)");
    }
    Py_INCREF(Py_None);
    return Py_None;
  }

  //  size() - returns n in vector(n)
  static PyObject* PhaseVector_size(PyObject *self, PyObject *args){
    pyORBIT_Object* pyPhaseVector = (pyORBIT_Object*) self;
		PhaseVector* cpp_PhaseVector = (PhaseVector*) pyPhaseVector->cpp_obj;
		return Py_BuildValue("i",cpp_PhaseVector->size());
  }

	//  copy() - returns the new PhaseVector that is equal to the parent one
  static PyObject* PhaseVector_copy(PyObject *self, PyObject *args){
    pyORBIT_Object* pyPhaseVector = (pyORBIT_Object*) self;
		PhaseVector* cpp_PhaseVector = (PhaseVector*) pyPhaseVector->cpp_obj;
		int n = cpp_PhaseVector->size();
		PyObject* mod = PyImport_ImportModule("orbit_utils");
		PyObject* pyVctr = PyObject_CallMethod(mod,"PhaseVector","i",n);
		pyORBIT_Object* pyPhaseVector_child = (pyORBIT_Object*) pyVctr;
		PhaseVector* cpp_PhaseVector_child = (PhaseVector*) pyPhaseVector_child->cpp_obj;
		if(!cpp_PhaseVector->copyTo(cpp_PhaseVector_child)){
			error("PyPhaseVector. Cannot copy the vector.");
		}
		Py_DECREF(mod);
    return pyVctr;	
  }
	
	//  copyTo() - copy values of the PhaseVector into another one
  static PyObject* PhaseVector_copyTo(PyObject *self, PyObject *args){
    pyORBIT_Object* pyPhaseVector = (pyORBIT_Object*) self;
		PhaseVector* cpp_PhaseVector = (PhaseVector*) pyPhaseVector->cpp_obj;
		PyObject *pyMatr;
		if(!PyArg_ParseTuple(args,"O:copyTo",&pyMatr)){
			error("PyPhaseVector - copyTo(matrix) - the target vector is needed.");
		}			
		pyORBIT_Object* pyPhaseVector_target = (pyORBIT_Object*) pyMatr;
		PhaseVector* cpp_PhaseVector_target = (PhaseVector*) pyPhaseVector_target->cpp_obj;
		int res = cpp_PhaseVector->copyTo(cpp_PhaseVector_target);
    return Py_BuildValue("i",res);	
  }
		
	//  zero() - sets all elements of the vector to 0
  static PyObject* PhaseVector_zero(PyObject *self, PyObject *args){
		((PhaseVector*) ((pyORBIT_Object*) self)->cpp_obj)->zero();
		Py_INCREF(Py_None);
    return Py_None;	
  }	
	
	//  norm() - returns norm of the vector
  static PyObject* PhaseVector_norm(PyObject *self, PyObject *args){
		double val = ((PhaseVector*) ((pyORBIT_Object*) self)->cpp_obj)->norm();
    return Py_BuildValue("d",val);;	
  }		
		
	//  add() - adds a number or vector to the this vector, returns a new vector
  static PyObject* PhaseVector_add(PyObject *self, PyObject *args){
    pyORBIT_Object* pyPhaseVector = (pyORBIT_Object*) self;
		PhaseVector* cpp_PhaseVector = (PhaseVector*) pyPhaseVector->cpp_obj;
		PyObject *pyIn;
		if(!PyArg_ParseTuple(args,"O:add",&pyIn)){
			error("PyPhaseVector - add(number or vector) - input parameter is needed.");
		}
		if(PyNumber_Check(pyIn)){
			double val;
			if(!PyArg_ParseTuple(args,"d:add",&val)){
				error("Py PhaseVector- add(number) - input parameter is needed.");
			}	
			PyObject* mod = PyImport_ImportModule("orbit_utils");
			PyObject* pyVctr = PyObject_CallMethod(mod,"PhaseVector","O",self);
			PhaseVector* cpp_vctr = (PhaseVector*) ((pyORBIT_Object*) pyVctr)->cpp_obj;
			cpp_vctr->add(val);
			Py_DECREF(mod);
			return pyVctr;
		}			
		PyObject* pyORBIT_PhaseVector_Type = getOrbitUtilsType("PhaseVector");
		if(PyObject_IsInstance(pyIn,pyORBIT_PhaseVector_Type)){
			PhaseVector* cpp_PhaseVector_In = (PhaseVector*) ((pyORBIT_Object*) pyIn)->cpp_obj;
			if(cpp_PhaseVector->size() != cpp_PhaseVector_In->size()){
				error("PyPhaseVector - add(vector) - unequal sizes of the vectors.");
			}
			PyObject* mod = PyImport_ImportModule("orbit_utils");
			PyObject* pyVctr = PyObject_CallMethod(mod,"PhaseVector","O",self);
			PhaseVector* cpp_vctr = (PhaseVector*) ((pyORBIT_Object*) pyVctr)->cpp_obj;
			cpp_vctr->add(cpp_PhaseVector_In);
			Py_DECREF(mod);
			return pyVctr;
		}	
		error("PyPhaseVector - add(number or vector) - input parameter is wrong.");	
		Py_INCREF(Py_None);
    return Py_None;			
  }
	
	//  mult() - mults the vector to number or a matrix. will return the new vector 
  static PyObject* PhaseVector_mult(PyObject *self, PyObject *args){
    pyORBIT_Object* pyPhaseVector = (pyORBIT_Object*) self;
		PhaseVector* cpp_PhaseVector = (PhaseVector*) pyPhaseVector->cpp_obj;
		PyObject *pyIn;
		if(!PyArg_ParseTuple(args,"O:mult",&pyIn)){
			error("PyPhaseVector - mult(number or vector or matrix) - input parameter is needed.");
		}			
		if(PyNumber_Check(pyIn)){
			double val;
			if(!PyArg_ParseTuple(args,"d:mult",&val)){
				error("Py PhaseVector- mult(number) - input parameter is needed.");
			}	
			PyObject* mod = PyImport_ImportModule("orbit_utils");
			PyObject* pyVctr = PyObject_CallMethod(mod,"PhaseVector","O",self);
			PhaseVector* cpp_vctr = (PhaseVector*) ((pyORBIT_Object*) pyVctr)->cpp_obj;
			cpp_vctr->mult(val);
			Py_DECREF(mod);
			return pyVctr;
		}		
		PyObject* pyORBIT_PhaseVector_Type = getOrbitUtilsType("PhaseVector");
		if(PyObject_IsInstance(pyIn,pyORBIT_PhaseVector_Type)){
			PhaseVector* cpp_PhaseVector_In = (PhaseVector*) ((pyORBIT_Object*) pyIn)->cpp_obj;
			if(cpp_PhaseVector->size() != cpp_PhaseVector_In->size()){
				error("PyPhaseVector - mult(vector) - unequal sizes of the vectors.");
			}
			return Py_BuildValue("d",cpp_PhaseVector->dot(cpp_PhaseVector_In));
		}
		PyObject* pyORBIT_Matrix_Type = getOrbitUtilsType("Matrix");
		if(PyObject_IsInstance(pyIn,pyORBIT_Matrix_Type)){
			Matrix* cpp_Matrix = (Matrix*) ((pyORBIT_Object*) pyIn)->cpp_obj;
			if(cpp_PhaseVector->size() != cpp_Matrix->rows()){
				error("PyPhaseVector - mult(matrix) - unequal sizes of vectors and rows.");
			}
			PyObject* mod = PyImport_ImportModule("orbit_utils");
			PyObject* pyVctr = PyObject_CallMethod(mod,"PhaseVector","i",cpp_Matrix->columns());
			MatrixOperations::mult(cpp_PhaseVector,cpp_Matrix,(PhaseVector*) ((pyORBIT_Object*)pyVctr)->cpp_obj);
			Py_DECREF(mod);
			return pyVctr;
		}		
		error("PyPhaseVector - mult(number or vector or matrix) - input parameter is wrong.");	
		Py_INCREF(Py_None);
    return Py_None;		
  }		
		
	// defenition of the methods of the python PhaseVector wrapper class
	// they will be vailable from python level
  static PyMethodDef PhaseVectorClassMethods[] = {
    { "size",       PhaseVector_size      ,METH_VARARGS,"Returns vector's number of elements"},
    { "get",        PhaseVector_get_set   ,METH_VARARGS,"Returns valuie(i)"},
    { "set",        PhaseVector_get_set   ,METH_VARARGS,"Sets the new value to (i) element - set (i,val)"},
		{ "copy",       PhaseVector_copy      ,METH_VARARGS,"Returns the copy of the vector"},
		{ "copyTo",     PhaseVector_copyTo    ,METH_VARARGS,"Copy the vector to the target vector"},
		{ "zero",       PhaseVector_zero      ,METH_VARARGS,"Sets all elements to 0."},
		{ "norm",       PhaseVector_norm      ,METH_VARARGS,"Returns the length of the vector"},
		{ "add",        PhaseVector_add       ,METH_VARARGS,"Adds a vector or number to the vector"},
		{ "mult",       PhaseVector_mult      ,METH_VARARGS,"Multiples a vector by a number, vector or matrix"},
    {NULL}
  };

	// defenition of the memebers of the python PhaseVector wrapper class
	// they will be vailable from python level
	static PyMemberDef PhaseVectorClassMembers [] = {
		{NULL}
	};

	//new python PhaseVector wrapper type definition
	static PyTypeObject pyORBIT_PhaseVector_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"PhaseVector", /*tp_name*/
		sizeof(pyORBIT_Object), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) PhaseVector_del , /*tp_dealloc*/
		0, /*tp_print*/
		0, /*tp_getattr*/
		0, /*tp_setattr*/
		0, /*tp_compare*/
		0, /*tp_repr*/
		0, /*tp_as_number*/
		0, /*tp_as_sequence*/
		0, /*tp_as_mapping*/
		0, /*tp_hash */
		0, /*tp_call*/
		0, /*tp_str*/
		0, /*tp_getattro*/
		0, /*tp_setattro*/
		0, /*tp_as_buffer*/
		Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
		"The PhaseVector python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		PhaseVectorClassMethods, /* tp_methods */
		PhaseVectorClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) PhaseVector_init, /* tp_init */
		0, /* tp_alloc */
		PhaseVector_new, /* tp_new */
	};	

	//--------------------------------------------------
	//Initialization function of the pyPhaseVector class
	//It will be called from Bunch wrapper initialization
	//--------------------------------------------------
  void initPhaseVector(PyObject* module){
		if (PyType_Ready(&pyORBIT_PhaseVector_Type) < 0) return;
		Py_INCREF(&pyORBIT_PhaseVector_Type);
		PyModule_AddObject(module, "PhaseVector", (PyObject *)&pyORBIT_PhaseVector_Type);
	}

#ifdef __cplusplus
}
#endif

//end of namespace wrap_utils_phase_vector
}
