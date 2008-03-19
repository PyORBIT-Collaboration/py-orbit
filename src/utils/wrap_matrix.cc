#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_utils.hh"
#include "wrap_matrix.hh"
#include "wrap_bunch.hh"

#include <iostream>

#include "Matrix.hh"
#include "PhaseVector.hh"
#include "MatrixOperations.hh"

using namespace OrbitUtils;
using namespace wrap_orbit_utils;

namespace wrap_utils_martix{

  void error(const char* msg){ ORBIT_MPI_Finalize(msg); }

#ifdef __cplusplus
extern "C" {
#endif
	//---------------------------------------------------------
	//Python Matrix class definition
	//---------------------------------------------------------

	//constructor for python class wrapping Matrix instance
	//It never will be called directly
	static PyObject* Matrix_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
	{
		pyORBIT_Object* self;
		self = (pyORBIT_Object *) type->tp_alloc(type, 0);
		self->cpp_obj = NULL;
		return (PyObject *) self;
	}

  //initializator for python  Matrix class
  //this is implementation of the __init__ method
  static int Matrix_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){
	  int nArgs = PyTuple_Size(args);
    if(nArgs == 1){
			PyObject* pyIn;
			if(!PyArg_ParseTuple(	args,"O:__init__",&pyIn)){
				error("PyMatrix - Matrix(matrix_parent) - constructor needs a parent matrix.");
			}		
			PyObject* pyORBIT_Matrix_Type = getOrbitUtilsType("Matrix");
			if(!PyObject_IsInstance(pyIn,pyORBIT_Matrix_Type)){
				error("PyMatrix - Matrix(matrix_parent) - constructor needs a parent matrix.");
			}
			Matrix* mtrx = (Matrix*)(((pyORBIT_Object*) pyIn)->cpp_obj);
			self->cpp_obj = new Matrix(mtrx->rows(),mtrx->columns());
			mtrx->copyTo((Matrix*) self->cpp_obj);
			return 0;
    }
		if(nArgs == 2){
			int n,m;
			if(!PyArg_ParseTuple(	args,"ii:__init__",&n,&m)){
				error("PyMatrix - Matrix(n,m) - a maririx size is needed.");
			}		
			self->cpp_obj = new Matrix(n,m);
		}
    return 0;
  }

  //-----------------------------------------------------
  //destructor for python Matrix class (__del__ method).
  //-----------------------------------------------------
  static void Matrix_del(pyORBIT_Object* self){
		//std::cerr<<"The Matrix __del__ has been called!"<<std::endl;
		delete ((Matrix*)self->cpp_obj);
		self->ob_type->tp_free((PyObject*)self);
  }

	//Sets or returns the value of the matrix element with (i,j) raw and column. 
  static PyObject* Matrix_get_set(PyObject *self, PyObject *args){
    //if nVars == 2 this is get value
    //if nVars == 3 this is set value
    int nVars = PyTuple_Size(args);
    pyORBIT_Object* pyMatrix = (pyORBIT_Object*) self;
    double val = 0.;
		int i,j;
		Matrix* cpp_Matrix = (Matrix*) pyMatrix->cpp_obj;
		if(!PyArg_ParseTuple(	args,"ii|d:set",&i,&j,&val)){
			error("PyMatrix - get/set(i,j[,value]) - something is missing.");
		}	
		if(i >= cpp_Matrix->rows() || j >= cpp_Matrix->columns() || i < 0 || j < 0){
			error("PyMatrix - get/set(i,j[,value]) - wrong i or j.");
		}				
    if(nVars == 2 ||  nVars == 3){
      if(nVars == 2){			
        val = cpp_Matrix->value(i,j);
      }
      else{
				cpp_Matrix->value(i,j) = val;
      }
			return Py_BuildValue("d",val);
    }
    else{
      error("PyMatrix. You should call get(i,j) or set(i,j,value)");
    }
    Py_INCREF(Py_None);
    return Py_None;
  }

  //  size() - returns (n,m) tuple
  static PyObject* Matrix_size(PyObject *self, PyObject *args){
    pyORBIT_Object* pyMatrix = (pyORBIT_Object*) self;
		Matrix* cpp_Matrix = (Matrix*) pyMatrix->cpp_obj;
		return Py_BuildValue("(i,i)",cpp_Matrix->rows(),cpp_Matrix->columns());
  }

	//  copy() - returns the new Matrix that is equal to the parent one
  static PyObject* Matrix_copy(PyObject *self, PyObject *args){
    pyORBIT_Object* pyMatrix = (pyORBIT_Object*) self;
		Matrix* cpp_Matrix = (Matrix*) pyMatrix->cpp_obj;
		int n = cpp_Matrix->rows();
		int m = cpp_Matrix->columns();
		PyObject* mod = PyImport_ImportModule("orbit_utils");
		PyObject* pyMatr = PyObject_CallMethod(mod,"Matrix","ii",n,m);
		pyORBIT_Object* pyMatrix_child = (pyORBIT_Object*) pyMatr;
		Matrix* cpp_Matrix_child = (Matrix*) pyMatrix_child->cpp_obj;
		if(!cpp_Matrix->copyTo(cpp_Matrix_child)){
			error("PyMatrix. Cannot copy the matrix.");
		}
		Py_DECREF(mod);
    return pyMatr;	
  }
	
	//  copyTo() - copy values of the Matrix into another one
  static PyObject* Matrix_copyTo(PyObject *self, PyObject *args){
    pyORBIT_Object* pyMatrix = (pyORBIT_Object*) self;
		Matrix* cpp_Matrix = (Matrix*) pyMatrix->cpp_obj;
		PyObject *pyMatr;
		if(!PyArg_ParseTuple(args,"O:copyTo",&pyMatr)){
			error("PyMatrix - copyTo(matrix) - the target matrix is needed.");
		}			
		pyORBIT_Object* pyMatrix_target = (pyORBIT_Object*) pyMatr;
		Matrix* cpp_Matrix_target = (Matrix*) pyMatrix_target->cpp_obj;
		int res = cpp_Matrix->copyTo(cpp_Matrix_target);
    return Py_BuildValue("i",res);	
  }
	
	//  transpose() - transpose the matrix on the place
  static PyObject* Matrix_transpose(PyObject *self, PyObject *args){
		((Matrix*) ((pyORBIT_Object*) self)->cpp_obj)->transpose();
		Py_INCREF(Py_None);
    return Py_None;	
  }	
	
	//  zero() - sets all elements of the matrix to 0
  static PyObject* Matrix_zero(PyObject *self, PyObject *args){
		((Matrix*) ((pyORBIT_Object*) self)->cpp_obj)->zero();
		Py_INCREF(Py_None);
    return Py_None;	
  }	
	
	//  unit() - makes the matrix equals to a unit matrix
  static PyObject* Matrix_unit(PyObject *self, PyObject *args){
		int res = ((Matrix*) ((pyORBIT_Object*) self)->cpp_obj)->unit();
    return Py_BuildValue("i",res);;	
  }		
	
	//  invert() - returns the inverted matrix, or None
  static PyObject* Matrix_invert(PyObject *self, PyObject *args){
		PyObject* mod = PyImport_ImportModule("orbit_utils");
		PyObject* pyMtrx = PyObject_CallMethod(mod,"Matrix","O",self);
		Matrix* cpp_mtrx = (Matrix*) ((pyORBIT_Object*) pyMtrx)->cpp_obj;
		Py_DECREF(mod);
		int res = MatrixOperations::invert(cpp_mtrx);
		if(res == 0){
			Py_DECREF(pyMtrx);
			Py_INCREF(Py_None);
			return Py_None;				
		}
		return pyMtrx;
  }	
	
	//  add() - adds a number or matrix to the this matrix. Returns a new matrix 
  static PyObject* Matrix_add(PyObject *self, PyObject *args){
		PyObject *pyIn;
		if(!PyArg_ParseTuple(args,"O:add",&pyIn)){
			error("PyMatrix - add(number or matrix) - input parameter is needed.");
		}			
		if(PyNumber_Check(pyIn)){
			double val;
			if(!PyArg_ParseTuple(args,"d:add",&val)){
				error("PyMatrix - add(number) - input parameter is needed.");
			}		
			PyObject* mod = PyImport_ImportModule("orbit_utils");
			PyObject* pyMtrx = PyObject_CallMethod(mod,"Matrix","O",self);
			Matrix* cpp_mtrx = (Matrix*) ((pyORBIT_Object*) pyMtrx)->cpp_obj;
			cpp_mtrx->add(val);
			Py_DECREF(mod);
			return pyMtrx;			
		}
		PyObject* pyORBIT_Matrix_Type = getOrbitUtilsType("Matrix");
		if(PyObject_IsInstance(pyIn,pyORBIT_Matrix_Type)){
			PyObject* mod = PyImport_ImportModule("orbit_utils");
			PyObject* pyMtrx = PyObject_CallMethod(mod,"Matrix","O",self);
			Matrix* cpp_mtrx = (Matrix*) ((pyORBIT_Object*) pyMtrx)->cpp_obj;
			Py_DECREF(mod);
			cpp_mtrx->add((Matrix*)(((pyORBIT_Object*) pyIn)->cpp_obj));
			return pyMtrx;
		}				
		error("PyPhaseMaatrix - add(number or matrix) - input parameter is wrong.");	
		Py_INCREF(Py_None);
    return Py_None;		
  }
	
	//  mult() - multiplies a number,matrix or vector to the this matrix. Returns a new matrix. 
  static PyObject* Matrix_mult(PyObject *self, PyObject *args){
    pyORBIT_Object* pyMatrix = (pyORBIT_Object*) self;
		Matrix* cpp_Matrix = (Matrix*) pyMatrix->cpp_obj;
		PyObject *pyIn;
		if(!PyArg_ParseTuple(args,"O:mult",&pyIn)){
			error("PyMatrix - mult(number,matrix, or vector) - input parameter is needed.");
		}			
		if(PyNumber_Check(pyIn)){
			double val;
			if(!PyArg_ParseTuple(args,"d:mult",&val)){
				error("PyMatrix - mult(number) - input parameter is needed.");
			}		
			PyObject* mod = PyImport_ImportModule("orbit_utils");
			PyObject* pyMtrx = PyObject_CallMethod(mod,"Matrix","O",self);
			Matrix* cpp_mtrx = (Matrix*) ((pyORBIT_Object*) pyMtrx)->cpp_obj;
			cpp_mtrx->mult(val);
			Py_DECREF(mod);
			return pyMtrx;
		}
		PyObject* pyORBIT_Matrix_Type = getOrbitUtilsType("Matrix");
		if(PyObject_IsInstance(pyIn,pyORBIT_Matrix_Type)){
			PyObject* mod = PyImport_ImportModule("orbit_utils");
			PyObject* pyMtrx = PyObject_CallMethod(mod,"Matrix","O",self);
			Matrix* cpp_mtrx = (Matrix*) ((pyORBIT_Object*) pyMtrx)->cpp_obj;
			Py_DECREF(mod);
			cpp_mtrx->mult((Matrix*)(((pyORBIT_Object*) pyIn)->cpp_obj));
			return pyMtrx;
		}		
		PyObject* pyORBIT_PhaseVector_Type = getOrbitUtilsType("PhaseVector");
		if(PyObject_IsInstance(pyIn,pyORBIT_PhaseVector_Type)){
			PhaseVector* cpp_PhaseVector = (PhaseVector*) ((pyORBIT_Object*) pyIn)->cpp_obj;
			if(cpp_Matrix->columns() != cpp_PhaseVector->size()){
				error("PyMatrix - mult(vector) - unequal sizes of columns and vectors.");
			}
			PyObject* mod = PyImport_ImportModule("orbit_utils");
			PyObject* pyVctr = PyObject_CallMethod(mod,"PhaseVector","i",cpp_Matrix->rows());
			MatrixOperations::mult(cpp_Matrix,cpp_PhaseVector,(PhaseVector*) ((pyORBIT_Object*)pyVctr)->cpp_obj);
			Py_DECREF(mod);
			return pyVctr;
		}
		error("PyPhaseMaatrix - mult(number, matrix, or vector) - input parameter is wrong.");	
		Py_INCREF(Py_None);
    return Py_None;	
  }	
	
	//  track() - trackiplies a number,matrix or vector to the this matrix. Returns a new matrix. 
  static PyObject* Matrix_track(PyObject *self, PyObject *args){
    pyORBIT_Object* pyMatrix = (pyORBIT_Object*) self;
		Matrix* cpp_Matrix = (Matrix*) pyMatrix->cpp_obj;
		PyObject *pyIn;
		if(!PyArg_ParseTuple(args,"O:track",&pyIn)){
			error("PyMatrix - track(Bunch) - Bunch is needed.");
		}			
		PyObject* pyBunchType = wrap_orbit_bunch::getBunchType("Bunch");
		if((!PyObject_IsInstance(pyIn,pyBunchType))){
			error("PyMatrix - track(Bunch) - input parameter is not Bunch");
		}		
		Bunch* bunch = (Bunch*) ((pyORBIT_Object*) pyIn)->cpp_obj;
		MatrixOperations::track(bunch,cpp_Matrix);	
		Py_INCREF(Py_None);
    return Py_None;	
  }	
	
	// defenition of the methods of the python Matrix wrapper class
	// they will be vailable from python level
  static PyMethodDef MatrixClassMethods[] = {
    { "size",       Matrix_size      ,METH_VARARGS,"Returns tuple with (n,m)"},
    { "get",        Matrix_get_set   ,METH_VARARGS,"Returns valuie(i,j)"},
    { "set",        Matrix_get_set   ,METH_VARARGS,"Sets the new value to (i,j) element - set (i,j,val)"},
		{ "copy",       Matrix_copy      ,METH_VARARGS,"Returns the copy of the matrix"},
		{ "copyTo",     Matrix_copyTo    ,METH_VARARGS,"Copy the matrix to the target matrix"},
		{ "transpose",  Matrix_transpose ,METH_VARARGS,"Transposes the matrix on the place."},
		{ "invert",     Matrix_invert    ,METH_VARARGS,"Returns the inverted the matrix, or None."},
		{ "zero",       Matrix_zero      ,METH_VARARGS,"Sets all elements to 0."},
		{ "unit",       Matrix_unit      ,METH_VARARGS,"Sets the matrix to a unit matrix"},
		{ "add",        Matrix_add       ,METH_VARARGS,"Adds a matrix or number to the matrix"},
		{ "mult",       Matrix_mult      ,METH_VARARGS,"Multiples a matrix by number, vector or matrix. Returns a new matrix"},
		{ "track",      Matrix_track     ,METH_VARARGS,"Tracks Bunch through the matrix"},
    {NULL}
  };

	// defenition of the memebers of the python Matrix wrapper class
	// they will be vailable from python level
	static PyMemberDef MatrixClassMembers [] = {
		{NULL}
	};

	//new python Matrix wrapper type definition
	static PyTypeObject pyORBIT_Matrix_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"Matrix", /*tp_name*/
		sizeof(pyORBIT_Object), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) Matrix_del , /*tp_dealloc*/
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
		"The Matrix python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		MatrixClassMethods, /* tp_methods */
		MatrixClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) Matrix_init, /* tp_init */
		0, /* tp_alloc */
		Matrix_new, /* tp_new */
	};

	//--------------------------------------------------
	//Initialization function of the pyMatrix class
	//It will be called from Bunch wrapper initialization
	//--------------------------------------------------
  void initMatrix(PyObject* module){
		if (PyType_Ready(&pyORBIT_Matrix_Type) < 0) return;
		Py_INCREF(&pyORBIT_Matrix_Type);
		PyModule_AddObject(module, "Matrix", (PyObject *)&pyORBIT_Matrix_Type);
	}

#ifdef __cplusplus
}
#endif

//end of namespace wrap_utils_martix
}
