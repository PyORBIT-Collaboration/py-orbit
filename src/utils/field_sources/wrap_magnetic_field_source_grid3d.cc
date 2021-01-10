#include <iostream>

#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"
#include "wrap_spacecharge.hh"

#include "wrap_utils.hh"
#include "MagnetFieldSourceGrid3D.hh"

using namespace OrbitUtils;
using namespace wrap_orbit_utils;

namespace wrap_field_source_grid3d{

  void error(const char* msg){ ORBIT_MPI_Finalize(msg); }

#ifdef __cplusplus
extern "C" {
#endif
	//---------------------------------------------------------
	//Python MagnetFieldSourceGrid3D class definition
	//---------------------------------------------------------

	//constructor for python class wrapping MagnetFieldSourceGrid3D instance
	//It never will be called directly
	static PyObject* MagnetFieldSourceGrid3D_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
	{
		pyORBIT_Object* self;
		self = (pyORBIT_Object *) type->tp_alloc(type, 0);
		self->cpp_obj = NULL;
		return (PyObject *) self;
	}

  //initializator for python  MagnetFieldSourceGrid3D class
  //this is implementation of the __init__ method
  static int MagnetFieldSourceGrid3D_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){
 		PyObject* pyBxGrid3D;
 		PyObject* pyByGrid3D;
 		PyObject* pyBzGrid3D;
		if(!PyArg_ParseTuple(	args,"OOO:__init__",&pyBxGrid3D,&pyByGrid3D,&pyBzGrid3D)){
			error("MagnetFieldSourceGrid3D(grid3D,grid3D,grid3D) - constructor needs 3x Grid3D with Bx,By,Bz fields.");
		}
		// check type
		PyObject* pyORBIT_Grid3D_Type = getSpaceChargeType("Grid3D");		
		if(!PyObject_IsInstance(pyBxGrid3D,pyORBIT_Grid3D_Type) 
			 || !PyObject_IsInstance(pyByGrid3D,pyORBIT_Grid3D_Type) 
		   || !PyObject_IsInstance(pyBzGrid3D,pyORBIT_Grid3D_Type)){
			error("MagnetFieldSourceGrid3D(grid3D,grid3D,grid3D) - constructor needs 3x Grid3D");
		}
		//make new C++ instance of MagnetFieldSourceGrid3D
		Grid3D* BxGrid3D = (Grid3D*) ((pyORBIT_Object*) pyBxGrid3D)->cpp_obj;
		Grid3D* ByGrid3D = (Grid3D*) ((pyORBIT_Object*) pyByGrid3D)->cpp_obj;
		Grid3D* BzGrid3D = (Grid3D*) ((pyORBIT_Object*) pyBzGrid3D)->cpp_obj;
		self->cpp_obj = new MagnetFieldSourceGrid3D(BxGrid3D,ByGrid3D,BzGrid3D);
		Py_INCREF(pyBxGrid3D);
		Py_INCREF(pyByGrid3D);
		Py_INCREF(pyBzGrid3D);	
		((MagnetFieldSourceGrid3D*) self->cpp_obj)->setPyWrapper((PyObject*) self);
    return 0;
  }

  /** Sets or returns X,Y,Z axis symmetries */
  static PyObject* MagnetFieldSourceGrid3D_symmetry(PyObject *self, PyObject *args){
	  MagnetFieldSourceGrid3D* cpp_fieldSource = (MagnetFieldSourceGrid3D*)((pyORBIT_Object*) self)->cpp_obj;
	  int nArgs = PyTuple_Size(args);
	  int symmetry_x;
	  int symmetry_y;
	  int symmetry_z;	  
	  if(nArgs == 3){
	  	if(!PyArg_ParseTuple(args,"iii:symmetry",&symmetry_x,&symmetry_y,&symmetry_z)){
	  		error("MagnetFieldSourceGrid3D.symmetry(symmetry_x,symmetry_y,symmetry_z) - params needed.");
	  	}
	  	cpp_fieldSource->setSymmetry(symmetry_x,symmetry_y,symmetry_z);
	  	Py_INCREF(Py_None);
	  	return Py_None;
	  }
	  if(nArgs == 0){
	  	cpp_fieldSource->getSymmetry(symmetry_x,symmetry_y,symmetry_z);
	  	return Py_BuildValue("(iii)",symmetry_x,symmetry_y,symmetry_z);
	  }
	  error("MagnetFieldSourceGrid3D.symmetry(...) - 0 or 3 int params are needed.");
	  Py_INCREF(Py_None);
	  return Py_None;		
  }	 
  
  /** Sets - Returns signs for fields in different quadrants that defined by signs of  signX, signY, signZ */
  static PyObject* MagnetFieldSourceGrid3D_fieldSignsForQuadrants(PyObject *self, PyObject *args){
  	MagnetFieldSourceGrid3D* cpp_fieldSource = (MagnetFieldSourceGrid3D*)((pyORBIT_Object*) self)->cpp_obj;
  	int nArgs = PyTuple_Size(args);
  	int signX;
  	int signY;
  	int signZ;
  	int signBx;
  	int signBy;
  	int signBz;
  	if(nArgs == 3){
  		if(!PyArg_ParseTuple(args,"iii:fieldSignForQuadrant",&signX,&signY,&signZ)){
	  		error("MagnetFieldSourceGrid3D.fieldSignForQuadrant(signX,signY,signZ) - params for quadrant indexes are needed.");
	  	}
	  	cpp_fieldSource->getFieldSignsForQuadrants(signX,signY,signZ,signBx,signBy,signBz);
	  	return Py_BuildValue("(iii)",signBx,signBy,signBz);
  	}
  	if(nArgs == 6){
  		if(!PyArg_ParseTuple(args,"iiiiii:fieldSignForQuadrant",&signX,&signY,&signZ,&signBx,&signBy,&signBz)){
	  		error("MagnetFieldSourceGrid3D.fieldSignForQuadrant(signX,signY,signZ,signBx,signBy,signBz) - params for quadrant indexes and field signs are needed.");
	  	} 
	  	cpp_fieldSource->setFieldSignsForQuadrants(signX,signY,signZ,signBx,signBy,signBz);
	  	cpp_fieldSource->getFieldSignsForQuadrants(signX,signY,signZ,signBx,signBy,signBz);
	  	return Py_BuildValue("(iii)",signBx,signBy,signBz);
  	}
	  error("MagnetFieldSourceGrid3D.fieldSignForQuadrant(...) - 3 or 6 int params are needed.");
	  Py_INCREF(Py_None);
	  return Py_None;	  	
  }

  /** Sets or Returns the scaling coefficient for inner fields in Grid3D instances */
  static PyObject* MagnetFieldSourceGrid3D_fieldCoeff(PyObject *self, PyObject *args){
	  MagnetFieldSourceGrid3D* cpp_fieldSource = (MagnetFieldSourceGrid3D*)((pyORBIT_Object*) self)->cpp_obj;
	  int nArgs = PyTuple_Size(args);
	  double field_coeff;  
	  if(nArgs == 1){
	  	if(!PyArg_ParseTuple(args,"d:fieldCoeff",&field_coeff)){
	  		error("MagnetFieldSourceGrid3D.fieldCoeff(field_coeff) - param is needed.");
	  	}
	  	cpp_fieldSource->setFieldCoeff(field_coeff);
	  }
	  field_coeff = cpp_fieldSource->getFieldCoeff();
	  return Py_BuildValue("d",field_coeff);
  }	  

   /** Returns the Ex,Ey,Ez, Bx,By,Bz fields at (x,y,z) point */
  static PyObject* MagnetFieldSourceGrid3D_getFields(PyObject *self, PyObject *args){
  	MagnetFieldSourceGrid3D* cpp_fieldSource = (MagnetFieldSourceGrid3D*)((pyORBIT_Object*) self)->cpp_obj;
  	double x,y,z;
  	if(!PyArg_ParseTuple(args,"ddd:getFields",&x,&y,&z)){
  		ORBIT_MPI_Finalize("MagnetFieldSourceGrid3D.getFields(x,y,z) - params needed.");
  	}
  	double fe_x; double fe_y; double fe_z;
  	double fm_x; double fm_y; double fm_z;
  	double t = 0.;
  	cpp_fieldSource->getElectricMagneticField(x,y,z,t,fe_x,fe_y,fe_z,fm_x,fm_y,fm_z);
  	return Py_BuildValue("(dddddd)",fe_x,fe_y,fe_z,fm_x,fm_y,fm_z);
  }
  
  /** Sets / Returns the coordinates transformation matrix 4x4 from external to inner system */
  static PyObject* MagnetFieldSourceGrid3D_transormfMatrix(PyObject *self, PyObject *args){
	  MagnetFieldSourceGrid3D* cpp_fieldSource = (MagnetFieldSourceGrid3D*)((pyORBIT_Object*) self)->cpp_obj;
	  int nArgs = PyTuple_Size(args);
	  PyObject* pyMatrix;
	  Matrix* cpp_matrix;
	  if(nArgs == 1){
	  	if(!PyArg_ParseTuple(args,"O:transormfMatrix",&pyMatrix)){
	  		error("MagnetFieldSourceGrid3D.transormfMatrix(Matrix) - parameter is needed.");
	  	}
	  	PyObject* pyORBIT_Matrix_Type = getOrbitUtilsType("Matrix");
	  	if(!PyObject_IsInstance(pyMatrix,pyORBIT_Matrix_Type)){
	  		error("MagnetFieldSourceGrid3D.transormfMatrix(Matrix) - parameter is not Matrix.");
	  	}
	  	cpp_matrix = (Matrix*) ((pyORBIT_Object*) pyMatrix)->cpp_obj;
	  	if(cpp_matrix->rows() != 4 || cpp_matrix->columns() != 4){
	  		error("MagnetFieldSourceGrid3D.transormfMatrix(Matrix) - Matrix is not 4x4.");
	  	}
	  	// the Py_INCREF(pyMatrix) call will be performed inside setCoordsTransformMatrix(...) method 	  	
	  	cpp_fieldSource->setCoordsTransformMatrix(cpp_matrix);
	  	Py_INCREF(Py_None);
	  	return Py_None;			  	
	  }
	  cpp_matrix = cpp_fieldSource->getCoordsTransformMatrix();
	  pyMatrix = (PyObject*) ((pyORBIT_Object*) cpp_matrix->getPyWrapper());
	  if(pyMatrix == NULL){
	  	error("MagnetFieldSourceGrid3D.transormfMatrix() - cannot return Matrix 4x4. You have to assign it first.");
	  }
	  Py_INCREF(pyMatrix);
	  return pyMatrix;
  }    
  
  //-----------------------------------------------------
  //destructor for python MagnetFieldSourceGrid3D class (__del__ method).
  //-----------------------------------------------------
  static void MagnetFieldSourceGrid3D_del(pyORBIT_Object* self){
		delete ((MagnetFieldSourceGrid3D*)self->cpp_obj);
		self->ob_type->tp_free((PyObject*)self);
  }
	
	// defenition of the methods of the python MagnetFieldSourceGrid3D wrapper class
	// they will be vailable from python level
  static PyMethodDef MagnetFieldSourceGrid3DClassMethods[] = {
    { "symmetry",              MagnetFieldSourceGrid3D_symmetry              ,METH_VARARGS,"Sets or returns X,Y,Z axis symmetries as 0 or 1 numbers - eg. (0,0,0)"},
    { "getFields",             MagnetFieldSourceGrid3D_getFields             ,METH_VARARGS,"Returns E and B fields (Ex,Ey,Ez,Bx,By,Bz) for (x,y,z) point"},
    { "fieldSignsForQuadrants",MagnetFieldSourceGrid3D_fieldSignsForQuadrants,METH_VARARGS, "Sets or returns field signs for quadrants x-y-z"},
    { "transormfMatrix",       MagnetFieldSourceGrid3D_transormfMatrix,       METH_VARARGS, "Sets or returns the coordinates transformation matrix 4x4 from external to inner system"},
    { "fieldCoeff",            MagnetFieldSourceGrid3D_fieldCoeff,            METH_VARARGS, "Sets or Returns the scaling coefficient for inner fields in Grid3D instances"},
     {NULL}
  };

	// defenition of the memebers of the python MagnetFieldSourceGrid3D wrapper class
	// they will be vailable from python level
	static PyMemberDef MagnetFieldSourceGrid3DClassMembers [] = {
		{NULL}
	};

	//new python MagnetFieldSourceGrid3D wrapper type definition
	static PyTypeObject pyORBIT_MagnetFieldSourceGrid3D_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"MagnetFieldSourceGrid3D", /*tp_name*/
		sizeof(pyORBIT_Object), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) MagnetFieldSourceGrid3D_del , /*tp_dealloc*/
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
		"The MagnetFieldSourceGrid3D python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		MagnetFieldSourceGrid3DClassMethods, /* tp_methods */
		MagnetFieldSourceGrid3DClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) MagnetFieldSourceGrid3D_init, /* tp_init */
		0, /* tp_alloc */
		MagnetFieldSourceGrid3D_new, /* tp_new */
	};

	//--------------------------------------------------
	//Initialization function of the pyMagnetFieldSourceGrid3D class
	//It will be called from wrap_field_sources_module
	//--------------------------------------------------
  void initMagnetFieldSourceGrid3D(PyObject* module){
		if (PyType_Ready(&pyORBIT_MagnetFieldSourceGrid3D_Type) < 0) return;
		Py_INCREF(&pyORBIT_MagnetFieldSourceGrid3D_Type);
		PyModule_AddObject(module, const_cast<char*>("MagnetFieldSourceGrid3D"), (PyObject *)&pyORBIT_MagnetFieldSourceGrid3D_Type);
	}

#ifdef __cplusplus
}
#endif

//end of namespace wrap_field_source_grid3d
}
