#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_matrix_generator.hh"
#include "wrap_teapotbase.hh"
#include "wrap_bunch.hh"
#include "wrap_utils.hh"

#include <iostream>

#include "MatrixGenerator.hh"
#include "MatrixOperations.hh"

using namespace OrbitUtils;
using namespace teapot_base;
using namespace wrap_teapotbase;

namespace wrap_teapotbase_matrix_generator
{
    void error(const char* msg){ ORBIT_MPI_Finalize(msg); }

#ifdef __cplusplus
extern "C"
{
#endif

    //---------------------------------------------------------
    //Python MatrixGenerator class definition
    //---------------------------------------------------------

    //Constructor for python class wrapping MatrixGenerator instance
    //It never will be called directly
    static PyObject* MatrixGenerator_new(PyTypeObject *type,
                                         PyObject *args, PyObject *kwds)
    {
        pyORBIT_Object* self;
        self = (pyORBIT_Object *) type->tp_alloc(type, 0);
        self->cpp_obj = NULL;
        return (PyObject *) self;
    }

    //Initializator for python  MatrixGenerator class (implementation of the __init__ )
    static int MatrixGenerator_init(pyORBIT_Object *self,
                                    PyObject *args, PyObject *kwds)
    {
        self->cpp_obj = new MatrixGenerator();
        return 0;
    }

    //-----------------------------------------------------
    //Destructor for python MatrixGenerator class (__del__ method).
    //-----------------------------------------------------
    static void MatrixGenerator_del(pyORBIT_Object* self)
    {
        delete ((MatrixGenerator*)self->cpp_obj);
        self->ob_type->tp_free((PyObject*)self);
    }

    //Sets or returns the value of the phase vector element with particular index
    static PyObject* MatrixGenerator_step_get_set(PyObject *self,
                                                  PyObject *args)
    {
        //if nVars == 1 this is get value
        //if nVars == 2 this is set value
        int nVars = PyTuple_Size(args);
        pyORBIT_Object* pyMatrixGenerator = (pyORBIT_Object*) self;
        double val = 0.;
        int i;
        MatrixGenerator* cpp_MatrixGenerator = (MatrixGenerator*) pyMatrixGenerator->cpp_obj;
        if(!PyArg_ParseTuple(	args, "i|d:set", &i, &val))
        {
            error("PyMatrixGenerator - getStep/setStep(i[,value]) - something is missing.");
        }
        if(i > 5 || i < 0)
        {
            error("PyMatrixGenerator - getStep/setStep(i[,value]) - wrong i:0-5.");
        }
        if(nVars == 1 ||  nVars == 2)
        {
            if(nVars == 1)
            {			
                val = cpp_MatrixGenerator->step(i);
            }
            else
            {
                cpp_MatrixGenerator->step(i) = val;
            }
            return Py_BuildValue("d",val);
        }
        else
        {
            error("PyMatrixGenerator. You should call getStep(i) or setStep(i,value)");
        }
        Py_INCREF(Py_None);
        return Py_None;
    }

    //  initBunch() - initialize bunch
    static PyObject* MatrixGenerator_initBunch(PyObject *self,
                                               PyObject *args)
    {
        pyORBIT_Object* pyMatrixGenerator = (pyORBIT_Object*) self;
        MatrixGenerator* cpp_MatrixGenerator = (MatrixGenerator*) pyMatrixGenerator->cpp_obj;
        PyObject *pyIn;
        if(!PyArg_ParseTuple(args,"O:initBunch",&pyIn))
        {
            error("PyMatrixGenerator - initBunch(Bunch) - input parameter is needed.");
        }
        PyObject* pyBunchType = wrap_orbit_bunch::getBunchType("Bunch");
        if(PyObject_IsInstance(pyIn,pyBunchType))
        {
            cpp_MatrixGenerator->initBunch((Bunch*) ((pyORBIT_Object*) pyIn)->cpp_obj);
            Py_INCREF(Py_None);
            return Py_None;
        }
        error("PyMatrixGenerator - initBunch(Bunch) - input parameter is wrong.");
        Py_INCREF(Py_None);
        return Py_None;
    }

    //  calcMatrix() - initialize bunch
    static PyObject* MatrixGenerator_calcMatrix(PyObject *self,
                                                PyObject *args)
    {
        pyORBIT_Object* pyMatrixGenerator = (pyORBIT_Object*) self;
        MatrixGenerator* cpp_MatrixGenerator = (MatrixGenerator*) pyMatrixGenerator->cpp_obj;
        PyObject *pyInB;
        PyObject *pyInM;
        if(!PyArg_ParseTuple(args,"OO:calcMatrix",&pyInB,&pyInM))
        {
            error("PyMatrixGenerator - calcMatrix(Bunch,Matrix) - input parameters are needed.");
        }
        PyObject* pyBunchType = wrap_orbit_bunch::getBunchType("Bunch");
        PyObject* pyORBIT_Matrix_Type = wrap_orbit_utils::getOrbitUtilsType("Matrix");
        if((!PyObject_IsInstance(pyInB,pyBunchType)) ||
           (!PyObject_IsInstance(pyInM,pyORBIT_Matrix_Type)))
        {
            error("PyMatrixGenerator - calcMatrix(Bunch,Matrix) - input parameters are wrong.");
        }
        Matrix* mtrx = (Matrix*) ((pyORBIT_Object*) pyInM)->cpp_obj;
        Bunch* bunch = (Bunch*) ((pyORBIT_Object*) pyInB)->cpp_obj;
        if(mtrx->rows() < 6 || mtrx->columns() < 6)
        {
            error("PyMatrixGenerator - calcMatrix(Bunch,Matrix) - Matrix should be more than 6x6.");
        }
        cpp_MatrixGenerator->calculateMatrix(bunch, mtrx);
        Py_INCREF(Py_None);
        return Py_None;
    }

    //  initBunchForChromaticityCoeff() - initialize bunch for Chromaticity calculation
    static PyObject* MatrixGenerator_initBunchChromCoeff(PyObject *self,
                                                         PyObject *args)
    {
        pyORBIT_Object* pyMatrixGenerator = (pyORBIT_Object*) self;
        MatrixGenerator* cpp_MatrixGenerator = (MatrixGenerator*) pyMatrixGenerator->cpp_obj;
        PyObject *pyIn;
        if(!PyArg_ParseTuple(args,"O:initBunchChromCoeff",&pyIn))
        {
            error("PyMatrixGenerator - initBunchChromCoeff(Bunch) - input parameter is needed.");
        }
        PyObject* pyBunchType = wrap_orbit_bunch::getBunchType("Bunch");
        if(PyObject_IsInstance(pyIn,pyBunchType))
        {
            cpp_MatrixGenerator->initBunchForChromaticityCoeff((Bunch*) ((pyORBIT_Object*) pyIn)->cpp_obj);
            Py_INCREF(Py_None);
            return Py_None;
        }
        error("PyMatrixGenerator - initBunchChromCoeff(Bunch) - input parameter is wrong.");
        Py_INCREF(Py_None);
        return Py_None;
    }

    //calculateChromaticityCoeff() - calculates the coeff.  coeff_x_dE, coeff_xp_dE, coeff_y_dE, coeff_yp_dE
    static PyObject* MatrixGenerator_calcChromCoeff(PyObject *self,
                                                    PyObject *args)
    {
        pyORBIT_Object* pyMatrixGenerator = (pyORBIT_Object*) self;
        MatrixGenerator* cpp_MatrixGenerator = (MatrixGenerator*) pyMatrixGenerator->cpp_obj;
        PyObject *pyIn;
        if(!PyArg_ParseTuple(args,"O:calcChromCoeff",&pyIn))
        {
            error("PyMatrixGenerator - calcChromCoeff(Bunch) - input parameter is needed.");
        }
        PyObject* pyBunchType = wrap_orbit_bunch::getBunchType("Bunch");
        if(PyObject_IsInstance(pyIn,pyBunchType))
        {
            double coeff_x_dE, coeff_xp_dE, coeff_y_dE, coeff_yp_dE;
            cpp_MatrixGenerator->calculateChromaticityCoeff((Bunch*) ((pyORBIT_Object*) pyIn)->cpp_obj,
                                                            coeff_x_dE, coeff_xp_dE, coeff_y_dE, coeff_yp_dE);
            return Py_BuildValue("(dddd)",coeff_x_dE, coeff_xp_dE,
                                 coeff_y_dE, coeff_yp_dE);
        }
        error("PyMatrixGenerator - calcChromCoeff(Bunch) - input parameter is wrong.");
        Py_INCREF(Py_None);
        return Py_None;
    }

    // Definition of methods of the python MatrixGenerator wrapper class
    // They will be vailable from python level
    static PyMethodDef MatrixGeneratorClassMethods[] =
    {
        { "getStep",             MatrixGenerator_step_get_set        ,METH_VARARGS, "Returns step(i)"},
        { "setStep",             MatrixGenerator_step_get_set        ,METH_VARARGS, "Sets the new value to step(i)"},
        { "initBunch",           MatrixGenerator_initBunch           ,METH_VARARGS, "Initializes the bunch"},
        { "calculateMatrix",     MatrixGenerator_calcMatrix          ,METH_VARARGS, "Calculates the transport matrix"},
        { "initBunchChromCoeff", MatrixGenerator_initBunchChromCoeff ,METH_VARARGS, "Initializes the bunch for chromaticity calc."},
        { "calcChromCoeff",      MatrixGenerator_calcChromCoeff      ,METH_VARARGS, "Calculates the cromaticity coeff."},
        {NULL}
    };

    // Definition of the memebers of the python MatrixGenerator wrapper class
    // They will be vailable from python level
    static PyMemberDef MatrixGeneratorClassMembers [] =
    {
        {NULL}
    };

    //New python MatrixGenerator wrapper type definition
    static PyTypeObject pyORBIT_MatrixGenerator_Type =
    {
        PyObject_HEAD_INIT(NULL)
        0, /*ob_size*/
        "MatrixGenerator", /*tp_name*/
        sizeof(pyORBIT_Object), /*tp_basicsize*/
        0, /*tp_itemsize*/
        (destructor) MatrixGenerator_del , /*tp_dealloc*/
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
        "The MatrixGenerator python wrapper", /* tp_doc */
        0, /* tp_traverse */
        0, /* tp_clear */
        0, /* tp_richcompare */
        0, /* tp_weaklistoffset */
        0, /* tp_iter */
        0, /* tp_iternext */
        MatrixGeneratorClassMethods, /* tp_methods */
        MatrixGeneratorClassMembers, /* tp_members */
        0, /* tp_getset */
        0, /* tp_base */
        0, /* tp_dict */
        0, /* tp_descr_get */
        0, /* tp_descr_set */
        0, /* tp_dictoffset */
        (initproc) MatrixGenerator_init, /* tp_init */
        0, /* tp_alloc */
        MatrixGenerator_new, /* tp_new */
    };

    //--------------------------------------------------
    //Initialization function of the pyMatrixGenerator class
    //It will be called from Bunch wrapper initialization
    //--------------------------------------------------
    void initMatrixGenerator(PyObject* module)
    {
        if (PyType_Ready(&pyORBIT_MatrixGenerator_Type) < 0) return;
        Py_INCREF(&pyORBIT_MatrixGenerator_Type);
        PyModule_AddObject(module, "MatrixGenerator",
                           (PyObject *)&pyORBIT_MatrixGenerator_Type);
    }

#ifdef __cplusplus
}
#endif

//end of namespace wrap_teapotbase_matrix_generator
}

