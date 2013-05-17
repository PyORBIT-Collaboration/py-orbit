#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_aperture.hh"
#include "wrap_bunch.hh"

#include <iostream>

#include "wrap_FieldTracker.hh"

namespace wrap_FieldTracker{
    
    void error(const char* msg){ ORBIT_MPI_Finalize(msg); }
    
#ifdef __cplusplus
    extern "C" {
#endif
        
        /**
         Constructor for python class wrapping c++ Aperture instance.
         It never will be called directly.
         */
        static PyObject* Aperture_new(PyTypeObject *type, PyObject *args, PyObject *kwds){
            pyORBIT_Object* self;
            self = (pyORBIT_Object *) type->tp_alloc(type, 0);
            self->cpp_obj = NULL;
            return (PyObject *) self;
        }
        
        /** This is implementation of the __init__ method */
        static int Aperture_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){
            
            int shape = 0.;
            double a = 0.;
            double b = 0.;
            double c = 0.;
            double d = 0.;
            
            //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
            if(!PyArg_ParseTuple(	args,"idddd:arguments",&shape,&a,&b,&c,&d)){
                error("PyBunch - addParticle - cannot parse arguments! It should be (shape,a,b,c,d)");
            }
            self->cpp_obj =  new Aperture(shape,a,b,c,d);
            ((Aperture*) self->cpp_obj)->setPyWrapper((PyObject*) self);
            return 0;
        }
        
        /** Performs the collimation tracking of the bunch */
        static PyObject* Aperture_checkBunch(PyObject *self, PyObject *args){
            Aperture* cpp_Aperture = (Aperture*)((pyORBIT_Object*) self)->cpp_obj;
            PyObject* pyBunch;
            PyObject* pyLostBunch;
            if(!PyArg_ParseTuple(args,"OO:checkBunch",&pyBunch, &pyLostBunch)){
                ORBIT_MPI_Finalize("Aperture - checkBunch(Bunch* bunch, Bunch* bunch) - parameter are needed.");
            }
            PyObject* pyORBIT_Bunch_Type = wrap_orbit_bunch::getBunchType("Bunch");
            if(!PyObject_IsInstance(pyBunch,pyORBIT_Bunch_Type) || !PyObject_IsInstance(pyLostBunch,pyORBIT_Bunch_Type)){
                ORBIT_MPI_Finalize("Aperture - checkBunch(Bunch* bunch, Bunch* bunch) - method needs a Bunch.");
            }
            
            Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object*)pyBunch)->cpp_obj;
            Bunch* cpp_lostbunch = (Bunch*) ((pyORBIT_Object*)pyLostBunch)->cpp_obj;
            cpp_Aperture->checkBunch(cpp_bunch, cpp_lostbunch);
            Py_INCREF(Py_None);
            return Py_None;
        }
		
        
        //-----------------------------------------------------
        //destructor for python Aperture class (__del__ method).
        //-----------------------------------------------------
        static void Aperture_del(pyORBIT_Object* self){
            //std::cerr<<"The Aperture __del__ has been called!"<<std::endl;
            delete ((Aperture*)self->cpp_obj);
            self->ob_type->tp_free((PyObject*)self);
        }
        
        // definition of the methods of the python Aperture wrapper class
        // they will be vailable from python level
        static PyMethodDef ApertureClassMethods[] = {
            { "checkBunch",				 Aperture_checkBunch,    	METH_VARARGS,"Performs the aperture check of the bunch."},
            {NULL}
        };
        
        static PyMemberDef ApertureClassMembers [] = {
            {NULL}
        };
        
        
        //new python Aperture wrapper type definition
        static PyTypeObject pyORBIT_Aperture_Type = {
            PyObject_HEAD_INIT(NULL)
            0, /*ob_size*/
            "Aperture", /*tp_name*/
            sizeof(pyORBIT_Object), /*tp_basicsize*/
            0, /*tp_itemsize*/
            (destructor) Aperture_del , /*tp_dealloc*/
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
            "The Aperture python wrapper", /* tp_doc */
            0, /* tp_traverse */
            0, /* tp_clear */
            0, /* tp_richcompare */
            0, /* tp_weaklistoffset */
            0, /* tp_iter */
            0, /* tp_iternext */
            ApertureClassMethods, /* tp_methods */
            ApertureClassMembers, /* tp_members */
            0, /* tp_getset */
            0, /* tp_base */
            0, /* tp_dict */
            0, /* tp_descr_get */
            0, /* tp_descr_set */
            0, /* tp_dictoffset */
            (initproc) Aperture_init, /* tp_init */
            0, /* tp_alloc */
            Aperture_new, /* tp_new */
        };	
        
        
        
        //--------------------------------------------------
        //Initialization Aperture of the pyBunchAperture class
        //--------------------------------------------------
        
        void initaperture(){
            //check that the Aperture wrapper is ready
            if (PyType_Ready(&pyORBIT_Aperture_Type) < 0) return;
            Py_INCREF(&pyORBIT_Aperture_Type);
            //create new module
            PyObject* module = Py_InitModule("aperture",ApertureClassMethods);
            PyModule_AddObject(module, "Aperture", (PyObject *)&pyORBIT_Aperture_Type);			
        }
        
#ifdef __cplusplus
    }
#endif
    
    
}
