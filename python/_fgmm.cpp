#include "fgmm++.hpp"
#include <Python.h>
#include <numpy/arrayobject.h>

/* the python object */
typedef struct {
  PyObject_HEAD
  Gmm * g;
} GMM;

static void 
Gmm_dealloc(GMM * self)
{
  self->ob_type->tp_free((PyObject*)self);
}

static PyObject *
Gmm_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    GMM *self;
    self = (GMM *)type->tp_alloc(type, 0);
    return (PyObject *) self;
}

static int
Gmm_init(GMM * self,PyObject *args,PyObject *kwds)
{
  int nstates;
  int dim;
  if(!PyArg_ParseTuple(args,"ii",&nstates,&dim))
    return NULL;
  self->g = new Gmm(nstates,dim);
  return 0;
}


static PyObject *
Gmm_init_random(GMM * self,PyObject *args, PyObject *kwds)
{
  PyObject * input_data;
  float * data;
  int data_size;
  if( ! PyArg_ParseTuple(args, "O", &input_data))
    return NULL;
  // had a check here .. 

  PyArrayObject* arr=NULL;
  bool new_arr = false;

  if( ! PyArray_Check(input_data))
    {
      return NULL;
    }

  if( !PyArray_ISCONTIGUOUS(input_data) )
    {
      arr = (PyArrayObject*) PyArray_ContiguousFromAny(input_data,NPY_FLOAT,0,0);
							
      new_arr = true;
    }
  else 
    arr = (PyArrayObject *) input_data;

  data = (float *) PyArray_DATA(arr);
  data_size = PyArray_DIM(arr,0);

  self->g->init(data,data_size);

  if(new_arr)
    Py_DECREF(arr);
  return PyInt_FromLong(1);
}


static PyObject *
Gmm_doEm(GMM * self,PyObject * args, PyObject *kwds)
{
  PyObject * input_data;
  float * data;
  int data_size;
  if( !PyArg_ParseTuple(args, "O", &input_data))
    return NULL;

  PyArrayObject* arr=NULL;
  bool new_arr = false;

  if( ! PyArray_Check(input_data))
    {
      return NULL;
    }

  if( !PyArray_ISCONTIGUOUS(input_data) )
    {
      arr = (PyArrayObject*) PyArray_ContiguousFromAny(input_data,NPY_FLOAT,0,0);
							
      new_arr = true;
    }
  else 
    arr = (PyArrayObject *) input_data;

  data = (float *) PyArray_DATA(arr);
  data_size = PyArray_DIM(arr,0);
  
  int steps = self->g->Em(data,data_size);
  if(new_arr)
    Py_DECREF(arr);
  return PyInt_FromLong(steps);
}


static PyObject *
Gmm_Dump(GMM * self)
{
  self->g->Dump();
  return PyInt_FromLong(1);
}

static PyObject *
Gmm_InitRegression(GMM * self, PyObject * args)
{ 

  int input_dim = 0;
  if( ! PyArg_ParseTuple(args, "i", &input_dim))
    return NULL;

  printf("init :: %d \n",input_dim);
  self->g->InitRegression(input_dim);
  
  return PyInt_FromLong(1);
}


static PyObject *
Gmm_DoRegression(GMM * self, PyObject * args)
{ 
  PyObject * input_data;
  float * output;

  PyObject * output_data;
  if( ! PyArg_ParseTuple(args, "O", &input_data))
    return NULL;

  if( !PyArray_Check(input_data))
    return NULL;

  PyArrayObject * arr=NULL;
  bool new_arr = false;
  if( !PyArray_ISCONTIGUOUS(input_data) )
    {
      arr = (PyArrayObject*) PyArray_ContiguousFromAny(input_data,NPY_FLOAT,0,0);
							
      new_arr = true;
    }
  else 
    arr = (PyArrayObject *) input_data;

  if( PyArray_DIM(input_data,0) != self->g->ninput) 
    {
      return NULL;
    }

  output = (float *) malloc( sizeof(float) * (self->g->dim - self->g->ninput));
  
  self->g->DoRegression((float *) PyArray_DATA(input_data), 
			output);

  npy_intp dims[1];
  dims[0] = self->g->dim - self->g->ninput;

  output_data = PyArray_SimpleNewFromData(1,dims,NPY_FLOAT,(void *) output);
  return output_data;
}

static PyObject *
Gmm_DoSamplingRegression(GMM * self, PyObject * args)
{ 
  PyObject * input_data;
  float * output;

  PyObject * output_data;
  if( ! PyArg_ParseTuple(args, "O", &input_data))
    return NULL;

  if( !PyArray_Check(input_data))
    return NULL;

  PyArrayObject * arr=NULL;
  bool new_arr = false;
  if( !PyArray_ISCONTIGUOUS(input_data) )
    {
      arr = (PyArrayObject*) PyArray_ContiguousFromAny(input_data,NPY_FLOAT,0,0);
							
      new_arr = true;
    }
  else 
    arr = (PyArrayObject *) input_data;

  if( PyArray_DIM(input_data,0) != self->g->ninput) 
    {
      return NULL;
    }

  output = (float *) malloc( sizeof(float) * (self->g->dim - self->g->ninput));
  
  self->g->DoSamplingRegression((float *) PyArray_DATA(input_data), 
			output);

  npy_intp dims[1];
  dims[0] = self->g->dim - self->g->ninput;

  output_data = PyArray_SimpleNewFromData(1,dims,NPY_FLOAT,(void *) output);
  return output_data;
}


static PyObject *
Gmm_Draw(GMM * self)
{
  float * output;
  PyObject * result;
  output = (float *) malloc(sizeof(float) * self->g->dim);
  self->g->Draw(output);

  npy_intp dims[1];
  dims[0] = self->g->dim;
  result = PyArray_SimpleNewFromData(1,dims,NPY_FLOAT,(void *) output);
  return result;
}
  
/*
static PyObject *
Gmm_getMeans(Gmm * self)
{
  
}

static PyObject *
Gmm_getCovar(Gmm * self)
{
  PyObject * covars;
  covars = PyList_New(0);
  for(int i=0;i<self->g->nState;i++)
    {
      PyList_Append(covars,matrix_to_list(self->g->sigma[i]));
    }
  return covars;
}
		 
static PyObject *
Gmm_gmr(Gmm * self, PyObject * args, PyObject * kwds)
{
  PyObject * result;
  PyObject * input_list;
  PyObject * like_py;
  PyObject * dimi_list; // input dimensions
  PyObject * dimo_list; // output dimension python list
  Vector input_vector;
  Vector input_dim;
  Vector output_dim;
  PyObject * result_sigma;
  if(! PyArg_ParseTuple(args, "OOO", &input_list,&dimi_list,&dimo_list))
    return NULL;
  input_vector = list_to_vector(input_list);
  input_dim = list_to_vector(dimi_list);
  output_dim = list_to_vector(dimo_list);
  Matrix sigma;
  double likelihood;
  Vector  result_m = self->g->doRegression(input_vector,sigma,
					   input_dim,output_dim,likelihood);
  result = vector_to_list(result_m);
  result_sigma = matrix_to_list(sigma);
  like_py = PyFloat_FromDouble(likelihood);
  return PyTuple_Pack(3,result,result_sigma,like_py);
}
  

static PyObject *
Gmm_getPriors(Gmm * self)
{
  PyObject * priors;
  priors = PyList_New(0);
  for(int i=0;i<self->g->nState;i++)
    PyList_Append(priors,PyFloat_FromDouble(self->g->priors[i]));
  return priors;
}
  
static PyObject *
Gmm_save(Gmm * self,PyObject * args, PyObject *kwds)
{
  const char * filename;
  if( !PyArg_ParseTuple(args, "s", &filename))
    return NULL;
  self->g->saveParams(filename);
  return PyInt_FromLong(1);
}

static PyObject *
Gmm_load(Gmm * self,PyObject * args, PyObject *kwds)
{
  const char * filename;
  if( !PyArg_ParseTuple(args, "s", &filename))
    return NULL;
  self->g->loadParams(filename);
  return PyInt_FromLong(1);
}

*/


static PyMethodDef Gmm_methods[] = {
  {"init",(PyCFunction) Gmm_init_random,METH_VARARGS | METH_KEYWORDS, "initialize EM process ."},
  {"Em",(PyCFunction) Gmm_doEm,METH_VARARGS | METH_KEYWORDS, "perform EM."},
  {"Dump",(PyCFunction) Gmm_Dump,METH_NOARGS,"dump gmm parameters on stdout"},
  {"InitRegression",(PyCFunction) Gmm_InitRegression,METH_VARARGS, "initialize Regression"},
  {"DoRegression",(PyCFunction) Gmm_DoRegression,METH_VARARGS, "Do Regression"},
  {"DoSamplingRegression",(PyCFunction) Gmm_DoSamplingRegression,METH_VARARGS, ""},
  {"Draw",(PyCFunction) Gmm_Draw,METH_NOARGS,"draw a random point from model"},
  /*  {"getMeans",(PyCFunction) Gmm_getMeans,METH_NOARGS, "means"},
  {"getCovar",(PyCFunction) Gmm_getCovar,METH_NOARGS, "covariances"},
  {"getPriors",(PyCFunction) Gmm_getPriors,METH_NOARGS, "priors .."},
  {"gmr",(PyCFunction) Gmm_gmr,METH_VARARGS | METH_KEYWORDS, "do regression"},
  {"save",(PyCFunction) Gmm_save,METH_VARARGS | METH_KEYWORDS, "save Model to file"},
  {"load",(PyCFunction) Gmm_load,METH_VARARGS | METH_KEYWORDS, "load model from file"},*/
  {NULL}
};

static PyTypeObject GmmType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "fgmm.GMM",                /*tp_name*/
    sizeof(GMM),               /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)Gmm_dealloc,   /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
    "Gaussian Mixture object",/* tp_doc */
    0,		               /* tp_traverse */
    0,		               /* tp_clear */
    0,		               /* tp_richcompare */
    0,		               /* tp_weaklistoffset */
    0,		               /* tp_iter */
    0,		               /* tp_iternext */
    Gmm_methods,               /* tp_methods */
    0,                         /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)Gmm_init,        /* tp_init */
    0,                         /* tp_alloc */
    Gmm_new,                   /* tp_new */
};


static PyMethodDef ModuleMethods[] = {
  {NULL, NULL, 0, NULL}        /* Sentinel */
};

PyMODINIT_FUNC
initfgmm(void)
{
    PyObject* m;

    /* this initializes the numpy lib, EXTREMELY IMPORTANT and confusing .. 
       intialize the &PyArray_Type for instance .. */

    import_array();  

    if (PyType_Ready(&GmmType) < 0)
        return;
    
   m = Py_InitModule("fgmm", ModuleMethods);

   if (m == NULL)
      return;
   
   Py_INCREF(&GmmType);
   PyModule_AddObject(m, "GMM", (PyObject *)&GmmType); 
    
}

