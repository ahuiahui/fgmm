#include "fgmm++.hpp"
#include <Python.h>
#include <numpy/arrayobject.h>

/* helper func to convert a PyObject (numpy array .. ) 
   to a _fgmm_real * array .. 
   
   /!\ the result must be DECREF after use .. 
*/

PyArrayObject * get_real_array(PyObject * input_data)
{

  PyArrayObject* arr=NULL;
  PyArrayObject* res=NULL;
  if(!PyArray_Check(input_data))
    {
      return NULL;
    }
  arr = (PyArrayObject*) PyArray_FromAny(input_data,NULL,0,2,NPY_CARRAY_RO,NULL);
#ifdef FGMM_USE_DOUBLE
  res = (PyArrayObject *) PyArray_Cast(arr, NPY_DOUBLE);
#else 
  res = (PyArrayObject *) PyArray_Cast(arr, NPY_FLOAT);
#endif // USE DOUBLE
  Py_DECREF(arr);
  return res;
};
 
PyObject * get_vector(int dim, _fgmm_real * dat)
{
  
  npy_intp dims[1];
  dims[0] = dim;
  PyObject * output_data;
#ifdef FGMM_USE_DOUBLE
  output_data = PyArray_SimpleNewFromData(1,dims,NPY_DOUBLE,(void *) dat);
#else
  output_data = PyArray_SimpleNewFromData(1,dims,NPY_FLOAT,(void *) dat);
#endif // USE DOUBLE
};

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
Gmm_init_kmeans(GMM * self, PyObject *args, PyObject *kwds)
{
  PyObject * input_data;
  _fgmm_real * data;
  int data_size;
  if( ! PyArg_ParseTuple(args, "O", &input_data))
    return NULL;
  // had a check here .. 

  PyArrayObject* arr= get_real_array(input_data);
  if(arr == NULL)
    return NULL;
  data = (_fgmm_real *) PyArray_DATA(arr);
  data_size = PyArray_DIM(arr,0);
  printf("K-means data :: %d\n", data_size);
  printf("K-Means data :: %f %f\n", data[2 * (data_size-1)], data[2 * (data_size-1) +1 ]);
  self->g->initKmeans(data,data_size);
  Py_DECREF(arr);
  return PyInt_FromLong(1);
}
  

static PyObject *
Gmm_init_random(GMM * self,PyObject *args, PyObject *kwds)
{
  PyObject * input_data;
  _fgmm_real * data;
  int data_size;
  if( ! PyArg_ParseTuple(args, "O", &input_data))
    return NULL;
  // had a check here .. 

  PyArrayObject* arr= get_real_array(input_data);
  if(arr == NULL)
    return NULL;
  data = (_fgmm_real *) PyArray_DATA(arr);
  data_size = PyArray_DIM(arr,0);

  self->g->init(data,data_size);

  Py_DECREF(arr);
  return PyInt_FromLong(1);
}


static PyObject *
Gmm_doEm(GMM * self,PyObject * args, PyObject *kwds)
{
  PyObject * input_data;
  _fgmm_real * data;
  int data_size;
  _fgmm_real epsilon = 1e-4;
  int covar_t = COVARIANCE_FULL;

  if( !PyArg_ParseTuple(args, "O|fb", &input_data,&epsilon,&covar_t))
    return NULL;

  printf("delta epsilon %e\n",epsilon);

  PyArrayObject* arr= get_real_array(input_data);
  if(arr == NULL)
    return NULL;
  data = (_fgmm_real *) PyArray_DATA(arr);

  printf("EM data :: %f %f %f\n", data[0], data[1], data[2]);

  data_size = PyArray_DIM(arr,0);
  
  int steps = self->g->Em(data,data_size,epsilon,COVARIANCE_TYPE(covar_t));
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
Gmm_Pdf(GMM * self, PyObject * args)
{
  PyObject * input_data;
  _fgmm_real * output;

  PyObject * output_data;
  if( ! PyArg_ParseTuple(args, "O", &input_data))
    return NULL;


  PyArrayObject * arr = get_real_array(input_data);
  if(arr == NULL)
    return NULL;

  output = (_fgmm_real *) malloc( sizeof(_fgmm_real) * (self->g->nstates));
  self->g->Pdf((_fgmm_real *) PyArray_DATA(arr), 
	       output);
  
  Py_DECREF(arr);
  output_data = get_vector( self->g->nstates, output);
  return output_data;
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
  _fgmm_real * output;
  PyObject * output_data;
  if( ! PyArg_ParseTuple(args, "O", &input_data))
    return NULL;

  PyArrayObject * arr= get_real_array(input_data);
  
  if( arr == NULL || PyArray_DIM(input_data,0) != self->g->ninput) 
    {
      return NULL;
    }

  output = (_fgmm_real *) malloc( sizeof(_fgmm_real) * (self->g->dim - self->g->ninput));
  
  self->g->DoRegression((_fgmm_real *) PyArray_DATA(input_data), 
			output);
  Py_DECREF(arr);
  output_data =  get_vector( self->g->dim - self->g->ninput, output);
  return output_data;
}

static PyObject *
Gmm_DoSamplingRegression(GMM * self, PyObject * args)
{ 
  PyObject * input_data;
  _fgmm_real * output;

  PyObject * output_data;
  if( ! PyArg_ParseTuple(args, "O", &input_data))
    return NULL;


  PyArrayObject * arr= get_real_array(input_data);

  if(arr==NULL)
    {
      return NULL;
    }

  output = (_fgmm_real *) malloc( sizeof(_fgmm_real) * (self->g->dim - self->g->ninput));
  
  self->g->DoSamplingRegression((_fgmm_real *) PyArray_DATA(arr), 
			output);

  Py_DECREF(arr);
  output_data = get_vector(self->g->dim - self->g->ninput, output);
  return output_data;
}


static PyObject *
Gmm_Draw(GMM * self)
{
  _fgmm_real * output;
  PyObject * result;
  output = (_fgmm_real *) malloc(sizeof(_fgmm_real) * self->g->dim);
  self->g->Draw(output);
  result = get_vector(self->g->dim, output);
  return result;
}
  

static PyObject *
Gmm_getMean(GMM * self,PyObject * args)
{
  _fgmm_real * output;
  PyObject * result;

  int state = 0;
  if( ! PyArg_ParseTuple(args, "i", &state))
    return NULL;

  output = (_fgmm_real *) malloc(sizeof(_fgmm_real) * self->g->dim);

  if(state > self->g->nstates) 
    {
      printf("Wrong index for GetMean, got %d out of %d\n",state,self->g->nstates);
      return NULL;
    }
  self->g->GetMean(state,output);
  printf("Mean :: %f %f \n", output[0], output[1]);
  result = get_vector(self->g->dim, output);
  return result;
}


static PyObject *
Gmm_getCovariance(GMM * self,PyObject * args)
{
  _fgmm_real * output;
  PyObject * result;

  int state = 0;
  if( ! PyArg_ParseTuple(args, "i", &state))
    return NULL;

  output = (_fgmm_real *) malloc(sizeof(_fgmm_real) * self->g->dim * self->g->dim );
  self->g->GetCovariance(state,output);

  npy_intp dims[2];
  dims[0] = self->g->dim;
  dims[1] = self->g->dim;
#ifdef FGMM_USE_DOUBLE
  result = PyArray_SimpleNewFromData(2,dims,NPY_DOUBLE,(void *) output);
#else
  result = PyArray_SimpleNewFromData(2,dims,NPY_FLOAT,(void *) output);
#endif // USE_DOUBLE
  return result;

}


static PyObject *
Gmm_getPrior(GMM * self, PyObject * args)
{
  int state=0;
  if( ! PyArg_ParseTuple(args, "i", &state))
    return NULL;
  float p = self->g->GetPrior(state);
  return PyFloat_FromDouble((double) p);
}

static PyObject *
Gmm_getstate(GMM * self, PyObject * args)
{
  int res=0;
  PyObject * obs=NULL;

  if( ! PyArg_ParseTuple(args,"O", &obs))
    return NULL;

  
  PyArrayObject * arr= get_real_array(obs);
  if(arr == NULL)
    return NULL;

  if( PyArray_DIM(arr,0) != self->g->dim) 
    {
      return NULL;
    }

  res = self->g->GetLikelyState((_fgmm_real *)PyArray_DATA(arr));
  Py_DECREF(arr);

  return PyInt_FromLong(res);
}

  

static PyMethodDef Gmm_methods[] = {
  {"init",(PyCFunction) Gmm_init_random,METH_VARARGS | METH_KEYWORDS, "initialize EM process ."},
  {"kmeans",(PyCFunction) Gmm_init_kmeans,METH_VARARGS | METH_KEYWORDS, "initialize EM with k-means ."},
  {"Em",(PyCFunction) Gmm_doEm,METH_VARARGS | METH_KEYWORDS, "perform EM."},
  {"Dump",(PyCFunction) Gmm_Dump,METH_NOARGS,"dump gmm parameters on stdout"},
  {"Pdf",(PyCFunction) Gmm_Pdf,METH_VARARGS,"pdf"},
  {"InitRegression",(PyCFunction) Gmm_InitRegression,METH_VARARGS, "initialize Regression"},
  {"DoRegression",(PyCFunction) Gmm_DoRegression,METH_VARARGS, "Do Regression"},
  {"DoSamplingRegression",(PyCFunction) Gmm_DoSamplingRegression,METH_VARARGS, ""},
  {"Draw",(PyCFunction) Gmm_Draw,METH_NOARGS,"draw a random point from model"},
  {"GetCovariance",(PyCFunction) Gmm_getCovariance,METH_VARARGS, "covariances"},
  {"GetMean",(PyCFunction) Gmm_getMean,METH_VARARGS, "covariances"},
  {"GetPrior",(PyCFunction) Gmm_getPrior,METH_NOARGS, "priors .."},
  {"GetLikelyState",(PyCFunction) Gmm_getstate,METH_VARARGS, "get the most likely state index"},
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
   PyModule_AddIntConstant(m,"COVARIANCE_FULL", COVARIANCE_FULL);
   PyModule_AddIntConstant(m,"COVARIANCE_DIAG", COVARIANCE_DIAG);
   PyModule_AddIntConstant(m,"COVARIANCE_SPHERE", COVARIANCE_SPHERE);
   
}

