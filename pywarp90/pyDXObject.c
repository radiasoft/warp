#include "Python.h"
#include <Numeric/arrayobject.h>
#include "dx/dx.h"
#include <unistd.h>

#define PyDXObject_Check(op) ((op)->ob_type == &PyDXObject_Type)

static PyObject *ErrorObject;

/* ----------------------------------------------------- */

/* Declarations for objects of type PyDXObject */

typedef struct {
  PyObject_HEAD
  Object dxobj;
} PyDXObject;

staticforward PyTypeObject PyDXObject_Type;


/* ---------- */
static PyDXObject *
newPyDXObject(Object dxobj)
{
  PyDXObject *self;
  self = PyObject_NEW(PyDXObject, &PyDXObject_Type);
  if (self == NULL) return NULL;
  if (dxobj != NULL) self->dxobj = dxobj;
  return self;
}


static void
PyDXObject_dealloc(PyObject *self)
{
  /* DXDelete(self->dxobj); */
  PyMem_DEL(self);
}

static int
PyDXObject_print(PyObject *self,FILE *fp,int flags)
{
  fprintf(fp,"<PyDXObject %d>",self);
  return 0;
}

static PyObject *
PyDXObject_repr(PyObject *self)
{
  PyObject *s;
  s = Py_BuildValue("s","<PyDXObject>");
  return s;
}

static PyObject *
PyDXObject_str(PyObject *self)
{
  PyObject *s;
  s = Py_BuildValue("s","<PyDXObject>");
  return s;
}

static char PyDXObject_isnull__doc__[] = "Returns true when DX object is null";

static PyObject *
PyDXObject_isnull(PyDXObject *self)
{
  if (self->dxobj == NULL) return Py_BuildValue("i",1);
  else                     return Py_BuildValue("i",0);
}

/* ---------------------------------------------------------------- */

static struct PyMethodDef PyDXObject_methods[] = {
 {"isnull",(PyCFunction)PyDXObject_isnull,METH_VARARGS,PyDXObject_isnull__doc__},
 {NULL,NULL} /* sentinel */
};

static PyObject *
PxDXObject_getattr(PyDXObject *self,char *name)
{
  return Py_FindMethod(PyDXObject_methods, (PyObject *)self, name);
}

/* ---------------------------------------------------------------- */
static char PyDXObject_Type__doc__[] = 
"OpenDX Object"
;

static PyTypeObject PyDXObject_Type = {
  PyObject_HEAD_INIT(&PyType_Type)
  0,                                /*ob_size*/
  "PyDXObject",                    /*tp_name*/
  sizeof(PyDXObject),              /*tp_basicsize*/
  0,                                /*tp_itemsize*/
  /* methods */
  (destructor)PyDXObject_dealloc,  /*tp_dealloc*/
  (printfunc)PyDXObject_print,     /*tp_print*/
  (getattrfunc)PxDXObject_getattr,  /*tp_getattr*/
  (setattrfunc)0,                   /*tp_setattr*/
  (cmpfunc)0,                       /*tp_compare*/
  (reprfunc)PyDXObject_repr,       /*tp_repr*/
  0,                                /*tp_as_number*/
  0,                                /*tp_as_sequence*/
  0,                                /*tp_as_mapping*/
  (hashfunc)0,                      /*tp_hash*/
  (ternaryfunc)0,                   /*tp_call*/
  (reprfunc)PyDXObject_str,        /*tp_str*/
  
  /* Space for future expansion */
  0L,0L,0L,0L,
  PyDXObject_Type__doc__ /* Documentation string */
};

/* End of code for PyDXObject objects */
/* -------------------------------------------------------- */

static char pyDXObject_DXNewArray__doc__[] = "DXNewArray";

static PyObject *
pyDXObject_DXNewArray(PyObject *self,PyObject *args)
{
  int t,c,rank;
  int shape[]={0,0,0,0,0,0,0};
  Array a=NULL;
  PyObject *result;
  if (!PyArg_ParseTuple(args, "iii|iiiiiii",&t,&c,&rank,
                        &shape[0],&shape[1],&shape[2],&shape[3],
                        &shape[4],&shape[5],&shape[6])) return NULL;
  a = DXNewArrayV(t,c,rank,shape);
  if (a != NULL) {
    result = (PyObject *)newPyDXObject((Object)a);
    return result;
    }
  else {
    PyErr_SetString(ErrorObject,"Error during DXNewArray");
    return NULL;
    }
}

static char pyDXObject_DXAddArrayData__doc__[] = "DXAddArrayData";

static PyObject *
pyDXObject_DXAddArrayData(PyObject *self,PyObject *args)
{
  PyObject *a,*data;
  PyArrayObject *pyarr;
  Array e;
  int start,n,type=PyArray_DOUBLE;
  if (!PyArg_ParseTuple(args, "O!iiO",&PyDXObject_Type,&a,&start,&n,&data))
    return NULL;
  pyarr = (PyArrayObject *)PyArray_ContiguousFromObject(data,
                                ((PyArrayObject *)data)->descr->type_num,0,0);
  e = DXAddArrayData((Array)((PyDXObject *)a)->dxobj,start,n,
                     (Pointer)(pyarr->data));
  if (e != NULL) {
    Py_INCREF(Py_None);
    return Py_None;
    }
  else {
    PyErr_SetString(ErrorObject,"Error during DXAddArrayData");
    return NULL;
    }
}

static char pyDXObject_DXNewField__doc__[] = "DXNewField";

static PyObject *
pyDXObject_DXNewField(PyObject *self,PyObject *args)
{
  Field f=NULL;
  PyObject *result;
  if (!PyArg_ParseTuple(args, "")) return NULL;
  f = DXNewField();
  result = (PyObject *)newPyDXObject((Object)f);
  if (result != NULL) {
    return result;
    }
  else {
    PyErr_SetString(ErrorObject,"Error during DXNewField");
    return NULL;
    }
}

static char pyDXObject_DXSetStringAttribute__doc__[] = "DXSetStringAttribute";

static PyObject *
pyDXObject_DXSetStringAttribute(PyObject *self,PyObject *args)
{
  PyObject *a;
  Object e;
  char *name,*x;
  if (!PyArg_ParseTuple(args, "O!ss",&PyDXObject_Type,&a,&name,&x)) return NULL;
  e = DXSetStringAttribute(((PyDXObject *)a)->dxobj,name,x);
  if (e != NULL) {
    Py_INCREF(Py_None);
    return Py_None;
    }
  else {
    PyErr_SetString(ErrorObject,"Error during DXSetStringAttribute");
    return NULL;
    }
}

static char pyDXObject_DXSetComponentValue__doc__[] = "DXSetComponentValue";

static PyObject *
pyDXObject_DXSetComponentValue(PyObject *self,PyObject *args)
{
  PyObject *f,*value;
  char *name;
  Field e;
  if (!PyArg_ParseTuple(args, "O!sO!",&PyDXObject_Type,&f,&name,
                                      &PyDXObject_Type,&value)) return NULL;
  e = DXSetComponentValue((Field)((PyDXObject *)f)->dxobj,name,
                          ((PyDXObject *)value)->dxobj);
  if (e != NULL) {
    Py_INCREF(Py_None);
    return Py_None;
    }
  else {
    PyErr_SetString(ErrorObject,"Error during DXSetComponentValue");
    return NULL;
    }
}

static char pyDXObject_DXEndField__doc__[] = "DXEndField";

static PyObject *
pyDXObject_DXEndField(PyObject *self,PyObject *args)
{
  PyObject *f;
  Field e;
  if (!PyArg_ParseTuple(args, "O!",&PyDXObject_Type,&f)) return NULL;
  e = DXEndField((Field)((PyDXObject *)f)->dxobj);
  if (e != NULL) {
    Py_INCREF(Py_None);
    return Py_None;
    }
  else {
    PyErr_SetString(ErrorObject,"Error during DXEndField");
    return NULL;
    }
}

static char pyDXObject_DXMakeStringList__doc__[] = "DXMakeStringList";

static PyObject *
pyDXObject_DXMakeStringList(PyObject *self,PyObject *args)
{
  PyObject *l;
  int n,i;
  char **s;
  Array sl;
  if (!PyArg_ParseTuple(args, "O!",&PyList_Type,&l)) return NULL;
  n = PyList_Size(l);
  s = (char **)malloc(sizeof(char*)*n);
  for (i=0;i<n;i++) *(s+i) = PyString_AsString(PyList_GetItem(l,i));
  sl = DXMakeStringListV(n,s);
  if (sl != NULL) {
    return (PyObject *)newPyDXObject((Object)sl);
    }
  else {
    PyErr_SetString(ErrorObject,"Error during DXMakeStringList");
    return NULL;
    }
}

/*
  a = DXMakeGridConnectionsV(dims, counts);
  a = DXMakeGridPositionsV(dims, counts, origins, deltas);
*/

static char pyDXObject_DXObject_fromarray__doc__[] =
"Returns PyDXObject"
;

static PyObject *
pyDXObject_DXObject_fromarray(PyObject *self,PyObject * args)
{
        PyObject *pyarr,*result;
        PyObject *originslist=NULL,*deltaslist=NULL;
        PyArrayObject *pyarrcont;
        double *data;
        int dims,*counts;

  Array a=NULL;
  Field f=NULL;
  int numelements, i, j;
  float ff,*deltas, *origins;

  if (!PyArg_ParseTuple(args, "O|OO",&pyarr,&originslist,&deltaslist))
     return NULL;
  pyarrcont = (PyArrayObject *)PyArray_ContiguousFromObject(
                                  (PyObject *)pyarr,PyArray_DOUBLE,0,0);

  data = (double *)(pyarrcont->data);
  dims = pyarrcont->nd;
  counts = pyarrcont->dimensions;

  origins = (float *) malloc(dims*sizeof(float));
  deltas  = (float *) malloc(dims*dims*sizeof(float));

  if (originslist != NULL) {
    if (PyList_Size(originslist) != dims) {
      PyErr_SetString(PyExc_AssertionError,"Len of data origins must be the same as the number of dims in the data");
      goto err;
      }
    for (i=0;i<dims;i++)
      origins[i] = (float) PyFloat_AsDouble(PyList_GetItem(originslist,i));
    }
  else {
    for (i=0;i<dims;i++) origins[i] = 0.;
    }

  if (deltaslist != NULL) {
    if (PyList_Size(deltaslist) != dims) {
      PyErr_SetString(PyExc_AssertionError,"Len of data deltas must be the same as the number of dims in the data");
      goto err;
      }
    for (i=0; i<dims; i++) {
      ff = (float) PyFloat_AsDouble(PyList_GetItem(deltaslist,i));
      for (j=0; j<dims; j++) {
        if (i==j) deltas[i*dims + j] = ff;
        else      deltas[i*dims + j] = 0.0;
        }
      }
    }
  else {
    for (i=0; i<dims; i++) {
      for (j=0; j<dims; j++) {
        if (i==j) deltas[i*dims + j] = 1.0;
        else      deltas[i*dims + j] = 0.0;
        }
      }
    }


  /* Convert data to DX object */
  /* make a new data array (scalar) */
  a = DXNewArray(TYPE_DOUBLE, CATEGORY_REAL, 0);

  /* figure out how many elements there are in the data array */
  for (i=0, numelements=1; i<dims; numelements *= counts[i], i++); 
  
  /* allocate space in the data array */
  DXAddArrayData(a, 0, numelements, data);

  /* get a pointer to the data array */
  data = (double *)DXGetArrayData(a);

  /* create a new field */
  f = DXNewField();

  /* set the dependency of the data to be on positions */
  DXSetStringAttribute((Object)a, "dep", "positions");

  /* set the data array as the data component of f */
  DXSetComponentValue(f, "data", (Object)a);
  a=NULL;

  /* create the connections array */
  /* DXMakeGridConnections will set the element type */
  a = DXMakeGridConnectionsV(dims, counts);
  DXSetComponentValue(f, "connections", (Object)a);
  a=NULL;

  a = DXMakeGridPositionsV(dims, counts, origins, deltas);
  DXSetComponentValue(f, "positions", (Object)a);
  a=NULL; 

  /* EndField sets default attributes (such as setting the connections */
  /* attribute to be "ref" positions), and creates the bounding box.  */
  DXEndField(f);

  result = (PyObject *)newPyDXObject((Object)f);
  return result;

err:
  free(origins);
  free(deltas);
  return NULL;
}

static char pyDXObject_DXCallModule__doc__[] =
""
;

static PyObject *
pyDXObject_DXCallModule(PyObject *self,PyObject * args)
{
  char *modname;
  PyObject *inputdict;
  PyObject *outputlist;
  PyObject *key,*value;
  ModuleInput  minput[30];
  ModuleOutput moutput[30];
  int ninput=0,noutput=0,i;
  int pos=0;
  char *varname;
  float fvalue;
  int ivalue;
  char *svalue;
  PyDXObject *ovalue;
  PyObject *result;
  Error e;

  if (!PyArg_ParseTuple(args, "sO!O!",&modname,
                                      &PyDict_Type,&inputdict,
                                      &PyList_Type,&outputlist)) return NULL;

  /* First, iterate over the input dictionary, setting values */
  while (PyDict_Next(inputdict,&pos,&key,&value)){
    varname = PyString_AsString(key);
    if (PyFloat_Check(value)) {
      fvalue = (float)PyFloat_AsDouble(value);
      DXModSetFloatInput(&minput[ninput],varname,fvalue);
      }
    else if (PyInt_Check(value)) {
      ivalue = (int)PyInt_AsLong(value);
      DXModSetIntegerInput(&minput[ninput],varname,ivalue);
      }
    else if (PyString_Check(value)) {
      svalue = PyString_AsString(value);
      DXModSetStringInput(&minput[ninput],varname,svalue);
      }
    else if (PyDXObject_Check(value)) {
      DXModSetObjectInput(&minput[ninput],varname,((PyDXObject *)value)->dxobj);
      }
    else {
      PyErr_Format(PyExc_TypeError,
        "The value for \"%s\" does not have a valid type",varname);
      }
    ninput++;
    }

  /* Create tuple of PyDXObjects to be returned */
  noutput = PyList_Size(outputlist);
  result = PyTuple_New(noutput);
  for (i=0;i<noutput;i++)
    PyTuple_SetItem(result,i,(PyObject *)newPyDXObject((Object)NULL));

  /* Loop over output list */
  for (i=0;i<noutput;i++) {
    key = PyList_GetItem(outputlist,i);
    varname = PyString_AsString(key);
    ovalue = (PyDXObject *)PyTuple_GetItem(result,i);
    DXModSetObjectOutput(&moutput[i],varname,&(ovalue->dxobj));
    }

  /* Now, make call to the module */
  e = DXCallModule(modname,ninput,minput,noutput,moutput);

  if (e == OK) {
    return result;
    }
  else {
    PyErr_SetString(ErrorObject,"Error during CallModule");
    return NULL;
    }
}

static char pyDXObject_DXReference__doc__[] =
""
;

static PyObject *
pyDXObject_DXReference(PyObject *self,PyObject * args)
{
  PyObject *pyvalue;
  Object e;
  if (!PyArg_ParseTuple(args, "O!",&PyDXObject_Type,&pyvalue)) return NULL;
  e = DXReference(((PyDXObject *)pyvalue)->dxobj);
  if (e != NULL) {
    Py_INCREF(Py_None);
    return Py_None;
    }
  else {
    PyErr_SetString(ErrorObject,"Error during DXReference");
    return NULL;
    }
}

static char pyDXObject_DXCheckRIH__doc__[] =
""
;

static PyObject *
pyDXObject_DXCheckRIH(PyObject *self,PyObject * args)
{
  int block;
  int e;

  if (!PyArg_ParseTuple(args, "i",&block)) return NULL;
  e = DXCheckRIH(block);
  return Py_BuildValue("(i)",e);
}

static char pyDXObject_DXDelete__doc__[] =
""
;

static PyObject *
pyDXObject_DXDelete(PyObject *self,PyObject * args)
{
  PyObject *pyvalue;
  Error e;
  if (!PyArg_ParseTuple(args, "O!",&PyDXObject_Type,&pyvalue)) return NULL;
  e = DXDelete(((PyDXObject *)pyvalue)->dxobj);
  if (e == OK) {
    Py_INCREF(Py_None);
    return Py_None;
    }
  else {
    PyErr_SetString(ErrorObject,"Error during DXDelete");
    return NULL;
    }
}

/* Handle events */
static PyObject *pyInputHandler;
pyDXObject_InputHandler(int fd, Pointer arg)
{
    /* getchar(); */
    PyEval_CallObject(pyInputHandler,Py_BuildValue("()"));
}

static char pyDXObject_DXRegisterInputHandler__doc__[] =
"Set function for input handler to call";

static PyObject *
pyDXObject_DXRegisterInputHandler(PyObject *self,PyObject * args)
{
  PyObject *result = NULL;
  PyObject *func;
  Error e;
  if (!PyArg_ParseTuple(args, "O", &func)) return NULL;
  if (!PyCallable_Check(func)) {
    PyErr_SetString(PyExc_TypeError, "parameter must be callable");
    return NULL;
    }
  Py_XINCREF(func);
  Py_XDECREF(pyInputHandler);
  pyInputHandler = func;
  e = DXRegisterInputHandler(pyDXObject_InputHandler,STDIN_FILENO,0);
  if (e == OK) {
    Py_INCREF(Py_None);
    return Py_None;
    }
  else {
    PyErr_SetString(ErrorObject,"Error during DXRegisterInputHandler");
    return NULL;
    }
}


/* List of methods defined in the module */

static struct PyMethodDef pyDXObject_methods[] = {
 {"DXObject_fromarray",(PyCFunction)pyDXObject_DXObject_fromarray,
         METH_VARARGS,pyDXObject_DXObject_fromarray__doc__},
 {"DXCallModule",(PyCFunction)pyDXObject_DXCallModule,
         METH_VARARGS,pyDXObject_DXCallModule__doc__},
 {"DXReference",(PyCFunction)pyDXObject_DXReference,
         METH_VARARGS,pyDXObject_DXReference__doc__},
 {"DXCheckRIH",(PyCFunction)pyDXObject_DXCheckRIH,
         METH_VARARGS,pyDXObject_DXCheckRIH__doc__},
 {"DXDelete",(PyCFunction)pyDXObject_DXDelete,
         METH_VARARGS,pyDXObject_DXDelete__doc__},
 {"DXRegisterInputHandler",(PyCFunction)pyDXObject_DXRegisterInputHandler,
         METH_VARARGS,pyDXObject_DXRegisterInputHandler__doc__},
 {"DXNewArray",(PyCFunction)pyDXObject_DXNewArray,
         METH_VARARGS,pyDXObject_DXNewArray__doc__},
 {"DXAddArrayData",(PyCFunction)pyDXObject_DXAddArrayData,
         METH_VARARGS,pyDXObject_DXAddArrayData__doc__},
 {"DXNewField",(PyCFunction)pyDXObject_DXNewField,
         METH_VARARGS,pyDXObject_DXNewField__doc__},
 {"DXSetStringAttribute",(PyCFunction)pyDXObject_DXSetStringAttribute,
         METH_VARARGS,pyDXObject_DXSetStringAttribute__doc__},
 {"DXSetComponentValue",(PyCFunction)pyDXObject_DXSetComponentValue,
         METH_VARARGS,pyDXObject_DXSetComponentValue__doc__},
 {"DXEndField",(PyCFunction)pyDXObject_DXEndField,
         METH_VARARGS,pyDXObject_DXEndField__doc__},
 {"DXMakeStringList",(PyCFunction)pyDXObject_DXMakeStringList,
         METH_VARARGS,pyDXObject_DXMakeStringList__doc__},
 {NULL,(PyCFunction)NULL, 0, NULL}                /* sentinel */
};


/* Initialization function for the module (*must* be called initpyDXObject) */

static char pyDXObject_module_documentation[] = "";

void
initpyDXObject()
{
  PyObject *m, *d;

  /* Create the module and add the functions */
  m = Py_InitModule4("pyDXObject", pyDXObject_methods,
                pyDXObject_module_documentation,
                (PyObject*)NULL,PYTHON_API_VERSION);

  /* Add some symbolic constants to the module */
  d = PyModule_GetDict(m);
  ErrorObject = PyString_FromString("pyDXObject.error");
  PyDict_SetItemString(d, "error", ErrorObject);

  /* XXXX Add constants here */
        
  PyDict_SetItemString(d,"TYPE_BYTE",PyInt_FromLong((long)TYPE_BYTE));
  PyDict_SetItemString(d,"TYPE_HYPER",PyInt_FromLong((long)TYPE_HYPER));
  PyDict_SetItemString(d,"TYPE_SHORT",PyInt_FromLong((long)TYPE_SHORT));
  PyDict_SetItemString(d,"TYPE_UBYTE",PyInt_FromLong((long)TYPE_UBYTE));
  PyDict_SetItemString(d,"TYPE_INT",PyInt_FromLong((long)TYPE_INT));
  PyDict_SetItemString(d,"TYPE_USHORT",PyInt_FromLong((long)TYPE_USHORT));
  PyDict_SetItemString(d,"TYPE_DOUBLE",PyInt_FromLong((long)TYPE_DOUBLE));
  PyDict_SetItemString(d,"TYPE_UINT",PyInt_FromLong((long)TYPE_UINT));
  PyDict_SetItemString(d,"TYPE_STRING",PyInt_FromLong((long)TYPE_STRING));
  PyDict_SetItemString(d,"TYPE_FLOAT",PyInt_FromLong((long)TYPE_FLOAT));
  PyDict_SetItemString(d,"CATEGORY_REAL",PyInt_FromLong((long)CATEGORY_REAL));
  PyDict_SetItemString(d,"CATEGORY_COMPLEX",PyInt_FromLong((long)CATEGORY_COMPLEX));

  /* Check for errors */
  if (PyErr_Occurred())
          Py_FatalError("can't initialize module pyDXObject");

  import_array();
  DXInitModules();
}




