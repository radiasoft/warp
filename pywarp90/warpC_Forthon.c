/* Created by David P. Grote */
/* $Id: warpC_Forthon.c,v 1.2 2004/02/10 01:09:55 dave Exp $ */
/* This is a stub module which calls the init functions of all of            */
/* the modules that are part of WARP. This is needed since the modules       */
/* depend on each other and so must be incorporated into one shared          */
/* object file.                                                              */
#include "Forthon.h"

/* ######################################################################### */
/* # Method list                                                             */
static struct PyMethodDef warpC_methods[] = {
  {NULL,NULL}};

/* ######################################################################### */
/* # The initialization function                                             */
void initwarpC()
{
  PyObject *m, *d;
  m = Py_InitModule("warpC", warpC_methods);
  d = PyModule_GetDict(m);
  ErrorObject = PyString_FromString("warpC.error");
  PyDict_SetItemString(d, "error", ErrorObject);
  if (PyErr_Occurred())
    Py_FatalError("can not initialize module warpC");

  printf("Forthon edition\n");

  import_array();

  inittoppy();
  initenvpy();
  initw3dpy();
  initf3dpy();
  initwxypy();
  initfxypy();
  initwrzpy();
  initfrzpy();
  initcirpy();
  initherpy();
  initchopy();
}


