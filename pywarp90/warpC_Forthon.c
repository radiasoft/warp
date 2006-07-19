/* Created by David P. Grote */
/* $Id: warpC_Forthon.c,v 1.5 2006/07/19 18:32:49 jlvay Exp $ */
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
void initwarpC(void)
{
  PyObject *m, *d;
  PyObject *pystdout;
  m = Py_InitModule("warpC", warpC_methods);
  d = PyModule_GetDict(m);
  ErrorObject = PyString_FromString("warpC.error");
  PyDict_SetItemString(d, "error", ErrorObject);
  if (PyErr_Occurred())
    Py_FatalError("can not initialize module warpC");

  pystdout = PySys_GetObject("stdout");
  PyFile_WriteString("Forthon edition\n",pystdout);

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
  initem2dpy();
}


