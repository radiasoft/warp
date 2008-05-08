/* Created by David P. Grote */
/* $Id: warpC_Forthon.c,v 1.8 2008/05/08 20:04:31 jlvay Exp $ */
/* This is a stub module which calls the init functions of all of            */
/* the modules that are part of WARP. This is needed since the modules       */
/* depend on each other and so must be incorporated into one shared          */
/* object file.                                                              */
#include "Forthon.h"

#ifndef PyMODINIT_FUNC
#define PyMODINIT_FUNC void
#endif

extern PyMODINIT_FUNC inittoppy(void);
extern PyMODINIT_FUNC initenvpy(void);
extern PyMODINIT_FUNC initw3dpy(void);
extern PyMODINIT_FUNC initf3dpy(void);
extern PyMODINIT_FUNC initwxypy(void);
extern PyMODINIT_FUNC initfxypy(void);
extern PyMODINIT_FUNC initwrzpy(void);
extern PyMODINIT_FUNC initfrzpy(void);
extern PyMODINIT_FUNC initcirpy(void);
extern PyMODINIT_FUNC initherpy(void);
extern PyMODINIT_FUNC initchopy(void);
extern PyMODINIT_FUNC initem2dpy(void);
extern PyMODINIT_FUNC initem3dpy(void);

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
  initem3dpy();
}


