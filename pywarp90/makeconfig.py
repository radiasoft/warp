#!/bin/env python
"""
Added the Warp init routines to python's config.c file.
"""
import sys
import os
import re

configfile = os.path.join(sys.exec_prefix,'lib','python'+sys.version[:3],'config','config.c')

ff = open(configfile,'r')
configlines = ff.readlines()
ff.close()

addmodule1 = """
extern void initpybasisC();
extern void inittoppy();
extern void initenvpy();
extern void initf3dpy();
extern void initw3dpy();
extern void initfxypy();
extern void initwxypy();
extern void initfrzpy();
extern void initwrzpy();
extern void initcirpy();
extern void initherpy();
extern void initchopy();
"""

addmodule2 = """
        {"pybasisC", initpybasisC},
        {"toppy", inittoppy},
        {"envpy", initenvpy},
        {"f3dpy", initf3dpy},
        {"w3dpy", initw3dpy},
        {"fxypy", initfxypy},
        {"wxypy", initwxypy},
        {"frzpy", initfrzpy},
        {"wrzpy", initwrzpy},
        {"cirpy", initcirpy},
        {"herpy", initherpy},
        {"chopy", initchopy},
"""

ff = open("config.c",'w')
for line in configlines:
  ff.write(line)
  if re.search("ADDMODULE MARKER 1",line): ff.write(addmodule1)
  if re.search("ADDMODULE MARKER 2",line): ff.write(addmodule2)

ff.close()
