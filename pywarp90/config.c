/* Generated automatically from ./config.c.in by makesetup. */
/* -*- C -*- ***********************************************
Copyright 1991-1995 by Stichting Mathematisch Centrum, Amsterdam,
The Netherlands.

                        All Rights Reserved

Permission to use, copy, modify, and distribute this software and its
documentation for any purpose and without fee is hereby granted,
provided that the above copyright notice appear in all copies and that
both that copyright notice and this permission notice appear in
supporting documentation, and that the names of Stichting Mathematisch
Centrum or CWI or Corporation for National Research Initiatives or
CNRI not be used in advertising or publicity pertaining to
distribution of the software without specific, written prior
permission.

While CWI is the initial source for this software, a modified version
is made available by the Corporation for National Research Initiatives
(CNRI) at the Internet address ftp://ftp.python.org.

STICHTING MATHEMATISCH CENTRUM AND CNRI DISCLAIM ALL WARRANTIES WITH
REGARD TO THIS SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS, IN NO EVENT SHALL STICHTING MATHEMATISCH
CENTRUM OR CNRI BE LIABLE FOR ANY SPECIAL, INDIRECT OR CONSEQUENTIAL
DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR
PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER
TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
PERFORMANCE OF THIS SOFTWARE.

******************************************************************/

/* Module configuration */

/* !!! !!! !!! This file is edited by the makesetup script !!! !!! !!! */

/* This file contains the table of built-in modules.
   See init_builtin() in import.c. */

#include "Python.h"


extern void initarrayfns();
extern void init_numpy();
extern void initmultiarray();
extern void initumath();
extern void initfftpack();
extern void initlapack_lite();
extern void initranlib();
extern void initgistC();
extern void initpypdb();
extern void initRNG();
extern void initpybasisC();
extern void inittoppy();
extern void initenvpy();
extern void initf3dpy();
extern void initw3dpy();
extern void initfxypy();
extern void initwxypy();
extern void initcirpy();
extern void initherpy();
extern void initregex();
extern void initpcre();
extern void initposix();
extern void initsignal();
extern void initarray();
extern void initcmath();
extern void initmath();
extern void initstrop();
extern void initstruct();
extern void inittime();
extern void initoperator();
extern void initfcntl();
extern void initpwd();
extern void initgrp();
extern void initselect();
extern void initsocket();
extern void initerrno();
extern void initmd5();
extern void initsha();
extern void init_tkinter();
extern void initrotor();
extern void initnew();
extern void initbinascii();
extern void initparser();
extern void initcStringIO();
extern void initcPickle();

/* -- ADDMODULE MARKER 1 -- */

extern void PyMarshal_Init();
extern void initimp();

struct _inittab _PyImport_Inittab[] = {

	{"arrayfns", initarrayfns},
	{"_numpy", init_numpy},
	{"multiarray", initmultiarray},
	{"umath", initumath},
	{"fftpack", initfftpack},
	{"lapack_lite", initlapack_lite},
	{"ranlib", initranlib},
	{"gistC", initgistC},
	{"pypdb", initpypdb},
	{"RNG", initRNG},
	{"pybasisC", initpybasisC},
	{"toppy", inittoppy},
	{"envpy", initenvpy},
	{"f3dpy", initf3dpy},
	{"w3dpy", initw3dpy},
	{"fxypy", initfxypy},
	{"wxypy", initwxypy},
	{"cirpy", initcirpy},
	{"herpy", initherpy},
	{"regex", initregex},
	{"pcre", initpcre},
	{"posix", initposix},
	{"signal", initsignal},
	{"array", initarray},
	{"cmath", initcmath},
	{"math", initmath},
	{"strop", initstrop},
	{"struct", initstruct},
	{"time", inittime},
	{"operator", initoperator},
	{"fcntl", initfcntl},
	{"pwd", initpwd},
	{"grp", initgrp},
	{"select", initselect},
	{"socket", initsocket},
	{"errno", initerrno},
	{"md5", initmd5},
	{"sha", initsha},
	{"_tkinter", init_tkinter},
	{"rotor", initrotor},
	{"new", initnew},
	{"binascii", initbinascii},
	{"parser", initparser},
	{"cStringIO", initcStringIO},
	{"cPickle", initcPickle},

/* -- ADDMODULE MARKER 2 -- */

	/* This module "lives in" with marshal.c */
	{"marshal", PyMarshal_Init},

	/* This lives it with import.c */
	{"imp", initimp},

	/* These entries are here for sys.builtin_module_names */
	{"__main__", NULL},
	{"__builtin__", NULL},
	{"sys", NULL},

	/* Sentinel */
	{0, 0}
};
