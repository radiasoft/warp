/* Generated automatically from ./config.c.in by makesetup. */
/* -*- C -*- ***********************************************
Copyright (c) 2000, BeOpen.com.
Copyright (c) 1995-2000, Corporation for National Research Initiatives.
Copyright (c) 1990-1995, Stichting Mathematisch Centrum.
All rights reserved.

See the file "Misc/COPYRIGHT" for information on usage and
redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
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
/* extern void initpypdb(); */
extern void initRNG();
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

extern void initthread(void);
extern void initgc(void);
extern void initbsddb(void);
extern void initregex(void);
extern void initpcre(void);
extern void initposix(void);
extern void initsignal(void);
extern void init_sre(void);
extern void initarray(void);
extern void initcmath(void);
extern void initmath(void);
extern void initstrop(void);
extern void initstruct(void);
extern void inittime(void);
extern void initoperator(void);
extern void init_codecs(void);
extern void initunicodedata(void);
extern void initucnhash(void);
extern void init_locale(void);
extern void initfcntl(void);
extern void initpwd(void);
extern void initgrp(void);
extern void initerrno(void);
extern void initselect(void);
extern void init_socket(void);
extern void initmmap(void);
extern void initmd5(void);
extern void initsha(void);
extern void init_tkinter(void);
extern void initrotor(void);
extern void initnew(void);
extern void initbinascii(void);
extern void initparser(void);
extern void initcStringIO(void);
extern void initcPickle(void);

/* -- ADDMODULE MARKER 1 -- */

extern void PyMarshal_Init(void);
extern void initimp(void);

struct _inittab _PyImport_Inittab[] = {

        {"arrayfns", initarrayfns},
        {"_numpy", init_numpy},
        {"multiarray", initmultiarray},
        {"umath", initumath},
        {"fftpack", initfftpack},
        {"lapack_lite", initlapack_lite},
        {"ranlib", initranlib},
        {"gistC", initgistC},
        /* {"pypdb", initpypdb}, */
        {"RNG", initRNG},
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

	{"thread", initthread},
	{"gc", initgc},
	{"bsddb", initbsddb},
	{"regex", initregex},
	{"pcre", initpcre},
	{"posix", initposix},
	{"signal", initsignal},
	{"_sre", init_sre},
	{"array", initarray},
	{"cmath", initcmath},
	{"math", initmath},
	{"strop", initstrop},
	{"struct", initstruct},
	{"time", inittime},
	{"operator", initoperator},
	{"_codecs", init_codecs},
	{"unicodedata", initunicodedata},
	{"ucnhash", initucnhash},
	{"_locale", init_locale},
	{"fcntl", initfcntl},
	{"pwd", initpwd},
	{"grp", initgrp},
	{"errno", initerrno},
	{"select", initselect},
	{"_socket", init_socket},
	{"mmap", initmmap},
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
	{"exceptions", NULL},

	/* Sentinel */
	{0, 0}
};
