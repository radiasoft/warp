ARCH = SP
MACHINE = SP
HOME = /usr/common/homes/g/grote
PYVERS = 1.5
PYTHON = $(HOME)/Python-$(PYVERS)
PYMAC = $(HOME)/pywrapper/pymac.py --f90 -t $(MACHINE)
PYCTOF = $(HOME)/pywrapper/pyctof.py
OPT = -O3 -qstrict -qarch=pwr3 -qtune=pwr3
POPT = -O
PYMODULEOPT = 
F90FREE  = mpxlf -c -qmaxmem=8192 -u -qdpc=e -qintsize=4 -WF,-DPARALLEL -qsuffix=f=f90:cpp=F90 -qfree=f90
F90FIXED = mpxlf -c -qmaxmem=8192 -u -qdpc=e -qintsize=4 -WF,-DPARALLEL -qfixed=132
CC = xlc -qmaxmem=8192
SO = sl
LD = mpcc
FORTLIB = -lxlf90 -lxlopt -lxlf -lxlomp_ser
LINKLIB =
DLLIB =
XLIB = 
PACTLIB = -L$(HOME)/pact/lib
YORICKLIB = -L$(HOME)/yorick/lib/yorick/1.4/lib
TCLLIB = 
LAPACKLIB = 
NUMPYINCLUDE = 
OLD = OLD1

WARPSCRIPTS = $(HOME)/scriptsnew
INCLUDE = -I$(PYTHON)/include/python$(PYVERS) -I$(HOME)/python
PYTHONLIB = -L$(PYTHON)/lib/python$(PYVERS)/config -lpython$(PYVERS)
PYTHONOBJ = $(PYTHON)/lib/python$(PYVERS)/config/python.o
OTHERLIBS = $(PACTLIB) -lpdb -lpml -lscore \
            $(YORICKLIB) -lgist \
            $(LAPACK) -lessl
RNGOBJ = ranf.o pmath_rng.o
RNGSRC = ranf.c pmath_rng.c
PYMPI = /usr/common/homes/g/grote/pyMPI
PYMPIOBJS = $(PYMPI)/initmpi.o $(PYMPI)/driver.o $(PYMPI)/patched152main.o

SHAREDOBJECTS = envpy.$(SO) frzpy.$(SO) wxypy.$(SO) f3dpy.$(SO) fxypy.$(SO) \
                w3dpy.$(SO) toppy.$(SO) wrzpy.$(SO) cirpy.$(SO) herpy.$(SO) \
		pybasisC.$(SO)

WARPOBJECTS = env.o \
              frz.o \
              wxy.o \
              f3d.o \
              fxy.o \
              util.o \
              w3d.o \
              dtop.o \
              f3d_mgrid.o \
              dw3d.o \
              dwrz.o \
              fft.o \
              top.o \
              wrz.o \
              cir.o \
              her.o

WARPSLAVEOBJECTS = topslave.o w3dslave.o f3dslave.o

MODULEOBJECTS = toppymodule.o envpymodule.o w3dpymodule.o f3dpymodule.o \
                wxypymodule.o fxypymodule.o wrzpymodule.o frzpymodule.o \
		cirpymodule.o herpymodule.o \
                top.p.o env.p.o w3d.p.o f3d.p.o wxy.p.o fxy.p.o \
		wrz.p.o frz.p.o \
		cir.p.o her.p.o pybasisC.o
MODULESOURCES = toppymodule.c envpymodule.c w3dpymodule.c f3dpymodule.c \
                wxypymodule.c fxypymodule.c wrzpymodule.c frzpymodule.c \
		cirpymodule.c herpymodule.c \
		pybasisC.c

static: $(WARPOBJECTS) $(WARPSLAVEOBJECTS) $(MODULEOBJECTS) config.o
	$(LD) -o pywarp90 $(PYMPIOBJS) $(WARPOBJECTS) $(WARPSLAVEOBJECTS) \
	$(MODULEOBJECTS) config.o \
	-lm -lX11 $(PYTHONLIB) $(LINKLIB) $(OTHERLIBS) $(FORTLIB)

warp: $(SHAREDOBJECTS)
	#cp *.$(SO) $(WARPSCRIPTS)
	#tar czf linuxpywarpso.tgz *.so
	#scp linuxpywarpso.tgz ife1:pub

config.o:config.c
	$(CC) -c config.c $(INCLUDE)

pybasisC.$(SO):pybasisC.o 
	#$(LD) -o pybasisC.$(SO) pybasisC.o $(RANFOBJS) $(INCLUDE) $(OTHERLIBS)
	$(LD) -o pybasisC.$(SO) pybasisC.o $(FORTLIB)
pybasisC.o:pybasisC.c
	$(CC) -g -c pybasisC.c $(INCLUDE) $(NUMPYINCLUDE)
pybasisC.c:$(HOME)/pywrapper/pybasisC.m pybasis.h
	cp $(HOME)/pywrapper/pybasisC.m .
	$(PYCTOF) --f90 -t$(ARCH) pybasisC
pybasis.h:$(HOME)/pywrapper/pybasis.h
	cp $(HOME)/pywrapper/pybasis.h pybasis.h

toppy.$(SO):toppymodule.o top.p.o top.o dtop.o fft.o util.o
	$(LD) -o toppy.$(SO) toppymodule.o top.p.o top.o dtop.o fft.o util.o -lm $(XLIB) -lX11 $(INCLUDE) $(OTHERLIBS) $(DLLIB) $(FORTLIB)

top.p.o:top.p.F90
	$(F90FREE) -c top.p.F90
	#$(F90FREE) $(POPT) -c top.p.F90
top.o:top.F top.p.o
	$(F90FIXED) $(OPT) -c top.F
dtop.o:dtop.F top.p.o
	$(F90FIXED) $(OPT) -c dtop.F
topslave.o:topslave.F top.p.o
	$(F90FIXED) $(OPT) -c topslave.F
fft.o:fft.F
	$(F90FIXED) $(OPT) -c fft.F
util.o:util.F
	$(F90FIXED) $(OPT) -c util.F

toppymodule.o:toppymodule.c pybasis.h
	$(CC) $(PYMODULEOPT) -c toppymodule.c $(INCLUDE) $(NUMPYINCLUDE)
toppymodule.c top.p.F90:top.v
	$(PYMAC) -a top.v

envpy.$(SO):envpymodule.o env.p.o env.o
	$(LD) -o envpy.$(SO) envpymodule.o env.p.o env.o -lm $(XLIB) -lX11 $(INCLUDE) $(OTHERLIBS) $(DLLIB) $(FORTLIB)

env.p.o:env.p.F90
	$(F90FREE) $(POPT) -c env.p.F90
env.o:env.F env.p.o top.p.o
	$(F90FIXED) $(OPT) -c env.F

envpymodule.o:envpymodule.c pybasis.h
	$(CC) $(PYMODULEOPT) -c envpymodule.c $(INCLUDE) $(NUMPYINCLUDE)
envpymodule.c env.p.F90:env.v
	$(PYMAC) -a env.v top.v


w3dpy.$(SO):w3dpymodule.o w3d.p.o w3d.o dw3d.o
	$(LD) -o w3dpy.$(SO) w3dpymodule.o w3d.p.o w3d.o dw3d.o -lm $(XLIB) -lX11 $(INCLUDE) $(OTHERLIBS) $(DLLIB) $(FORTLIB)

w3d.p.o:w3d.p.F90
	$(F90FREE) $(POPT) -c w3d.p.F90
w3d.o:w3d.F w3d.p.o top.p.o
	$(F90FIXED) $(OPT) -c w3d.F
dw3d.o:dw3d.F w3d.p.o top.p.o
	$(F90FIXED) $(OPT) -c dw3d.F
w3dslave.o:w3dslave.F w3d.p.o top.p.o
	$(F90FIXED) $(OPT) -c w3dslave.F

w3dpymodule.o:w3dpymodule.c pybasis.h
	$(CC) $(PYMODULEOPT) -c w3dpymodule.c $(INCLUDE) $(NUMPYINCLUDE)
w3dpymodule.c w3d.p.F90:w3d.v
	$(PYMAC) -a w3d.v top.v


f3dpy.$(SO):f3dpymodule.o f3d.p.o f3d.o f3d_mgrid.o
	$(LD) -o f3dpy.$(SO) f3dpymodule.o f3d.p.o f3d.o f3d_mgrid.o -lm $(XLIB) -lX11 $(INCLUDE) $(OTHERLIBS) $(DLLIB) $(FORTLIB)

f3d.p.o:f3d.p.F90
	$(F90FREE) $(POPT) -c f3d.p.F90
f3d.o:f3d.F f3d.p.o top.p.o
	$(F90FIXED) $(OPT) -c f3d.F
f3d_mgrid.o:f3d_mgrid.F f3d.p.o top.p.o
	$(F90FIXED) $(OPT) -c f3d_mgrid.F
f3dslave.o:f3dslave.F f3d.p.o top.p.o
	$(F90FIXED) $(OPT) -c f3dslave.F

f3dpymodule.o:f3dpymodule.c pybasis.h
	$(CC) $(PYMODULEOPT) -c f3dpymodule.c $(INCLUDE) $(NUMPYINCLUDE)
f3dpymodule.c f3d.p.F90:f3d.v
	$(PYMAC) -a f3d.v top.v

wrzpy.$(SO):wrzpymodule.o wrz.p.o wrz.o dwrz.o
	$(LD) -o wrzpy.$(SO) wrzpymodule.o wrz.p.o wrz.o dwrz.o -lm $(XLIB) -lX11 $(INCLUDE) $(OTHERLIBS) $(DLLIB) $(FORTLIB)

wrz.p.o:wrz.p.F90
	$(F90FREE) $(POPT) -c wrz.p.F90
wrz.o:wrz.F wrz.p.o top.p.o
	$(F90FIXED) $(OPT) -c wrz.F
dwrz.o:dwrz.F wrz.p.o top.p.o
	$(F90FIXED) $(OPT) -c dwrz.F

wrzpymodule.o:wrzpymodule.c pybasis.h
	$(CC) $(PYMODULEOPT) -c wrzpymodule.c $(INCLUDE) $(NUMPYINCLUDE)
wrzpymodule.c wrz.p.F90:wrz.v
	$(PYMAC) -a wrz.v top.v


frzpy.$(SO):frzpymodule.o frz.p.o frz.o
	$(LD) -o frzpy.$(SO) frzpymodule.o frz.p.o frz.o -lm $(XLIB) -lX11 $(INCLUDE) $(OTHERLIBS) $(DLLIB) $(FORTLIB)

frz.p.o:frz.p.F90
	$(F90FREE) $(POPT) -c frz.p.F90
frz.o:frz.F frz.p.o top.p.o
	$(F90FIXED) $(OPT) -c frz.F

frzpymodule.o:frzpymodule.c pybasis.h
	$(CC) $(PYMODULEOPT) -c frzpymodule.c $(INCLUDE) $(NUMPYINCLUDE)
frzpymodule.c frz.p.F90:frz.v
	$(PYMAC) -a frz.v top.v


wxypy.$(SO):wxypymodule.o wxy.p.o wxy.o
	$(LD) -o wxypy.$(SO) wxypymodule.o wxy.p.o wxy.o -lm $(XLIB) -lX11 $(INCLUDE) $(OTHERLIBS) $(DLLIB) $(FORTLIB)

wxy.p.o:wxy.p.F90
	$(F90FREE) $(POPT) -c wxy.p.F90
wxy.o:wxy.F wxy.p.o top.p.o w3d.p.o
	$(F90FIXED) -c wxy.F

wxypymodule.o:wxypymodule.c pybasis.h
	$(CC) $(PYMODULEOPT) -c wxypymodule.c $(INCLUDE) $(NUMPYINCLUDE)
wxypymodule.c wxy.p.F90:wxy.v
	$(PYMAC) -a wxy.v top.v


fxypy.$(SO):fxypymodule.o fxy.p.o fxy.o
	$(LD) -o fxypy.$(SO) fxypymodule.o fxy.p.o fxy.o -lm $(XLIB) -lX11 $(INCLUDE) $(OTHERLIBS) $(DLLIB) $(FORTLIB)

fxy.p.o:fxy.p.F90
	$(F90FREE) $(POPT) -c fxy.p.F90
fxy.o:fxy.F fxy.p.o top.p.o
	$(F90FIXED) $(OPT) -c fxy.F

fxypymodule.o:fxypymodule.c pybasis.h
	$(CC) $(PYMODULEOPT) -c fxypymodule.c $(INCLUDE) $(NUMPYINCLUDE)
fxypymodule.c fxy.p.F90:fxy.v
	$(PYMAC) -a fxy.v top.v


cirpy.$(SO):cirpymodule.o cir.p.o cir.o
	$(LD) -o cirpy.$(SO) cirpymodule.o cir.p.o cir.o -lm $(XLIB) -lX11 $(INCLUDE) $(OTHERLIBS) $(DLLIB) $(FORTLIB)

cir.p.o:cir.p.F90
	$(F90FREE) $(POPT) -c cir.p.F90
cir.o:cir.F cir.p.o top.p.o
	$(F90FIXED) $(OPT) -c cir.F

cirpymodule.o:cirpymodule.c pybasis.h
	$(CC) $(PYMODULEOPT) -c cirpymodule.c $(INCLUDE) $(NUMPYINCLUDE)
cirpymodule.c cir.p.F90:cir.v
	$(PYMAC) -a cir.v top.v


herpy.$(SO):herpymodule.o her.p.o her.o
	$(LD) -o herpy.$(SO) herpymodule.o her.p.o her.o -lm $(XLIB) -lX11 $(INCLUDE) $(OTHERLIBS) $(DLLIB) $(FORTLIB)

her.p.o:her.p.F90
	$(F90FREE) $(POPT) -c her.p.F90
her.o:her.F her.p.o top.p.o
	$(F90FIXED) $(OPT) -c her.F

herpymodule.o:herpymodule.c pybasis.h
	$(CC) $(PYMODULEOPT) -c herpymodule.c $(INCLUDE) $(NUMPYINCLUDE)
herpymodule.c her.p.F90:her.v
	$(PYMAC) -a her.v top.v


clean:
	rm -f *.mod *.o *.p.F90 *.p.F90 *module.c *.scalars *.$(SO) pybasis*
