ARCH = SP
MACHINE = SP
HOME = /usr/common/homes/g/grote
PYVERS = 1.5
PYTHON = $(HOME)/Python-$(PYVERS)
PYWRAPPERHOME = $(HOME)/pywrapper
PYWRAPPERARGS =
#VAMPIR = $(VTINC) $(VTLIB)
OPT = -O3 -qstrict -qarch=pwr3 -qtune=pwr3
POPT = -O
PYMODULEOPT = 
DATASIZE = -bmaxdata:0x20000000
F90FREE  = mpxlf -c -qmaxmem=8192 -u -qdpc=e -qintsize=4 -WF,-DPARALLEL -qsuffix=f=f90:cpp=F90 -qfree=f90 $(VAMPIR) $(DATASIZE)
F90FIXED = mpxlf -c -qmaxmem=8192 -u -qdpc=e -qintsize=4 -WF,-DPARALLEL -qfixed=132 $(VAMPIR) $(DATASIZE)
CC = mpcc -qmaxmem=8192 $(DATASIZE)
SO = so
LD = mpcc $(VAMPIR) $(DATASIZE)
LDSHARED = $(PYTHON)/lib/python$(PYVERS)/config/ld_so_aix mpcc -bI:$(PYTHON)/lib/python$(PYVERS)/config/python.exp
FORTLIB = -lxlf90 -lxlopt -lxlf -lxlomp_ser
LINKLIB = -lpthread
DLLIB =
XLIB = 
PACTLIB = -L$(HOME)/pact/lib
YORICKLIB = -L$(HOME)/yorick/lib/yorick/1.4/lib
TCLLIB = 
LAPACKLIB = 
NUMPYINCLUDE = 
OLD = OLD1

WARPSCRIPTS = $(HOME)/scriptsnew
PYTHONSOS = $(PYTHON)/lib/python$(PYVERS)/site-packages/*/*.so
OTHERLIBS = $(PACTLIB) -lpdb -lpml -lscore \
            $(YORICKLIB) -lgist \
            $(LAPACK) -lessl
PYMPI = /usr/common/homes/g/grote/pyMPI
PYMPIOBJS = $(PYMPI)/initmpi.o $(PYMPI)/driver.o $(PYMPI)/patched152main.o

warp: build_static

include Makefile.dependencies

build_static: $(WARPOBJECTS) $(WARPSLAVEOBJECTS) $(MODULEOBJECTS) config.o
	$(LD) -o pywarp90 $(PYMPIOBJS) $(WARPOBJECTS) $(WARPSLAVEOBJECTS) \
	$(MODULEOBJECTS) config.o \
	-lm -lX11 $(PYTHONLIB) $(LINKLIB) $(OTHERLIBS) $(FORTLIB)

public:
	adddate.py pywarp90
	mv pywarp902* $(HOME)/pub
	chmod go+rx $(HOME)/pub/pywarp90*

build_dynamic: warpC.$(SO)
	#cp warpC.$(SO) $(WARPSCRIPTS)
warpC.$(SO): warpC.c $(MODULEOBJECTS) $(WARPOBJECTS) $(WARPSLAVEOBJECTS)
	$(CC) -c warpC.c $(INCLUDE) $(NUMPYINCLUDE)
	$(LD) -o warpC.$(SO) warpC.o $(MODULEOBJECTS) $(WARPOBJECTS) \
	$(WARPSLAVEOBJECTS) $(OTHERLIBS) $(FORTLIB)

