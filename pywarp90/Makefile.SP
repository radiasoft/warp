ARCH = SP
MACHINE = SP
HOME = /usr/common/homes/g/grote
PYVERS = 2.2
PYTHON = $(HOME)/Python-$(PYVERS)
PYWRAPPERHOME = $(HOME)/pywrapper
PYWRAPPERARGS =
#VAMPIR = $(VTINC) $(VTLIB)
OPT = -O3 -qstrict -qarch=pwr3 -qtune=pwr3
#OPT = -O4 -qarch=pwr3 -qtune=pwr3 -qunroll
#OPT = -O3 -qarch=pwr3 -qtune=pwr3 -qhot=vector
#OPT = -g -C -qsigtrap
POPT = -O
PYMODULEOPT = 
DATASIZE = -bmaxdata:0x80000000 -bmaxstack:0x10000000	
#DATASIZE = -bmaxdata:0x70000000 -bmaxstack:0x10000000	# for debugging
ESSL = -WF,-DESSL
F90FREE  = mpxlf -c -qmaxmem=8192 -u -qdpc=e -qintsize=4 -qsave=defaultinit -WF,-DMPIPARALLEL -qsuffix=f=f90:cpp=F90 -qfree=f90 $(VAMPIR) $(DATASIZE) $(ESSL)
F90FIXED = mpxlf -c -qmaxmem=8192 -u -qdpc=e -qintsize=4 -qsave=defaultinit -WF,-DMPIPARALLEL -qfixed=132 $(VAMPIR) $(DATASIZE) $(ESSL)
CC = mpcc_r -qmaxmem=8192 $(DATASIZE)
SO = so
LD = mpcc_r $(VAMPIR) $(DATASIZE)
LD = mpxlf_r $(VAMPIR) $(DATASIZE)
LDSHARED = $(PYTHON)/lib/python$(PYVERS)/config/ld_so_aix mpcc_r -bI:$(PYTHON)/lib/python$(PYVERS)/config/python.exp
FINCLUDE = $(HPMTOOLKIT)
FORTLIB = 
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
PYTHONSOS = $(PYTHON)/lib/python$(PYVERS)/site-packages/*/*.so \
            $(PYTHON)/lib/python$(PYVERS)/site-packages/*/*/*.so
OTHERLIBS = $(PACTLIB) -lpdb -lpml -lscore \
            $(YORICKLIB) -lgist \
            $(LAPACK) -lessl \
            -lld -lreadline -ltermcap -L/usr/common/homes/g/grote/readline/lib -ldl
#PYMPI = /usr/common/homes/g/grote/pyMPI
#PYMPIOBJS = $(PYMPI)/initmpi.o $(PYMPI)/driver.o $(PYMPI)/patched152main.o
PYMPI = /usr/common/homes/g/grote/pyMPI-1.2a6
PYMPIOBJS = $(PYMPI)/driver.o $(PYMPI)/patchedmain.o -L$(PYMPI) -lpympi

warp: build_static

include Makefile.dependencies

build_static: $(WARPOBJECTS) $(WARPSLAVEOBJECTS) $(MODULEOBJECTS) config.o
	$(LD) -g -o pywarp90 $(PYMPIOBJS) \
	$(WARPOBJECTS) $(WARPSLAVEOBJECTS) $(MODULEOBJECTS) config.o \
	$(PYTHONSOS) \
	-bE:$(PYTHON)/lib/python2.2/config/python.exp \
	-lm -lX11 $(PYTHONLIB) $(LINKLIB) $(OTHERLIBS) $(FORTLIB)  \
	$(SP_TIMER)
	#$(HPMTOOLKIT) $(TRACE_MPIF) $(TRACE_MPI)

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

