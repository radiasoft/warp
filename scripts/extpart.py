"""Creates classes for handling extrapolated particle windows, ExtPart, and z-crossing locations, ZCrossingParticles.
Two functions are available for saving an ExtPart object in a file.
 - dumpExtPart(object,filename)
 - restoreExtPart(object,filename)
"""
__all__ = ['ExtPart','dumpExtPart','restoreExtPart','ZCrossingParticles']
from warp import *
from appendablearray import *
import cPickle
extpart_version = "$Id: extpart.py,v 1.81 2011/09/29 00:16:10 grote Exp $"

def extpartdoc():
    import extpart
    print extpart.__doc__

_extforcenorestore = 0
def extforcenorestore():
    global _extforcenorestore
    _extforcenorestore = 1
def extnoforcenorestore():
    global _extforcenorestore
    _extforcenorestore = 0

############################################################################
class ParticleAccumulator(object):
    """This virtual class defines a container to setup and keep track of
particle data. It can optionally accumulate the data over multiple time steps.
The creator options are:
 - iz: grid location where the extrapolated data is saved.
 - zz: lab location where data is saved.
 - nmax: max size of the arrays. Defaults to 3*top.pnumz[iz] if non-zero,
         otherwise 10000.
 - laccumulate=0: when true, particles are accumulated over multiple steps.
 - name=None: descriptive name for location. It must be unique for each
              instance of the diagnostic.
 - lautodump=0: when true, after the grid moves beyond the z location,
                automatically dump the data to a file, clear the arrays and
                disable itself. Also must have name set - name is appended to
                the name of the files where the data is stored. The object
                writes itself out to a pdb file.
 - dumptofile=0: when true, the particle data is always dumped to a file
                 and not saved in memory. Name must be set - name is appended
                 to the name of the files where the data is stored. Setting this
                 to true implies that the data is accumulated.
                 Use the restoredata method to read the data back in.

One of iz or zz must be specified.

Available methods:
 - setaccumulate(v=1): Turns accumulation on or off. If turned on, all
                       currently saved data is deleted.
 - clear(): Will clear the existing data.
 - disable(): Turns off collecting of data
 - enable(): Turns on collecting of data (only needed after disable)

The follow all take an optional argument to specify species number.
 - getn: Get number of particles
 - gett: Get time at which particle was saved
 - getx, y, ux, uy, uz, vx, vy, vz, xp, yp, r, theta, rp: Get the various
     coordinates or velocity of particles
 - getpid: Get the particle ID info

The following are available plot routines. All take an optional argument to
specify species number. Additional arguments are the same as the 'pp' plotting
routines (such as ppxxp).
 - pxy, pxvx, pyvy, pvxvy, pxxp, pyyp, pxpyp, prrp,
 - ptx, pty, ptxp, ptyp, ptux, ptuy, ptuz, ptvx
 - ptvy, ptvz, ptrace
    """

    """The following names must be setup in the inheriting class
self.type
self.topnwinname
self.topizwinname
self.topzzwinname
self.topwzwinname
self.topnmaxname
self.topnpidmaxname
self.topnname
self.toptname
self.topxname
self.topyname
self.topuxname
self.topuyname
self.topuzname
self.toppidname
self.topgroupname
    """
    name_cache = []

    def _topproperties(name):
        def fget(self):
            return getattr(top,getattr(self,name))
        def fset(self,value):
            return setattr(top,getattr(self,name),value)
        return fget,fset,None,None

    topnwin = property(*_topproperties('topnwinname'))
    topizwin = property(*_topproperties('topizwinname'))
    topzzwin = property(*_topproperties('topzzwinname'))
    topwzwin = property(*_topproperties('topwzwinname'))
    topnmax = property(*_topproperties('topnmaxname'))
    topnpidmax = property(*_topproperties('topnpidmaxname'))
    topn = property(*_topproperties('topnname'))
    topt = property(*_topproperties('toptname'))
    topx = property(*_topproperties('topxname'))
    topy = property(*_topproperties('topyname'))
    topux = property(*_topproperties('topuxname'))
    topuy = property(*_topproperties('topuyname'))
    topuz = property(*_topproperties('topuzname'))
    toppid = property(*_topproperties('toppidname'))

    def __init__(self,iz=-1,zz=0.,nmax=None,laccumulate=0,
                 name=None,lautodump=0,dumptofile=0):
        # --- Save input values, getting default values when needed
        assert type(iz) is IntType,"iz must be an integer"
        assert iz >= 0 or zz is not None,"Either iz or zz must be specified"
        self.iz = iz
        self.zz = zz
        self.laccumulate = laccumulate
        self.lautodump = lautodump
        if name is not None:
          if name in self.name_cache:
            raise ValueError('the name "%s" is not unique'%name)
          else:
            self.name_cache.append(name)
        self.name = name
        self.dumptofile = dumptofile
        self.dt = top.dt

        if nmax is None:
            self.nmax = 10000
            if top.allocated("pnumz") and 0 <= self.getiz() <= top.nzmmnt:
                if top.pnumz[self.getiz(),-1] > 0:
                    if top.nszmmnt > 1:
                        self.nmax = nint(max(top.pnumz[self.getiz(),:-1])*3)
                    else:
                        self.nmax = nint(top.pnumz[self.getiz(),-1]*3)
        else:
            self.nmax = nmax

        # --- Add this new window to the group in top
        self.enabled = 0
        self.enable()
        # --- Setup empty arrays for accumulation if laccumulate if true.
        # --- Otherwise, the arrays will just point to the data in the
        # --- top group.
        self.setuparrays(top.ns)

    def getiz(self):
        if self.iz >= 0:
            return self.iz
        else:
            if top.dzm != 0.:
                return int((self.zz - top.zmmntmin)*top.dzmi)
            else:
                return -1

    def getzz(self):
        if self.iz >= 0:
            return top.zmmntmin + self.iz*top.dzm
        else:
            return self.zz

    def setuparrays(self,ns,bump=None):
        self.restored = False
        if self.laccumulate and not self.dumptofile:
            if bump is None: bump = self.nmax
            self.t = []
            self.x = []
            self.y = []
            self.ux = []
            self.uy = []
            self.uz = []
            self.pid = []
            for js in range(ns):
                self.t.append(AppendableArray(self.nmax,typecode='d',
                                                autobump=bump))
                self.x.append(AppendableArray(self.nmax,typecode='d',
                                                autobump=bump))
                self.y.append(AppendableArray(self.nmax,typecode='d',
                                                autobump=bump))
                self.ux.append(AppendableArray(self.nmax,typecode='d',
                                                 autobump=bump))
                self.uy.append(AppendableArray(self.nmax,typecode='d',
                                                 autobump=bump))
                self.uz.append(AppendableArray(self.nmax,typecode='d',
                                                 autobump=bump))
                self.pid.append(AppendableArray(self.nmax,typecode='d',
                                                  autobump=bump,
                                                  unitshape=(self.topnpidmax,)))
        else:
            self.t = ns*[zeros(0,'d')]
            self.x = ns*[zeros(0,'d')]
            self.y = ns*[zeros(0,'d')]
            self.ux = ns*[zeros(0,'d')]
            self.uy = ns*[zeros(0,'d')]
            self.uz = ns*[zeros(0,'d')]
            self.pid = ns*[zeros((0,self.topnpidmax),'d')]

    def addspecies(self):
        if self.laccumulate and not self.dumptofile:
            for js in range(self.getns(),top.ns):
                bump = self.nmax
                self.t.append(AppendableArray(self.nmax,typecode='d',
                                                autobump=bump))
                self.x.append(AppendableArray(self.nmax,typecode='d',
                                                autobump=bump))
                self.y.append(AppendableArray(self.nmax,typecode='d',
                                                autobump=bump))
                self.ux.append(AppendableArray(self.nmax,typecode='d',
                                                 autobump=bump))
                self.uy.append(AppendableArray(self.nmax,typecode='d',
                                                 autobump=bump))
                self.uz.append(AppendableArray(self.nmax,typecode='d',
                                                 autobump=bump))
                self.pid.append(AppendableArray(self.nmax,typecode='d',
                                                  autobump=bump,
                                                  unitshape=(self.topnpidmax,)))
        else:
            self.t = top.ns*[zeros(0,'d')]
            self.x = top.ns*[zeros(0,'d')]
            self.y = top.ns*[zeros(0,'d')]
            self.ux = top.ns*[zeros(0,'d')]
            self.uy = top.ns*[zeros(0,'d')]
            self.uz = top.ns*[zeros(0,'d')]
            self.pid = top.ns*[zeros((0,self.topnpidmax),'d')]

    def clear(self):
        self.setuparrays(top.ns)

    def getid(self,safe=0):
        'If safe, then return None if id is not found rather than raising error'
        assert self.enabled,"This window is disabled and there is no associated id"
        for i in range(self.topnwin):
            if self.topizwin[i] == self.iz and self.iz >= 0: return i
            if self.topzzwin[i] == self.zz and self.iz == -1: return i
        if not safe:
            raise "Uh Ooh! Somehow the window was deleted! I can't continue! "+self.titleright()
        else:
            return None

    def setupid(self):
        self.topnwin = self.topnwin + 1
        if self.topnmax < self.nmax: self.topnmax = self.nmax
        err = gchange(self.topgroupname)
        self.topizwin[-1] = self.iz
        self.topzzwin[-1] = self.zz

    def updatenpidmax(self):
        # --- First check if npid has changed since self.topnpidmax was last set
        if self.topnpidmax != top.npid:
            self.topnpidmax = top.npid
            gchange(self.topgroupname)
        # --- Now make sure that the data arrays are changed appropriately.
        # --- This is here in case self.topnpidmax had been update elsewhere
        # --- but self.pid had not.
        if self.laccumulate and not self.dumptofile:
            for js in range(self.getns()):
                if self.pid[js].unitshape()[0] != self.topnpidmax:
                    self.pid[js].reshape((self.topnpidmax,))

    def enable(self):
        """
Enable this diagnostic. This happens automatically when the diagnostic is
created. This can be used to turn the diagnostic back on after it had been
disabled. If the data is being accumulated, any new data will be added after
the existing data.
        """
        if self.enabled: return
        self.setupid()
        # --- Set so accumulate method is called after time steps
        installafterstep(self.accumulate)
        self.enabled = 1

    def disable(self):
        """
Disable the diagnostic. Particle data will no longer be gathered nor
accumulated. If the data is being accumulated, any existing data is preserved.
        """
        if not self.enabled: return
        # --- Set so accumulate method is not called after time steps
        uninstallafterstep(self.accumulate)
        # --- Check if setupwz is installed, and if so, uninstall it
        if isinstalledafterfs(self.setupwz): uninstallafterfs(self.setupwz)
        # --- Remove this window from the list. Turn safe on when getting
        # --- the id, since it may for some reason not be consistent.
        id = self.getid(safe=1)
        if id is not None:
            for i in range(id,self.topnwin-1):
                self.topizwin[i] = self.topizwin[i+1]
                self.topzzwin[i] = self.topzzwin[i+1]
                self.topwzwin[i] = self.topwzwin[i+1]
            self.topnwin = self.topnwin - 1
            gchange(self.topgroupname)
        self.enabled = 0

    def __setstate__(self,dict):
        self.__dict__.update(dict)
        if not self.enabled: return
        id = self.getid(safe=1)
        if id is None: self.setupid()
        self.restoredata()
        if not isinstalledafterstep(self.accumulate):
            installafterstep(self.accumulate)

    def setaccumulate(self,v=1):
        self.laccumulate = v
        if self.laccumulate: self.setuparrays(top.ns)

    def accumulate(self):
        # --- If self.topnwin is 0 then something is really wrong - this routine
        # --- should never be called if self.topnwin is zero.
        if self.topnwin == 0: return
        # --- If the data is currently restored from files, which assumes that
        # --- dumptofiles is true, then reset the arrays so that they are
        # --- ready to accumulate more data and are prepared for writing
        # --- that data out to the files.
        if self.restored:
          self.setuparrays(self.getns())
        # --- Check if the number of species has changed. This is done to ensure
        # --- crashes don't happen.
        if top.ns > self.getns():
            self.addspecies()
        # --- If this windows is outside of the grid, then just return.
        if (self.iz == -1 and
            (self.zz+self.wz < w3d.zmmin+top.zbeam or
             self.zz-self.wz > w3d.zmmax+top.zbeam)):
            self.autodump()
            return
        id = self.getid()
        self.updatenpidmax()
        # --- Loop over species, collecting only ones where some particles
        # --- were saved.
        for js in range(top.ns):
            # --- Gather the data.
            # --- In serial, the arrays are just returned as is.
            # --- In parallel, if dumptofile is not being done, then the data
            # --- is gathered onto PE0 and empty arrays are returned on other
            # --- processors. If dumptofile is on, then the data is kept local
            # --- and each processor writes out its own file.

            # --- First, get the local data.
            nn = self.topn[id,js]
            t = self.topt[:nn,id,js]
            x = self.topx[:nn,id,js]
            y = self.topy[:nn,id,js]
            ux = self.topux[:nn,id,js]
            uy = self.topuy[:nn,id,js]
            uz = self.topuz[:nn,id,js]
            if self.topnpidmax > 0:
                pid = self.toppid[:nn,:,id,js]
            else:
                pid = zeros((nn,0),'d')

            if not self.dumptofile:
                # --- Gather the data onto PE0.
                ntot = globalsum(nn)
                if ntot == 0: continue

                t = gatherarray(t,othersempty=1)
                x = gatherarray(x,othersempty=1)
                y = gatherarray(y,othersempty=1)
                ux = gatherarray(ux,othersempty=1)
                uy = gatherarray(uy,othersempty=1)
                uz = gatherarray(uz,othersempty=1)
                if self.topnpidmax > 0:
                    pid = gatherarray(pid,othersempty=1)
                else:
                    pid = zeros((ntot,0),'d')

            if self.laccumulate and not self.dumptofile:

                # --- Only PE0 accumulates the data.
                if me == 0:
                    self.t[js].append(t.copy())
                    self.x[js].append(x.copy())
                    self.y[js].append(y.copy())
                    self.ux[js].append(ux.copy())
                    self.uy[js].append(uy.copy())
                    self.uz[js].append(uz.copy())
                    self.pid[js].append(pid.copy())

            else:

                self.t[js] = t
                self.x[js] = x
                self.y[js] = y
                self.ux[js] = ux
                self.uy[js] = uy
                self.uz[js] = uz
                self.pid[js] = pid

        if self.dumptofile: self.dodumptofile()
        # --- Force n to zero to ensure that particles are not saved twice.
        self.topn[id,:] = 0

    ############################################################################
    def autodump(self):
        if not self.lautodump or self.name is None: return
        if not self.laccumulate and not self.dumptofile: return
        if self.iz >= 0: return
        if self.zz+self.wz > w3d.zmmin+top.zbeam: return
        # --- Check if there is any data. If there is none, then don't make
        # --- a dump.
        ntot = 0
        for js in range(self.getns()):
            ntot = ntot + self.getn(js=js)
        if ntot > 0:
            ff = None
            #try:
            #    ff = PWpyt.PW(self.name+'_%sdump.pyt'%self.type)
            #    dumpsmode = 1
            #except:
            ff = PW.PW(self.name+'_%05d_%05d_%sdump.pdb'%(me,npes,self.type))
            dumpsmode = 0
            if ff is None:
                 print "%s: %s unable to dump data to file."%(self.topgroupname,self.name)
                 return
            ff.write(self.name+'@pickle',cPickle.dumps(self,dumpsmode))
            ff.close()
        self.nmax = 1
        self.clear()
        # --- Disable is done last so that the object written out to the
        # --- file is still enabled. That flag is used in restoredata to
        # --- determine whether or not to restore the data. The logic is set
        # --- so that the object in an autodump file will restore the data
        # --- but one is a generic dump file won't (unless it was not auto
        # --- dumped, in which case the object in the generic dump is the only
        # --- copy). There will of course be exceptions, so restoredata takes
        # --- and option argument to force restoration of data, and the
        # --- extforcenorestore function turns any restores off.
        self.disable()

    ############################################################################
    def dodumptofile(self):
        #self.dodumptofilePDB()
        self.dodumptofilePickle()

    def dodumptofilePDB(self):
        #if me != 0: return
        ff = None
        for js in range(self.getns()):
            if len(self.t[js][:]) > 0:
                if ff is None:
                    # --- Only create the file if there is data to write out.
                    ff = PW.PW(self.name+'_%s_%05d_%05d.pdb'%(self.type,me,npes),'a',verbose=0)
                if ff is None:
                    print "%s: %s unable to dump data to file."%(self.topgroupname,self.name)
                    return
                suffix = "_%d_%d"%(top.it,js)
                ff.write('n'+suffix,len(self.t[js][:]))
                ff.write('t'+suffix,self.t[js][:])
                ff.write('x'+suffix,self.x[js][:])
                ff.write('y'+suffix,self.y[js][:])
                ff.write('ux'+suffix,self.ux[js][:])
                ff.write('uy'+suffix,self.uy[js][:])
                ff.write('uz'+suffix,self.uz[js][:])
                ff.write('pid'+suffix,self.pid[js][...])
        if ff is not None:
            ff.close()

    def dodumptofilePickle(self):
        #if me != 0: return
        ff = None
        for js in range(self.getns()):
            if len(self.t[js][:]) > 0:
                if ff is None:
                    # --- Only create the file if there is data to write out.
                    if npes > 1:
                        ff = open(self.name+'_%s_%05d_%05d.pkl'%(self.type,me,npes),'a')
                    else:
                        ff = open(self.name+'_%s.pkl'%self.type,'a')
                if ff is None:
                    print "%s: %s unable to dump data to file."%(self.topgroupname,self.name)
                    return
                suffix = "_%d_%d"%(top.it,js)
                cPickle.dump(('n'+suffix,len(self.t[js][:])),ff,-1)
                cPickle.dump(('t'+suffix,self.t[js][:]),ff,-1)
                cPickle.dump(('x'+suffix,self.x[js][:]),ff,-1)
                cPickle.dump(('y'+suffix,self.y[js][:]),ff,-1)
                cPickle.dump(('ux'+suffix,self.ux[js][:]),ff,-1)
                cPickle.dump(('uy'+suffix,self.uy[js][:]),ff,-1)
                cPickle.dump(('uz'+suffix,self.uz[js][:]),ff,-1)
                cPickle.dump(('pid'+suffix,self.pid[js][...]),ff,-1)
        if ff is not None:
            ff.close()

    ############################################################################
    def restoredata(self,lforce=0,files=None,nprocs=None):
        """
Restores the data that was written out to a file. This is used when doing
post processing of the saved data when the flag dumptofile was turned on.  
        """
        if self.restored: return
        #self.restoredataPDB(0,files)
        self.restoredataPickle(0,files,nprocs=nprocs)
        self.restored = True

    def restoredataPickle(self,lforce=0,files=None,names=None,nprocs=None):
        """
Restores data dumped to a file. Note that this turns off the dumptofile
feature.
  - lforce=0: if true, force a restore, despite the value of enabled.
        """
        if files is None: files = []
        if names is None: names = []
        if not self.dumptofile: return
        if not lforce and (not self.enabled or _extforcenorestore): return

        # --- Check some input, converting to lists if needed.
        if not isinstance(files,ListType):
            files = [files]
        if not isinstance(names,ListType):
            names = [names]

        # --- If the list of files was not given, then it needs to be
        # --- generated.
        fnametries = []
        if len(files) == 0:
            # --- If a name was not given, use the instance's name
            if len(names) == 0:
                names = [self.name]
            for name in names:
                if nprocs is None and npes <= 1:
                    # --- Use the serial naming
                    files = [self.name+'_%s.pkl'%self.type]
                else:
                    if npes > 1:
                        # --- If currently running in parallel, only read in
                        # --- the data for this processor.
                        nplist = [me]
                        nprocs = npes
                    else:
                        # --- Read the data in from all processors.
                        nplist = range(nprocs)
                    for iproc in nplist:
                        fname = self.name+'_%s_%05d_%05d.pkl'%(self.type,iproc,nprocs)
                        fnametries.append(fname)
                        if os.path.exists(fname):
                            files.append(fname)

        if len(files) == 0:
            print "%s restoredata: warning, no files were found, nothing will be restored"%self.topgroupname
            if len(fnametries) > 0:
                print "Tried the filenames:",fnametries

        #datadict = self.getPDBdatadict(files)
        datadict = self.getPickledatadict(files)

        # --- Get total number of particles
        # --- jsmax is calculated since this object could be restored
        # --- into a python session that did not have all of the species
        # --- declared.
        ntot = []
        jsmax = top.ns - 1
        for var,val in datadict.items():
            if var[0] == 'n':
                ss = var.split('_')
                jsmax = max(jsmax,int(ss[2]))
                while jsmax >= len(ntot): ntot.append(0)
                ntot[jsmax] = ntot[jsmax] + val

        if npes > 1:
            # --- Make sure that all processors have the same number of
            # --- species.
            jsmax = globalmax(jsmax)

        # --- Get the size of the pid array. If self.topnpidmax is zero, that
        # --- probably means that this data is being read in for
        # --- postprocessing, since the pid arrays was not setup.
        # --- If self.topnpidmax is nonzero and different than the size
        # --- of pid, raise an exception since there is likely something
        # --- wrong.
        for var,val in datadict.items():
            if var[0:3] == 'pid':
                try:
                    npid = val.shape[1]
                except IndexError:
                    # --- This is needed so that old datasets can be read in.
                    # --- Old datasets were not written properly, leaving
                    # --- pid a 1-D array.
                    npid = 1
                if self.topnpidmax == 0:
                    self.topnpidmax = npid
                    top.npid = npid
                else:
                    assert self.topnpidmax == npid,\
                           'npid is different than in the run where the %s data was saved'%self.topgroupname
                break

        if len(ntot) == 0:
            # --- No data was read in, so this value won't be used.
            # --- This can happen in a parallel run since not all processors
            # --- will necessarily write out data.
            # --- This case is needed since max(ntot) would give an error.
            bump = 1
        else:
            # --- Set the bump size so that the arrays will be made just
            # --- large enough to hold the data.
            bump = max(ntot) + 1

        # --- Temporarily set dumptofile=0 so that setuparrays will make
        # --- the arrays appendable, and ready to take the data to be
        # --- read in.
        save_dumptofile = self.dumptofile
        save_laccumulate = self.laccumulate
        self.dumptofile = 0
        self.laccumulate = 1
        self.setuparrays(jsmax+1,bump=bump)
        self.dumptofile = save_dumptofile
        self.laccumulate = save_laccumulate

        # --- This loop must be ordered because of the append
        varlist = datadict.keys()
        varlist.sort()
        for var in varlist:
            if var[0] == 'n':
                ss = var.split('_')
                js = int(ss[2])
                suffix = var[1:]
                self.t[js].append(datadict['t%s'%suffix])
                self.x[js].append(datadict['x%s'%suffix])
                self.y[js].append(datadict['y%s'%suffix])
                self.ux[js].append(datadict['ux%s'%suffix])
                self.uy[js].append(datadict['uy%s'%suffix])
                self.uz[js].append(datadict['uz%s'%suffix])
                pid = datadict['pid%s'%suffix]
                if len(pid.shape) == 1:
                    pid.shape = (pid.shape[0],1)
                self.pid[js].append(pid)

    def getPickledatadict(self,files):
        # --- Read in all of the data into a dictionary.
        # --- The file name is appended to the end of the variable
        # --- name to avoid duplicate names. For example, when reading
        # --- in the data from a parallel run, on a time step, data may have
        # --- been written from multiple processors, duplicating the names
        # --- in the files from each of the processors.
        datadict = {}
        for file in files:
            ff = open(file,'r')
            while 1:
                try:
                    data = cPickle.load(ff)
                except:
                    break
                datadict[data[0]+'_'+file] = data[1]
            ff.close()
        return datadict

    def getPDBdatadict(self,files):
        # --- Read in all of the data into a dictionary.
        datadict = {}
        for file in files:
            try:
                ff = PRpyt.PR(file,verbose=0)
            except:
                ff = PR.PR(file,verbose=0)
            for var in ff.inquire_names():
                datadict[var] = ff.read(var)
            ff.close()
        return datadict

    def restoredataPDBold(self,lforce=0,files=[]):
        """
This is kept around for compatibility with old data files.
Restores data dumped to a file. Note that this turns off the dumptofile
feature.
  - lforce=0: if true, force a restore, despite the value of enabled.
        """
        if not self.dumptofile: return
        if not lforce and (not self.enabled or _extforcenorestore): return
        self.dumptofile = 0
        self.laccumulate = 1
        try:
            ff = PRpyt.PR(self.name+'_%s.pyt'%self.type,'a',verbose=0)
        except:
            ff = PR.PR(self.name+'_%s.pdb'%self.type,'a',verbose=0)
        # --- Get total number of particles
        ntot = []
        jsmax = 0
        varlist = list(ff.inquire_names())
        varlist.sort()
        for var in varlist:
            if var[0] == 'n':
                name,ii,js = var.split('_')
                jsmax = max(jsmax,eval(js))
                while jsmax >= len(ntot): ntot.append(0)
                ntot[jsmax] = ntot[jsmax] + ff.read(var)
        self.setuparrays(jsmax+1,bump=max(array(ntot))+1)
        for var in varlist:
            if var[0] == 'n':
                name,iis,jss = var.split('_')
                nn = ff.read(var)
                ii = eval(iis)
                js = eval(jss)
                self.t[js].append(ff.read('t_%d_%d'%(ii,js)))
                self.x[js].append(ff.read('x_%d_%d'%(ii,js)))
                self.y[js].append(ff.read('y_%d_%d'%(ii,js)))
                self.ux[js].append(ff.read('ux_%d_%d'%(ii,js)))
                self.uy[js].append(ff.read('uy_%d_%d'%(ii,js)))
                self.uz[js].append(ff.read('uz_%d_%d'%(ii,js)))
                self.pid[js].append(ff.read('pid_%d_%d'%(ii,js)))
        ff.close()

    ############################################################################
    def selectparticles(self,val,js=0,tc=None,wt=None,tp=None,z=None,v=None,
                        gather=1,bcast=1):
        """
This method shouldn't be called directly, but only through one of the "get"
methods.
 - val: The data to be selected from.
 - js=0: Species number to gather from.
         If js=='all', then the quantity is gathered from all species into
         a single array.
 - tc=None: time at which to gather particles from. When not given, returns
            all particles.
 - wt=top.dt: Width of region around tc from which to select particles.
 - tp=gett(js=js): Time value to use for the particle selection
 - z=None: when specified, projects the data to the given z location
 - gather=1: in parallel, when true, the particles are all gathered onto
             one processor
 - bcast=1: in parallel, when true, the gathered particles are broadcast
            to all processors
        """
        if js == 'all':
            nn = sum(map(len,val))
            rr = AppendableArray(initlen=nn,typecode='d')
            for js in range(len(val)):
                rr.append(self.selectparticles(val,js,tc,wt,tp,gather=0))
            result = rr[...]
        elif tc is None:
            result = val[js][...]
        else:
            if wt is None: wt = self.dt
            if tp is None: tp = self.t[js][:]
            ii = compress((tc-wt<tp)&(tp<tc+wt),arange(len(tp)))
            result = take(val[js][...],ii)

        if z is not None and v is not None:
            # --- Project to the new z value given the current velocity
            try:
                # --- v is either one of the 'get' methods
                v = v(js,tc,wt,tp)
            except TypeError:
                # --- or a constant
                pass
            zz = self.getzz()
            vz = self.getvz(js,tc,wt,tp)
            delt = (z - zz)/vz
            result = result + v*delt

        if lparallel and gather: result = gatherarray(result,bcast=bcast)
        return result

    def getns(self):
        """Get the number of species for which data is saved."""
        return len(self.t)
    def gett(self,js=0,tc=None,wt=None,tp=None,z=None,gather=1,bcast=1):
        """Get the time that each particle was saved. The same options as for :py:func:`selectparticles` apply."""
        return self.selectparticles(self.t,js,tc,wt,tp,z,1.,
                                    gather=gather,bcast=bcast)
    def getx(self,js=0,tc=None,wt=None,tp=None,z=None,gather=1,bcast=1):
        """Get the x position of the saved particles. The same options as for :py:func:`selectparticles` apply."""
        return self.selectparticles(self.x,js,tc,wt,tp,z,self.getux,
                                    gather=gather,bcast=bcast)
    def gety(self,js=0,tc=None,wt=None,tp=None,z=None,gather=1,bcast=1):
        """Get the y position of the saved particles. The same options as for :py:func:`selectparticles` apply."""
        return self.selectparticles(self.y,js,tc,wt,tp,z,self.getuy,
                                    gather=gather,bcast=bcast)
    def getux(self,js=0,tc=None,wt=None,tp=None,z=None,gather=1,bcast=1):
        """Get the x massless momentum of the saved particles. The same options as for :py:func:`selectparticles` apply."""
        return self.selectparticles(self.ux,js,tc,wt,tp,
                                    gather=gather,bcast=bcast)
    def getuy(self,js=0,tc=None,wt=None,tp=None,z=None,gather=1,bcast=1):
        """Get the y massless momentum of the saved particles. The same options as for :py:func:`selectparticles` apply."""
        return self.selectparticles(self.uy,js,tc,wt,tp,
                                    gather=gather,bcast=bcast)
    def getuz(self,js=0,tc=None,wt=None,tp=None,z=None,gather=1,bcast=1):
        """Get the z massless momentum of the saved particles. The same options as for :py:func:`selectparticles` apply."""
        return self.selectparticles(self.uz,js,tc,wt,tp,
                                    gather=gather,bcast=bcast)
    def getvx(self,js=0,tc=None,wt=None,tp=None,z=None,gather=1,bcast=1):
        """Get the x velocity of the saved particles. The same options as for :py:func:`selectparticles` apply."""
        return self.selectparticles(self.ux,js,tc,wt,tp,
                                    gather=gather,bcast=bcast)
    def getvy(self,js=0,tc=None,wt=None,tp=None,z=None,gather=1,bcast=1):
        """Get the y velocity of the saved particles. The same options as for :py:func:`selectparticles` apply."""
        return self.selectparticles(self.uy,js,tc,wt,tp,
                                    gather=gather,bcast=bcast)
    def getvz(self,js=0,tc=None,wt=None,tp=None,z=None,gather=1,bcast=1):
        """Get the z velocity of the saved particles. The same options as for :py:func:`selectparticles` apply."""
        return self.selectparticles(self.uz,js,tc,wt,tp,
                                    gather=gather,bcast=bcast)
    def getpid(self,js=0,tc=None,wt=None,tp=None,z=None,id=0,gather=1,bcast=1):
        """Get the pid of the saved particles.
  - id=0: which pid to return

The same options as for :py:func:`selectparticles` apply."""
        self.updatenpidmax()
        if self.topnpidmax > 0:
            return self.selectparticles(self.pid,js,tc,wt,tp,
                                        gather=gather,bcast=bcast)[:,id]
        else:
            return zeros((self.getn(js,tc,wt,tp,
                                    gather=gather,bcast=bcast),0),'d')

    def getxp(self,js=0,tc=None,wt=None,tp=None,z=None,gather=1,bcast=1):
        """Get the x' of the saved particles. The same options as for :py:func:`selectparticles` apply."""
        return (self.getux(js,tc,wt,tp,gather=gather,bcast=bcast)/
                self.getuz(js,tc,wt,tp,gather=gather,bcast=bcast))
    def getyp(self,js=0,tc=None,wt=None,tp=None,z=None,gather=1,bcast=1):
        """Get the y' of the saved particles. The same options as for :py:func:`selectparticles` apply."""
        return (self.getuy(js,tc,wt,tp,gather=gather,bcast=bcast)/
                self.getuz(js,tc,wt,tp,gather=gather,bcast=bcast))
    def getr(self,js=0,tc=None,wt=None,tp=None,z=None,gather=1,bcast=1):
        """Get the radial position of the saved particles. The same options as for :py:func:`selectparticles` apply."""
        return sqrt(self.getx(js,tc,wt,tp,z,gather=gather,bcast=bcast)**2 +
                    self.gety(js,tc,wt,tp,z,gather=gather,bcast=bcast)**2)
    def gettheta(self,js=0,tc=None,wt=None,tp=None,z=None,gather=1,bcast=1):
        """Get the angular position of the saved particles. The same options as for :py:func:`selectparticles` apply."""
        return arctan2(self.gety(js,tc,wt,tp,z,gather=gather,bcast=bcast),
                       self.getx(js,tc,wt,tp,z,gather=gather,bcast=bcast))
    def getvr(self,js=0,tc=None,wt=None,tp=None,z=None,gather=1,bcast=1):
        """Get the radial velocity of the saved particles. The same options as for :py:func:`selectparticles` apply."""
        theta = self.gettheta(js,tc,wt,tp,gather=gather,bcast=bcast)
        vx = self.getvx(js,tc,wt,tp,gather=gather,bcast=bcast)
        vy = self.getvy(js,tc,wt,tp,gather=gather,bcast=bcast)
        return (vx*cos(theta) + vy*sin(theta))
    def getrp(self,js=0,tc=None,wt=None,tp=None,z=None,gather=1,bcast=1):
        """Get the r' of the saved particles. The same options as for :py:func:`selectparticles` apply."""
        theta = self.gettheta(js,tc,wt,tp,gather=gather,bcast=bcast)
        xp = self.getxp(js,tc,wt,tp,gather=gather,bcast=bcast)
        yp = self.getyp(js,tc,wt,tp,gather=gather,bcast=bcast)
        return (xp*cos(theta) + yp*sin(theta))
    def getn(self,js=0,tc=None,wt=None,tp=None,z=None,gather=1,bcast=1):
        """Get the number of saved particles. The same options as for :py:func:`selectparticles` apply."""
        return len(self.gett(js,tc,wt,tp,gather=gather,bcast=bcast))

    def xxpslope(self,js=0,tc=None,wt=None,tp=None,z=None):
        """Get the x-x' slope of the saved particles. The same options as for :py:func:`selectparticles` apply."""
        if self.getn(js,tc,wt,tp) == 0:
            return 0.
        else:
            return (
                (ave(self.getx(js,tc,wt,tp,z)*self.getxp(js,tc,wt,tp)) -
                 ave(self.getx(js,tc,wt,tp,z))*ave(self.getxp(js,tc,wt,tp)))/
                (ave(self.getx(js,tc,wt,tp,z)*self.getx(js,tc,wt,tp,z)) -
                 ave(self.getx(js,tc,wt,tp,z))*ave(self.getx(js,tc,wt,tp,z))))
    def yypslope(self,js=0,tc=None,wt=None,tp=None,z=None):
        """Get the y-y' slope of the saved particles. The same options as for :py:func:`selectparticles` apply."""
        if self.getn(js,tc,wt,tp) == 0:
            return 0.
        else:
            return (
                (ave(self.gety(js,tc,wt,tp,z)*self.getyp(js,tc,wt,tp)) -
                 ave(self.gety(js,tc,wt,tp,z))*ave(self.getyp(js,tc,wt,tp)))/
                (ave(self.gety(js,tc,wt,tp,z)*self.gety(js,tc,wt,tp,z)) -
                 ave(self.gety(js,tc,wt,tp,z))*ave(self.gety(js,tc,wt,tp,z))))
    def rrpslope(self,js=0,tc=None,wt=None,tp=None,z=None):
        """Get the r-r' slope of the saved particles. The same options as for :py:func:`selectparticles` apply."""
        if self.getn(js,tc,wt,tp) == 0:
            return 0.
        else:
            return (ave(self.getr(js,tc,wt,tp,z)*self.getrp(js,tc,wt,tp))/
                    ave(self.getr(js,tc,wt,tp,z)**2))

    #getx.__doc__ += selectparticles.__doc__
    #getns.__doc__ += selectparticles.__doc__
    #gett.__doc__ += selectparticles.__doc__
    #getx.__doc__ += selectparticles.__doc__
    #gety.__doc__ += selectparticles.__doc__
    #getux.__doc__ += selectparticles.__doc__
    #getuy.__doc__ += selectparticles.__doc__
    #getuz.__doc__ += selectparticles.__doc__
    #getvx.__doc__ += selectparticles.__doc__
    #getvy.__doc__ += selectparticles.__doc__
    #getvz.__doc__ += selectparticles.__doc__
    #getxp.__doc__ += selectparticles.__doc__
    #getyp.__doc__ += selectparticles.__doc__
    #getr.__doc__ += selectparticles.__doc__
    #gettheta.__doc__ += selectparticles.__doc__
    #getvr.__doc__ += selectparticles.__doc__
    #getrp.__doc__ += selectparticles.__doc__
    #getn.__doc__ += selectparticles.__doc__
    #xxpslope.__doc__ += selectparticles.__doc__
    #yypslope.__doc__ += selectparticles.__doc__
    #rrpslope.__doc__ += selectparticles.__doc__
    #getpid.__doc__ += selectparticles.__doc__[:-4]+' - id=0: which pid to return\n'

    ############################################################################
    ############################################################################
    # --- Define plotting routines for the extrapolated particles.

    def checkplotargs(self,kw):
        """Convenience routine to check arguments of particle plot routines.
Warning: this has the side affect of adding the arguement allowbadargs to
the kw dictionary. This is done since the calls to these functions here to
make the plots may have unused arguements since the entire kw list passed
into each of the pp plotting routines is passed into each of these
functions.
        """
        badargs = ppgeneric(checkargs=1,kwdict=kw)
        kw['allowbadargs'] = 1
        if badargs: raise 'bad arguments ',' '.join(badargs.keys())

    def titleright(self,tc=None,wt=None,z=None,slope=None):
        if tc is None:
            ttext = ''
        else:
            if wt is None: wt = self.dt
            ttext = '  time = %e ^+_-%e'%(tc,wt)
        if self.iz >= 0:
            ztext =  'iz = %d (z = %f m)'%(self.iz,w3d.zmmin+self.iz*w3d.dz)
        else:
            ztext =  'z = %f m'%self.zz
        if z is not None:
            ztext = ztext + ' projected to z = %f m'%z
        if slope is None:
            slopetext = ''
        else:
            slopetext = '  slope=%7.4f'%slope
        return ztext + ttext + slopetext

    def ppmultispecies(self,pp,args,kw):
        """Check if js is -1 or a list, in which case call the pp function for
each species and each one in the list. Also assign colors accordingly
        """
        args = list(args)
        js = args[0]
        if js != -1 and type(js) != ListType:
            return false
        else:
            if js == -1: js = range(self.getns())
            ncolor = kw.get('ncolor',240)
            color = kw.get('color',range(0,ncolor,ncolor/len(js)))
            for i in xrange(len(js)):
                args[0] = js[i]
                kw['color'] = color[i]
                apply(pp,args,kw)
            return true

    ############################################################################
    def pxy(self,js=0,tc=None,wt=None,tp=None,z=None,**kw):
        """
Plots X-Y for extraploated particles.
The same arguments for :py:func:`selectparticles` and :py:func:`~warpplots.ppgeneric` apply.
        """
        self.checkplotargs(kw)
        if self.ppmultispecies(self.pxy,(js,tc,wt,tp,z),kw): return
        x = self.getx(js,tc,wt,tp,z)
        y = self.gety(js,tc,wt,tp,z)
        if 'pplimits' in kw.keys():
            kw['lframe'] = 1
        else:
            kw['pplimits'] = (top.xplmin,top.xplmax,top.yplmin,top.yplmax)
        settitles("Y vs X","X","Y",self.titleright(tc,wt,z))
        return ppgeneric(y,x,kwdict=kw)

    ############################################################################
    def pxvx(self,js=0,tc=None,wt=None,tp=None,z=None,slope=0.,offset=0.,
             **kw):
        """
Plots X-Vx for extraploated particles
 - slope=0.: slope subtracted from vx, it is calculated automatically
 - offset=0.: offset in x

The same arguments for :py:func:`selectparticles` and :py:func:`~warpplots.ppgeneric` apply.
        """
        self.checkplotargs(kw)
        if self.ppmultispecies(self.pxvx,(js,tc,wt,tp,z),kw): return
        x = self.getx(js,tc,wt,tp,z)
        vx = self.getvx(js,tc,wt,tp)
        if type(slope) == type(''):
            if len(x) > 0:
                slope = (ave(x*vx)-ave(x)*ave(vx))/(ave(x*x) - ave(x)**2)
                offset = ave(vx)-slope*ave(x)
            else:
                slope = 0.
                offset = 0.
        kw['slope'] = slope
        kw['offset'] = offset
        if 'pplimits' in kw.keys():
            kw['lframe'] = 1
        settitles("Vx vs X","X","Vx",self.titleright(tc,wt,z,slope))
        return ppgeneric(vx,x,kwdict=kw)

    ############################################################################
    def pyvy(self,js=0,tc=None,wt=None,tp=None,z=None,slope=0.,offset=0.,
             **kw):
        """
Plots Y-Vy for extraploated particles
 - slope=0.: slope subtracted from vy, it is calculated automatically
 - offset=0.: offset in y

The same arguments for :py:func:`selectparticles` and :py:func:`~warpplots.ppgeneric` apply.
        """
        self.checkplotargs(kw)
        if self.ppmultispecies(self.pyvy,(js,tc,wt,tp,z,slope,offset),kw):
            return
        y = self.gety(js,tc,wt,tp,z)
        vy = self.getvy(js,tc,wt,tp)
        if type(slope) == type(''):
            if len(y) > 0:
                slope = (ave(y*vy)-ave(y)*ave(vy))/(ave(y*y) - ave(y)**2)
                offset = ave(vy)-slope*ave(y)
            else:
                slope = 0.
                offset = 0.
        kw['slope'] = slope
        kw['offset'] = offset
        if 'pplimits' in kw.keys():
            kw['lframe'] = 1
        settitles("Vy vs Y","Y","Vy",self.titleright(tc,wt,z,slope))
        return ppgeneric(vy,y,kwdict=kw)

    ############################################################################
    def pvxvy(self,js=0,tc=None,wt=None,tp=None,z=None,**kw):
        """
Plots Vx-Vy for extraploated particles.
The same arguments for :py:func:`selectparticles` and :py:func:`~warpplots.ppgeneric` apply.
        """
        self.checkplotargs(kw)
        if self.ppmultispecies(self.pvxvy,(js,tc,wt,tp),kw): return
        vx = self.getvx(js,tc,wt,tp)
        vy = self.getvy(js,tc,wt,tp)
        if 'pplimits' in kw.keys():
            kw['lframe'] = 1
        settitles("Vy vs Vx","Vx","Vy",self.titleright(tc,wt,z))
        return ppgeneric(vy,vx,kwdict=kw)

    ############################################################################
    def prvr(self,js=0,tc=None,wt=None,tp=None,z=None,scale=0.,slope=0.,offset=0.,
             **kw):
        """
Plots R-Vr for extraploated particles
 - slope=0.: slope subtracted from the vr, it is calculated automatically
 - offset=0.: offset in r

The same arguments for :py:func:`selectparticles` and :py:func:`~warpplots.ppgeneric` apply.
        """
        self.checkplotargs(kw)
        if self.ppmultispecies(self.prvr,(js,tc,wt,tp,z,scale,slope,offset),kw):
            return
        x = self.getx(js,tc,wt,tp,z)
        y = self.gety(js,tc,wt,tp,z)
        vx = self.getvx(js,tc,wt,tp)
        vy = self.getvy(js,tc,wt,tp)
        xscale = 1.
        yscale = 1.
        vxscale = 1.
        vyscale = 1.
        if scale:
            xscale = 2.*sqrt(ave(x*x) - ave(x)**2)
            yscale = 2.*sqrt(ave(y*y) - ave(y)**2)
            vxscale = 2.*sqrt(ave(vx*vx) - ave(vx)**2)
            vyscale = 2.*sqrt(ave(vy*vy) - ave(vy)**2)
        x = x/xscale
        y = y/yscale
        vx = vx/vxscale
        vy = vy/vyscale
        r = sqrt(x**2 + y**2)
        t = arctan2(y,x)
        vr = vx*cos(t) + vy*sin(t)
        if type(slope) == type(''):
            if len(r) > 0:
              slope = ave(r*vr)/ave(r*r)
            else:
                slope = 0.
        kw['slope'] = slope
        if 'pplimits' in kw.keys():
            kw['lframe'] = 1
        settitles("Vr vs R","R","Vr",self.titleright(tc,wt,z,slope))
        return ppgeneric(vr,r,kwdict=kw)

    ############################################################################
    def pxxp(self,js=0,tc=None,wt=None,tp=None,z=None,slope=0.,offset=0.,
             **kw):
        """
Plots X-X' for extraploated particles
 - slope=0.: slope subtracted from xp, it is calculated automatically
 - offset=0.: offset in x

The same arguments for :py:func:`selectparticles` and :py:func:`~warpplots.ppgeneric` apply.
        """
        self.checkplotargs(kw)
        if self.ppmultispecies(self.pxxp,(js,tc,wt,tp,z),kw): return
        x = self.getx(js,tc,wt,tp,z)
        xp = self.getxp(js,tc,wt,tp)
        if type(slope) == type(''):
            if len(x) > 0:
                slope = (ave(x*xp)-ave(x)*ave(xp))/(ave(x*x) - ave(x)**2)
                offset = ave(xp)-slope*ave(x)
            else:
                slope = 0.
                offset = 0.
        kw['slope'] = slope
        kw['offset'] = offset
        if 'pplimits' in kw.keys():
            kw['lframe'] = 1
        else:
            kw['pplimits'] = (top.xplmin,top.xplmax,top.xpplmin,top.xpplmax)
        settitles("X' vs X","X","X'",self.titleright(tc,wt,z,slope))
        return ppgeneric(xp,x,kwdict=kw)

    ############################################################################
    def pyyp(self,js=0,tc=None,wt=None,tp=None,z=None,slope=0.,offset=0.,
             **kw):
        """
Plots Y-Y' for extraploated particles
 - slope=0.: slope subtracted from yp, it is calculated automatically
 - offset=0.: offset in y

The same arguments for :py:func:`selectparticles` and :py:func:`~warpplots.ppgeneric` apply.
        """
        self.checkplotargs(kw)
        if self.ppmultispecies(self.pyyp,(js,tc,wt,tp,z,slope,offset),kw):
            return
        y = self.gety(js,tc,wt,tp,z)
        yp = self.getyp(js,tc,wt,tp)
        if type(slope) == type(''):
            if len(y) > 0:
                slope = (ave(y*yp)-ave(y)*ave(yp))/(ave(y*y) - ave(y)**2)
                offset = ave(yp)-slope*ave(y)
            else:
                slope = 0.
                offset = 0.
        kw['slope'] = slope
        kw['offset'] = offset
        if 'pplimits' in kw.keys():
            kw['lframe'] = 1
        else:
            kw['pplimits'] = (top.yplmin,top.yplmax,top.ypplmin,top.ypplmax)
        settitles("Y' vs Y","Y","Y'",self.titleright(tc,wt,z,slope))
        return ppgeneric(yp,y,kwdict=kw)

    ############################################################################
    def pxpyp(self,js=0,tc=None,wt=None,tp=None,z=None,**kw):
        """
Plots X'-Y' for extraploated particles.
The same arguments for :py:func:`selectparticles` and :py:func:`~warpplots.ppgeneric` apply.
        """
        self.checkplotargs(kw)
        if self.ppmultispecies(self.pxpyp,(js,tc,wt,tp),kw): return
        xp = self.getxp(js,tc,wt,tp)
        yp = self.getyp(js,tc,wt,tp)
        if 'pplimits' in kw.keys():
            kw['lframe'] = 1
        else:
            kw['pplimits'] = (top.xpplmin,top.xpplmax,top.ypplmin,top.ypplmax)
        settitles("Y' vs X'","X'","Y'",self.titleright(tc,wt,z))
        return ppgeneric(yp,xp,kwdict=kw)

    ############################################################################
    def prrp(self,js=0,tc=None,wt=None,tp=None,z=None,scale=0.,slope=0.,offset=0.,
             **kw):
        """
Plots R-R' for extraploated particles
 - slope=0.: slope subtracted from the rp, it is calculated automatically
 - offset=0.: offset in r

The same arguments for :py:func:`selectparticles` and :py:func:`~warpplots.ppgeneric` apply.
        """
        self.checkplotargs(kw)
        if self.ppmultispecies(self.prrp,(js,tc,wt,tp,z,scale,slope,offset),kw):
            return
        x = self.getx(js,tc,wt,tp,z)
        y = self.gety(js,tc,wt,tp,z)
        xp = self.getxp(js,tc,wt,tp)
        yp = self.getyp(js,tc,wt,tp)
        xscale = 1.
        yscale = 1.
        xpscale = 1.
        ypscale = 1.
        if scale:
            xscale = 2.*sqrt(ave(x*x) - ave(x)**2)
            yscale = 2.*sqrt(ave(y*y) - ave(y)**2)
            xpscale = 2.*sqrt(ave(xp*xp) - ave(xp)**2)
            ypscale = 2.*sqrt(ave(yp*yp) - ave(yp)**2)
        x = x/xscale
        y = y/yscale
        xp = xp/xpscale
        yp = yp/ypscale
        r = sqrt(x**2 + y**2)
        t = arctan2(y,x)
        rp = xp*cos(t) + yp*sin(t)
        if type(slope) == type(''):
            if len(r) > 0:
              slope = ave(r*rp)/ave(r*r)
            else:
                slope = 0.
        kw['slope'] = slope
        if 'pplimits' in kw.keys():
            kw['lframe'] = 1
        else:
            kw['pplimits'] = (0.,max(top.xplmax/xscale,top.yplmax/yscale),
                              top.xpplmin/xpscale,top.xpplmax/ypscale)
        settitles("R' vs R","R","R'",self.titleright(tc,wt,z,slope))
        return ppgeneric(rp,r,kwdict=kw)

    ############################################################################
    def ptx(self,js=0,tc=None,wt=None,tp=None,z=None,**kw):
        """
Plots time-X for extraploated particles.
The same arguments for :py:func:`selectparticles` and :py:func:`~warpplots.ppgeneric` apply.
        """
        self.checkplotargs(kw)
        if self.ppmultispecies(self.ptx,(js,tc,wt,tp,z),kw): return
        t = self.gett(js,tc,wt,tp,z)
        x = self.getx(js,tc,wt,tp,z)
        if 'pplimits' in kw.keys():
            kw['lframe'] = 1
        else:
            kw['pplimits'] = ('e','e',top.xplmin,top.xplmax)
        settitles("X vs time","time","X",self.titleright(tc,wt,z))
        return ppgeneric(x,t,kwdict=kw)

    ############################################################################
    def pty(self,js=0,tc=None,wt=None,tp=None,z=None,**kw):
        """
Plots time-Y for extraploated particles.
The same arguments for :py:func:`selectparticles` and :py:func:`~warpplots.ppgeneric` apply.
        """
        self.checkplotargs(kw)
        if self.ppmultispecies(self.pty,(js,tc,wt,tp,z),kw): return
        t = self.gett(js,tc,wt,tp,z)
        y = self.gety(js,tc,wt,tp,z)
        if 'pplimits' in kw.keys():
            kw['lframe'] = 1
        else:
            kw['pplimits'] = ('e','e',top.yplmin,top.yplmax)
        settitles("Y vs time","time","Y",self.titleright(tc,wt,z))
        return ppgeneric(y,t,kwdict=kw)

    ############################################################################
    def ptxp(self,js=0,tc=None,wt=None,tp=None,z=None,**kw):
        """
Plots time-X' for extraploated particles.
The same arguments for :py:func:`selectparticles` and :py:func:`~warpplots.ppgeneric` apply.
        """
        self.checkplotargs(kw)
        if self.ppmultispecies(self.ptxp,(js,tc,wt,tp,z),kw): return
        t = self.gett(js,tc,wt,tp,z)
        xp = self.getxp(js,tc,wt,tp)
        if 'pplimits' in kw.keys():
            kw['lframe'] = 1
        else:
            kw['pplimits'] = ('e','e',top.xpplmin,top.xpplmax)
        settitles("X' vs time","time","X'",self.titleright(tc,wt,z))
        return ppgeneric(xp,t,kwdict=kw)

    ############################################################################
    def ptyp(self,js=0,tc=None,wt=None,tp=None,z=None,**kw):
        """
Plots time-Y' for extraploated particles.
The same arguments for :py:func:`selectparticles` and :py:func:`~warpplots.ppgeneric` apply.
        """
        self.checkplotargs(kw)
        if self.ppmultispecies(self.ptyp,(js,tc,wt,tp,z),kw): return
        t = self.gett(js,tc,wt,tp,z)
        yp = self.getyp(js,tc,wt,tp)
        if 'pplimits' in kw.keys():
            kw['lframe'] = 1
        else:
            kw['pplimits'] = ('e','e',top.ypplmin,top.ypplmax)
        settitles("Y' vs time","time","Y'",self.titleright(tc,wt,z))
        return ppgeneric(yp,t,kwdict=kw)

    ############################################################################
    def ptux(self,js=0,tc=None,wt=None,tp=None,z=None,**kw):
        """
Plots time-ux for extraploated particles.
The same arguments for :py:func:`selectparticles` and :py:func:`~warpplots.ppgeneric` apply.
        """
        self.checkplotargs(kw)
        if self.ppmultispecies(self.ptux,(js,tc,wt,tp,z),kw): return
        t = self.gett(js,tc,wt,tp,z)
        ux = self.getux(js,tc,wt,tp)
        if 'pplimits' in kw.keys():
            kw['lframe'] = 1
        settitles("ux vs time","time","ux",self.titleright(tc,wt,z))
        return ppgeneric(ux,t,kwdict=kw)

    ############################################################################
    def ptuy(self,js=0,tc=None,wt=None,tp=None,z=None,**kw):
        """
Plots time-uy for extraploated particles.
The same arguments for :py:func:`selectparticles` and :py:func:`~warpplots.ppgeneric` apply.
        """
        self.checkplotargs(kw)
        if self.ppmultispecies(self.ptuy,(js,tc,wt,tp,z),kw): return
        t = self.gett(js,tc,wt,tp,z)
        uy = self.getuy(js,tc,wt,tp)
        if 'pplimits' in kw.keys():
            kw['lframe'] = 1
        settitles("uy vs time","time","uy",self.titleright(tc,wt,z))
        return ppgeneric(uy,t,kwdict=kw)

    ############################################################################
    def ptuz(self,js=0,tc=None,wt=None,tp=None,z=None,**kw):
        """
Plots time-uz for extraploated particles.
The same arguments for :py:func:`selectparticles` and :py:func:`~warpplots.ppgeneric` apply.
        """
        self.checkplotargs(kw)
        if self.ppmultispecies(self.ptuz,(js,tc,wt,tp,z),kw): return
        t = self.gett(js,tc,wt,tp,z)
        uz = self.getuz(js,tc,wt,tp)
        if 'pplimits' in kw.keys():
            kw['lframe'] = 1
        settitles("uz vs time","time","uz",self.titleright(tc,wt,z))
        return ppgeneric(uz,t,kwdict=kw)

    ############################################################################
    def ptvx(self,js=0,tc=None,wt=None,tp=None,z=None,**kw):
        """
Plots time-Vx for extraploated particles.
The same arguments for :py:func:`selectparticles` and :py:func:`~warpplots.ppgeneric` apply.
        """
        self.checkplotargs(kw)
        if self.ppmultispecies(self.ptvx,(js,tc,wt,tp,z),kw): return
        t = self.gett(js,tc,wt,tp,z)
        vx = self.getvx(js,tc,wt,tp)
        if 'pplimits' in kw.keys():
            kw['lframe'] = 1
        settitles("Vx vs time","time","Vx",self.titleright(tc,wt,z))
        return ppgeneric(vx,t,kwdict=kw)

    ############################################################################
    def ptvy(self,js=0,tc=None,wt=None,tp=None,z=None,**kw):
        """
Plots time-Vy for extraploated particles.
The same arguments for :py:func:`selectparticles` and :py:func:`~warpplots.ppgeneric` apply.
        """
        self.checkplotargs(kw)
        if self.ppmultispecies(self.ptvy,(js,tc,wt,tp,z),kw): return
        t = self.gett(js,tc,wt,tp,z)
        vy = self.getvy(js,tc,wt,tp)
        if 'pplimits' in kw.keys():
            kw['lframe'] = 1
        settitles("Vy vs time","time","Vy",self.titleright(tc,wt,z))
        return ppgeneric(vy,t,kwdict=kw)

    ############################################################################
    def ptvz(self,js=0,tc=None,wt=None,tp=None,z=None,**kw):
        """
Plots time-Vz for extraploated particles.
The same arguments for :py:func:`selectparticles` and :py:func:`~warpplots.ppgeneric` apply.
        """
        self.checkplotargs(kw)
        if self.ppmultispecies(self.ptvz,(js,tc,wt,tp,z),kw): return
        t = self.gett(js,tc,wt,tp,z)
        vz = self.getvz(js,tc,wt,tp)
        if 'pplimits' in kw.keys():
            kw['lframe'] = 1
        settitles("Vz vs time","time","Vz",self.titleright(tc,wt,z))
        return ppgeneric(vz,t,kwdict=kw)

    ############################################################################
    def ptkez(self,js=0,tc=None,wt=None,tp=None,z=None,**kw):
        """
Plots time-kinetic energy for extraploated particles.
The same arguments for :py:func:`selectparticles` and :py:func:`~warpplots.ppgeneric` apply.
        """
        self.checkplotargs(kw)
        if self.ppmultispecies(self.ptkez,(js,tc,wt,tp,z),kw): return
        t = self.gett(js,tc,wt,tp,z)
        vz = self.getvz(js,tc,wt,tp)
        kez = 0.5*top.pgroup.sm[js]*vz**2/jperev
        if 'pplimits' in kw.keys():
            kw['lframe'] = 1
        settitles("KEz vs time","time","KE (volts)",self.titleright(tc,wt,z))
        return ppgeneric(kez,t,kwdict=kw)

    ############################################################################
    def ptrace(self,js=0,tc=None,wt=None,tp=None,z=None,slope=0.,
               pplimits=None,**kw):
        """
Plots X-Y, X-X', Y-Y', Y'-X' in single page
 - slope=0.: slope subtracted from the angle. If 'auto', it is calculated
             automatically for the X-X' and Y-Y' plots.
 - pplimits=None: An optional list of up to four tuples, one for each phase
                  space plot. If any of the tuples are empty, the limits used
                  will be the usual ones for that plot.

The same arguments for :py:func:`selectparticles` and :py:func:`~warpplots.ppgeneric` apply.
        """
        self.checkplotargs(kw)
        if self.ppmultispecies(self.ptrace,(js,tc,wt,tp,z,slope,pplimits),kw):
            return
        x = self.getx(js,tc,wt,tp,z)
        y = self.gety(js,tc,wt,tp,z)
        xp = self.getxp(js,tc,wt,tp)
        yp = self.getyp(js,tc,wt,tp)
        titler = self.titleright(tc,wt,z)
        defaultpplimits = [(top.xplmin,top.xplmax,top.yplmin,top.yplmax),
                           (top.yplmin,top.yplmax,top.ypplmin,top.ypplmax),
                           (top.xplmin,top.xplmax,top.xpplmin,top.xpplmax),
                           (top.ypplmin,top.ypplmax,top.xpplmin,top.xpplmax)]
        if pplimits is None:
            pplimits = defaultpplimits
        else:
            kw['lframe'] = 1
            if type(pplimits[0]) != type(()):
                pplimits = 4*[pplimits]
            else:
                for i in range(4):
                    if i == len(pplimits): pplimits.append(defaultpplimits[i])
                    if not pplimits[i]: pplimits[i] = defaultpplimits[i]

        kw['view'] = 3
        kw['pplimits'] = pplimits[0]
        if type(slope)==type(''):
            kw['slope'] = 0.
        settitles("Y vs X","X","Y",titler)
        ppgeneric(y,x,kwdict=kw)

        kw['view'] = 4
        kw['pplimits'] = pplimits[1]
        if type(slope)==type(''):
            kw['slope'] = (ave(y*yp)-ave(y)*ave(yp))/dvnz(ave(y*y) - ave(y)**2)
        settitles("Y' vs Y","Y","Y'",titler)
        ppgeneric(yp,y,kwdict=kw)

        kw['view'] = 5
        kw['pplimits'] = pplimits[2]
        if type(slope)==type(''):
            kw['slope'] = (ave(x*xp)-ave(x)*ave(xp))/dvnz(ave(x*x) - ave(x)**2)
        settitles("X' vs X","X","X'",titler)
        ppgeneric(xp,x,kwdict=kw)

        kw['view'] = 6
        kw['pplimits'] = pplimits[3]
        if type(slope)==type(''):
            kw['slope'] = 0.
        settitles("X' vs Y'","Y'","X'",titler)
        ppgeneric(xp,yp,kwdict=kw)

    #pxy.__doc__ = pxy.__doc__ + selectparticles.__doc__[:-4]+'plus all ppgeneric options\n'
    #pxxp.__doc__ = pxxp.__doc__ + selectparticles.__doc__[:-4]+'plus all ppgeneric options\n'
    #pyyp.__doc__ = pyyp.__doc__ + selectparticles.__doc__[:-4]+'plus all ppgeneric options\n'
    #pxpyp.__doc__ = pxpyp.__doc__ + selectparticles.__doc__[:-4]+'plus all ppgeneric options\n'
    #prrp.__doc__ = prrp.__doc__ + selectparticles.__doc__[:-4]+'plus all ppgeneric options\n'
    #ptx.__doc__ = ptx.__doc__ + selectparticles.__doc__[:-4]+'plus all ppgeneric options\n'
    #pty.__doc__ = pty.__doc__ + selectparticles.__doc__[:-4]+'plus all ppgeneric options\n'
    #ptxp.__doc__ = ptxp.__doc__ + selectparticles.__doc__[:-4]+'plus all ppgeneric options\n'
    #ptyp.__doc__ = ptyp.__doc__ + selectparticles.__doc__[:-4]+'plus all ppgeneric options\n'
    #ptux.__doc__ = ptux.__doc__ + selectparticles.__doc__[:-4]+'plus all ppgeneric options\n'
    #ptuy.__doc__ = ptuy.__doc__ + selectparticles.__doc__[:-4]+'plus all ppgeneric options\n'
    #ptuz.__doc__ = ptuz.__doc__ + selectparticles.__doc__[:-4]+'plus all ppgeneric options\n'
    #ptvx.__doc__ = ptvx.__doc__ + selectparticles.__doc__[:-4]+'plus all ppgeneric options\n'
    #ptvy.__doc__ = ptvy.__doc__ + selectparticles.__doc__[:-4]+'plus all ppgeneric options\n'
    #ptvz.__doc__ = ptvz.__doc__ + selectparticles.__doc__[:-4]+'plus all ppgeneric options\n'
    #ptrace.__doc__ = ptrace.__doc__ + selectparticles.__doc__[:-4]+'plus all ppgeneric options\n'

############################################################################
class ExtPart(ParticleAccumulator):
    """This class defines a container to setup and keep track of extropolated
particle data. It can optionally accumulate the data over multiple time steps.
Each time step, all particles that are within +-wz of the specified z location
are extrapolated to the z location and saved. The extrapolation is done by
first estimating the forces on the particle using a backward difference of
the velocity, the current velocity minus the previous velocity. This force is
applied to the current velocity to get the extrapolated velocity. The
extrapolated position is given by ze = zi + vi*dt + 1/2 a*dt**2, where
zi and vi are the current position and velocity, a is the estimated force
over m.

Note that because the extrapolation is inconsistent with the particle advance
used during the time steps, there can be gaps in the resulting data at time
step boundaries. It is recommended that :py:func:`ZCrossingParticles` be used
instead.  

The creator options are:

 - iz: grid location where the extrapolated data is saved.
 - zz: lab location where data is saved.
 - wz: width of lab window
 - nepmax: max size of the arrays. Defaults to 3*top.pnumz[iz] if non-zero,
           otherwise 10000.
 - laccumulate=0: when true, particles are accumulated over multiple steps.
 - lepsaveonce=0: when true, particles are saved only once, when they are
                  closest to the z location. Otherwise, extrapolated
                  information from a particle will be saved multiple times
                  if the particle is near the plane over multiple time steps.
                  Note that this is a global option - the value must be the
                  same in all instances; it sets flag top.lepsaveonce.
 - name=None: descriptive name for location. It must be unique for each
              instance of the diagnostic.
 - lautodump=0: when true, after the grid moves beyond the z location,
                automatically dump the data to a file, clear the arrays and
                disable itself. Also must have name set - name is appended to
                the name of the files where the data is stored. The object
                writes itself out to a pdb file.
 - dumptofile=0: when true, the particle data is always dumped to a file
                 and not saved in memory. Name must be set - name is appended
                 to the name of the files where the data is stored. Setting this
                 to true implies that the data is accumulated.
                 Use the restoredata method to read the data back in.

One of iz or zz must be specified.

Available methods:

 - setaccumulate(v=1): Turns accumulation on or off. If turned on, all
                       currently saved data is deleted.
 - clear(): Will clear the existing data.
 - disable(): Turns off collecting of data
 - enable(): Turns on collecting of data (only needed after disable)

The follow all take an optional argument to specify species number.

 - getn: Get number of particles
 - gett: Get time at which particle was saved
 - getx, y, ux, uy, uz, vx, vy, vz, xp, yp, r, theta, rp: Get the various
     coordinates or velocity of particles
 - getpid: Get the particle ID info

The following are available plot routines. All take an optional argument to
specify species number. Additional arguments are the same as the 'pp' plotting
routines (such as ppxxp).

 - pxy, pxxp, pyyp, pxpyp, prrp, ptx, pty, ptxp, ptyp, ptux, ptuy, ptuz, ptvx
 - ptvy, ptvz, ptrace
    """
    name_cache = []

    def __init__(self,iz=-1,zz=0.,wz=None,nepmax=None,laccumulate=0,
                 lepsaveonce=None,name=None,lautodump=0,dumptofile=0):

        if wz is None and w3d.dz != 0.: self.wz = w3d.dz
        else:                           self.wz = wz
        if lepsaveonce is not None: top.lepsaveonce = lepsaveonce

        self.type = 'ep'
        self.topnwinname  = 'nepwin'
        self.topizwinname  = 'izepwin'
        self.topzzwinname  = 'zzepwin'
        self.topwzwinname  = 'wzepwin'
        self.topnmaxname  = 'nepmax'
        self.topnpidmaxname  = 'npidepmax'
        self.topnname  = 'nep'
        self.toptname  = 'tep'
        self.topxname  = 'xep'
        self.topyname  = 'yep'
        self.topuxname  = 'uxep'
        self.topuyname  = 'uyep'
        self.topuzname  = 'uzep'
        self.toppidname  = 'pidep'
        self.topgroupname = 'ExtPart'

        ParticleAccumulator.__init__(self,iz=iz,zz=zz,nmax=nepmax,
                                     laccumulate=laccumulate,
                                     name=name,
                                     lautodump=lautodump,dumptofile=dumptofile)

    def setupid(self):
        ParticleAccumulator.setupid(self)
        if self.wz is not None:
            top.wzepwin[-1] = self.wz
        else:
            installafterfs(self.setupwz)

    def setupwz(self):
        """This is needed in case the instance is created before the generate
           and wz is not explicitly specified. Presumably, w3d.dz should be
           set after a field solve is done."""
        if w3d.dz > 0. or self.wz is not None:
            if self.wz is None: self.wz = w3d.dz
            id = self.getid(safe=1)
            if id is not None: top.wzepwin[id] = self.wz
            uninstallafterfs(self.setupwz)

############################################################################
class ZCrossingParticles(ParticleAccumulator):
    """This class defines a container to setup and keep track of particles
crossing a Z location. It can optionally accumulate the data over multiple
time steps. Each time step, during the particle advance, a check is made if
the particles will cross the z location during the advance, and if so, the
particle data is saved. The current velocity is saved, as well as the position
which is advanced the fraction of the time step the brings the particle z to
the diagnostic z location.

This gives nice results since the calculation of the particle data at the
diagnostic location is consistent with the particle advance used during the
simulation.

The creator options are:

 - iz: grid location where the extrapolated data is saved.
 - zz: lab location where data is saved.
 - nmax: max size of the arrays. Defaults to 3*top.pnumz[iz] if non-zero,
         otherwise 10000.
 - laccumulate=0: when true, particles are accumulated over multiple steps.
 - name=None: descriptive name for location. It must be unique for each
              instance of the diagnostic.
 - lautodump=0: when true, after the grid moves beyond the z location,
                automatically dump the data to a file, clear the arrays and
                disable itself. Also must have name set - name is appended to
                the name of the files where the data is stored. The object
                writes itself out to a pdb file.
 - dumptofile=0: when true, the particle data is always dumped to a file
                 and not saved in memory. Name must be set - name is appended
                 to the name of the files where the data is stored.
                 Setting this to true implies that the data is accumulated.
                 Use the restoredata method to read the data back in.

One of iz or zz must be specified.

Available methods:

 - setaccumulate(v=1): Turns accumulation on or off. If turned on, all
                       currently saved data is deleted.
 - clear(): Will clear the existing data.
 - disable(): Turns off collecting of data
 - enable(): Turns on collecting of data (only needed after disable)

The follow all take an optional argument to specify species number.

 - getn: Get number of particles
 - gett: Get time at which particle was saved
 - getx, y, ux, uy, uz, vx, vy, vz, xp, yp, r, theta, rp: Get the various
     coordinates or velocity of particles
 - getpid: Get the particle ID info

The following are available plot routines. All take an optional argument to
specify species number. Additional arguments are the same as the 'pp' plotting
routines (such as ppxxp).

 - pxy, pxxp, pyyp, pxpyp, prrp, ptx, pty, ptxp, ptyp, ptux, ptuy, ptuz, ptvx
 - ptvy, ptvz, ptrace
    """
    name_cache = []

    def __init__(self,iz=-1,zz=0.,nmax=None,laccumulate=0,
                 name=None,lautodump=0,dumptofile=0):

        if w3d.dz != 0.:
            self.wz = w3d.dz
        else:
            self.wz = None

        self.type = 'zc'
        self.topnwinname  = 'nzcwin'
        self.topizwinname  = 'izzcwin'
        self.topzzwinname  = 'zzzcwin'
        self.topwzwinname  = 'wzzcwin'
        self.topnmaxname  = 'nzcmax'
        self.topnpidmaxname  = 'npidzcmax'
        self.topnname  = 'nzc'
        self.toptname  = 'tzc'
        self.topxname  = 'xzc'
        self.topyname  = 'yzc'
        self.topuxname  = 'uxzc'
        self.topuyname  = 'uyzc'
        self.topuzname  = 'uzzc'
        self.toppidname  = 'pidzc'
        self.topgroupname = 'ZCrossingParticles'

        ParticleAccumulator.__init__(self,iz=iz,zz=zz,nmax=nmax,
                                     laccumulate=laccumulate,
                                     name=name,
                                     lautodump=lautodump,dumptofile=dumptofile)

    def setupid(self):
        ParticleAccumulator.setupid(self)
        if self.wz is None:
            installaftergenerate(self.setupwz)

    def setupwz(self):
        """This is needed in case the instance is created before the generate
           and wz is not explicitly specified. Presumably, w3d.dz should be
           set after the generate."""
        if w3d.dz > 0. or self.wz is not None:
            if self.wz is None: self.wz = w3d.dz
            uninstallaftergenerate(self.setupwz)

##############################################################################
def dumpExtPart(object,filename):
    """Dump the saved extrapolated data to a file
   - filename: The name of the file to save the data in"""
    if me == 0:
        # --- Only PE0 writes the object to the file since it is the processor
        # --- where the data is gathered.
        ff = open(filename,'w')
        cPickle.dump(object,ff,1)
        ff.close()

def restoreExtPart(object,filename):
    """Restore extrapolated data from the given file"""
    if me == 0:
        # --- Only PE0 wrote the object to the file since it is the processor
        # --- where the data was gathered.
        ff = open(filename,'r')
        result = cPickle.load(ff)
        ff.close()
        result.enable()
        # --- Get the value of iz
        iz = result.iz
    else:
        # --- Create temp iz
        iz = 0
    # --- PE0 broadcasts its value of iz to all of the other processors
    # --- which create new instances of the ExtPart class.
    iz = parallel.broadcast(iz)
    if me > 0: result = ExtPart(iz)
    return result


##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
############################################################################
# --- This is only here to allow reading in old data files.
class ExtPartDeprecated:
    def __init__(self,iz=-1,zz=0.,wz=None,nepmax=None,laccumulate=0,
                 lepsaveonce=None,name=None,lautodump=0,dumptofile=0):
        # --- Save input values, getting default values when needed
        assert type(iz) is IntType,"iz must be an integer"
        assert iz >= 0 or zz is not None,"Either iz or zz must be specified"
        self.iz = iz
        self.zz = zz
        if wz is None and w3d.dz != 0.: self.wz = w3d.dz
        else:                           self.wz = wz
        self.laccumulate = laccumulate
        if lepsaveonce is not None: top.lepsaveonce = lepsaveonce
        self.lautodump = lautodump
        self.name = name
        self.dumptofile = dumptofile
        self.dt = top.dt
        if nepmax is None:
            self.nepmax = 10000
            if top.allocated("pnumz") and 0 <= self.getiz() <= top.nzmmnt:
                if top.pnumz[self.getiz(),-1] > 0:
                    if top.nszmmnt > 1:
                        self.nepmax = nint(max(top.pnumz[self.getiz(),:-1])*3)
                    else:
                        self.nepmax = nint(top.pnumz[self.getiz(),-1]*3)
        else:
            self.nepmax = nepmax
        # --- Add this new window to the ExtPart group in top
        self.enabled = 0
        self.enable()
        # --- Setup empty arrays for accumulation if laccumulate if true.
        # --- Otherwise, the arrays will just point to the data in ExtPart.
        self.setuparrays(top.ns)

    def getiz(self):
        if self.iz >= 0:
            return self.iz
        else:
            if top.dzm != 0.:
                return int((self.zz - top.zmmntmin)*top.dzmi)
            else:
                return -1

    def getzz(self):
        if self.iz >= 0:
            return top.zmmntmin + self.iz*top.dzm
        else:
            return self.zz

    def setuparrays(self,ns,bump=None):
        if self.laccumulate and not self.dumptofile:
            if bump is None: bump = self.nepmax
            self.tep = []
            self.xep = []
            self.yep = []
            self.uxep = []
            self.uyep = []
            self.uzep = []
            self.pidep = []
            for js in range(ns):
                self.tep.append(AppendableArray(self.nepmax,typecode='d',
                                                autobump=bump))
                self.xep.append(AppendableArray(self.nepmax,typecode='d',
                                                autobump=bump))
                self.yep.append(AppendableArray(self.nepmax,typecode='d',
                                                autobump=bump))
                self.uxep.append(AppendableArray(self.nepmax,typecode='d',
                                                 autobump=bump))
                self.uyep.append(AppendableArray(self.nepmax,typecode='d',
                                                 autobump=bump))
                self.uzep.append(AppendableArray(self.nepmax,typecode='d',
                                                 autobump=bump))
                self.pidep.append(AppendableArray(self.nepmax,typecode='d',
                                                  autobump=bump,
                                                  unitshape=(top.npidepmax,)))
        else:
            self.tep = ns*[zeros(0,'d')]
            self.xep = ns*[zeros(0,'d')]
            self.yep = ns*[zeros(0,'d')]
            self.uxep = ns*[zeros(0,'d')]
            self.uyep = ns*[zeros(0,'d')]
            self.uzep = ns*[zeros(0,'d')]
            self.pidep = ns*[zeros((0,top.npidepmax),'d')]

    def addspecies(self):
        if self.laccumulate and not self.dumptofile:
            for js in range(len(self.tep),top.ns):
                bump = self.nepmax
                self.tep.append(AppendableArray(self.nepmax,typecode='d',
                                                autobump=bump))
                self.xep.append(AppendableArray(self.nepmax,typecode='d',
                                                autobump=bump))
                self.yep.append(AppendableArray(self.nepmax,typecode='d',
                                                autobump=bump))
                self.uxep.append(AppendableArray(self.nepmax,typecode='d',
                                                 autobump=bump))
                self.uyep.append(AppendableArray(self.nepmax,typecode='d',
                                                 autobump=bump))
                self.uzep.append(AppendableArray(self.nepmax,typecode='d',
                                                 autobump=bump))
                self.pidep.append(AppendableArray(self.nepmax,typecode='d',
                                                  autobump=bump,
                                                  unitshape=(top.npidepmax,)))
        else:
            self.tep = top.ns*[zeros(0,'d')]
            self.xep = top.ns*[zeros(0,'d')]
            self.yep = top.ns*[zeros(0,'d')]
            self.uxep = top.ns*[zeros(0,'d')]
            self.uyep = top.ns*[zeros(0,'d')]
            self.uzep = top.ns*[zeros(0,'d')]
            self.pidep = top.ns*[zeros((0,top.npidepmax),'d')]

    def clear(self):
        self.setuparrays(top.ns)

    def getid(self,safe=0):
        'If safe, then return None if id is not found rather than raising error'
        assert self.enabled,"This window is disabled and there is no associated id"
        for i in range(top.nepwin):
            if top.izepwin[i] == self.iz and self.iz >= 0: return i
            if top.zzepwin[i] == self.zz and self.iz == -1: return i
        if not safe:
            raise "Uh Ooh! Somehow the window was deleted! I can't continue! "+self.titleright()
        else:
            return None

    def setupid(self):
        top.nepwin = top.nepwin + 1
        if top.nepmax < self.nepmax: top.nepmax = self.nepmax
        err = gchange("ExtPart")
        top.izepwin[-1] = self.iz
        top.zzepwin[-1] = self.zz
        if self.wz is not None:
            top.wzepwin[-1] = self.wz
        else:
            installafterfs(self.setupwz)

    def setupwz(self):
        """This is needed in case the instance is created before the generate
           and wz is not explicitly specified. Presumably, w3d.dz should be
           set after a field solve is done."""
        if w3d.dz > 0. or self.wz is not None:
            if self.wz is None: self.wz = w3d.dz
            id = self.getid(safe=1)
            if id is not None: top.wzepwin[id] = self.wz
            uninstallafterfs(self.setupwz)

    def updatenpidepmax(self):
        # --- First check if npid has changed since top.npidepmax was last set
        if top.npidepmax != top.npid:
            top.npidepmax = top.npid
            gchange('ExtPart')
        # --- Now make sure that the data arrays are changed appropriately.
        # --- This is here in case top.npidepmax had been update elsewhere
        # --- but self.pidep had not.
        if self.laccumulate and not self.dumptofile:
            for js in range(top.ns):
                if self.pidep[js].unitshape()[0] != top.npidepmax:
                    self.pidep[js].reshape((top.npidepmax,))

    def enable(self):
        # --- Add this window to the list
        # --- Only add this location to the list if it is not already there.
        # --- Note that it is not an error to have more than one instance
        # --- have the same location. For example one could be accumulating
        # --- while another isn't or the widths could be different.
        if self.enabled: return
        self.setupid()
        # --- Set so accumulate method is called after time steps
        installafterstep(self.accumulate)
        self.enabled = 1

    def disable(self):
        if not self.enabled: return
        # --- Set so accumulate method is not called after time steps
        uninstallafterstep(self.accumulate)
        # --- Check if setupwz is installed, and if so, uninstall it
        if isinstalledafterfs(self.setupwz): uninstallafterfs(self.setupwz)
        # --- Remove this window from the list. Turn safe on when gettin
        # --- the id, since it may for some reason not be consistent.
        id = self.getid(safe=1)
        if id is not None:
            for i in range(id,top.nepwin-1):
                top.izepwin[i] = top.izepwin[i+1]
                top.zzepwin[i] = top.zzepwin[i+1]
                top.wzepwin[i] = top.wzepwin[i+1]
            top.nepwin = top.nepwin - 1
            gchange("ExtPart")
        self.enabled = 0

    def __setstate__(self,dict):
        self.__dict__.update(dict)
        if not self.enabled: return
        id = self.getid(safe=1)
        if id is None: self.setupid()
        self.restoredata()
        if not isinstalledafterstep(self.accumulate):
            installafterstep(self.accumulate)

    def setaccumulate(self,v=1):
        self.laccumulate = v
        if self.laccumulate: self.setuparrays(top.ns)

    def accumulate(self):
        # --- If top.nepwin is 0 then something is really wrong - this routine
        # --- should never be called if top.nepwin is zero.
        if top.nepwin == 0: return
        # --- Check if the number of species has changed. This is done to ensure
        # --- crashes don't happen.
        if top.ns > self.getns():
            self.addspecies()
        # --- If this windows is outside of the grid, then just return.
        if (self.iz == -1 and
            (self.zz+self.wz < w3d.zmmin+top.zbeam or
             self.zz-self.wz > w3d.zmmax+top.zbeam)):
            self.autodump()
            return
        id = self.getid()
        self.updatenpidepmax()
        # --- Loop over species, collecting only ones where some particles
        # --- were saved.
        for js in range(top.ns):
            # --- Gather the data.
            # --- In serial, the arrays are just returned as is.
            # --- In parallel, if dumptofile is not being done, then the data
            # --- is gathered onto PE0 and empty arrays are returned on other
            # --- processors. If dumptofile is on, then the data is kept local
            # --- and each processor writes out its own file.

            # --- First, get the local data.
            nn = top.nep[id,js]
            t = top.tep[:nn,id,js]
            x = top.xep[:nn,id,js]
            y = top.yep[:nn,id,js]
            ux = top.uxep[:nn,id,js]
            uy = top.uyep[:nn,id,js]
            uz = top.uzep[:nn,id,js]
            if top.npidepmax > 0:
                pid = top.pidep[:nn,:,id,js]
            else:
                pid = zeros((nn,0),'d')

            if not self.dumptofile:
                # --- Gather the data onto PE0.
                ntot = globalsum(nn)
                if ntot == 0: continue

                t = gatherarray(t,othersempty=1)
                x = gatherarray(x,othersempty=1)
                y = gatherarray(y,othersempty=1)
                ux = gatherarray(ux,othersempty=1)
                uy = gatherarray(uy,othersempty=1)
                uz = gatherarray(uz,othersempty=1)
                if top.npidepmax > 0:
                    pid = gatherarray(pid,othersempty=1)
                else:
                    pid = zeros((ntot,0),'d')

            if self.laccumulate and not self.dumptofile:

                # --- Only PE0 accumulates the data.
                if me == 0:
                    self.tep[js].append(t.copy())
                    self.xep[js].append(x.copy())
                    self.yep[js].append(y.copy())
                    self.uxep[js].append(ux.copy())
                    self.uyep[js].append(uy.copy())
                    self.uzep[js].append(uz.copy())
                    self.pidep[js].append(pid.copy())

            else:

                self.tep[js] = t
                self.xep[js] = x
                self.yep[js] = y
                self.uxep[js] = ux
                self.uyep[js] = uy
                self.uzep[js] = uz
                self.pidep[js] = pid

        if self.dumptofile: self.dodumptofile()
        # --- Force nep to zero to ensure that particles are not saved twice.
        top.nep[id,:] = 0

    ############################################################################
    def autodump(self):
        if not self.lautodump or self.name is None: return
        if not self.laccumulate and not self.dumptofile: return
        if self.iz >= 0: return
        if self.zz+self.wz > w3d.zmmin+top.zbeam: return
        # --- Check if there is any data. If there is none, then don't make
        # --- a dump.
        ntot = 0
        for js in range(self.getns()):
            ntot = ntot + self.getn(js=js)
        if ntot > 0:
            ff = None
            #try:
            #    ff = PWpyt.PW(self.name+'_epdump.pyt')
            #    dumpsmode = 1
            #except:
            ff = PW.PW(self.name+'_%05d_%05d_epdump.pdb'%(me,npes))
            dumpsmode = 0
            if ff is None:
                 print "ExtPart: %s unable to dump data to file."%self.name
                 return
            ff.write(self.name+'@pickle',cPickle.dumps(self,dumpsmode))
            ff.close()
        self.nepmax = 1
        self.clear()
        # --- Disable is done last so that the object written out to the
        # --- file is still enabled. That flag is used in restoredata to
        # --- determine whether or not to restore the data. The logic is set
        # --- so that the object in an autodump file will restore the data
        # --- but one is a generic dump file won't (unless it was not auto
        # --- dumped, in which case the object in the generic dump is the only
        # --- copy). There will of course be exceptions, so restoredata takes
        # --- and option argument to force restoration of data, and the
        # --- extforcenorestore function turns any restores off.
        self.disable()

    ############################################################################
    def dodumptofile(self):
        #self.dodumptofilePDB()
        self.dodumptofilePickle()

    def dodumptofilePDB(self):
        #if me != 0: return
        ff = None
        for js in range(top.ns):
            if len(self.tep[js][:]) > 0:
                if ff is None:
                    # --- Only create the file if there is data to write out.
                    ff = PW.PW(self.name+'_ep_%05d_%05d.pdb'%(me,npes),'a',verbose=0)
                if ff is None:
                    print "ExtPart: %s unable to dump data to file."%self.name
                    return
                suffix = "_%d_%d"%(top.it,js)
                ff.write('n'+suffix,len(self.tep[js][:]))
                ff.write('t'+suffix,self.tep[js][:])
                ff.write('x'+suffix,self.xep[js][:])
                ff.write('y'+suffix,self.yep[js][:])
                ff.write('ux'+suffix,self.uxep[js][:])
                ff.write('uy'+suffix,self.uyep[js][:])
                ff.write('uz'+suffix,self.uzep[js][:])
                ff.write('pid'+suffix,self.pidep[js][...])
        if ff is not None:
            ff.close()

    def dodumptofilePickle(self):
        #if me != 0: return
        ff = None
        for js in range(top.ns):
            if len(self.tep[js][:]) > 0:
                if ff is None:
                    # --- Only create the file if there is data to write out.
                    if npes > 1:
                        ff = open(self.name+'_ep_%05d_%05d.pkl'%(me,npes),'a')
                    else:
                        ff = open(self.name+'_ep.pkl','a')
                if ff is None:
                    print "ExtPart: %s unable to dump data to file."%self.name
                    return
                suffix = "_%d_%d"%(top.it,js)
                cPickle.dump(('n'+suffix,len(self.tep[js][:])),ff,-1)
                cPickle.dump(('t'+suffix,self.tep[js][:]),ff,-1)
                cPickle.dump(('x'+suffix,self.xep[js][:]),ff,-1)
                cPickle.dump(('y'+suffix,self.yep[js][:]),ff,-1)
                cPickle.dump(('ux'+suffix,self.uxep[js][:]),ff,-1)
                cPickle.dump(('uy'+suffix,self.uyep[js][:]),ff,-1)
                cPickle.dump(('uz'+suffix,self.uzep[js][:]),ff,-1)
                cPickle.dump(('pid'+suffix,self.pidep[js][...]),ff,-1)
        if ff is not None:
            ff.close()

    ############################################################################
    def restoredata(self,lforce=0,files=None):
        #self.restoredataPDB(0,files)
        self.restoredataPickle(0,files)

    def restoredataPickle(self,lforce=0,files=None,names=None,nprocs=None):
        """
Restores data dumped to a file. Note that this turns off the dumptofile
feature.
  - lforce=0: if true, force a restore, despite the value of enabled.
        """
        if files is None: files = []
        if names is None: names = []
        if not self.dumptofile: return
        if not lforce and (not self.enabled or _extforcenorestore): return
        self.dumptofile = 0
        self.laccumulate = 1

        # --- Check some input, converting to lists if needed.
        if not isinstance(files,ListType):
            files = [files]
        if not isinstance(names,ListType):
            names = [names]

        # --- If the list of files was not given, then it needs to be
        # --- generated.
        fnametries = []
        if len(files) == 0:
            # --- If a name was not given, use the instance's name
            if len(names) == 0:
                names = [self.name]
            for name in names:
                if nprocs is None and npes <= 1:
                    # --- Use the serial naming
                    files = [self.name+'_ep.pkl']
                else:
                    if npes > 1:
                        # --- If currently running in parallel, only read in
                        # --- the data for this processor.
                        nplist = [me]
                        nprocs = npes
                    else:
                        # --- Read the data in from all processors.
                        nplist = range(nprocs)
                    for iproc in nplist:
                        fname = self.name+'_ep_%05d_%05d.pkl'%(iproc,nprocs)
                        fnametries.append(fname)
                        if os.path.exists(fname):
                            files.append(fname)

        if len(files) == 0:
            print "ExtPart restoredata: warning, no files were found, nothing will be restored"
            if len(fnametries) > 0:
                print "Tried the filenames:",fnametries

        #datadict = self.getPDBdatadict(files)
        datadict = self.getPickledatadict(files)

        # --- Get total number of particles
        ntot = []
        jsmax = 0
        for var,val in datadict.items():
            if var[0] == 'n':
                name,ii,js = var.split('_')
                jsmax = max(jsmax,eval(js))
                while jsmax >= len(ntot): ntot.append(0)
                ntot[jsmax] = ntot[jsmax] + val

        # --- Get the size of the pid array. If top.npidepmax is zero, that
        # --- probably means that this data is being read in for
        # --- postprocessing, since the pid arrays was not setup.
        # --- If top.npidepmax is nonzero and different than the size
        # --- of pid, raise an exception since there is likely something
        # --- wrong.
        for var,val in datadict.items():
            if var[0:3] == 'pid':
                try:
                    npid = val.shape[1]
                except IndexError:
                    # --- This is needed so that old datasets can be read in.
                    # --- Old datasets were not written properly, leaving
                    # --- pid a 1-D array.
                    npid = 1
                if top.npidepmax == 0:
                    top.npidepmax = npid
                    top.npid = npid
                else:
                    assert top.npidepmax == npid,\
                           'npid is different than in the run where the ExtPart data was saved'
                break


        if len(ntot) == 0:
            # --- No data was read in, so this value won't be used.
            # --- This can happen in a parallel run since not all processors
            # --- will necessarily write out data.
            # --- This case is needed since max(ntot) would give an error.
            bump = 1
        else:
            # --- Set the bump size so that the arrays will be made just
            # --- large enough to hold the data.
            bump = max(ntot) + 1
        self.setuparrays(jsmax+1,bump=bump)

        # --- This loop must be ordered because of the append
        varlist = datadict.keys()
        varlist.sort()
        for var in varlist:
            if var[0] == 'n':
                name,iis,jss = var.split('_')
                nn = datadict[var]
                ii = eval(iis)
                js = eval(jss)
                self.tep[js].append(datadict['t_%d_%d'%(ii,js)])
                self.xep[js].append(datadict['x_%d_%d'%(ii,js)])
                self.yep[js].append(datadict['y_%d_%d'%(ii,js)])
                self.uxep[js].append(datadict['ux_%d_%d'%(ii,js)])
                self.uyep[js].append(datadict['uy_%d_%d'%(ii,js)])
                self.uzep[js].append(datadict['uz_%d_%d'%(ii,js)])
                pid = datadict['pid_%d_%d'%(ii,js)]
                if len(pid.shape) == 1:
                    pid.shape = (pid.shape[0],1)
                self.pidep[js].append(pid)

    def getPickledatadict(self,files):
        # --- Read in all of the data into a dictionary.
        datadict = {}
        for file in files:
            ff = open(file,'r')
            while 1:
                try:
                    data = cPickle.load(ff)
                except:
                    break
                datadict[data[0]] = data[1]
            ff.close()
        return datadict

    def getPDBdatadict(self,files):
        # --- Read in all of the data into a dictionary.
        datadict = {}
        for file in files:
            try:
                ff = PRpyt.PR(file,verbose=0)
            except:
                ff = PR.PR(file,verbose=0)
            for var in ff.inquire_names():
                datadict[var] = ff.read(var)
            ff.close()
        return datadict

    def restoredataPDBold(self,lforce=0,files=[]):
        """
This is kept around for compatibility with old data files.
Restores data dumped to a file. Note that this turns off the dumptofile
feature.
  - lforce=0: if true, force a restore, despite the value of enabled.
        """
        if not self.dumptofile: return
        if not lforce and (not self.enabled or _extforcenorestore): return
        self.dumptofile = 0
        self.laccumulate = 1
        try:
            ff = PRpyt.PR(self.name+'_ep.pyt','a',verbose=0)
        except:
            ff = PR.PR(self.name+'_ep.pdb','a',verbose=0)
        # --- Get total number of particles
        ntot = []
        jsmax = 0
        varlist = list(ff.inquire_names())
        varlist.sort()
        for var in varlist:
            if var[0] == 'n':
                name,ii,js = var.split('_')
                jsmax = max(jsmax,eval(js))
                while jsmax >= len(ntot): ntot.append(0)
                ntot[jsmax] = ntot[jsmax] + ff.read(var)
        self.setuparrays(jsmax+1,bump=max(array(ntot))+1)
        for var in varlist:
            if var[0] == 'n':
                name,iis,jss = var.split('_')
                nn = ff.read(var)
                ii = eval(iis)
                js = eval(jss)
                self.tep[js].append(ff.read('t_%d_%d'%(ii,js)))
                self.xep[js].append(ff.read('x_%d_%d'%(ii,js)))
                self.yep[js].append(ff.read('y_%d_%d'%(ii,js)))
                self.uxep[js].append(ff.read('ux_%d_%d'%(ii,js)))
                self.uyep[js].append(ff.read('uy_%d_%d'%(ii,js)))
                self.uzep[js].append(ff.read('uz_%d_%d'%(ii,js)))
                self.pidep[js].append(ff.read('pid_%d_%d'%(ii,js)))
        ff.close()

    ############################################################################
    def selectparticles(self,val,js=0,tc=None,wt=None,tp=None,z=None,v=None,
                        gather=1,bcast=1):
        """
 - js=0: Species number to gather from.
         If js=='all', then the quantity is gathered from all species into
         a single array.
 - tc=None: time at which to gather particles from. When not given, returns
            all particles.
 - wt=top.dt: Width of region around tc from which to select particles.
 - tp=gett(js=js): Time value to use for the particle selection
 - z=None: when specified, projects the data to the given z location
 - gather=1: in parallel, when true, the particles are all gathered onto
             one processor
 - bcast=1: in parallel, when true, the gathered particles are broadcast
            to all processors
        """
        if js == 'all':
            nn = sum(map(len,val))
            rr = AppendableArray(initlen=nn,typecode='d')
            for js in range(len(val)):
                rr.append(self.selectparticles(val,js,tc,wt,tp))
            result = rr[...]
        elif tc is None:
            result = val[js][...]
        else:
            if wt is None: wt = self.dt
            if tp is None: tp = self.tep[js][:]
            ii = compress((tc-wt<tp)&(tp<tc+wt),arange(len(tp)))
            result = take(val[js][...],ii)

        if z is not None and v is not None:
            # --- Project to the new z value given the current velocity
            try:
                # --- v is either one of the 'get' methods
                v = v(js,tc,wt,tp)
            except TypeError:
                # --- or a constant
                pass
            zep = self.getzz()
            vzep = self.getvz(js,tc,wt,tp)
            delt = (z - zep)/vzep
            result = result + v*delt

        if lparallel and gather: result = gatherarray(result,bcast=bcast)
        return result

    def getns(self):
        return len(self.tep)
    def gett(self,js=0,tc=None,wt=None,tp=None,z=None,gather=1,bcast=1):
        return self.selectparticles(self.tep,js,tc,wt,tp,z,1.,
                                    gather=gather,bcast=bcast)
    def getx(self,js=0,tc=None,wt=None,tp=None,z=None,gather=1,bcast=1):
        return self.selectparticles(self.xep,js,tc,wt,tp,z,self.getux,
                                    gather=gather,bcast=bcast)
    def gety(self,js=0,tc=None,wt=None,tp=None,z=None,gather=1,bcast=1):
        return self.selectparticles(self.yep,js,tc,wt,tp,z,self.getuy,
                                    gather=gather,bcast=bcast)
    def getux(self,js=0,tc=None,wt=None,tp=None,z=None,gather=1,bcast=1):
        return self.selectparticles(self.uxep,js,tc,wt,tp,
                                    gather=gather,bcast=bcast)
    def getuy(self,js=0,tc=None,wt=None,tp=None,z=None,gather=1,bcast=1):
        return self.selectparticles(self.uyep,js,tc,wt,tp,
                                    gather=gather,bcast=bcast)
    def getuz(self,js=0,tc=None,wt=None,tp=None,z=None,gather=1,bcast=1):
        return self.selectparticles(self.uzep,js,tc,wt,tp,
                                    gather=gather,bcast=bcast)
    def getvx(self,js=0,tc=None,wt=None,tp=None,z=None,gather=1,bcast=1):
        return self.selectparticles(self.uxep,js,tc,wt,tp,
                                    gather=gather,bcast=bcast)
    def getvy(self,js=0,tc=None,wt=None,tp=None,z=None,gather=1,bcast=1):
        return self.selectparticles(self.uyep,js,tc,wt,tp,
                                    gather=gather,bcast=bcast)
    def getvz(self,js=0,tc=None,wt=None,tp=None,z=None,gather=1,bcast=1):
        return self.selectparticles(self.uzep,js,tc,wt,tp,
                                    gather=gather,bcast=bcast)
    def getpid(self,js=0,tc=None,wt=None,tp=None,z=None,id=0,gather=1,bcast=1):
        self.updatenpidepmax()
        if top.npidepmax > 0:
            return self.selectparticles(self.pidep,js,tc,wt,tp,
                                        gather=gather,bcast=bcast)[:,id]
        else:
            return zeros((self.getn(js,tc,wt,tp,
                                    gather=gather,bcast=bcast),0),'d')

    def getxp(self,js=0,tc=None,wt=None,tp=None,z=None,gather=1,bcast=1):
        return (self.getux(js,tc,wt,tp,gather=gather,bcast=bcast)/
                self.getuz(js,tc,wt,tp,gather=gather,bcast=bcast))
    def getyp(self,js=0,tc=None,wt=None,tp=None,z=None,gather=1,bcast=1):
        return (self.getuy(js,tc,wt,tp,gather=gather,bcast=bcast)/
                self.getuz(js,tc,wt,tp,gather=gather,bcast=bcast))
    def getr(self,js=0,tc=None,wt=None,tp=None,z=None,gather=1,bcast=1):
        return sqrt(self.getx(js,tc,wt,tp,z,gather=gather,bcast=bcast)**2 +
                    self.gety(js,tc,wt,tp,z,gather=gather,bcast=bcast)**2)
    def gettheta(self,js=0,tc=None,wt=None,tp=None,z=None,gather=1,bcast=1):
        return arctan2(self.gety(js,tc,wt,tp,z,gather=gather,bcast=bcast),
                       self.getx(js,tc,wt,tp,z,gather=gather,bcast=bcast))
    def getrp(self,js=0,tc=None,wt=None,tp=None,z=None,gather=1,bcast=1):
        return (self.getxp(js,tc,wt,tp,gather=gather,bcast=bcast)*
                    cos(self.gettheta(js,tc,wt,tp,gather=gather,bcast=bcast)) +
                self.getyp(js,tc,wt,tp,gather=gather,bcast=bcast)*
                    sin(self.gettheta(js,tc,wt,tp,gather=gather,bcast=bcast)))
    def getn(self,js=0,tc=None,wt=None,tp=None,z=None,gather=1,bcast=1):
        return len(self.gett(js,tc,wt,tp,gather=gather,bcast=bcast))

    def xxpslope(self,js=0,tc=None,wt=None,tp=None,z=None):
        if self.getn(js,tc,wt,tp) == 0:
            return 0.
        else:
            return (
                (ave(self.getx(js,tc,wt,tp,z)*self.getxp(js,tc,wt,tp)) -
                 ave(self.getx(js,tc,wt,tp,z))*ave(self.getxp(js,tc,wt,tp)))/
                (ave(self.getx(js,tc,wt,tp,z)*self.getx(js,tc,wt,tp,z)) -
                 ave(self.getx(js,tc,wt,tp,z))*ave(self.getx(js,tc,wt,tp,z))))
    def yypslope(self,js=0,tc=None,wt=None,tp=None,z=None):
        if self.getn(js,tc,wt,tp) == 0:
            return 0.
        else:
            return (
                (ave(self.gety(js,tc,wt,tp,z)*self.getyp(js,tc,wt,tp)) -
                 ave(self.gety(js,tc,wt,tp,z))*ave(self.getyp(js,tc,wt,tp)))/
                (ave(self.gety(js,tc,wt,tp,z)*self.gety(js,tc,wt,tp,z)) -
                 ave(self.gety(js,tc,wt,tp,z))*ave(self.gety(js,tc,wt,tp,z))))
    def rrpslope(self,js=0,tc=None,wt=None,tp=None,z=None):
        if self.getn(js,tc,wt,tp) == 0:
            return 0.
        else:
            return (ave(self.getr(js,tc,wt,tp,z)*self.getrp(js,tc,wt,tp))/
                    ave(self.getr(js,tc,wt,tp,z)**2))

    ############################################################################
    ############################################################################
    # --- Define plotting routines for the extrapolated particles.

    def checkplotargs(self,kw):
        """Convenience routine to check arguments of particle plot routines.
Warning: this has the side affect of adding the arguement allowbadargs to
the kw dictionary. This is done since the calls to these functions here to
make the plots may have unused arguements since the entire kw list passed
into each of the pp plotting routines is passed into each of these
functions.
        """
        badargs = ppgeneric(checkargs=1,kwdict=kw)
        kw['allowbadargs'] = 1
        if badargs: raise 'bad arguments ',' '.join(badargs.keys())

    def titleright(self,tc=None,wt=None,z=None,slope=None):
        if tc is None:
            ttext = ''
        else:
            if wt is None: wt = self.dt
            ttext = '  time = %e ^+_-%e'%(tc,wt)
        if self.iz >= 0:
            ztext =  'iz = %d (z = %f m)'%(self.iz,w3d.zmmin+self.iz*w3d.dz)
        else:
            ztext =  'z = %f m'%self.zz
        if z is not None:
            ztext = ztext + ' projected to z = %f m'%z
        if slope is None:
            slopetext = ''
        else:
            slopetext = '  slope=%7.4f'%slope
        return ztext + ttext + slopetext

    def ppmultispecies(self,pp,args,kw):
        """Check if js is -1 or a list, in which case call the pp function for
each species and each one in the list. Also assign colors accordingly
        """
        args = list(args)
        js = args[0]
        if js != -1 and type(js) != ListType:
            return false
        else:
            if js == -1: js = range(self.getns())
            ncolor = kw.get('ncolor',240)
            color = kw.get('color',range(0,ncolor,ncolor/len(js)))
            for i in xrange(len(js)):
                args[0] = js[i]
                kw['color'] = color[i]
                apply(pp,args,kw)
            return true

    ############################################################################
    def pxy(self,js=0,tc=None,wt=None,tp=None,z=None,**kw):
        """Plots X-Y for extraploated particles"""
        self.checkplotargs(kw)
        if self.ppmultispecies(self.pxy,(js,tc,wt,tp,z),kw): return
        x = self.getx(js,tc,wt,tp,z)
        y = self.gety(js,tc,wt,tp,z)
        if 'pplimits' in kw.keys():
            kw['lframe'] = 1
        else:
            kw['pplimits'] = (top.xplmin,top.xplmax,top.yplmin,top.yplmax)
        settitles("Y vs X","X","Y",self.titleright(tc,wt,z))
        return ppgeneric(y,x,kwdict=kw)

    ############################################################################
    def pxxp(self,js=0,tc=None,wt=None,tp=None,z=None,slope=0.,offset=0.,
             **kw):
        """Plots X-X' for extraploated particles
 - slope=0.: slope subtracted from xp, it is calculated automatically
 - offset=0.: offset in x
        """
        self.checkplotargs(kw)
        if self.ppmultispecies(self.pxxp,(js,tc,wt,tp,z),kw): return
        x = self.getx(js,tc,wt,tp,z)
        xp = self.getxp(js,tc,wt,tp)
        if type(slope) == type(''):
            if len(x) > 0:
                slope = (ave(x*xp)-ave(x)*ave(xp))/(ave(x*x) - ave(x)**2)
                offset = ave(xp)-slope*ave(x)
            else:
                slope = 0.
                offset = 0.
        kw['slope'] = slope
        kw['offset'] = offset
        if 'pplimits' in kw.keys():
            kw['lframe'] = 1
        else:
            kw['pplimits'] = (top.xplmin,top.xplmax,top.xpplmin,top.xpplmax)
        settitles("X' vs X","X","X'",self.titleright(tc,wt,z,slope))
        return ppgeneric(xp,x,kwdict=kw)

    ############################################################################
    def pyyp(self,js=0,tc=None,wt=None,tp=None,z=None,slope=0.,offset=0.,
             **kw):
        """Plots Y-Y' for extraploated particles
 - slope=0.: slope subtracted from yp, it is calculated automatically
 - offset=0.: offset in y
        """
        self.checkplotargs(kw)
        if self.ppmultispecies(self.pyyp,(js,tc,wt,tp,z,slope,offset),kw):
            return
        y = self.gety(js,tc,wt,tp,z)
        yp = self.getyp(js,tc,wt,tp)
        if type(slope) == type(''):
            if len(y) > 0:
                slope = (ave(y*yp)-ave(y)*ave(yp))/(ave(y*y) - ave(y)**2)
                offset = ave(yp)-slope*ave(y)
            else:
                slope = 0.
                offset = 0.
        kw['slope'] = slope
        kw['offset'] = offset
        if 'pplimits' in kw.keys():
            kw['lframe'] = 1
        else:
            kw['pplimits'] = (top.yplmin,top.yplmax,top.ypplmin,top.ypplmax)
        settitles("Y' vs Y","Y","Y'",self.titleright(tc,wt,z,slope))
        return ppgeneric(yp,y,kwdict=kw)

    ############################################################################
    def pxpyp(self,js=0,tc=None,wt=None,tp=None,z=None,**kw):
        """Plots X'-Y' for extraploated particles"""
        self.checkplotargs(kw)
        if self.ppmultispecies(self.pxpyp,(js,tc,wt,tp),kw): return
        xp = self.getxp(js,tc,wt,tp)
        yp = self.getyp(js,tc,wt,tp)
        if 'pplimits' in kw.keys():
            kw['lframe'] = 1
        else:
            kw['pplimits'] = (top.xpplmin,top.xpplmax,top.ypplmin,top.ypplmax)
        settitles("Y' vs X'","X'","Y'",self.titleright(tc,wt,z))
        return ppgeneric(yp,xp,kwdict=kw)

    ############################################################################
    def prrp(self,js=0,tc=None,wt=None,tp=None,z=None,scale=0.,slope=0.,offset=0.,
             **kw):
        """Plots R-R' for extraploated particles
 - slope=0.: slope subtracted from the rp, it is calculated automatically
 - offset=0.: offset in r
        """
        self.checkplotargs(kw)
        if self.ppmultispecies(self.prrp,(js,tc,wt,tp,z,scale,slope,offset),kw):
            return
        x = self.getx(js,tc,wt,tp,z)
        y = self.gety(js,tc,wt,tp,z)
        xp = self.getxp(js,tc,wt,tp)
        yp = self.getyp(js,tc,wt,tp)
        xscale = 1.
        yscale = 1.
        xpscale = 1.
        ypscale = 1.
        if scale:
            xscale = 2.*sqrt(ave(x*x) - ave(x)**2)
            yscale = 2.*sqrt(ave(y*y) - ave(y)**2)
            xpscale = 2.*sqrt(ave(xp*xp) - ave(xp)**2)
            ypscale = 2.*sqrt(ave(yp*yp) - ave(yp)**2)
        x = x/xscale
        y = y/yscale
        xp = xp/xpscale
        yp = yp/ypscale
        r = sqrt(x**2 + y**2)
        t = arctan2(y,x)
        rp = xp*cos(t) + yp*sin(t)
        if type(slope) == type(''):
            if len(r) > 0:
              slope = ave(r*rp)/ave(r*r)
            else:
                slope = 0.
        kw['slope'] = slope
        if 'pplimits' in kw.keys():
            kw['lframe'] = 1
        else:
            kw['pplimits'] = (0.,max(top.xplmax/xscale,top.yplmax/yscale),
                              top.xpplmin/xpscale,top.xpplmax/ypscale)
        settitles("R' vs R","R","R'",self.titleright(tc,wt,z,slope))
        return ppgeneric(rp,r,kwdict=kw)

    ############################################################################
    def ptx(self,js=0,tc=None,wt=None,tp=None,z=None,**kw):
        """Plots time-X for extraploated particles"""
        self.checkplotargs(kw)
        if self.ppmultispecies(self.ptx,(js,tc,wt,tp,z),kw): return
        t = self.gett(js,tc,wt,tp,z)
        x = self.getx(js,tc,wt,tp,z)
        if 'pplimits' in kw.keys():
            kw['lframe'] = 1
        else:
            kw['pplimits'] = ('e','e',top.xplmin,top.xplmax)
        settitles("X vs time","time","X",self.titleright(tc,wt,z))
        return ppgeneric(x,t,kwdict=kw)

    ############################################################################
    def pty(self,js=0,tc=None,wt=None,tp=None,z=None,**kw):
        """Plots time-Y for extraploated particles"""
        self.checkplotargs(kw)
        if self.ppmultispecies(self.pty,(js,tc,wt,tp,z),kw): return
        t = self.gett(js,tc,wt,tp,z)
        y = self.gety(js,tc,wt,tp,z)
        if 'pplimits' in kw.keys():
            kw['lframe'] = 1
        else:
            kw['pplimits'] = ('e','e',top.yplmin,top.yplmax)
        settitles("Y vs time","time","Y",self.titleright(tc,wt,z))
        return ppgeneric(y,t,kwdict=kw)

    ############################################################################
    def ptxp(self,js=0,tc=None,wt=None,tp=None,z=None,**kw):
        """Plots time-X' for extraploated particles"""
        self.checkplotargs(kw)
        if self.ppmultispecies(self.ptxp,(js,tc,wt,tp,z),kw): return
        t = self.gett(js,tc,wt,tp,z)
        xp = self.getxp(js,tc,wt,tp)
        if 'pplimits' in kw.keys():
            kw['lframe'] = 1
        else:
            kw['pplimits'] = ('e','e',top.xpplmin,top.xpplmax)
        settitles("X' vs time","time","X'",self.titleright(tc,wt,z))
        return ppgeneric(xp,t,kwdict=kw)

    ############################################################################
    def ptyp(self,js=0,tc=None,wt=None,tp=None,z=None,**kw):
        """Plots time-Y' for extraploated particles"""
        self.checkplotargs(kw)
        if self.ppmultispecies(self.ptyp,(js,tc,wt,tp,z),kw): return
        t = self.gett(js,tc,wt,tp,z)
        yp = self.getyp(js,tc,wt,tp)
        if 'pplimits' in kw.keys():
            kw['lframe'] = 1
        else:
            kw['pplimits'] = ('e','e',top.ypplmin,top.ypplmax)
        settitles("Y' vs time","time","Y'",self.titleright(tc,wt,z))
        return ppgeneric(yp,t,kwdict=kw)

    ############################################################################
    def ptux(self,js=0,tc=None,wt=None,tp=None,z=None,**kw):
        """Plots time-ux for extraploated particles"""
        self.checkplotargs(kw)
        if self.ppmultispecies(self.ptux,(js,tc,wt,tp,z),kw): return
        t = self.gett(js,tc,wt,tp,z)
        ux = self.getux(js,tc,wt,tp)
        if 'pplimits' in kw.keys():
            kw['lframe'] = 1
        settitles("ux vs time","time","ux",self.titleright(tc,wt,z))
        return ppgeneric(ux,t,kwdict=kw)

    ############################################################################
    def ptuy(self,js=0,tc=None,wt=None,tp=None,z=None,**kw):
        """Plots time-uy for extraploated particles"""
        self.checkplotargs(kw)
        if self.ppmultispecies(self.ptuy,(js,tc,wt,tp,z),kw): return
        t = self.gett(js,tc,wt,tp,z)
        uy = self.getuy(js,tc,wt,tp)
        if 'pplimits' in kw.keys():
            kw['lframe'] = 1
        settitles("uy vs time","time","uy",self.titleright(tc,wt,z))
        return ppgeneric(uy,t,kwdict=kw)

    ############################################################################
    def ptuz(self,js=0,tc=None,wt=None,tp=None,z=None,**kw):
        """Plots time-uz for extraploated particles"""
        self.checkplotargs(kw)
        if self.ppmultispecies(self.ptuz,(js,tc,wt,tp,z),kw): return
        t = self.gett(js,tc,wt,tp,z)
        uz = self.getuz(js,tc,wt,tp)
        if 'pplimits' in kw.keys():
            kw['lframe'] = 1
        settitles("uz vs time","time","uz",self.titleright(tc,wt,z))
        return ppgeneric(uz,t,kwdict=kw)

    ############################################################################
    def ptvx(self,js=0,tc=None,wt=None,tp=None,z=None,**kw):
        """Plots time-Vx for extraploated particles"""
        self.checkplotargs(kw)
        if self.ppmultispecies(self.ptvx,(js,tc,wt,tp,z),kw): return
        t = self.gett(js,tc,wt,tp,z)
        vx = self.getvx(js,tc,wt,tp)
        if 'pplimits' in kw.keys():
            kw['lframe'] = 1
        settitles("Vx vs time","time","Vx",self.titleright(tc,wt,z))
        return ppgeneric(vx,t,kwdict=kw)

    ############################################################################
    def ptvy(self,js=0,tc=None,wt=None,tp=None,z=None,**kw):
        """Plots time-Vy for extraploated particles"""
        self.checkplotargs(kw)
        if self.ppmultispecies(self.ptvy,(js,tc,wt,tp,z),kw): return
        t = self.gett(js,tc,wt,tp,z)
        vy = self.getvy(js,tc,wt,tp)
        if 'pplimits' in kw.keys():
            kw['lframe'] = 1
        settitles("Vy vs time","time","Vy",self.titleright(tc,wt,z))
        return ppgeneric(vy,t,kwdict=kw)

    ############################################################################
    def ptvz(self,js=0,tc=None,wt=None,tp=None,z=None,**kw):
        """Plots time-Vz for extraploated particles"""
        self.checkplotargs(kw)
        if self.ppmultispecies(self.ptvz,(js,tc,wt,tp,z),kw): return
        t = self.gett(js,tc,wt,tp,z)
        vz = self.getvz(js,tc,wt,tp)
        if 'pplimits' in kw.keys():
            kw['lframe'] = 1
        settitles("Vz vs time","time","Vz",self.titleright(tc,wt,z))
        return ppgeneric(vz,t,kwdict=kw)

    ############################################################################
    def ptkez(self,js=0,tc=None,wt=None,tp=None,z=None,**kw):
        """Plots time-kinetic energy for extraploated particles"""
        self.checkplotargs(kw)
        if self.ppmultispecies(self.ptkez,(js,tc,wt,tp,z),kw): return
        t = self.gett(js,tc,wt,tp,z)
        vz = self.getvz(js,tc,wt,tp)
        kez = 0.5*top.pgroup.sm[js]*vz**2/jperev
        if 'pplimits' in kw.keys():
            kw['lframe'] = 1
        settitles("KEz vs time","time","KE (volts)",self.titleright(tc,wt,z))
        return ppgeneric(kez,t,kwdict=kw)

    ############################################################################
    def ptrace(self,js=0,tc=None,wt=None,tp=None,z=None,slope=0.,
               pplimits=None,**kw):
        """
Plots X-Y, X-X', Y-Y', Y'-X' in single page
 - slope=0.: slope subtracted from the angle. If 'auto', it is calculated
             automatically for the X-X' and Y-Y' plots.
 - pplimits=None: An optional list of up to four tuples, one for each phase
                  space plot. If any of the tuples are empty, the limits used
                  will be the usual ones for that plot.
        """
        self.checkplotargs(kw)
        if self.ppmultispecies(self.ptrace,(js,tc,wt,tp,z,slope,pplimits),kw):
            return
        x = self.getx(js,tc,wt,tp,z)
        y = self.gety(js,tc,wt,tp,z)
        xp = self.getxp(js,tc,wt,tp)
        yp = self.getyp(js,tc,wt,tp)
        titler = self.titleright(tc,wt,z)
        defaultpplimits = [(top.xplmin,top.xplmax,top.yplmin,top.yplmax),
                           (top.yplmin,top.yplmax,top.ypplmin,top.ypplmax),
                           (top.xplmin,top.xplmax,top.xpplmin,top.xpplmax),
                           (top.ypplmin,top.ypplmax,top.xpplmin,top.xpplmax)]
        if pplimits is None:
            pplimits = defaultpplimits
        else:
            kw['lframe'] = 1
            if type(pplimits[0]) != type(()):
                pplimits = 4*[pplimits]
            else:
                for i in range(4):
                    if i == len(pplimits): pplimits.append(defaultpplimits[i])
                    if not pplimits[i]: pplimits[i] = defaultpplimits[i]

        kw['view'] = 3
        kw['pplimits'] = pplimits[0]
        if type(slope)==type(''):
            kw['slope'] = 0.
        settitles("Y vs X","X","Y",titler)
        ppgeneric(y,x,kwdict=kw)

        kw['view'] = 4
        kw['pplimits'] = pplimits[1]
        if type(slope)==type(''):
            kw['slope'] = (ave(y*yp)-ave(y)*ave(yp))/dvnz(ave(y*y) - ave(y)**2)
        settitles("Y' vs Y","Y","Y'",titler)
        ppgeneric(yp,y,kwdict=kw)

        kw['view'] = 5
        kw['pplimits'] = pplimits[2]
        if type(slope)==type(''):
            kw['slope'] = (ave(x*xp)-ave(x)*ave(xp))/dvnz(ave(x*x) - ave(x)**2)
        settitles("X' vs X","X","X'",titler)
        ppgeneric(xp,x,kwdict=kw)

        kw['view'] = 6
        kw['pplimits'] = pplimits[3]
        if type(slope)==type(''):
            kw['slope'] = 0.
        settitles("X' vs Y'","Y'","X'",titler)
        ppgeneric(xp,yp,kwdict=kw)

