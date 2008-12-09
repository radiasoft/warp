"""Diagnostics captured as particles cross grid cells.
"""
__all__ = ['GridCrossingDiags']
from warp import *

class GridCrossingDiags(object):
    """

Sets up diagnostics at z grid cells that are gathered from particles that
cross the cell.
  - js: species of particles to include. Currently can handle only a single
        species.
  - zmmin,zmmax,dz,nz: grid parameters, all taken from w3d if not supplied.
  - nzscale=1: multiplier on nz - makes it easy to use the w3d grid parameters
               but with a differing number of grid points.
  - nr,rmax: radial extent of radial profile diagnostic. If not given,
             the radial diagnostic is not done.
  - ztarget: Z location with the radial profile is calculated.
  - dumptofile=None: When given, the data will be written to a file with
                     the given name as a prefix.
  - starttime,endtime=None: If given, the data will be collected for times
                            only between the two values.

The following quantities are calculated:
count: count of the number of particles that cross the cell each time
current: current, same as count but scaled by particle charge and 1/top.dt.
rrms:
rprms:

    """

    def __init__(self,js,zmmin=None,zmmax=None,dz=None,nz=None,nzscale=1,
                 nhist=None,nr=None,rmax=None,ztarget=None,dumptofile=None,
                 starttime=None,endtime=None):
        self.js = js
        self.zmmin = zmmin
        self.zmmax = zmmax
        self.dz = dz
        self.nz = nz
        self.nzscale = nzscale
        self.nhist = nhist
        self.nr = nr
        self.rmax = rmax
        self.ztarget = ztarget
        self.dumptofile = dumptofile
        self.starttime = starttime
        self.endtime = endtime

        self.zoldpid = nextpid()

        self.count = None
        self.current = None
        self.rrms = None
        self.rprms = None

        installafterstep(self.getdiagnostics)

    def disable(self):
        uninstallafterstep(self.getdiagnostics)

    def initializegrid(self):
        # --- Initialize grid parameters if needed. This is done here
        # --- in case the grid had not been setup yet when init was called.
        if self.zmmin is None: self.zmmin = w3d.zmmin
        if self.zmmax is None: self.zmmax = w3d.zmmax
        if self.dz is None: self.dz = w3d.dz/self.nzscale
        if self.nz is None:
            self.nz = nint((self.zmmax - self.zmmin)/self.dz)
            self.dz = (self.zmmax - self.zmmin)/self.nz

        self.ldoradialdiag = (self.nr is not None)

    def getdiagnostics(self):

        # --- Check if particle was advanced
        if not top.pgroup.ldts[self.js]: return

        # --- Check the start and end times
        if self.starttime is not None:
            if top.time < self.starttime: return
        if self.endtime is not None:
            if top.time > self.endtime: return

        self.initializegrid()

        # --- Create handy locals.
        js = self.js
        zmmin = self.zmmin
        zmmax = self.zmmax
        dz = self.dz
        nz = self.nz

        # --- Do some error checking
        if zmmax < zmmin: return

        rmax = self.rmax
        nr = self.nr
        ldoradialdiag = self.ldoradialdiag
        if ldoradialdiag:
            dr = rmax/nr

        zbeam = top.zbeam
        zoldpid = self.zoldpid

        # --- Initialize the data lists.
        if self.count is None:
            self.time = []
            self.count = [zeros(1+nz,'d')]
            self.current = [zeros(1+nz,'d')]
            self.rrms = [zeros(1+nz,'d')]
            self.rprms = [zeros(1+nz,'d')]
            self.dummy = zeros(1+nz,'d')
            if ldoradialdiag:
                self.rprofile = [zeros((1+nr,1+nz),'d')]
                self.rprofilecount = zeros((1+nr,1+nz),'d')
                #self.rprofilemesh = iota(0,nr)*dr

        if self.nhist is None: nhist = top.nhist
        else:                  nhist = self.nhist

        # --- The data is gathered from top.it-nhist/2 to top.it+nhist/2-1.
        # --- This resets things at the half way point.
        if top.it%nhist == int(nhist/2):
            co = self.count[-1]
            cu = self.current[-1]
            rr = self.rrms[-1]
            rp = self.rprms[-1]
            if ldoradialdiag:
                rprof = self.rprofile[-1]

            # --- Finish the calculation, gathering data from all processors and
            # --- dividing out the count.
            co[...] = parallelsum(co)
            rr[...] = parallelsum(rr)
            rr[...] = sqrt(rr/where(co==0.,1.,co))
            rp[...] = parallelsum(rp)
            rp[...] = sqrt(rp/where(co==0.,1.,co))

            if ldoradialdiag:
                rprof[...] = parallelsum(rprof)

            # --- Scale the current appropriately.
            cu[...] = co*(top.pgroup.sq[js]*top.pgroup.sw[js]/(top.dt*nhist))

            if self.dumptofile: self.dodumptofile()

            # --- Create space for the next set of data.
            if me == 0 and not self.dumptofile:
                self.time.append(top.time)
                self.count.append(zeros(1+nz,'d'))
                self.current.append(zeros(1+nz,'d'))
                self.rrms.append(zeros(1+nz,'d'))
                self.rprms.append(zeros(1+nz,'d'))
                if ldoradialdiag:
                    self.rprofile.append(zeros((1+nr,1+nz),'d'))
            else:
                # --- On other processors or if the data is being dumped to a
                # --- file, just zero out the existing arrays.
                # --- There's no reason to keep the history on all processors.
                self.count[0][:] = 0.
                self.current[0][:] = 0.
                self.rrms[0][:] = 0.
                self.rprms[0][:] = 0.
                if ldoradialdiag:
                    self.rprofile[0][...] = 0.

        # --- The code below is all local and can be skipped if there are no
        # --- particles locally.
        if top.pgroup.nps[js] == 0: return

        rnew = getr(js=js,gather=0)
        rpnew = getrp(js=js,gather=0)
        znew = getz(js=js,gather=0)
        zold = getpid(js=js,id=zoldpid-1,gather=0)

        iznew = int((znew - (zbeam + zmmin))/dz)
        izold = int((zold - (zbeam + zmmin))/dz)

        icrossed = (iznew > izold)

        izc = compress(icrossed,iznew)
        rc = compress(icrossed,rnew)
        rpc = compress(icrossed,rpnew)
        np = len(izc)

        if top.wpid > 0:
            weight = getpid(js=js,id=top.wpid-1,gather=0)
            ww = compress(icrossed,weight)
        else:
            ww = ones(np,'d')

        zc = izc.astype('d')
        deposgrid1d(1,np,zc,ww,nz,self.count[-1],self.dummy,0.,nz)
        deposgrid1d(1,np,zc,ww*rc**2,nz,self.rrms[-1],self.dummy,0.,nz)
        deposgrid1d(1,np,zc,ww*rpc**2,nz,self.rprms[-1],self.dummy,0.,nz)

        if ldoradialdiag:
            vz = compress(icrossed,getvz(gather=0))
            ke = 0.5*top.pgroup.sm[js]*vz**2
            if top.wpid > 0:
                weight = getpid(js=js,id=top.wpid-1,gather=0)
                ww = compress(icrossed,weight)*top.pgroup.sw[js]
            else:
                ww = top.pgroup.sw[js]
            deposgrid2d(1,np,zc,rc,ke*ww,nz,nr,transpose(self.rprofile[-1]),
                        transpose(self.rprofilecount),0.,nz,0.,rmax)

        # --- Save particle z positions.
        i1 = top.pgroup.ins[js] - 1
        i2 = i1 + top.pgroup.nps[js]
        top.pgroup.pid[i1:i2,zoldpid-1] = top.pgroup.zp[i1:i2]

    def dodumptofile(self):
        #self.dodumptofilePDB()
        self.dodumptofilePickle()

    def dodumptofilePDB(self):
        if me != 0: return
        ff = PW.PW(self.dumptofile+'_gridcrossing.pdb','a',verbose=0)
        suffix = "_%08d"%(top.it)
        ff.write('time'+suffix,top.time)
        ff.write('count'+suffix,self.count[0])
        ff.write('current'+suffix,self.current[0])
        ff.write('rrms'+suffix,self.rrms[0])
        ff.write('rprms'+suffix,self.rprms[0])
        if self.ldoradialdiag:
            ff.write('rprofile'+suffix,self.rprofile[0])
        ff.close()

    def dodumptofilePickle(self):
        if me != 0: return
        import cPickle
        ff = open(self.dumptofile+'_gridcrossing.pkl','a')
        suffix = "_%08d"%(top.it)
        cPickle.dump(('time'+suffix,top.time),ff,-1)
        cPickle.dump(('count'+suffix,self.count[0]),ff,-1)
        cPickle.dump(('current'+suffix,self.current[0]),ff,-1)
        cPickle.dump(('rrms'+suffix,self.rrms[0]),ff,-1)
        cPickle.dump(('rprms'+suffix,self.rprms[0]),ff,-1)
        if self.ldoradialdiag:
            cPickle.dump(('rprofile'+suffix,self.rprofile[0]),ff,-1)
        ff.close()

    def restorefromfile(self,files=[]):
        #self.restorefromfilePDB(files)
        self.restorefromfilePickle(files)

    def restorefromfilePDB(self,files=[]):
        if me != 0: return
        ff = PR.PR(self.dumptofile+'_gridcrossing.pdb')

        self.time = []
        self.count = []
        self.current = []
        self.rrms = []
        self.rprms = []
        # --- At this point, getdiagnostics may not have been executed, so
        # --- self.ldoradialdiag may not be set. So assume that it is and
        # --- create the rprofile list.
        self.rprofile = []

        varlist = list(ff.inquire_names())
        varlist.sort()
        for var in varlist:
            if var[0] == 't':
                name,it = string.split(var,'_')
                suffix = "_%d"%(it)
                self.time.append(ff.read('time'+suffix))
                self.count.append(ff.read('count'+suffix))
                self.current.append(ff.read('current'+suffix))
                self.rrms.append(ff.read('rrms'+suffix))
                self.rprms.append(ff.read('rprms'+suffix))
                try:
                    self.rprofile.append(ff.read('rprofile'+suffix))
                except:
                    # --- This just means that there is no rprofile data
                    pass

        ff.close()

        # --- If there is no rprofile data, then delete the attribute
        if len(self.rprofile) == 0:
            del self.rprofile

    def restorefromfilePickle(self,files=[],
                              starttime=-largepos,endtime=+largepos):
        if me != 0: return

        if not isinstance(files,ListType):
            files = list([files])
        if len(files) == 0:
            files = [self.dumptofile+'_gridcrossing.pkl']

        # --- Read all of the data in. Only save the data if the time is
        # --- between start and endtime.
        import cPickle
        savedata = 0
        datadict = {}
        for file in files:
            ff = open(file,'r')
            while 1:
                try:
                    data = cPickle.load(ff)
                except:
                    break
                if data[0][:4] == 'time':
                    if starttime <= data[1] <= endtime:
                        savedata = 1
                if savedata:
                    datadict[data[0]] = data[1]
            ff.close()

        # --- Fix old bad naming
        varlist = datadict.keys()
        for var in varlist:
            name,it = string.split(var,'_')
            if len(it) < 8:
                 newname = name + '_' + (8-len(it))*'0' + it
                 datadict[newname] = datadict[var]
                 del datadict[var]

        self.time = []
        self.count = []
        self.current = []
        self.rrms = []
        self.rprms = []
        # --- At this point, getdiagnostics may not have been executed, so
        # --- self.ldoradialdiag may not be set. So assume that it is and
        # --- create the rprofile list.
        self.rprofile = []

        varlist = datadict.keys()
        varlist.sort()
        for var in varlist:
            if var[0] == 't':
                name,it = string.split(var,'_')
                suffix = "_%s"%(it)
                self.time.append(datadict['time'+suffix])
                self.count.append(datadict['count'+suffix])
                self.current.append(datadict['current'+suffix])
                self.rrms.append(datadict['rrms'+suffix])
                self.rprms.append(datadict['rprms'+suffix])
                try:
                    self.rprofile.append(datadict['rprofile'+suffix])
                except:
                    # --- This just means that there is no rprofile data
                    pass

        # --- If there is no rprofile data, then delete the attribute
        if len(self.rprofile) == 0:
            del self.rprofile

