"""Diagnostics captured as particles cross grid cells.
"""
__all__ = ['GridCrossingDiags']
from warp import *
import cPickle

class GridCrossingDiags(object):
    """

Sets up diagnostics at z grid cells that are gathered from particles that
cross the cell.
  - js: species of particles to include. Currently can handle only a single
        species. Can be either the species index number, or an instance of the
        Species class.
  - zmmin,zmmax,dz,nz: grid parameters
        zmmin and zmmax default to w3d.zmmin and w3d.zmmax.
        dz defaults to w3d.dz/nzscale.
        nz defaults to nint((zmmax-zmmin)/dz) and dz will then be adjusted so
        that dz = (zmmax-zmmin)/nz.
        Note that the defaults are only calculated the first time that the
        diagnostic is done.
  - nzscale=1: multiplier on nz - makes it easy to use the w3d grid parameters
               but with a differing number of grid points.
  - nr,rmax: radial extent of radial profile diagnostic. Both must be given
             for the radial diagnostic to be done.
  - scintxmin,scintxmax,scintymin,scintymax,scintzmin,scintzmax,scintnx,scintny:
      specifies a volume where scintillator planes will be. All parameters
      must be specified. Notice that this will save the time history of a
      three-dimensional array and so can get very large. The size of the
      volume should be of minimal size. The resulting data will be stored
      in the scintillator attribute.
  - dumptofile=None: When given, the data will be written to a file with
                     the given name as a prefix.
  - starttime,endtime=None: If given, the data will be collected for times
                            only between the two values.
  - lmoving_frame=false: When true, the diagnostic moves with the beam frame.

The following quantities are calculated:
count: count of the number of particles that cross the cell each time
current: current, same as count but scaled by particle charge and 1/top.dt.
vzbar:
xbar, ybar:
xsqbar, ysqbar:
xrms, yrms:
rrms:
rprms:

Note that on the first time step, there is no old z data so the determination
if particles have crossed a grid cell can not be done so the results will
be unreliable.

    """

    def __init__(self,js,zmmin=None,zmmax=None,dz=None,nz=None,nzscale=1,
                 nhist=None,nr=None,rmax=None,ztarget=None,
                 scintxmin=None,scintxmax=None,
                 scintymin=None,scintymax=None,
                 scintzmin=None,scintzmax=None,
                 scintnx=None,
                 scintny=None,
                 dumptofile=None,
                 starttime=None,endtime=None,lmoving_frame=0):
        if isinstance(js,Species):
            self.js = js.jslist[0]
        else:
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

        self.scintxmin = scintxmin
        self.scintxmax = scintxmax
        self.scintymin = scintymin
        self.scintymax = scintymax
        self.scintzmin = scintzmin
        self.scintzmax = scintzmax
        self.scintnx = scintnx
        self.scintny = scintny

        self.dumptofile = dumptofile
        self.starttime = starttime
        self.endtime = endtime
        self.lmoving_frame = lmoving_frame

        self.zoldpid = nextpid()

        self.initializedata()
        self.enable()

    def initializedata(self):
        # --- Note to users: when retrieving results, the following attributes
        # --- should be accessed through their associated properties,
        # --- for example self.current, without the underscore.
        self._time = []
        self._zbeam = []
        self._count = []
        self._current = []
        self._vzbar = []
        self._xbar = []
        self._ybar = []
        self._xsqbar = []
        self._ysqbar = []
        self._xrms = []
        self._yrms = []
        self._rrms = []
        self._rprms = []

        self.ldoradialdiag = ((self.nr is not None) and
                              (self.rmax is not None))
        if self.ldoradialdiag:
            self._rprofile = []

        self.ldoscintillator = (self.scintzmin is not None)
        if self.ldoscintillator:
            self._scinttime = []
            self._scintillator = []

    def enable(self):
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
        assert abs((self.zmmax-self.zmmin)-self.nz*self.dz) < 1.e-5*self.dz,\
            "zmmin, zmmax, dz, and nz are not consistent with each other"

        if self.ldoscintillator:
            self.scintzmin = nint((self.scintzmin-self.zmmin)/self.dz)*self.dz + self.zmmin
            self.scintzmax = nint((self.scintzmax-self.zmmin)/self.dz)*self.dz + self.zmmin
            self.scintnz = nint((self.scintzmax-self.scintzmin)/self.dz)
            if self.scintxmin is None: self.scintxmin = -self.scintxmax
            if self.scintymin is None: self.scintymin = -self.scintymax
            self.scintdx = (self.scintxmax - self.scintxmin)/self.scintnx
            self.scintdy = (self.scintymax - self.scintymin)/self.scintny

    def appendnextarrays(self,zbeam):
        nz = self.nz
        self._time.append(top.time)
        self._zbeam.append(zbeam)
        self._count.append(zeros(1+nz,'d'))
        self._current.append(zeros(1+nz,'d'))
        self._vzbar.append(zeros(1+nz,'d'))
        self._xbar.append(zeros(1+nz,'d'))
        self._ybar.append(zeros(1+nz,'d'))
        self._xsqbar.append(zeros(1+nz,'d'))
        self._ysqbar.append(zeros(1+nz,'d'))
        self._xrms.append(None)
        self._yrms.append(None)
        self._rrms.append(None)
        self._rprms.append(zeros(1+nz,'d'))
        if self.ldoradialdiag:
            nr = self.nr
            self._rprofile.append(zeros((1+nr,1+nz),'d'))
        if self.ldoscintillator:
            if (len(self._scintillator) == 0 or
                maxnd(self._scintillator[-1]) > 0.):
                # --- Note that the data is only saved if it is nonzero
                scintnx = self.scintnx
                scintny = self.scintny
                scintnz = self.scintnz
                self._scintillator.append(zeros((1+scintnx,1+scintny,1+scintnz)))
                self._scinttime.append(top.time)

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
        if self.ldoradialdiag:
            dr = rmax/nr

        if self.ldoscintillator:
            scintnx = self.scintnx
            scintny = self.scintny
            scintnz = self.scintnz

        if self.lmoving_frame:
            zbeam = top.zbeam
        else:
            zbeam = 0.
        zoldpid = self.zoldpid

        # --- Create temporary work space
        gridcount = zeros(1+nz,'d')
        if self.ldoradialdiag:
            rprofilecount = zeros((1+nr,1+nz),'d')
            #self.rprofilemesh = iota(0,nr)*dr
        if self.ldoscintillator:
            scintillatorcount = zeros((1+scintnx,1+scintny,1+scintnz))

        if self.nhist is None: nhist = top.nhist
        else:                  nhist = self.nhist

        # --- The data is gathered from top.it-nhist/2 to top.it+nhist/2-1.
        # --- At the half way point, create space for the next set of data.
        if len(self._time) == 0 or (top.it-1)%nhist == int(nhist/2):

            if me == 0 and not self.dumptofile:
                self.appendnextarrays(zbeam)
            else:
                # --- On other processors or if the data is being dumped to a
                # --- file, just zero out the existing arrays.
                # --- There's no reason to keep the history on all processors.
                # --- The arrays are created the first time the diagnostic
                # --- is done.
                if len(self._time) == 0:
                    self.appendnextarrays(zbeam)
                else:
                    self._count[0].fill(0.)
                    self._current[0].fill(0.)
                    self._vzbar[0].fill(0.)
                    self._xbar[0].fill(0.)
                    self._ybar[0].fill(0.)
                    self._xsqbar[0].fill(0.)
                    self._ysqbar[0].fill(0.)
                    self._rprms[0].fill(0.)
                    if self.ldoradialdiag:
                        self._rprofile[0].fill(0.)
                    if self.ldoscintillator:
                        self._scintillator[0].fill(0.)

        # --- The code below is all local and can be skipped if there are no
        # --- particles locally.
        if top.pgroup.nps[js] > 0:

            xnew = getx(js=js,gather=0)
            ynew = gety(js=js,gather=0)
            rpnew = getrp(js=js,gather=0)
            znew = getz(js=js,gather=0)
            vznew = getvz(js=js,gather=0)
            zold = getpid(js=js,id=zoldpid-1,gather=0)

            iznew = int((znew - (zbeam + zmmin))/dz)
            izold = int((zold - (zbeam + zmmin))/dz)

            icrossed = (iznew > izold)

            izc = iznew[icrossed]
            xc = xnew[icrossed]
            yc = ynew[icrossed]
            rpc = rpnew[icrossed]
            vzc = vznew[icrossed]
            np = len(izc)

            if top.wpid > 0:
                weight = getpid(js=js,id=top.wpid-1,gather=0)
                ww = weight[icrossed]
            else:
                ww = ones(np,'d')

            zc = izc.astype('d')
            deposgrid1d(1,np,zc,ww,nz,self._count[-1],gridcount,0.,nz)
            deposgrid1d(1,np,zc,ww*vzc,nz,self._vzbar[-1],gridcount,0.,nz)
            deposgrid1d(1,np,zc,ww*xc,nz,self._xbar[-1],gridcount,0.,nz)
            deposgrid1d(1,np,zc,ww*yc,nz,self._ybar[-1],gridcount,0.,nz)
            deposgrid1d(1,np,zc,ww*xc**2,nz,self._xsqbar[-1],gridcount,0.,nz)
            deposgrid1d(1,np,zc,ww*yc**2,nz,self._ysqbar[-1],gridcount,0.,nz)
            deposgrid1d(1,np,zc,ww*rpc**2,nz,self._rprms[-1],gridcount,0.,nz)

            if self.ldoradialdiag or self.ldoscintillator:
                vz = getvz(js=js,gather=0)[icrossed]
                ke = 0.5*top.pgroup.sm[js]*vz**2
                ww *= top.pgroup.sw[js]

            if self.ldoradialdiag:
                rc = sqrt(xc**2 + yc**2)
                deposgrid2d(1,np,zc,rc,ke*ww,nz,nr,transpose(self._rprofile[-1]),
                            transpose(rprofilecount),0.,nz,0.,rmax)

            if self.ldoscintillator:
                izmin = (self.scintzmin - (zbeam + zmmin))/dz
                izmax = (self.scintzmax - (zbeam + zmmin))/dz
                deposgrid3d(1,np,zc,yc,xc,ke*ww,scintnz,scintny,scintnx,
                            transpose(self._scintillator[-1]),
                            transpose(scintillatorcount),
                            izmin,izmax,
                            self.scintymin,self.scintymax,
                            self.scintxmin,self.scintxmax)

            # --- Save particle z positions.
            i1 = top.pgroup.ins[js] - 1
            i2 = i1 + top.pgroup.nps[js]
            top.pgroup.pid[i1:i2,zoldpid-1] = top.pgroup.zp[i1:i2]

        # --- The data is gathered from top.it-nhist/2 to top.it+nhist/2-1.
        # --- At the half way point, finish the calculation by summing over
        # --- processors, dividing by the counts to get the averages, and
        # --- calculating the rms quantities.
        if top.it%nhist == int(nhist/2):
            co = self._count[-1]
            cu = self._current[-1]
            vzbar = self._vzbar[-1]
            xbar = self._xbar[-1]
            ybar = self._ybar[-1]
            xsqbar = self._xsqbar[-1]
            ysqbar = self._ysqbar[-1]
            rp = self._rprms[-1]

            # --- Finish the calculation, gathering data from all processors and
            # --- dividing out the count.
            co[...] = parallelsum(co)
            vzbar[...] = parallelsum(vzbar)
            xbar[...] = parallelsum(xbar)
            ybar[...] = parallelsum(ybar)
            xsqbar[...] = parallelsum(xsqbar)
            ysqbar[...] = parallelsum(ysqbar)
            rp[...] = parallelsum(rp)

            cotemp = where(co==0.,1.,co)
            vzbar[...] = vzbar/cotemp
            xbar[...] = xbar/cotemp
            ybar[...] = ybar/cotemp
            xsqbar[...] = xsqbar/cotemp
            ysqbar[...] = ysqbar/cotemp
            rp[...] = sqrt(rp/cotemp)

            self._xrms[-1] = sqrt(abs(xsqbar - xbar**2))
            self._yrms[-1] = sqrt(abs(ysqbar - ybar**2))
            self._rrms[-1] = sqrt(abs(xsqbar + ysqbar - xbar**2 - ybar**2))

            # --- Scale the current appropriately.
            cu[...] = co*(top.pgroup.sq[js]*top.pgroup.sw[js]/(top.dt*nhist))

            if self.ldoradialdiag:
                rprof = self._rprofile[-1]
                rprof[...] = parallelsum(rprof)

            if self.ldoscintillator:
                scint = self._scintillator[-1]
                scint[...] = parallelsum(scint)

            if self.dumptofile: self.dodumptofile(zbeam)

    # ----------------------------------------------------------------------
    def dodumptofile(self,zbeam):
        #self.dodumptofilePDB(zbeam)
        self.dodumptofilePickle(zbeam)

    def dodumptofilePDB(self,zbeam):
        if me != 0: return
        ff = PW.PW(self.dumptofile+'_gridcrossing.pdb','a',verbose=0)
        suffix = "_%08d"%(top.it)
        ff.write('time'+suffix,top.time)
        ff.write('zbeam'+suffix,zbeam)
        ff.write('count'+suffix,self._count[0])
        ff.write('current'+suffix,self._current[0])
        ff.write('vzbar'+suffix,self._vzbar[0])
        ff.write('xbar'+suffix,self._xbar[0])
        ff.write('ybar'+suffix,self._ybar[0])
        ff.write('xsqbar'+suffix,self._xsqbar[0])
        ff.write('ysqbar'+suffix,self._ysqbar[0])
        ff.write('xrms'+suffix,self._xrms[0])
        ff.write('yrms'+suffix,self._yrms[0])
        ff.write('rrms'+suffix,self._rrms[0])
        ff.write('rprms'+suffix,self._rprms[0])
        if self.ldoradialdiag:
            ff.write('rprofile'+suffix,self._rprofile[0])
        if self.ldoscintillator:
            if maxnd(self._scintillator[0]) > 0.:
                # --- Note that the data is only saved if it is nonzero
                ff.write('scinttime'+suffix,top.time)
                ff.write('scintillator'+suffix,self._scintillator[0])
        ff.close()

    def dodumptofilePickle(self,zbeam):
        if me != 0: return
        if not os.path.exists(self.dumptofile+'_gridcrossing.pkl'):
            ff = open(self.dumptofile+'_gridcrossing.pkl','w')
            # --- Save the input parameters to the file.
            cPickle.dump(('js',self.js),ff,-1)
            cPickle.dump(('zmmin',self.zmmin),ff,-1)
            cPickle.dump(('zmmax',self.zmmax),ff,-1)
            cPickle.dump(('dz',self.dz),ff,-1)
            cPickle.dump(('nz',self.nz),ff,-1)
            cPickle.dump(('nzscale',self.nzscale),ff,-1)
            cPickle.dump(('nhist',self.nhist),ff,-1)
            cPickle.dump(('nr',self.nr),ff,-1)
            cPickle.dump(('rmax',self.rmax),ff,-1)
            cPickle.dump(('ztarget',self.ztarget),ff,-1)
            cPickle.dump(('dumptofile',self.dumptofile),ff,-1)
            cPickle.dump(('starttime',self.starttime),ff,-1)
            cPickle.dump(('endtime',self.endtime),ff,-1)
            cPickle.dump(('ldoradialdiag',self.ldoradialdiag),ff,-1)
            cPickle.dump(('ldoscintillator',self.ldoscintillator),ff,-1)
            if self.ldoscintillator:
                cPickle.dump(('scintxmin',self.scintxmin),ff,-1)
                cPickle.dump(('scintxmax',self.scintxmax),ff,-1)
                cPickle.dump(('scintymin',self.scintymin),ff,-1)
                cPickle.dump(('scintymax',self.scintymax),ff,-1)
                cPickle.dump(('scintzmin',self.scintzmin),ff,-1)
                cPickle.dump(('scintzmax',self.scintzmax),ff,-1)
                cPickle.dump(('scintnx',self.scintnx),ff,-1)
                cPickle.dump(('scintny',self.scintny),ff,-1)
                cPickle.dump(('scintnz',self.scintnz),ff,-1)
                cPickle.dump(('scintdx',self.scintdx),ff,-1)
                cPickle.dump(('scintdy',self.scintdy),ff,-1)
        else:
            ff = open(self.dumptofile+'_gridcrossing.pkl','a')
        suffix = "_%08d"%(top.it)
        cPickle.dump(('time'+suffix,top.time),ff,-1)
        cPickle.dump(('zbeam'+suffix,zbeam),ff,-1)
        cPickle.dump(('count'+suffix,self._count[0]),ff,-1)
        cPickle.dump(('current'+suffix,self._current[0]),ff,-1)
        cPickle.dump(('vzbar'+suffix,self._vzbar[0]),ff,-1)
        cPickle.dump(('xbar'+suffix,self._xbar[0]),ff,-1)
        cPickle.dump(('ybar'+suffix,self._ybar[0]),ff,-1)
        cPickle.dump(('xsqbar'+suffix,self._xsqbar[0]),ff,-1)
        cPickle.dump(('ysqbar'+suffix,self._ysqbar[0]),ff,-1)
        cPickle.dump(('xrms'+suffix,self._xrms[0]),ff,-1)
        cPickle.dump(('yrms'+suffix,self._yrms[0]),ff,-1)
        cPickle.dump(('rrms'+suffix,self._rrms[0]),ff,-1)
        cPickle.dump(('rprms'+suffix,self._rprms[0]),ff,-1)
        if self.ldoradialdiag:
            cPickle.dump(('rprofile'+suffix,self._rprofile[0]),ff,-1)
        if self.ldoscintillator:
            if maxnd(self._scintillator[0]) > 0.:
                # --- Note that the data is only saved if it is nonzero
                cPickle.dump(('scinttime'+suffix,top.time),ff,-1)
                cPickle.dump(('scintillator'+suffix,self._scintillator[0]),ff,-1)
        ff.close()

    def restorefromfile(self,files=[],readscintillator=1):
        #self.restorefromfilePDB(files,readscintillator)
        self.restorefromfilePickle(files,readscintillator=readscintillator)

    def restorefromfilePDB(self,files=[],readscintillator=1):
        if me != 0: return
        ff = PR.PR(self.dumptofile+'_gridcrossing.pdb')

        self._time = []
        self._zbeam = []
        self._count = []
        self._current = []
        self._vzbar = []
        self._xbar = []
        self._ybar = []
        self._xsqbar = []
        self._ysqbar = []
        self._xrms = []
        self._yrms = []
        self._rrms = []
        self._rprms = []
        # --- At this point, getdiagnostics may not have been executed, so
        # --- self.ldoradialdiag may not be set. So assume that it is and
        # --- create the rprofile list.
        self._rprofile = []
        self._scintillator = []

        varlist = list(ff.inquire_names())
        varlist.sort()
        for var in varlist:
            if var[0] == 't':
                name,it = string.split(var,'_')
                suffix = "_%d"%(it)
                self._time.append(ff.read('time'+suffix))
                self._zbeam.append(ff.read('zbeam'+suffix))
                self._count.append(ff.read('count'+suffix))
                self._current.append(ff.read('current'+suffix))
                self._vzbar.append(ff.read('vzbar'+suffix))
                self._xbar.append(ff.read('xbar'+suffix))
                self._ybar.append(ff.read('ybar'+suffix))
                self._xsqbar.append(ff.read('xsqbar'+suffix))
                self._ysqbar.append(ff.read('ysqbar'+suffix))
                self._xrms.append(ff.read('xrms'+suffix))
                self._yrms.append(ff.read('yrms'+suffix))
                self._rrms.append(ff.read('rrms'+suffix))
                self._rprms.append(ff.read('rprms'+suffix))
                try:
                    self._rprofile.append(ff.read('rprofile'+suffix))
                except:
                    # --- This just means that there is no rprofile data
                    pass
                try:
                    self._scinttime.append(ff.read('scinttime'+suffix))
                    self._scintillator.append(ff.read('scintillator'+suffix))
                except:
                    # --- This just means that there is no scintillator data
                    pass

        ff.close()

        # --- If there is no rprofile data, then delete the attribute
        if len(self._rprofile) == 0:
            del self._rprofile
        if len(self._scintillator) == 0:
            del self._scintillator

    def restorefromfilePickle(self,files=[],
                              starttime=-largepos,endtime=+largepos,
                              readscintillator=1):
        if me != 0: return

        if not isinstance(files,ListType):
            files = list([files])
        if len(files) == 0:
            files = [self.dumptofile+'_gridcrossing.pkl']

        # --- First, read in the input parameters, if they were saved.
        # --- This reads in everything at the beginning of the file until
        # --- the time data is found, which starts the data section of the
        # --- file.
        ff = open(files[0],'r')
        data = cPickle.load(ff)
        while data[0][0:4] != 'time':
            setattr(self,data[0],data[1])
            data = cPickle.load(ff)
        ff.close()

        # --- Read all of the data in. Only save the data if the time is
        # --- between start and endtime.
        savedata = 0
        datadict = {}
        for file in files:
            ff = open(file,'r')
            while 1:
                try:
                    tell = ff.tell()
                    data = cPickle.load(ff)
                except:
                    break
                if data[0][:4] == 'time':
                    if starttime <= data[1] <= endtime:
                        savedata = 1
                if not readscintillator and data[0][:12] == 'scintillator':
                    data = (data[0],tell)
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

        self._time = []
        self._zbeam = []
        self._count = []
        self._current = []
        self._vzbar = []
        self._xbar = []
        self._ybar = []
        self._xsqbar = []
        self._ysqbar = []
        self._xrms = []
        self._yrms = []
        self._rrms = []
        self._rprms = []
        # --- At this point, getdiagnostics may not have been executed, so
        # --- self.ldoradialdiag may not be set. So assume that it is and
        # --- create the rprofile list.
        self._rprofile = []
        self._scintillator = []

        varlist = datadict.keys()
        varlist.sort()
        for var in varlist:
            if var[0:4] == 'time':
                name,it = string.split(var,'_')
                suffix = "_%s"%(it)
                self._time.append(datadict['time'+suffix])
                self._zbeam.append(datadict['zbeam'+suffix])
                self._count.append(datadict['count'+suffix])
                self._current.append(datadict['current'+suffix])
                self._vzbar.append(datadict['vzbar'+suffix])
                self._xbar.append(datadict['xbar'+suffix])
                self._ybar.append(datadict['ybar'+suffix])
                self._xsqbar.append(datadict['xsqbar'+suffix])
                self._ysqbar.append(datadict['ysqbar'+suffix])
                self._xrms.append(datadict['xrms'+suffix])
                self._yrms.append(datadict['yrms'+suffix])
                self._rrms.append(datadict['rrms'+suffix])
                self._rprms.append(datadict['rprms'+suffix])
                try:
                    self._rprofile.append(datadict['rprofile'+suffix])
                except:
                    # --- This just means that there is no rprofile data
                    pass
                try:
                    self._scinttime.append(datadict['scinttime'+suffix])
                    self._scintillator.append(datadict['scintillator'+suffix])
                except:
                    # --- This just means that there is no scintillator data
                    pass

        # --- If there is no rprofile data, then delete the attribute
        if len(self._rprofile) == 0:
            del self._rprofile
        if len(self._scintillator) == 0:
            del self._scintillator

    def readscintillator(self,i,file=None):
        if file is None:
            file = self.dumptofile+'_gridcrossing.pkl'

        ff = open(file,'r')
        ff.seek(self._scintillator[i])
        data = cPickle.load(ff)
        ff.close()
        return data[1]

    # ----------------------------------------------------------------------
    def setupanalysis(self):
        self.arraytime = array(self.time)
        self.arraycurrent = array(self.current)
        self.arrayradius = array(self.rrms)
        self.zmesh = self.zmmin + arange(0,self.nz+1,dtype='l')*self.dz

        self.currentmax = zeros(self.nz+1,'d')
        self.ratcurrentmax = zeros(self.nz+1,'d')
        for iz in range(self.nz+1):
            # --- Find the max current over time at the location iz
            ii = argmax(self.arraycurrent[:,iz])
            # --- Save the current and beam radius at that time
            self.currentmax[iz] = self.arraycurrent[ii,iz]
            self.ratcurrentmax[iz] = self.arrayradius[ii,iz]*100.

        if self.ldoradialdiag:
            self.arrayrprofile = array(self.rprofile)
            dr = self.rmax/self.nr
            self.rprofilemesh = iota(0,self.nr)*dr
            aa = pi*2.*self.rprofilemesh*dr # --- Is this correct???
            aa[0] = pi*0.25*dr**2
            aa *= 10000.
            self.aa = aa

    def saveresults(self,filename):
        ff = PW.PW(filename)
        ff.zmesh = self.zmesh
        ff.currentmax = self.currentmax
        ff.ratcurrentmax = self.ratcurrentmax
        if self.ldoradialdiag:
            ff.aa = self.aa
            ff.Esum = self.Esum
            ff.Etot = self.Etot
            ff.rprofilemesh = self.rprofilemesh
        ff.close()

    def ppcurrmax(self):
        plp(self.currentmax,self.zmesh,msize=3)
        plp(self.ratcurrentmax,self.zmesh,color=blue,msize=3)
        ptitles('spot size','Z (m)','Current (Amps)',
                'Black is peak current, Blue is corresponding radius')

    def ppfluence(self,Esum):
        """Plots the fluence as a function of radius"""
        Etot = sum(Esum)
        self.Esum = Esum
        self.Etot = Etot
        # --- Plot the energy density versus radius,
        # --- summed over the time window.
        plg(Esum/self.aa,self.rprofilemesh*100)
        ptitles('%d KeV'%ee,'R (cm)','joules/sq-cm','Energy deposition on target, summed over 5 ns')
        plt("Etot = %7.2f mJ"%(Etot*1000.),.45,.82)

    def ppfluenceattarget(self,ztarget,deltat):
        """Plot the fluence on the target, integrating over the time +/- deltat
around the peak current."""
        iztarget = int((ztarget - self.zmmin)/self.dz)
        ii = argmax(self.arraycurrent[:,iztarget])
        di = int(deltat/top.dt/self.nhist)
        Esum = sum(self.arrayrprofile[ii-di:ii+di,:,iztarget],0)
        self.ppfluence(Esum)

    def ppfluenceatspot(self,deltat=None,currmin=None,tslice=slice(None)):
        iztarget = argmin(ratcurrentmax[tslice])
        if deltat is not None:
            ii = argmax(self.arraycurrent[:,iztarget])
            di = int(deltat/top.dt/self.nhist)
            Esum = sum(self.arrayrprofile[ii,:,iztarget],0)
        elif currmin is not None:
            ii = (gridcurrent[:,iztarget] > currmin)
            Esum = sum(self.arrayrprofile[ii-di:ii+di,:,iztarget],0)
        self.ppfluence(Esum)

    # ----------------------------------------------------------------------
    def _pp2d(self,data,lbeamframe=1,**kw):
        zmesh = self.zmmin + arange(0,self.nz+1,dtype='l')*self.dz
        if lbeamframe:
            zz = zmesh[:,newaxis]*ones(data.shape[0])[newaxis,:]
        else:
            zz = zmesh[:,newaxis] + self.zbeam[newaxis,:]
        tt = self.time[newaxis,:]*ones(self.nz+1)[:,newaxis]
        ppgeneric(gridt=data,xmesh=zz,ymesh=tt,**kw)

    def pp2dcount(self,**kw):
        self._pp2d(self.count,**kw)
    def pp2dcurrent(self,**kw):
        self._pp2d(self.current,**kw)
    def pp2dvzbar(self,**kw):
        self._pp2d(self.vzbar,**kw)
    def pp2dxbar(self,**kw):
        self._pp2d(self.xbar,**kw)
    def pp2dybar(self,**kw):
        self._pp2d(self.ybar,**kw)
    def pp2dxsqbar(self,**kw):
        self._pp2d(self.xsqbar,**kw)
    def pp2dysqbar(self,**kw):
        self._pp2d(self.ysqbar,**kw)
    def pp2dxrms(self,**kw):
        self._pp2d(self.xrms,**kw)
    def pp2dyrms(self,**kw):
        self._pp2d(self.yrms,**kw)
    def pp2drrms(self,**kw):
        self._pp2d(self.rrms,**kw)
    def pp2drprms(self,**kw):
        self._pp2d(self.rprms,**kw)

    # ----------------------------------------------------------------------
    def _gettimehistory(self,data,z):
        d = []
        t = []
        for i in range(data.shape[0]):
            if self.zmmin <= z-self.zbeam[i] <= self.zmmax:
                t.append(self.time[i])
                zz = (z - self.zbeam[i] - self.zmmin)/self.dz
                iz = int(zz)
                wz = zz - iz
                if iz < self.nz:
                    d.append(data[i,iz]*(1. - wz) + data[i,iz+1]*wz)
                else:
                    d.append(data[i,iz])
        return array(d),array(t)
        
    def hcount(self,**kw):
        return self._gettimehistory(self.count,**kw)
    def hcurrent(self,**kw):
        return self._gettimehistory(self.current,**kw)
    def hvzbar(self,**kw):
        return self._gettimehistory(self.vzbar,**kw)
    def hxbar(self,**kw):
        return self._gettimehistory(self.xbar,**kw)
    def hybar(self,**kw):
        return self._gettimehistory(self.ybar,**kw)
    def hxsqbar(self,**kw):
        return self._gettimehistory(self.xsqbar,**kw)
    def hysqbar(self,**kw):
        return self._gettimehistory(self.ysqbar,**kw)
    def hxrms(self,**kw):
        return self._gettimehistory(self.xrms,**kw)
    def hyrms(self,**kw):
        return self._gettimehistory(self.yrms,**kw)
    def hrrms(self,**kw):
        return self._gettimehistory(self.rrms,**kw)
    def hrprms(self,**kw):
        return self._gettimehistory(self.rprms,**kw)

    # ----------------------------------------------------------------------
    def _timeintegrate(self,data,laverage):

        try:
            self.zmesh
        except AttributeError:
            self.zmesh = self.zmmin + arange(0,self.nz+1)*self.dz

        zmin = self.zmesh[0] + self.zbeam.min()
        zmax = self.zmesh[-1] + self.zbeam.max()
        nz = nint((zmax - zmin)/self.dz)
        dz = (zmax - zmin)/nz

        grid = zeros(1+nz,'d')
        gridcount = zeros(1+nz,'d')
        gridmesh = zmin + arange(nz+1)*dz

        if laverage:
            count = self.count

        for i in range(data.shape[0]):
            if laverage:
                deposgrid1dw(1,data.shape[1],
                            self.zmesh+self.zbeam[i],
                            data[i,:],
                            count[i,:],
                            nz,grid,gridcount,zmin,zmax)
            else:
                deposgrid1d(1,data.shape[1],
                            self.zmesh+self.zbeam[i],
                            data[i,:],
                            nz,grid,gridcount,zmin,zmax)

        if laverage:
            result = grid/where(gridcount > 0.,gridcount,1.)
        else:
            result = grid

        return result,gridmesh

    def timeintegratedcount(self):
        return self._timeintegrate(self.count,laverage=0)

    def timeintegratedcurrent(self):
        return self._timeintegrate(self.current,laverage=0)

    def timeintegratedvzbar(self):
        return self._timeintegrate(self.vzbar,laverage=1)

    def timeintegratedxbar(self):
        return self._timeintegrate(self.xbar,laverage=1)

    def timeintegratedybar(self):
        return self._timeintegrate(self.ybar,laverage=1)

    def timeintegratedxsqbar(self):
        return self._timeintegrate(self.xsqbar,laverage=1)

    def timeintegratedysqbar(self):
        return self._timeintegrate(self.ysqbar,laverage=1)

    def timeintegratedxrms(self):
        data = self.xrms**2
        result,gridmesh = self._timeintegrate(data,laverage=1)
        result = sqrt(maximum(0.,result))
        return result,gridmesh

    def timeintegratedyrms(self):
        data = self.yrms**2
        result,gridmesh = self._timeintegrate(data,laverage=1)
        result = sqrt(maximum(0.,result))
        return result,gridmesh

    def timeintegratedxprms(self):
        data = self.xprms**2
        result,gridmesh = self._timeintegrate(data,laverage=1)
        result = sqrt(maximum(0.,result))
        return result,gridmesh

    def timeintegratedyprms(self):
        data = self.yprms**2
        result,gridmesh = self._timeintegrate(data,laverage=1)
        result = sqrt(maximum(0.,result))
        return result,gridmesh

    def timeintegratedcorkscrew(self):
        xbarint,gridmesh = self.timeintegratedxbar()
        ybarint,gridmesh = self.timeintegratedybar()
        xbarsqint,gridmesh = self._timeintegrate(self.xbar**2,laverage=1)
        ybarsqint,gridmesh = self._timeintegrate(self.ybar**2,laverage=1)
        corkscrew = sqrt(maximum(0.,xbarsqint - xbarint**2 + ybarsqint - ybarint**2))
        return corkscrew,gridmesh

    # ----------------------------------------------------------------------
    # --- Setup the properties so that the last set of data which is
    # --- still being accumulated is not returned, and so that the
    # --- data is converted to an array.
    def _setupproperty(name,doc=None):
        def fget(self):
            if self.nhist is None: nhist = top.nhist
            else:                  nhist = self.nhist
            # --- Get the data, removing the last element if the accumulation
            # --- of the data is not complete.
            result = getattr(self,'_'+name)
            if top.it%nhist != int(nhist/2):
                result = result[:-1]

            # --- Check if there is a cached array.
            # --- If so, and if it is the same size as reult, then return it,
            # --- otherwise convert result to an array and return it.
            cache = getattr(self,'_cache'+name,None)
            if cache is not None and len(cache) == len(result):
                result = cache
            else:
                try:
                    result = array(result)
                except ValueError:
                    # --- This can happen if self.nz changed at some point,
                    # --- which changed the length of the new data so that
                    # --- all of the elements do not have the same length.
                    pass
                setattr(self,'_cache'+name,result)
            return result
        return fget,None,None,doc

    time = property(*_setupproperty('time'))
    zbeam = property(*_setupproperty('zbeam'))
    count = property(*_setupproperty('count'))
    current = property(*_setupproperty('current'))
    vzbar = property(*_setupproperty('vzbar'))
    xbar = property(*_setupproperty('xbar'))
    ybar = property(*_setupproperty('ybar'))
    xsqbar = property(*_setupproperty('xsqbar'))
    ysqbar = property(*_setupproperty('ysqbar'))
    xrms = property(*_setupproperty('xrms'))
    yrms = property(*_setupproperty('yrms'))
    rrms = property(*_setupproperty('rrms'))
    rprms = property(*_setupproperty('rprms'))
    rprofile = property(*_setupproperty('rprofile'))
    scinttime = property(*_setupproperty('scinttime'))
    scintillator = property(*_setupproperty('scintillator'))
    del _setupproperty

