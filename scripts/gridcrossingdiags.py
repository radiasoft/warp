"""Diagnostics captured as particles cross grid cells.
"""
__all__ = ['GridCrossingDiags','GridCrossingDiagsOld']
from warp import *
import cPickle

gridcrossingdiags_version = "$Id: gridcrossingdiags.py,v 1.34 2010/04/21 00:15:22 dave Exp $"

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
  - nhist=top.nhist: Specifies how often the data is collected. Note that
                     nhist can be < 1., in which case the data is subdivided
                     in a time step based on the time that the particles
                     crossed the grid cell. The calculation assumes constant
                     velocity over the time step. The internal value of nhist
                     may differ slightly from the input value since it is
                     adjusted so that 1/nhist is an integer.
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
  - laccumulatedata=true: When true, the data is accumulated, otherwise the
                          data will be from only the most recent time step.
                          If false, this assumes that nhist <= 1.
  - starttime,endtime=None: If given, the data will be collected for times
                            only between the two values.
  - lmoving_frame=false: When true, the diagnostic moves with the beam frame.

The following quantities are calculated:
count: count of the number of particles that cross the cell each time
current: current, same as count but scaled by particle charge and 1/top.dt.
vzbar:
xbar, ybar:
xsqbar, ysqbar:
vxbar, vybar, vzbar:
vxsqbar, vysqbar, vzsqbar:
xrms, yrms:
vxrms, vyrms, vzrms:
epsnx, epsny:
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
                 laccumulatedata=true,
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
        self.laccumulatedata = laccumulatedata
        self.starttime = starttime
        self.endtime = endtime
        self.lmoving_frame = lmoving_frame

        self.zoldpid = nextpid()

        # --- Set these to None so that the routines below can check whether
        # --- they have been set up yet.
        self.gcindex = None
        self.gcmoments = None

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
        self._xbar = []
        self._ybar = []
        self._xsqbar = []
        self._ysqbar = []
        self._vxbar = []
        self._vybar = []
        self._vzbar = []
        self._vxsqbar = []
        self._vysqbar = []
        self._vzsqbar = []
        self._xvxbar = []
        self._yvybar = []
        self._xrms = []
        self._yrms = []
        self._vxrms = []
        self._vyrms = []
        self._vzrms = []
        self._epsnx = []
        self._epsny = []
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
        if self.gcmoments is None:
            # --- Only create the instance here when the diagnostic is
            # --- enabled.
            self.gcmoments = GridCrossing_MomentsType()

        # --- Check if one of the built in gcmoments is available for use.
        # --- If not, complain.
        if top.getpyobject('gcmoments1') is None:
            self.gcindex = 1
            top.gcmoments1 = self.gcmoments
        elif top.getpyobject('gcmoments2') is None:
            self.gcindex = 2
            top.gcmoments2 = self.gcmoments
        else:
            raise RunTimeError('Only two grid crossing moments can be enabled at a time')

        # --- This flag must be set for the diagnostics to be done.
        top.lgcmoments = true

        installbeforestep(self.initializegrid)
        installafterstep(self.getdiagnostics)

    def disable(self):
        # --- If not already enabled, then do nothing.
        if self.gcmoments is None: return
        if self.gcindex is None: return

        # --- Free up the gcmoments that was being used.
        if self.gcindex == 1:
            del top.gcmoments1
        elif self.gcindex == 2:
            del top.gcmoments2

        self.gcindex = None

        # --- If there are no other gcmoments active, then turn the flag off.
        if (top.getpyobject('gcmoments1') is None and
            top.getpyobject('gcmoments2') is None):
            top.lgcmoments = false
        
        uninstallbeforestep(self.initializegrid)
        uninstallafterstep(self.getdiagnostics)

    def initializegrid(self):
        # --- Initialize grid parameters if needed. This is done here
        # --- in case the grid had not been setup yet when init was called.
        # --- This is done beforestep so that the GridCrossing_Moments
        # --- group can be setup properly before gridcrossingmoments is
        # --- called.
        if self.zmmin is None: self.zmmin = w3d.zmmin
        if self.zmmax is None: self.zmmax = w3d.zmmax
        if self.dz is None: self.dz = w3d.dz/self.nzscale
        if self.nz is None:
            self.nz = nint((self.zmmax - self.zmmin)/self.dz)
            assert self.nz > 0,"GridCrossingDiags: zmmax - zmmin must be > dz"
            self.dz = (self.zmmax - self.zmmin)/self.nz
        assert abs((self.zmmax-self.zmmin)-self.nz*self.dz) < 1.e-5*self.dz,\
            "zmmin, zmmax, dz, and nz are not consistent with each other"

        self.zmesh = self.zmmin + arange(0,self.nz+1,dtype='l')*self.dz

        if self.nhist is not None and self.nhist < 1:
            self.nt = nint(1./self.nhist)
            # --- Make sure the nt and nhist are consistent.
            self.nhist = 1./self.nt
        else:
            self.nt = 1

        # --- Setup the GridCrossing_Moments group
        self.gcmoments.zmmingc = self.zmmin
        self.gcmoments.zmmaxgc = self.zmmax
        self.gcmoments.ntgc = self.nt
        self.gcmoments.nzgc = self.nz
        self.gcmoments.nszgc = top.ns
        self.gcmoments.dzgc = self.dz
        self.gcmoments.lmoving_framegc = self.lmoving_frame
        self.gcmoments.gchange()

        if self.ldoscintillator:
            self.scintzmin = nint((self.scintzmin-self.zmmin)/self.dz)*self.dz + self.zmmin
            self.scintzmax = nint((self.scintzmax-self.zmmin)/self.dz)*self.dz + self.zmmin
            self.scintnz = nint((self.scintzmax-self.scintzmin)/self.dz)
            if self.scintxmin is None: self.scintxmin = -self.scintxmax
            if self.scintymin is None: self.scintymin = -self.scintymax
            self.scintdx = (self.scintxmax - self.scintxmin)/self.scintnx
            self.scintdy = (self.scintymax - self.scintymin)/self.scintny

    def gettimezbeam(self,zbeam):
        # --- The time extends from just after the last time up to and
        # --- including the current time.
        tt = top.time - top.dt + arange(1,self.nt+1)*top.dt/self.nt
        zb = ones(self.nt)*zbeam
        return tt,zb

    def appendnextarrays(self,zbeam):
        tt,zb = self.gettimezbeam(zbeam)
        self._time.append(tt)
        self._zbeam.append(zb)
        self._count.append(None)
        self._xbar.append(None)
        self._ybar.append(None)
        self._xsqbar.append(None)
        self._ysqbar.append(None)
        self._rprms.append(None)
        self._vxbar.append(None)
        self._vybar.append(None)
        self._vzbar.append(None)
        self._vxsqbar.append(None)
        self._vysqbar.append(None)
        self._vzsqbar.append(None)
        self._xvxbar.append(None)
        self._yvybar.append(None)
        self._current.append(None)
        self._xrms.append(None)
        self._yrms.append(None)
        self._vxrms.append(None)
        self._vyrms.append(None)
        self._vzrms.append(None)
        self._epsnx.append(None)
        self._epsny.append(None)
        self._rrms.append(None)

        if self.ldoradialdiag:
            nr = self.nr
            nz = self.nz
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
        if getcurrpkg() != 'w3d': return

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
        nt = self.nt

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

        # --- Create temporary work space
        if self.ldoradialdiag:
            rprofilecount = zeros((1+nr,1+nz),'d')
            #self.rprofilemesh = iota(0,nr)*dr
        if self.ldoscintillator:
            scintillatorcount = zeros((1+scintnx,1+scintny,1+scintnz))

        if self.nhist is None: nhist = top.nhist
        else:                  nhist = self.nhist

        # --- The data is gathered from top.it-nhist/2 to top.it+nhist/2-1.
        # --- At the half way point, create space for the next set of data.
        # --- If nhist < 1, then collect data every time step.
        if (len(self._time) == 0 or nhist < 1. or
            (top.it-1)%nhist == int(nhist/2)):

            if self.laccumulatedata and me == 0 and not self.dumptofile:
                self.appendnextarrays(zbeam)
            else:
                # --- On other processors or if the data is being dumped to a
                # --- file or if the data is not being accumulated, just zero
                # --- out the existing arrays.
                # --- There's no reason to keep the history on all processors.
                # --- The arrays are created the first time the diagnostic
                # --- is done.
                if len(self._time) == 0:
                    self.appendnextarrays(zbeam)
                else:
                    tt,zb = self.gettimezbeam(zbeam)
                    self._time[0] = tt
                    self._zbeam[0] = zb
                    if self.ldoradialdiag:
                        self._rprofile[0].fill(0.)
                    if self.ldoscintillator:
                        self._scintillator[0].fill(0.)

        # --- The code below is all local and can be skipped if there are no
        # --- particles locally.
        if top.pgroup.nps[js] > 0:

            if nt == 1 or self.ldoradialdiag or self.ldoscintillator:

                xnew = getx(js=js,gather=0)
                ynew = gety(js=js,gather=0)
                rpnew = getrp(js=js,gather=0)
                znew = getz(js=js,gather=0)
                vznew = getvz(js=js,gather=0)
                zold = getpid(js=js,id=self.zoldpid-1,gather=0)

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

            if self.ldoradialdiag or self.ldoscintillator:
                np = len(izc)
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
            top.pgroup.pid[i1:i2,self.zoldpid-1] = top.pgroup.zp[i1:i2]

        # --- The data is gathered from top.it-nhist/2 to top.it+nhist/2-1.
        # --- At the half way point, finish the calculation by summing over
        # --- processors, dividing by the counts to get the averages, and
        # --- calculating the rms quantities.
        if nhist < 1. or top.it%nhist == int(nhist/2):
            self._count[-1] = self.gcmoments.pnumgc[1:,:,self.js].copy()
            self._xbar[-1] = self.gcmoments.xbargc[1:,:,self.js].copy()
            self._ybar[-1] = self.gcmoments.ybargc[1:,:,self.js].copy()
            self._xsqbar[-1] = self.gcmoments.xsqbargc[1:,:,self.js].copy()
            self._ysqbar[-1] = self.gcmoments.ysqbargc[1:,:,self.js].copy()
            self._rprms[-1] = self.gcmoments.rprmsgc[1:,:,self.js].copy()
            self._vxbar[-1] = self.gcmoments.vxbargc[1:,:,self.js].copy()
            self._vybar[-1] = self.gcmoments.vybargc[1:,:,self.js].copy()
            self._vzbar[-1] = self.gcmoments.vzbargc[1:,:,self.js].copy()
            self._vxsqbar[-1] = self.gcmoments.vxsqbargc[1:,:,self.js].copy()
            self._vysqbar[-1] = self.gcmoments.vysqbargc[1:,:,self.js].copy()
            self._vzsqbar[-1] = self.gcmoments.vzsqbargc[1:,:,self.js].copy()
            self._xvxbar[-1] = self.gcmoments.xvxbargc[1:,:,self.js].copy()
            self._yvybar[-1] = self.gcmoments.yvybargc[1:,:,self.js].copy()

            co = self._count[-1]
            xbar = self._xbar[-1]
            ybar = self._ybar[-1]
            xsqbar = self._xsqbar[-1]
            ysqbar = self._ysqbar[-1]
            rp = self._rprms[-1]
            vxbar = self._vxbar[-1]
            vybar = self._vybar[-1]
            vzbar = self._vzbar[-1]
            vxsqbar = self._vxsqbar[-1]
            vysqbar = self._vysqbar[-1]
            vzsqbar = self._vzsqbar[-1]
            xvxbar = self._xvxbar[-1]
            yvybar = self._yvybar[-1]

            # --- Finish the calculation, gathering data from all processors and
            # --- dividing out the count.
            co[...] = parallelsum(co)
            xbar[...] = parallelsum(xbar)
            ybar[...] = parallelsum(ybar)
            xsqbar[...] = parallelsum(xsqbar)
            ysqbar[...] = parallelsum(ysqbar)
            rp[...] = parallelsum(rp)
            vxbar[...] = parallelsum(vxbar)
            vybar[...] = parallelsum(vybar)
            vzbar[...] = parallelsum(vzbar)
            vxsqbar[...] = parallelsum(vxsqbar)
            vysqbar[...] = parallelsum(vysqbar)
            vzsqbar[...] = parallelsum(vzsqbar)
            xvxbar[...] = parallelsum(xvxbar)
            yvybar[...] = parallelsum(yvybar)

            cotemp = where(co==0.,1.,co)
            xbar[...] = xbar/cotemp
            ybar[...] = ybar/cotemp
            xsqbar[...] = xsqbar/cotemp
            ysqbar[...] = ysqbar/cotemp
            rp[...] = sqrt(rp/cotemp)
            vxbar[...] = vxbar/cotemp
            vybar[...] = vybar/cotemp
            vzbar[...] = vzbar/cotemp
            vxsqbar[...] = vxsqbar/cotemp
            vysqbar[...] = vysqbar/cotemp
            vzsqbar[...] = vzsqbar/cotemp
            xvxbar[...] = xvxbar/cotemp
            yvybar[...] = yvybar/cotemp

            self._xrms[-1] = sqrt(abs(xsqbar - xbar**2))
            self._yrms[-1] = sqrt(abs(ysqbar - ybar**2))
            self._rrms[-1] = sqrt(abs(xsqbar + ysqbar - xbar**2 - ybar**2))
            self._vxrms[-1] = sqrt(abs(vxsqbar - vxbar**2))
            self._vyrms[-1] = sqrt(abs(vysqbar - vybar**2))
            self._vzrms[-1] = sqrt(abs(vzsqbar - vzbar**2))
            self._epsnx[-1] = 4.*sqrt(abs((xsqbar-xbar**2)*(vxsqbar-vxbar**2)
                                          - (xvxbar - xbar*vxbar)**2))
            self._epsny[-1] = 4.*sqrt(abs((ysqbar-ybar**2)*(vysqbar-vybar**2)
                                          - (yvybar - ybar*vybar)**2))

            # --- Scale the current appropriately.
            self._current[-1] = co*(top.pgroup.sq[js]/(top.dt*nhist))

            # --- Now that the data was copied out, zero the top arrays
            # --- so that the data doesn't accumulate.
            self.gcmoments.pnumgc.fill(0.)
            self.gcmoments.xbargc.fill(0.)
            self.gcmoments.ybargc.fill(0.)
            self.gcmoments.xsqbargc.fill(0.)
            self.gcmoments.ysqbargc.fill(0.)
            self.gcmoments.rprmsgc.fill(0.)
            self.gcmoments.vxbargc.fill(0.)
            self.gcmoments.vybargc.fill(0.)
            self.gcmoments.vzbargc.fill(0.)
            self.gcmoments.vxsqbargc.fill(0.)
            self.gcmoments.vysqbargc.fill(0.)
            self.gcmoments.vzsqbargc.fill(0.)
            self.gcmoments.xvxbargc.fill(0.)
            self.gcmoments.yvybargc.fill(0.)


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
        ff.write('time'+suffix,self._time[0])
        ff.write('zbeam'+suffix,self._zbeam[0])
        ff.write('count'+suffix,self._count[0])
        ff.write('current'+suffix,self._current[0])
        ff.write('xbar'+suffix,self._xbar[0])
        ff.write('ybar'+suffix,self._ybar[0])
        ff.write('xsqbar'+suffix,self._xsqbar[0])
        ff.write('ysqbar'+suffix,self._ysqbar[0])
        ff.write('vxbar'+suffix,self._vxbar[0])
        ff.write('vybar'+suffix,self._vybar[0])
        ff.write('vzbar'+suffix,self._vzbar[0])
        ff.write('vxsqbar'+suffix,self._vxsqbar[0])
        ff.write('vysqbar'+suffix,self._vysqbar[0])
        ff.write('vzsqbar'+suffix,self._vzsqbar[0])
        ff.write('xvxbar'+suffix,self._xvxbar[0])
        ff.write('yvybar'+suffix,self._yvybar[0])
        ff.write('xrms'+suffix,self._xrms[0])
        ff.write('yrms'+suffix,self._yrms[0])
        ff.write('vxrms'+suffix,self._vxrms[0])
        ff.write('vyrms'+suffix,self._vyrms[0])
        ff.write('vzrms'+suffix,self._vzrms[0])
        ff.write('epsnx'+suffix,self._epsnx[0])
        ff.write('epsny'+suffix,self._epsny[0])
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
            cPickle.dump(('nt',self.nt),ff,-1)
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
        cPickle.dump(('time'+suffix,self._time[0]),ff,-1)
        cPickle.dump(('zbeam'+suffix,self._zbeam[0]),ff,-1)
        cPickle.dump(('count'+suffix,self._count[0]),ff,-1)
        cPickle.dump(('current'+suffix,self._current[0]),ff,-1)
        cPickle.dump(('xbar'+suffix,self._xbar[0]),ff,-1)
        cPickle.dump(('ybar'+suffix,self._ybar[0]),ff,-1)
        cPickle.dump(('xsqbar'+suffix,self._xsqbar[0]),ff,-1)
        cPickle.dump(('ysqbar'+suffix,self._ysqbar[0]),ff,-1)
        cPickle.dump(('vxbar'+suffix,self._vxbar[0]),ff,-1)
        cPickle.dump(('vybar'+suffix,self._vybar[0]),ff,-1)
        cPickle.dump(('vzbar'+suffix,self._vzbar[0]),ff,-1)
        cPickle.dump(('vxsqbar'+suffix,self._vxsqbar[0]),ff,-1)
        cPickle.dump(('vysqbar'+suffix,self._vysqbar[0]),ff,-1)
        cPickle.dump(('vzsqbar'+suffix,self._vzsqbar[0]),ff,-1)
        cPickle.dump(('xvxbar'+suffix,self._xvxbar[0]),ff,-1)
        cPickle.dump(('yvybar'+suffix,self._yvybar[0]),ff,-1)
        cPickle.dump(('xrms'+suffix,self._xrms[0]),ff,-1)
        cPickle.dump(('yrms'+suffix,self._yrms[0]),ff,-1)
        cPickle.dump(('vxrms'+suffix,self._vxrms[0]),ff,-1)
        cPickle.dump(('vyrms'+suffix,self._vyrms[0]),ff,-1)
        cPickle.dump(('vzrms'+suffix,self._vzrms[0]),ff,-1)
        cPickle.dump(('epsnx'+suffix,self._epsnx[0]),ff,-1)
        cPickle.dump(('epsny'+suffix,self._epsny[0]),ff,-1)
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
        self._xbar = []
        self._ybar = []
        self._xsqbar = []
        self._ysqbar = []
        self._vxbar = []
        self._vybar = []
        self._vzbar = []
        self._vxsqbar = []
        self._vysqbar = []
        self._vzsqbar = []
        self._xvxbar = []
        self._yvybar = []
        self._xrms = []
        self._yrms = []
        self._vxrms = []
        self._vyrms = []
        self._vzrms = []
        self._epsnx = []
        self._epsny = []
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
                self._xbar.append(ff.read('xbar'+suffix))
                self._ybar.append(ff.read('ybar'+suffix))
                self._xsqbar.append(ff.read('xsqbar'+suffix))
                self._ysqbar.append(ff.read('ysqbar'+suffix))
                self._vxbar.append(ff.read('vxbar'+suffix))
                self._vybar.append(ff.read('vybar'+suffix))
                self._vzbar.append(ff.read('vzbar'+suffix))
                self._vxsqbar.append(ff.read('vxsqbar'+suffix))
                self._vysqbar.append(ff.read('vysqbar'+suffix))
                self._vzsqbar.append(ff.read('vzsqbar'+suffix))
                self._xvxbar.append(ff.read('xvxbar'+suffix))
                self._yvybar.append(ff.read('yvybar'+suffix))
                self._xrms.append(ff.read('xrms'+suffix))
                self._yrms.append(ff.read('yrms'+suffix))
                self._vxrms.append(ff.read('vxrms'+suffix))
                self._vyrms.append(ff.read('vyrms'+suffix))
                self._vzrms.append(ff.read('vzrms'+suffix))
                self._epsnx.append(ff.read('epsnx'+suffix))
                self._epsny.append(ff.read('epsny'+suffix))
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

        # --- Read all of the data in. Only keep the data if the time is
        # --- between start and endtime.
        keepdata = 0
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
                    if self.nhist < 1.:
                        # --- Keep the data if any portion of it is with in
                        # --- the stand and end time.
                        t1 = data[1][0]
                        t2 = data[1][-1]
                    else:
                        t1 = t2 = data[1]
                    keepdata = (starttime <= t2 and t1 <= endtime)
                if not readscintillator and data[0][:12] == 'scintillator':
                    data = (data[0],tell)
                if keepdata:
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
        self._xbar = []
        self._ybar = []
        self._xsqbar = []
        self._ysqbar = []
        self._vxbar = []
        self._vybar = []
        self._vzbar = []
        self._vxsqbar = []
        self._vysqbar = []
        self._vzsqbar = []
        self._xvxbar = []
        self._yvybar = []
        self._xrms = []
        self._yrms = []
        self._vxrms = []
        self._vyrms = []
        self._vzrms = []
        self._epsnx = []
        self._epsny = []
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
                self._xbar.append(datadict['xbar'+suffix])
                self._ybar.append(datadict['ybar'+suffix])
                self._xsqbar.append(datadict['xsqbar'+suffix])
                self._ysqbar.append(datadict['ysqbar'+suffix])
                self._vxbar.append(datadict['vxbar'+suffix])
                self._vybar.append(datadict['vybar'+suffix])
                self._vzbar.append(datadict['vzbar'+suffix])
                self._vxsqbar.append(datadict['vxsqbar'+suffix])
                self._vysqbar.append(datadict['vysqbar'+suffix])
                self._vzsqbar.append(datadict['vzsqbar'+suffix])
                self._xvxbar.append(datadict['xvxbar'+suffix])
                self._yvybar.append(datadict['yvybar'+suffix])
                self._xrms.append(datadict['xrms'+suffix])
                self._yrms.append(datadict['yrms'+suffix])
                self._vxrms.append(datadict['vxrms'+suffix])
                self._vyrms.append(datadict['vyrms'+suffix])
                self._vzrms.append(datadict['vzrms'+suffix])
                self._epsnx.append(datadict['epsnx'+suffix])
                self._epsny.append(datadict['epsny'+suffix])
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
        self.arraytime = self.time
        self.arraycurrent = self.current
        self.arrayradius = self.rrms
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
    def pp2dvxbar(self,**kw):
        self._pp2d(self.vxbar,**kw)
    def pp2dvybar(self,**kw):
        self._pp2d(self.vybar,**kw)
    def pp2dvzbar(self,**kw):
        self._pp2d(self.vzbar,**kw)
    def pp2dvxsqbar(self,**kw):
        self._pp2d(self.vxsqbar,**kw)
    def pp2dvysqbar(self,**kw):
        self._pp2d(self.vysqbar,**kw)
    def pp2dvzsqbar(self,**kw):
        self._pp2d(self.vzsqbar,**kw)
    def pp2dxvxbar(self,**kw):
        self._pp2d(self.xvxbar,**kw)
    def pp2dyvybar(self,**kw):
        self._pp2d(self.yvybar,**kw)
    def pp2dxrms(self,**kw):
        self._pp2d(self.xrms,**kw)
    def pp2dyrms(self,**kw):
        self._pp2d(self.yrms,**kw)
    def pp2dvxrms(self,**kw):
        self._pp2d(self.vxrms,**kw)
    def pp2dvyrms(self,**kw):
        self._pp2d(self.vyrms,**kw)
    def pp2dvzrms(self,**kw):
        self._pp2d(self.vzrms,**kw)
    def pp2depsnx(self,**kw):
        self._pp2d(self.epsnx,**kw)
    def pp2depsny(self,**kw):
        self._pp2d(self.epsny,**kw)
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
        
    def hcount(self,z):
        return self._gettimehistory(self.count,z)
    def hcurrent(self,z):
        return self._gettimehistory(self.current,z)
    def hxbar(self,z):
        return self._gettimehistory(self.xbar,z)
    def hybar(self,z):
        return self._gettimehistory(self.ybar,z)
    def hxsqbar(self,z):
        return self._gettimehistory(self.xsqbar,z)
    def hysqbar(self,z):
        return self._gettimehistory(self.ysqbar,z)
    def hvxbar(self,z):
        return self._gettimehistory(self.vxbar,z)
    def hvybar(self,z):
        return self._gettimehistory(self.vybar,z)
    def hvzbar(self,z):
        return self._gettimehistory(self.vzbar,z)
    def hvxsqbar(self,z):
        return self._gettimehistory(self.vxsqbar,z)
    def hvysqbar(self,z):
        return self._gettimehistory(self.vysqbar,z)
    def hvzsqbar(self,z):
        return self._gettimehistory(self.vzsqbar,z)
    def hxvxbar(self,z):
        return self._gettimehistory(self.xvxbar,z)
    def hyvybar(self,z):
        return self._gettimehistory(self.yvybar,z)
    def hxrms(self,z):
        return self._gettimehistory(self.xrms,z)
    def hyrms(self,z):
        return self._gettimehistory(self.yrms,z)
    def hvxrms(self,z):
        return self._gettimehistory(self.vxrms,z)
    def hvyrms(self,z):
        return self._gettimehistory(self.vyrms,z)
    def hvzrms(self,z):
        return self._gettimehistory(self.vzrms,z)
    def hepsnx(self,z):
        return self._gettimehistory(self.epsnx,z)
    def hepsny(self,z):
        return self._gettimehistory(self.epsny,z)
    def hrrms(self,z):
        return self._gettimehistory(self.rrms,z)
    def hrprms(self,z):
        return self._gettimehistory(self.rprms,z)


    # ----------------------------------------------------------------------
    def _timeintegrate(self,data,laverage,weight=None):

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

        if weight is not None:
            # --- Use the given weight, with or without averaging
            count = weight
        elif laverage:
            # --- If no weight is given, average using the particle count
            count = self.count

        for i in range(data.shape[0]):
            if laverage or weight is not None:
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

    def timeintegratedcount(self,laverage=0,weight=None):
        return self._timeintegrate(self.count,laverage,weight)

    def timeintegratedcurrent(self,laverage=0,weight=None):
        return self._timeintegrate(self.current,laverage,weight)

    def timeintegratedvzbar(self,laverage=1,weight=None):
        return self._timeintegrate(self.vzbar,laverage,weight)

    def timeintegratedxbar(self,laverage=1,weight=None):
        return self._timeintegrate(self.xbar,laverage,weight)

    def timeintegratedybar(self,laverage=1,weight=None):
        return self._timeintegrate(self.ybar,laverage,weight)

    def timeintegratedxsqbar(self,laverage=1,weight=None):
        return self._timeintegrate(self.xsqbar,laverage,weight)

    def timeintegratedysqbar(self,laverage=1,weight=None):
        return self._timeintegrate(self.ysqbar,laverage,weight)

    def timeintegratedxrms(self,laverage=1,weight=None):
        data = self.xrms**2
        result,gridmesh = self._timeintegrate(data,laverage,weight)
        result = sqrt(maximum(0.,result))
        return result,gridmesh

    def timeintegratedyrms(self,laverage=1,weight=None):
        data = self.yrms**2
        result,gridmesh = self._timeintegrate(data,laverage,weight)
        result = sqrt(maximum(0.,result))
        return result,gridmesh

    def timeintegratedxprms(self,laverage=1,weight=None):
        data = self.xprms**2
        result,gridmesh = self._timeintegrate(data,laverage,weight)
        result = sqrt(maximum(0.,result))
        return result,gridmesh

    def timeintegratedyprms(self,laverage=1,weight=None):
        data = self.yprms**2
        result,gridmesh = self._timeintegrate(data,laverage,weight)
        result = sqrt(maximum(0.,result))
        return result,gridmesh

    def timeintegratedcorkscrew(self,laverage=1,weight=None):
        xbarint,gridmesh = self.timeintegratedxbar()
        ybarint,gridmesh = self.timeintegratedybar()
        xbarsqint,gridmesh = self._timeintegrate(self.xbar**2,laverage,weight)
        ybarsqint,gridmesh = self._timeintegrate(self.ybar**2,laverage,weight)
        corkscrew = sqrt(maximum(0.,xbarsqint - xbarint**2 + ybarsqint - ybarint**2))
        return corkscrew,gridmesh

    def timeintegratedvzrms(self,laverage=1,weight=None):
        '''If laverage is false, then sum up vzrms directly, otherwise sum
        vzrms**2 and return the sqrt of that'''
        if laverage:
            data = self.vzrms**2
        else:
            data = self.vzrms
        result,gridmesh = self._timeintegrate(data,laverage,weight)
        if laverage:
            result = sqrt(maximum(0.,result))
        return result,gridmesh

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
            if nhist > 1 and top.it%nhist != int(nhist/2):
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
                    if name not in ['rprofile','scinttime','scintillator']:
                        # --- Reshape, putting the time blocks into one
                        # --- dimension.
                        ss = result.shape
                        if len(ss) == 2:
                            result.shape = (ss[0]*ss[1],)
                        else:
                            result.shape = (ss[0]*ss[1],ss[2])
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
    xbar = property(*_setupproperty('xbar'))
    ybar = property(*_setupproperty('ybar'))
    xsqbar = property(*_setupproperty('xsqbar'))
    ysqbar = property(*_setupproperty('ysqbar'))
    vxbar = property(*_setupproperty('vxbar'))
    vybar = property(*_setupproperty('vybar'))
    vzbar = property(*_setupproperty('vzbar'))
    vxsqbar = property(*_setupproperty('vxsqbar'))
    vysqbar = property(*_setupproperty('vysqbar'))
    vzsqbar = property(*_setupproperty('vzsqbar'))
    xvxbar = property(*_setupproperty('xvxbar'))
    yvybar = property(*_setupproperty('yvybar'))
    xrms = property(*_setupproperty('xrms'))
    yrms = property(*_setupproperty('yrms'))
    vxrms = property(*_setupproperty('vxrms'))
    vyrms = property(*_setupproperty('vyrms'))
    vzrms = property(*_setupproperty('vzrms'))
    epsnx = property(*_setupproperty('epsnx'))
    epsny = property(*_setupproperty('epsny'))
    rrms = property(*_setupproperty('rrms'))
    rprms = property(*_setupproperty('rprms'))
    rprofile = property(*_setupproperty('rprofile'))
    scinttime = property(*_setupproperty('scinttime'))
    scintillator = property(*_setupproperty('scintillator'))
    del _setupproperty

class GridCrossingDiagsOld(object):
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
  - nhist=top.nhist: Specifies how often the data is collected. Note that
                     nhist can be < 1., in which case the data is subdivided
                     in a time step based on the time that the particles
                     crossed the grid cell. The calculation assumes constant
                     velocity over the time step. The internal value of nhist
                     may differ slightly from the input value since it is
                     adjusted so that 1/nhist is an integer.
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
        if nhist < 1.:
          # --- Also save the old velocity to give a better extrapolation.
          self.vxoldpid = nextpid()
          self.vyoldpid = nextpid()
          self.vzoldpid = nextpid()


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
            assert self.nz > 0,"GridCrossingDiags: zmmax - zmmin must be > dz"
            self.dz = (self.zmmax - self.zmmin)/self.nz
        assert abs((self.zmmax-self.zmmin)-self.nz*self.dz) < 1.e-5*self.dz,\
            "zmmin, zmmax, dz, and nz are not consistent with each other"

        self.zmesh = self.zmmin + arange(0,self.nz+1,dtype='l')*self.dz

        if self.nhist is not None and self.nhist < 1:
            self.nt = nint(1./self.nhist)
            # --- Make sure the nt and nhist are consistent.
            self.nhist = 1./self.nt
        else:
            self.nt = 1

        if self.ldoscintillator:
            self.scintzmin = nint((self.scintzmin-self.zmmin)/self.dz)*self.dz + self.zmmin
            self.scintzmax = nint((self.scintzmax-self.zmmin)/self.dz)*self.dz + self.zmmin
            self.scintnz = nint((self.scintzmax-self.scintzmin)/self.dz)
            if self.scintxmin is None: self.scintxmin = -self.scintxmax
            if self.scintymin is None: self.scintymin = -self.scintymax
            self.scintdx = (self.scintxmax - self.scintxmin)/self.scintnx
            self.scintdy = (self.scintymax - self.scintymin)/self.scintny

    def gettimezbeam(self,zbeam):
        if self.nt == 1:
            tt = top.time
            zb = zbeam
        else:
            # --- The time extends from just after the last time up to and
            # --- including the current time.
            tt = top.time - top.dt + arange(1,self.nt+1)*top.dt/self.nt
            zb = ones(self.nt)*zbeam
        return tt,zb

    def appendnextarrays(self,zbeam):
        nz = self.nz
        if self.nt == 1:
            ss = (1+nz,)
        else:
            ss = (self.nt,1+nz)

        tt,zb = self.gettimezbeam(zbeam)
        self._time.append(tt)
        self._zbeam.append(zb)
        self._count.append(zeros(ss,'d'))
        self._current.append(zeros(ss,'d'))
        self._vzbar.append(zeros(ss,'d'))
        self._xbar.append(zeros(ss,'d'))
        self._ybar.append(zeros(ss,'d'))
        self._xsqbar.append(zeros(ss,'d'))
        self._ysqbar.append(zeros(ss,'d'))
        self._xrms.append(None)
        self._yrms.append(None)
        self._rrms.append(None)
        self._rprms.append(zeros(ss,'d'))

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
        nt = self.nt

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
        if self.ldoradialdiag:
            rprofilecount = zeros((1+nr,1+nz),'d')
            #self.rprofilemesh = iota(0,nr)*dr
        if self.ldoscintillator:
            scintillatorcount = zeros((1+scintnx,1+scintny,1+scintnz))

        if self.nhist is None: nhist = top.nhist
        else:                  nhist = self.nhist

        # --- The data is gathered from top.it-nhist/2 to top.it+nhist/2-1.
        # --- At the half way point, create space for the next set of data.
        # --- If nhist < 1, then collect data every time step.
        if (len(self._time) == 0 or nhist < 1. or
            (top.it-1)%nhist == int(nhist/2)):

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
                    tt,zb = self.gettimezbeam(zbeam)
                    self._time[0] = tt
                    self._zbeam[0] = zb
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

            if nt == 1 or self.ldoradialdiag or self.ldoscintillator:

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

            if nt == 1:
                gcount = zeros_like(self._count[-1])

                deposgrid1d(1,np,zc,ww,nz,self._count[-1],gcount,0.,nz)
                deposgrid1d(1,np,zc,ww*vzc,nz,self._vzbar[-1],gcount,0.,nz)
                deposgrid1d(1,np,zc,ww*xc,nz,self._xbar[-1],gcount,0.,nz)
                deposgrid1d(1,np,zc,ww*yc,nz,self._ybar[-1],gcount,0.,nz)
                deposgrid1d(1,np,zc,ww*xc**2,nz,self._xsqbar[-1],gcount,0.,nz)
                deposgrid1d(1,np,zc,ww*yc**2,nz,self._ysqbar[-1],gcount,0.,nz)
                deposgrid1d(1,np,zc,ww*rpc**2,nz,self._rprms[-1],gcount,0.,nz)

            else:

                xnew = getx(js=js,gather=0)
                ynew = gety(js=js,gather=0)
                znew = getz(js=js,gather=0)
                vxnew = getvx(js=js,gather=0)
                vynew = getvy(js=js,gather=0)
                vznew = getvz(js=js,gather=0)
                zold = getpid(js=js,id=zoldpid-1,gather=0)
                vxold = getpid(js=js,id=self.vxoldpid-1,gather=0)
                vyold = getpid(js=js,id=self.vyoldpid-1,gather=0)
                vzold = getpid(js=js,id=self.vzoldpid-1,gather=0)
                np = len(xnew)

                if top.wpid > 0:
                    wwnew = getpid(js=js,id=top.wpid-1,gather=0)
                else:
                    wwnew = ones(np,'d')

                # --- The computation is too complicated for Python.
                # --- Do all the work in fortran.
                # --- This includes cases where particles can cross multiple
                # --- grid cells in a time step.
                gridcrossingmomentsold(np,wwnew,xnew,ynew,znew,vxnew,vynew,vznew,
                                    zold,vxold,vyold,vzold,
                                    top.dt,zmmin+zbeam,dz,
                                    nt,nz,
                                    self._count[-1],
                                    self._vzbar[-1],
                                    self._xbar[-1],
                                    self._ybar[-1],
                                    self._xsqbar[-1],
                                    self._ysqbar[-1],
                                    self._rprms[-1])

            if self.ldoradialdiag or self.ldoscintillator:
                np = len(izc)
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
            if self.nhist < 1.:
                top.pgroup.pid[i1:i2,self.vxoldpid-1] = top.pgroup.uxp[i1:i2]
                top.pgroup.pid[i1:i2,self.vyoldpid-1] = top.pgroup.uyp[i1:i2]
                top.pgroup.pid[i1:i2,self.vzoldpid-1] = top.pgroup.uzp[i1:i2]

        # --- The data is gathered from top.it-nhist/2 to top.it+nhist/2-1.
        # --- At the half way point, finish the calculation by summing over
        # --- processors, dividing by the counts to get the averages, and
        # --- calculating the rms quantities.
        if nhist < 1. or top.it%nhist == int(nhist/2):
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
        ff.write('time'+suffix,self._time[0])
        ff.write('zbeam'+suffix,self._zbeam[0])
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
            cPickle.dump(('nt',self.nt),ff,-1)
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
        cPickle.dump(('time'+suffix,self._time[0]),ff,-1)
        cPickle.dump(('zbeam'+suffix,self._zbeam[0]),ff,-1)
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

        # --- Read all of the data in. Only keep the data if the time is
        # --- between start and endtime.
        keepdata = 0
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
                    if self.nhist < 1.:
                        # --- Keep the data if any portion of it is with in
                        # --- the stand and end time.
                        t1 = data[1][0]
                        t2 = data[1][-1]
                    else:
                        t1 = t2 = data[1]
                    keepdata = (starttime <= t2 and t1 <= endtime)
                if not readscintillator and data[0][:12] == 'scintillator':
                    data = (data[0],tell)
                if keepdata:
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
        self.arraytime = self.time
        self.arraycurrent = self.current
        self.arrayradius = self.rrms
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
    def _timeintegrate(self,data,laverage,weight=None):

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

        if weight is not None:
            # --- Use the given weight, with or without averaging
            count = weight
        elif laverage:
            # --- If no weight is given, average using the particle count
            count = self.count

        for i in range(data.shape[0]):
            if laverage or weight is not None:
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

    def timeintegratedcount(self,laverage=0,weight=None):
        return self._timeintegrate(self.count,laverage,weight)

    def timeintegratedcurrent(self,laverage=0,weight=None):
        return self._timeintegrate(self.current,laverage,weight)

    def timeintegratedvzbar(self,laverage=1,weight=None):
        return self._timeintegrate(self.vzbar,laverage,weight)

    def timeintegratedxbar(self,laverage=1,weight=None):
        return self._timeintegrate(self.xbar,laverage,weight)

    def timeintegratedybar(self,laverage=1,weight=None):
        return self._timeintegrate(self.ybar,laverage,weight)

    def timeintegratedxsqbar(self,laverage=1,weight=None):
        return self._timeintegrate(self.xsqbar,laverage,weight)

    def timeintegratedysqbar(self,laverage=1,weight=None):
        return self._timeintegrate(self.ysqbar,laverage,weight)

    def timeintegratedxrms(self,laverage=1,weight=None):
        data = self.xrms**2
        result,gridmesh = self._timeintegrate(data,laverage,weight)
        result = sqrt(maximum(0.,result))
        return result,gridmesh

    def timeintegratedyrms(self,laverage=1,weight=None):
        data = self.yrms**2
        result,gridmesh = self._timeintegrate(data,laverage,weight)
        result = sqrt(maximum(0.,result))
        return result,gridmesh

    def timeintegratedxprms(self,laverage=1,weight=None):
        data = self.xprms**2
        result,gridmesh = self._timeintegrate(data,laverage,weight)
        result = sqrt(maximum(0.,result))
        return result,gridmesh

    def timeintegratedyprms(self,laverage=1,weight=None):
        data = self.yprms**2
        result,gridmesh = self._timeintegrate(data,laverage,weight)
        result = sqrt(maximum(0.,result))
        return result,gridmesh

    def timeintegratedcorkscrew(self,laverage=1,weight=None):
        xbarint,gridmesh = self.timeintegratedxbar()
        ybarint,gridmesh = self.timeintegratedybar()
        xbarsqint,gridmesh = self._timeintegrate(self.xbar**2,laverage,weight)
        ybarsqint,gridmesh = self._timeintegrate(self.ybar**2,laverage,weight)
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
            if nhist > 1 and top.it%nhist != int(nhist/2):
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
                    if (self.nt > 1 and
                        name not in ['rprofile','scinttime','scintillator']):
                        # --- Reshape, putting the time blocks into one
                        # --- dimension.
                        ss = result.shape
                        if len(ss) == 2:
                            result.shape = (ss[0]*ss[1],)
                        else:
                            result.shape = (ss[0]*ss[1],ss[2])
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

