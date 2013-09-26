"""PlaneRestore class that is used to restore particle and field data
at a specified z plane, that was saved by PlaneSave.
The two simulations are linked together.
"""

__all__ = ['PlaneRestore']

from warp import *
import cPickle

class PlaneRestore:
    """
Restores the particle data and phi from a file, incorporating them into the
current simulation. The saved phi is used as a boundary condition. The particles
are injected from the saved plane.
Input:
  - filename=runid.plane: filename where data is stored
  - zplane: location where simulations are linked together. Units of meters
            relative to the lab frame. Defaults to w3d.zmmin.
  - js: species which are saved. Defaults to all species. Can be single
        integer or a list of integers.
  - l_restore_phi=True: flag for restoring phi or not.
  - lrestoreparticles=True: flag for restoring particles
  - starttime=None: If specified, the time at which to start the simulation.
                    This can be used to skip part of the saved data, or to
                    start at an earlier time before saved data is available.
  - verbose=False: When True, prints messages about what is happening
    """

    def __init__(self,filename,zplane=None,js=None,
                 l_restore_phi=True,lrestoreparticles=True,starttime=None,verbose=False):

        # --- Save some input values
        self.filename = filename
        self.zplane = zplane
        self.js = js
        self.l_restore_phi = l_restore_phi
        self.lrestoreparticles = lrestoreparticles
        self.starttime = starttime
        self.verbose = verbose

        self.initted = False

        self.f = open(self.filename,'r')
        self.readinitdata()

        # --- Install the routines that do the work.
        self.lsavesynchronized = self.initdata['lsavesynchronized']

        installbeforefs(self.restoreplane_bfs)
        installafterfs(self.restoreplane_afs)

        if self.lsavesynchronized:
            # --- Turn on "all special" so that the particles are synchronized
            # --- when new synchronized particles are added.
            top.allspecl = true
            installafterstep(self.restoreparticles)
        else:
            installuserinjection(self.restoreparticles)

    def disable(self):
        "Turn off restoration of data"
        uninstallbeforefs(self.restoreplane_bfs)
        uninstallafterfs(self.restoreplane_afs)
        if self.lsavesynchronized:
            uninstallafterstep(self.restoreparticles)
        else:
            uninstalluserinjection(self.restoreparticles)

    def read(self):
        return cPickle.load(self.f)

    def readinitdata(self):
        "Read in the initial data"
        self.initdata = {}
        while True:
            name,val = self.read()
            self.initdata[name] = val
            if name == 'solvergeom': break

    def readnextstep(self):
        self.data = {}
        while True:
            try:
                name,val = self.read()
            except EOFError:
                # --- There is no more data. Turn the restore off.
                self.disable()
                # --- Put a dummy value in None to flag that the restore for
                # --- this step should be skipped.
                self.data['it'] = None
                if self.verbose: print "PlaneRestore: no more data, restore is ending"
                return
            self.data[name] = val
            if name.startswith('phiplane'): break

        # --- name will be 'phiplane%09d'%it
        self.data['it'] = int(name[8:])

        if self.verbose: print "PlaneRestore: read in data from step %d"%self.data['it']

    def initrestoreplane(self):
        self.initted = True

        if self.zplane is None: self.zplane = w3d.zmmin

        self.readnextstep()

        self.zshift = self.zplane - self.initdata['zplane']

        self.lsavephi = self.initdata['lsavephi']
        self.lsaveparticles = self.initdata['lsaveparticles']
        if not self.lsavephi: self.l_restore_phi = False
        if not self.lsaveparticles: self.lrestoreparticles = False

        # --- get time level of first plane and subtract 1
        self.it_restore = -1

        # --- get time step, tmin, tmax
        self.dt = self.initdata['dt']
        self.deltat = self.initdata['deltat']
        self.tmin = self.initdata['tmin']
        top.time = self.tmin
        self.time_restore = self.tmin

        if self.lrestoreparticles:
            # --- initializes list of species
            if self.js is None:
                self.jslist = range(top.ns)
            else:
                try:
                    list(self.js)
                    self.jslist = self.js
                except TypeError:
                    self.jslist= [self.js]

            # --- restore particle charge, mass, weight
            for js in self.jslist:
                top.pgroup.sq[js] = self.initdata['sq_%d'%js]
                top.pgroup.sm[js] = self.initdata['sm_%d'%js]
                top.pgroup.sw[js] = self.initdata['sw_%d'%js]

            # --- make sure that pid will be allocated
            #top.npid = self.initdata['npid']
            #setuppgroup(top.pgroup)

        if self.l_restore_phi:
            # --- restore solver geometry of the saved data
            try:
                self.solvergeom = self.initdata['solvergeom']
            except:
                self.solvergeom = w3d.XYZgeom

            # set up indices which specify transverse extent of saved and restored phi
            # _r for restored phi array, _s for saved phi array
            # '0' is minimum index, 'm' is maximum index

            self.sym_plane = self.initdata['sym_plane']
            ixa_plane = self.initdata['ixa_plane']
            iya_plane = self.initdata['iya_plane']
            self.nx_plane = self.initdata['nx_plane']
            self.ny_plane = self.initdata['ny_plane']
            self.xmmin = self.initdata['xmmin']
            self.xmmax = self.initdata['xmmax']
            self.ymmin = self.initdata['ymmin']
            self.ymmax = self.initdata['ymmax']

            #self.nx0_r = max(0,int(floor((self.xmmin - w3d.xmmin)/w3d.dx)))
            #self.ny0_r = max(0,int(floor((self.ymmin - w3d.ymmin)/w3d.dy)))
            #self.nxm_r = min(w3d.nx,int(floor((self.xmmax - w3d.xmmin)/w3d.dx)))
            #self.nym_r = min(w3d.ny,int(floor((self.ymmax - w3d.ymmin)/w3d.dy)))

            self.nx0_r = max(0, 0 - ixa_plane + w3d.ix_axis)
            self.ny0_r = max(0, 0 - iya_plane + w3d.iy_axis)
            self.nxm_r = min(w3d.nx, self.nx_plane - ixa_plane + w3d.ix_axis)
            self.nym_r = min(w3d.ny, self.ny_plane - iya_plane + w3d.iy_axis)
            self.nx0_s = self.nx0_r - w3d.ix_axis + ixa_plane
            self.ny0_s = self.ny0_r - w3d.iy_axis + iya_plane
            self.nxm_s = self.nxm_r - w3d.ix_axis + ixa_plane
            self.nym_s = self.nym_r - w3d.iy_axis + iya_plane

            # --- deal with symmetries
            # --- if saved is 2 or 4 fold symmetric and restored isn't,
            # --- lower half of restored is filled with inverted saved phi
            if ((self.sym_plane == 2 and (not w3d.l2symtry and not w3d.l4symtry)) or
                (self.sym_plane == 4 and (not w3d.l2symtry and not w3d.l4symtry))):
                self.ny0_r2 = max(0, - self.ny_plane - iya_plane + w3d.iy_axis)
                self.nym_r2 = min(w3d.ny, 0 - iya_plane + w3d.iy_axis)
                self.ny0_s2 = - self.ny0_r + w3d.iy_axis + iya_plane
                self.nym_s2 =   self.nym_r - w3d.iy_axis + iya_plane
            if ((self.sym_plane == 4 and (not w3d.l2symtry and not w3d.l4symtry)) or
                (self.sym_plane == 4 and (    w3d.l2symtry and not w3d.l4symtry))):
                self.nx0_r2 = max(0, - self.nx_plane - ixa_plane + w3d.ix_axis)
                self.nxm_r2 = min(w3d.nx, 0 - ixa_plane + w3d.ix_axis)
                self.nx0_s2 = self.nxm_r - w3d.ix_axis + ixa_plane
                self.nxm_s2 = - self.nx0_r + w3d.ix_axis + ixa_plane

        # --- Reset the time to the start time if specified.
        if self.starttime is not None:
            top.time = self.starttime
            # --- Advance time_restore until it gets to top.time.
            # --- Use self.deltat/2. to avoid round off problems that
            # --- could occur when comparing time_restore+deltat to top.time.
            while self.time_restore+self.deltat/2. < top.time:
                # --- increment the timelevel of the plane
                self.it_restore += 1
                self.time_restore += self.deltat
            # --- Setup phi at the start time
            self.restoreplane_bfs()

        if self.verbose:
            print "PlaneRestore: initial data"
            print "  File",self.filename
            print "  Restoring phi",self.lsavephi
            print "  Restoring particles",self.lsaveparticles
            print "  Start time",self.tmin
            print "  Time step",self.deltat

    ###########################################################################
    def disable_plane_restore(self):
        # for some reason, does not work!
        uninstalluserinjection(self.restoreparticles)
        uninstallbeforefs(self.restoreplane_bfs)
        uninstallafterfs(self.restoreplane_afs)

    def jumptotime(self,time):
        """Jump to the specified time and set the phi boundary condition.
    No particles are loaded."""

        # --- Do the initialization if it hasn't been done yet.
        if not self.initted:
            self.initrestoreplane()

        # --- Set the time to the desired time.
        top.time = time

        # --- Advance time_restore until it gets to top.time.
        # --- Use self.deltat/2. to avoid round off problems that
        # --- could occur when comparing time_restore+deltat to top.time.
        while self.time_restore+self.deltat/2. < top.time:
            # --- increment the timelevel of the plane
            self.it_restore += 1
            self.time_restore += self.deltat

        # --- restore phi only if between grid bounds
        if (self.zplane < w3d.zmmin+top.zbeam or
            self.zplane+top.zbeam >= w3d.zmmax): return

        # --- calculate grid location of new_plane
        iz = nint((self.zplane - top.zbeam - w3d.zmmin)/w3d.dz)

        # --- load saved phi into the phi array
        self.restore_phi(iz,self.it_restore)

    ###########################################################################
    def restoreparticles(self):
        "Restore the particles"
        # --- Do the initialization if it hasn't been done yet.
        if not self.initted:
            self.initrestoreplane()

        # --- restore only if between grid bounds
        if (self.zplane < w3d.zmmin+top.zbeam or
            self.zplane >= w3d.zmmax+top.zbeam): return

        if self.data['it'] is None: return

        # --- Loop over restored data, restoring the data up to and including
        # --- the current time level of the simulation. This allows the stored
        # --- data deltat to be different than the current simulation dt.
        # --- Use time_restore-deltat/2 to avoid round off problems that
        # --- could occur when comparing time_restore to top.time.
        while self.time_restore-self.deltat/2. < top.time:

            # --- increment the timelevel of the plane
            self.it_restore += 1
            self.time_restore += self.deltat

            while self.data['it'] < self.it_restore:
                self.readnextstep()
                if self.data['it'] is None: return

            # --- Apparently, no data was written for this step, so do nothing
            if self.data['it'] > self.it_restore:
                if self.verbose: print "PlaneRestore: no data for step",self.it_restore
                continue

            # --- load particles for each species
            for js in self.jslist:
                self.restoreparticlespecies(js,self.it_restore)

    def restoreparticlespecies(self,js=0,it=0):
        if not self.lrestoreparticles: return

        # --- put restored data into particle arrays, adjusting the z location
        suffix = '%09d_%d'%(it,js)

        # --- Check if data was written for this step.
        if 'xp'+suffix not in self.data:
            if self.verbose: print "PlaneRestore: no particle data for step",it
            return

        xx = self.data['xp'+suffix]
        yy = self.data['yp'+suffix]
        zz = self.data['zp'+suffix] + self.zshift
        ux = self.data['uxp'+suffix]
        uy = self.data['uyp'+suffix]
        uz = self.data['uzp'+suffix]
        gi = self.data['gaminv'+suffix]
        pid = self.data['pid'+suffix]

        # --- Do some fudging to get the shape of pid correct. This is not
        # --- perfect, since it will may munge the data in pid if things
        # --- are arranged differently.
        if pid.shape[1] < top.npid:
            newid = fzeros((pid.shape[0],top.npid),'d')
            newid[:,:pid.shape[1]] = pid
            pid = newid
        elif pid.shape[1] > top.npid:
            pid = pid[:,:top.npid]

        # --- If the current time step is larger than the saved time step, the particles may need
        # --- to be advanced to the current time level. If the particles were not synchronized, then
        # --- the velocity will be half a saved time step size back and will always need advancing.
        # --- If the particles were being saved every time step (i.e. deltat == dt), then particles saved
        # --- during the early part of the current time step will need to be advanced to the current
        # --- time level.
        # --- If the saved particles were synchronized and deltat == top.dt, then no advance is needed.
        # --- Note that in older data files, time was not saved, so this advance will be skipped.
        tname = 'time%09d'%it
        if tname in self.data:
            oldtime = self.data[tname]
            if (nint(top.dt/self.dt) > 1 and (not self.lsavesynchronized or
                                              (nint(self.deltat/self.dt) == 1 and oldtime+self.dt/2. < top.time))):
                # --- Make sure to only advance the local particles.
                xx,yy,zz,ux,uy,uz,gi,pid = self.localizeparticles(xx,yy,zz,ux,uy,uz,gi,pid)
                if len(xx) == 0: return
                self.advanceparticles(oldtime,xx,yy,zz,ux,uy,uz,gi,js)

        # --- Check if particles are being added out of bounds
        if zz.min() < top.zpmin+top.zbeam or zz.max() > top.zpmax+top.zbeam:
            print "PlaneRestore: restored particles are out of bounds."
            print "\nThe extent of the simulation is %f to %f"%(top.zpmin+top.zbeam,top.zpmax+top.zbeam)
            print "The extent of the restored particles is %f to %f\n"%(zz.min(),zz.max())
            raise Exception("PlaneRestore: restored particles are out of bounds.")

        # --- Note that all processors read in the data, but only particles
        # --- within the processors domain are added.
        if self.verbose: print "PlaneRestore: Restoring %d particles on step %d"%(len(xx),it)
        addparticles(xx,yy,zz,ux,uy,uz,gi,pid,
                     js=js,
                     lallindomain=false,
                     lmomentum=true,
                     resetrho=false,
                     lnewparticles=false)

    def localizeparticles(self,xx,yy,zz,ux,uy,uz,gi,pid):
        "Down select the particles, saving only particles within the local domain."
        # --- Note that this is not correct if domain decomposition is done transversely.
        xmmin,xmmax,ymmin,ymmax = getparticleextent()
        zmmin = top.zpminlocal + top.zgrid
        zmmax = top.zpmaxlocal + top.zgrid

        ii = logical_and(logical_and(logical_and(xmmin<=xx,xx<xmmax),
                                     logical_and(ymmin<=yy,yy<ymmax)),
                                     logical_and(zmmin<=zz,zz<zmmax))

        return xx[ii],yy[ii],zz[ii],ux[ii],uy[ii],uz[ii],gi[ii],pid[ii,:]

    def advanceparticles(self,oldtime,xx,yy,zz,ux,uy,uz,gi,js):
        "Advance particles to catch them up to the current time level"
        nn = len(xx)
        if nn == 0: return

        jsid = top.pgroup.sid[js]
        ndts = top.pgroup.ndts[js]
        q = top.pgroup.sq[js]
        m = top.pgroup.sm[js]

        # --- Gather the self and applied fields
        ex,ey,ez,bx,by,bz = zeros((6,nn))
        fetche3dfrompositions(jsid,ndts,nn,xx,yy,zz,ex,ey,ez,bx,by,bz)
        exap,eyap,ezap,bxap,byap,bzap = getappliedfields(xx,yy,zz,time=top.time,js=js)
        ex += exap
        ey += eyap
        ez += ezap
        bx += bxap
        by += byap
        bz += bzap

        if self.lsavesynchronized:
            if self.verbose: print "RestorePlane: synchornizing particles from time",oldtime,"to time",top.time

            # --- Do a split leap frog advance to get the particles to the
            # --- current time level, with positions and velocity synchronized.
            deltime = top.time - oldtime
            bpush3d(nn,ux,uy,uz,gi,bx,by,bz,q,m,deltime/2.,top.ibpush)
            epush3d(nn,ux,uy,uz,ex,ey,ez,q,m,deltime/2.)
            gammaadv(nn,gi,ux,uy,uz,top.gamadv,top.lrelativ)
            xpush3d(nn,xx,yy,zz,ux,uy,uz,gi,deltime)
            epush3d(nn,ux,uy,uz,ex,ey,ez,q,m,deltime/2.)
            gammaadv(nn,gi,ux,uy,uz,top.gamadv,top.lrelativ)
            bpush3d(nn,ux,uy,uz,gi,bx,by,bz,q,m,deltime/2.,top.ibpush)

        else:
            if self.verbose: print "RestorePlane: advancing particles from time",oldtime,"to time",top.time

            # --- This is messy and not really recommended. It would be much better
            # --- to save the data synchronized.

            # --- Do a split leap frog advance to get the particles to the
            # --- current time level, with positions at the current time level
            # --- and the velocities half a current time step size back.

            # --- First, synchronized the position and velocity, advancing velocities
            # --- a half of the saved time step size.
            epush3d(nn,ux,uy,uz,ex,ey,ez,q,m,self.dt/2.)
            gammaadv(nn,gi,ux,uy,uz,top.gamadv,top.lrelativ)
            bpush3d(nn,ux,uy,uz,gi,bx,by,bz,q,m,self.dt/2.,top.ibpush)

            # --- Do split lead frog to get the positions and velocities to the current time
            deltime = top.time - oldtime
            bpush3d(nn,ux,uy,uz,gi,bx,by,bz,q,m,deltime/2.,top.ibpush)
            epush3d(nn,ux,uy,uz,ex,ey,ez,q,m,deltime/2.)
            gammaadv(nn,gi,ux,uy,uz,top.gamadv,top.lrelativ)
            xpush3d(nn,xx,yy,zz,ux,uy,uz,gi,deltime)
            epush3d(nn,ux,uy,uz,ex,ey,ez,q,m,deltime/2.)
            gammaadv(nn,gi,ux,uy,uz,top.gamadv,top.lrelativ)
            bpush3d(nn,ux,uy,uz,gi,bx,by,bz,q,m,deltime/2.,top.ibpush)

            # --- Push velocities backward half a step - this is only approximate
            epush3d(nn,ux,uy,uz,ex,ey,ez,q,m,-top.dt/2.)
            gammaadv(nn,gi,ux,uy,uz,top.gamadv,top.lrelativ)
            bpush3d(nn,ux,uy,uz,gi,bx,by,bz,q,m,-top.dt/2.,top.ibpush)

    ###########################################################################
    # --- restore the next plane of data
    def restoreplane_bfs(self):

        # --- Do the initialization if it hasn't been done yet.
        if not self.initted:
            self.initrestoreplane()

        # --- restore only if between grid bounds
        if (self.zplane < w3d.zmmin+top.zbeam or
            self.zplane+top.zbeam >= w3d.zmmax): return

        if self.data['it'] is None: return

        # --- calculate grid location of new_plane
        iz = nint((self.zplane - top.zbeam - w3d.zmmin)/w3d.dz)

        # --- load saved phi into the phi array
        self.restore_phi(iz,self.it_restore)

        if self.verbose: print "PlaneRestore: Restoring phi on step %d"%(self.it_restore)

    def restoreplane_afs(self):
        # --- this routine resets the potential at the plane iz=-1 after the
        # --- field solve if this is needed

        # --- restore only if between grid bounds
        if (self.zplane < w3d.zmmin+top.zbeam or
            self.zplane+top.zbeam >= w3d.zmmax): return

        if self.data['it'] is None: return

        # --- calculate grid location of new_plane
        iz = nint((self.zplane - top.zbeam - w3d.zmmin)/w3d.dz)

        # --- reset phi at plane iz=-1 if zplane is at iz=0
        if (iz == 0):
            self.restore_phi(iz,self.it_restore)

    #######################################################################
    def restore_phi(self,iz,it):
        # --- return if flag indicates phi not to be restored
        if not self.l_restore_phi: return

        while self.data['it'] < it:
            self.readnextstep()
            if self.data['it'] is None: return

        # --- Apparently, no data was written for this step, so do nothing
        if self.data['it'] > it: return

        # --- Read in the phi data if it is available.
        if 'phiplane%09d'%it not in self.data: return

        savedphi = self.data['phiplane%09d'%it]

        if savedphi is None: return

        solver = getregisteredsolver()
        if solver is None: solver = w3d

        if self.solvergeom == solver.solvergeom and solver.solvergeom == w3d.XYZgeom:
            self.restore_phi_3d_to_3d(iz-top.izfsslave[me],it,savedphi,solver.phi,
                                      solver)
            if top.izpslave[me] != top.izfsslave[me]:
                # --- This is not really correct, since phip will have a different
                # --- shape and phi, so the dimensions should be passed in too.
                self.restore_phi_3d_to_3d(iz-top.izpslave[me],it,savedphi,solver.phip,
                                          solver)
        elif self.solvergeom == solver.solvergeom and solver.solvergeom == w3d.RZgeom:
            self.restore_phi_rz_to_rz(iz-top.izfsslave[me],it,savedphi,solver.phi,
                                      solver)
        elif self.solvergeom == w3d.RZgeom and solver.solvergeom == w3d.XYZgeom:
            self.restore_phi_rz_to_3d(iz,it,savedphi,solver.phi)
            if top.izpslave[me] != top.izfsslave[me]:
                self.restore_phi_rz_to_3d(iz,it,savedphi,solver.phip)

    #######################################################################
    # This routine copies the saved phi plane into the current phi array
    # making use of different numbers of grid cells and differing symmetries.
    # Both saved and restored phi are 3-D.
    def restore_phi_3d_to_3d(self,iz,it,savedphi,phi,solver):
        if iz < 0 or iz > solver.nzlocal: return
        for i in range(2):
            grid2grid(phi[self.nx0_r:self.nxm_r+1,self.ny0_r:self.nym_r+1,iz+i],
                      solver.nx,solver.ny,
                      solver.xmmin,solver.xmmax,solver.ymmin,solver.ymmax,
                      savedphi[...,i],self.nx_plane,self.ny_plane,
                      self.xmmin,self.xmmax,self.ymmin,self.ymmax)

        if ((self.sym_plane == 2 and (not solver.l2symtry and not solver.l4symtry)) or
            (self.sym_plane == 4 and (not solver.l2symtry and not solver.l4symtry))):
        #     phi(self.nx0_r:self.nxm_r,self.ny0_r2:self.nym_r2,iz-1:iz)=
        #      phi_plane(nx0_s:nxm_s,nym_s2:ny0_s2:-1,)
            for i in range(2):
                grid2grid(phi[self.nx0_r:self.nxm_r,self.ny0_r2:self.nym_r2,iz+i],
                          solver.nx,solver.ny,
                          solver.xmmin,solver.xmmax,solver.ymmin,solver.ymmax,
                          savedphi[...,i],self.nx_plane,self.ny_plane,
                          self.xmmin,self.xmmax,self.ymmin,self.ymmax)

        if ((self.sym_plane == 4 and ( solver.l2symtry and not solver.l4symtry))):
        #  phi(nx0_r2:nxm_r2,ny0_r:nym_r,iz-1:iz)=
        #    phi_plane(nx0_s2:nxm_s2:-1,ny0_s:nym_s,)
            for i in range(2):
                grid2grid(phi[self.nx0_r:self.nxm_r,self.ny0_r2:self.nym_r2,iz+i],
                          solver.nx,solver.ny,
                          solver.xmmin,solver.xmmax,solver.ymmin,solver.ymmax,
                          savedphi[...,i],self.nx_plane,self.ny_plane,
                          self.xmmin,self.xmmax,self.ymmin,self.ymmax)

        if ((self.sym_plane == 4 and (not solver.l2symtry and not solver.l4symtry))):
        #  phi(nx0_r2:nxm_r2,ny0_r2:nym_r,iz-1:iz)=
        #    phi_plane(nx0_s2:nxm_s2:-1,ny0_s2:nym_s,)
            for i in range(2):
                grid2grid(phi[self.nx0_r:self.nxm_r,self.ny0_r2:self.nym_r2,iz+i],
                          solver.nx,solver.ny,
                          solver.xmmin,solver.xmmax,solver.ymmin,solver.ymmax,
                          savedphi[...,i],self.nx_plane,self.ny_plane,
                          self.xmmin,self.xmmax,self.ymmin,self.ymmax)

    def restore_phi_rz_to_rz(self,iz,it,savedphi,phi,solver):
        if iz < 0 or iz > solver.nzlocal: return
        # --- For now, this assumes that the arrays are the same shape.
        phi[1:-1,0,iz] = savedphi[:,0]
        phi[1:-1,0,iz+1] = savedphi[:,1]

    #######################################################################
    # This routine copies the saved phi plane into the current phi array
    # making use of different numbers of grid cells and differing symmetries.
    # Saved phi is rz and restored phi is 3-D.
    def restore_phi_rz_to_3d(self,iz,it,savedphi,phi):
        if iz < 0 or iz > w3d.nzlocal: return
        try:
            self._rz_to_3d_inited
        except:
            self._rz_to_3d_inited = 1
            xmmin = w3d.xmmin + w3d.dx*self.nx0_r
            nx = self.nxm_r - self.nx0_r
            ymmin = w3d.ymmin + w3d.dy*self.ny0_r
            ny = self.nym_r - self.ny0_r
            xmesh,ymesh = getmesh2d(xmmin,w3d.dx,nx,ymmin,w3d.dy,ny)
            rmesh = sqrt(xmesh**2 + ymesh**2)
            dr = (self.xmmax - self.xmmin)/self.nx_plane
            self.irmesh = aint(rmesh/dr)
            self.wrmesh =     rmesh/dr  - self.irmesh
            self.wrmesh = where(self.irmesh >= self.nx_plane,1,self.wrmesh)
            self.irmesh = where(self.irmesh >= self.nx_plane,self.nx_plane-1,self.irmesh)

        i1 = self.nx0_r
        i2 = self.nxm_r+1
        j1 = self.ny0_r
        j2 = self.nym_r+1
        for i in range(2):
            phi[i1:i2,j1:j2,iz+i] = (
                take(savedphi[:,0,i],self.irmesh  )*(1.-self.wrmesh) +
                take(savedphi[:,0,i],self.irmesh+1)*self.wrmesh)
