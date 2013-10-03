"""General injection algorithms
"""
__all__ = ['InjectionGaussLaw',
           'particleinjection_version','particleinjection_doc']
from warp import *
import generateconductors
import copy

particleinjection_version = "$Id: particleinjection.py,v 1.17 2011/11/11 21:52:26 rcohen Exp $"
def particleinjection_doc():
    import particleinjection
    print particleinjection.__doc__

class InjectionGaussLaw(object):
    """Sets up injection using Gauss's law to determine the amount of charge to
inject.
  Qnew = eps0 sum(Enorm) - Qold
where Enorm is the normal E field on the surface of the dual cell, that
extends from i-1/2 to i+1/2.

 - conductors: a conductor or list of conductors which act as particle scrapers
               Note that each conductor MUST have a unique id.
 - inj_d: the location of the virtual surface away from the conductor surface. It is in units of
          the geometric mean of the grid cell sizes. It defaults to top.inj_d[0].
 - rnnmax: if set, is an upper bound to number of particles to inject per timestep per cell.
              Note in cylindrical geometry, it is the number injected per timestep per cell divided
              by r (which for a constant flux density is constant, since the cell volume is proportional to r)
 - relax: if set, it is the relaxation paramter for the injected number,
              rnn =  relax*rnn + (1.-relax)*rnn_old
 - includebadparticles=False: Sometimes, particles cannot be interpolated to the emitting surface.
                              Set this to true to include those particles (not recommended).
 - solvers=None: optional list of solvers to use when obtaining the E fields. When not given,
                 the solvers returned by registeredsolvers will be used. If multiple solvers
                 are given, they must all have the same grid.

After an instance is created, additional conductors can be added by calling
the method registerconductors which takes either a conductor or a list of
conductors are an argument.
    """
    def __init__(self,js=None,conductors=None,vthermal=0.,
                 lcorrectede=None,l_inj_addtempz_abs=None,lsmooth121=0,
                 grid=None,inj_d=None,rnnmax=None,relax=None,includebadparticles=False,
                 solvers=None):
        self.vthermal = vthermal
        self.lcorrectede = lcorrectede
        self.l_inj_addtempz_abs = l_inj_addtempz_abs
        self.lsmooth121 = lsmooth121
        self.inj_d = inj_d
        self.relax = relax
        self.includebadparticles = includebadparticles
        self.solvers = solvers

        self.inj_np = 0.   # initial "old" value of number of particles to inject

        #if js is None: js = range(top.ns)
        ## --- Make sure that js is a list
        #try:              js[0]
        #except TypeError: js = [js]
        if js is None: js = 0
        self.js = js

        self.usergrid = (grid is not None)
        # --- Don't create the grid until it is needed.
        self.grid = grid

        # --- register any initial conductors
        self.conductors = []
        if conductors is None:
            # --- Grab a copy of the list of all conductors created so far.
            conductors = copy.copy(generateconductors.listofallconductors)
        self.registerconductors(conductors)

        # --- If the user specified the grid, then add the conductors
        if self.usergrid: self.updateconductors()

        self.rnnmax = rnnmax

        self.injectedparticlesid = nextpid() - 1

        self.enable()

    def enable(self):
        installuserinjection(self.doinjection)
        installuserinjection2(self.finishinjection)

    def disable(self):
        if isinstalleduserinjection(self.doinjection):
            uninstalluserinjection(self.doinjection)
        if isinstalleduserinjection2(self.finishinjection):
            uninstalluserinjection2(self.finishinjection)

    def getlcorrectede(self):
        if self.lcorrectede is not None:
            return self.lcorrectede
        else:
            return f3d.lcorrectede

    def getinj_d(self):
        if self.inj_d is not None:
            return self.inj_d
        else:
            return top.inj_d[0]

    def getrelax(self):
        if self.relax is not None:
            return self.relax
        else:
            return top.inj_param

    def getEfields(self,solvers):
        Ex,Ey,Ez = self.getEfieldsfromsolver(solvers[0])
        for solver in solvers[1:]:
            Ex1,Ey1,Ez1 = self.getEfieldsfromsolver(solver)
            Ex += Ex1
            Ey += Ey1
            Ez += Ez1
        return Ex,Ey,Ez

    def getEfieldsfromsolver(self,solver):
        """Get the E fields from the active field solver.
This gets the E fields at the face centers of the dual cells. With the
grid cells at ijk, the fields are obtained at the following locations:
  Ex(i+-1/2,j     ,k     )
  Ey(i     ,j+-1/2,k     )
  Ez(i     ,j     ,k+-1/2)
The grid points extend inclusively from (0,0,0) to (nx,ny,nz).
The sizes of the E arrays will be:
  Ex(nx+2,ny+1,nz+1)
  Ey(nx+1,ny+2,nz+1)
  Ez(nx+1,ny+1,nz+2)
        """
        # --- This routine could do the appropriate sums if mesh refinement
        # --- is being used. It can also call a routine which includes
        # --- the conductor when calculating the field.

        # --- Test to see if solver is EM
        # --- If so, just copy the fields and return
        if isinstance(solver,EM3D):
            # --- We need the fields 1/2 cell displaced from each node on
            # --- either side of node.
            # --- So we need values in ghost cells.
            # --- The EM solver has solver.nxguard guard cells, uppper and lower, in x,
            # ---  solver.nyguard guard cells in y, etc.
            # --- All fields have nj+2*ngj+1 values in direction j
            # --- So for Ex there are values starting 1/2 cell above the first computational guard cell and ending
            # ---  1/2 cell above the upper end of the compuational grid (last guard-cell node).
            # --- This last row of Ex is not used.
            # ---  Hence the index of Ex corresponding to the one just below the first PHYSICAL cell
            # ---  is nxguard-1, and we use Ex values running from index nxguard-1 up to but not including
            # ---  index -nxguard, hence in Python slicing language, slicing  nxguard-1:-nxguard
            # ---  In the z direction we need Ex from index nzguard to but not including -nzguard, i.e.
            # ---  slicing nzguard:-nzguard.   Similar permutations for other field components.
            Exraw = solver.fields.Ex
            Eyraw = solver.fields.Ey
            Ezraw = solver.fields.Ez
            # --- Test to see if fields are staggered; if yes do nothing;
            # --- if no, stagger them, do calculation, and then convert back.
            convertback = 0
            if solver.fields.l_nodecentered:
                solver.node2yee3d()
                convertback = 1
            nxguard = solver.nxguard
            nyguard = solver.nyguard
            nzguard = solver.nzguard
            if shape(Exraw)[1]==1:
                # --- x-z or r-z
                Ex = Exraw[nxguard-1:-nxguard,:,nzguard:-nzguard]
                Ez = Ezraw[nxguard:-nxguard,:,nzguard-1:-nzguard]
                Ey = zeros((shape(Ez)[0],2,shape(Ex)[2]),'d')
                # --- formerly assumed node centered so did averaging, e.g.
                #          Ex = .5*(Exraw[2:-3,:,3:-3]+Exraw[3:-2,:,3:-3])
            else:
                # --- 3D
                Ex = Exraw[nxguard-1:-nxguard,nyguard:-nyguard,nzguard:-nzguard]
                Ey = Eyraw[nxguard:-nxguard,nyguard-1:-nyguard,nzguard:-nzguard]
                Ez = Ezraw[nxguard:-nxguard,nyguard:-nyguard,nzguard-1:-nzguard]
            if convertback:
                solver.yee2node3d()
            return Ex,Ey,Ez

        # --- Electrostatic field solvers
        if solver is w3d: phip = solver.phip
        else:             phip = solver.potentialp

        if self.getlcorrectede():

            # --- The shape includes a guard cell in the axis parallel
            # --- to the E field. The calculation of s assumes that there
            # --- is one guard cell on each boundary.
            ng = array([solver.nxguardphi,solver.nyguardphi,solver.nzguardphi])
            s = array(phip.shape) - 2*ng
            Ex = zeros((s[0]+1,s[1],s[2]),'d')
            Ey = zeros((s[0],s[1]+1,s[2]),'d')
            Ez = zeros((s[0],s[1],s[2]+1),'d')

            if solver is w3d:
                conductorobject = f3d.conductors
                solvertop = top
                solverf3d = f3d
            else:
                conductorobject = solver.getconductorobject('p')
                solvertop = solver
                solverf3d = solver

            setupconductorfielddata(solver.nx,solver.ny,solver.nz,
                                    solver.nxp,solver.nyp,solver.nzp,
                                    solver.dx,solver.dy,solver.dz,
                                    conductorobject,solvertop.ppdecomp)
            sete3dongridwithconductor(conductorobject,phip,
                                      solver.dx,solver.dy,solver.dz,
                                      solver.nxp,solver.nyp,solver.nzp,
                                      Ex,Ey,Ez,
                                      solver.nxguardphi,solver.nyguardphi,solver.nzguardphi,
                                      solverf3d.bounds)

        else:

            # --- Calculate E's directly from grid.
            dx = solver.dx
            dy = solver.dy
            dz = solver.dz

            xslice = slice(1,-1)
            yslice = slice(1,-1)
            zslice = slice(1,-1)

            # --- When using a 2D geometry, use all of the y dimension,
            # --- which should be of length 1.
            if solver.solvergeom in [w3d.RZgeom,w3d.XZgeom]:
                yslice = Ellipsis

            Ex = ((phip[:-1,yslice,zslice] - phip[1:,yslice,zslice])/dx)
            if solver.solvergeom in [w3d.RZgeom,w3d.XZgeom]:
                Ey = zeros((phip.shape[0]-2,2,phip.shape[2]-2),'d')
            else:
                Ey = ((phip[xslice,:-1,zslice] - phip[xslice,1:,zslice])/dy)
            Ez = ((phip[xslice,yslice,:-1] - phip[xslice,yslice,1:])/dz)

        #self.Ex = Ex
        #self.Ey = Ey
        #self.Ez = Ez
        return Ex,Ey,Ez

    def getintegratedcharge(self,solvers):
        Qold = self.getintegratedchargefromsolver(solvers[0])
        for solver in solvers[1:]:
            Qold1 = self.getintegratedchargefromsolver(solver)
            Qold += Qold1
        return Qold

    def getintegratedchargefromsolver(self,solver):
        """Get the charge, integrated over the dual cell. This is simply
the charge density at the grid point at the center of the dual cell times
the area of the dual cell.
        """
        if isinstance(solver,EM3D):
            if not solver.l_deposit_rho:
                raise Exception("The l_deposit_rho flag must be true for the EM solver")
            nxguard = solver.nxguard
            nyguard = solver.nyguard
            nzguard = solver.nzguard
            rhop = solver.fields.Rho[nxguard:-nxguard,nyguard:-nyguard,nzguard:-nzguard]
        else:
            if solver is w3d: rhop = solver.rhop
            else:             rhop = solver.sourcep
        dx = solver.dx
        dy = solver.dy
        dz = solver.dz
        Irho = rhop*dx*dy*dz

        return Irho

    def getsolvers(self):
        if self.solvers is None:
            solvers = getregisteredsolvers()
        else:
            solvers = self.solvers

        if len(solvers) == 0:
            solvers = [w3d]
        return solvers

    def doinjection(self):
        # --- This is true when the egun model is being used
        if top.inject == 100: return

        self.updategrid()
        solvers = self.getsolvers()
        solver = solvers[0]
        self.l_2d = (solver.solvergeom in [w3d.XYgeom,w3d.RZgeom])
        self.lcylindrical = (solver.solvergeom==w3d.RZgeom)

        dx = solver.dx
        dy = solver.dy
        dz = solver.dz
        dxyz = sqrt(dx**2 + dy**2 + dz**2)

        # --- Get the E fields on the face centers of the dual cell
        Ex,Ey,Ez = self.getEfields(solvers)

        # --- Get the charge integrated over the dual cell
        Qold = self.getintegratedcharge(solvers)

        # --- Do the integrals of E normal over the sides of the dual cell
        Enorm  = Ex[1:,:,:]*dy*dz
        Enorm -= Ex[:-1,:,:]*dy*dz
        if not self.l_2d:
            Enorm += Ey[:,1:,:]*dx*dz
            Enorm -= Ey[:,:-1,:]*dx*dz
        Enorm += Ez[:,:,1:]*dx*dy
        Enorm -= Ez[:,:,:-1]*dx*dy
        Enorm *= eps0

        Qnew = Enorm - Qold

        # --- Only inject particle for cells in or near conductors.
        Qnew = where(self.grid.isinside == 0.,0.,Qnew)
        if self.lsmooth121: self.smooth121(Qnew)
        #if self.lsmooth121: Qnew = where(self.isdeepinside == 1.,0.,Qnew)

        # --- Calculate the number of new particles to add at each grid cell.
        rnn = Qnew/(top.pgroup.sq[self.js]*top.pgroup.sw[self.js])

        # --- Make sure it is positive or zero
        rnn = maximum(rnn,0.)

        # --- If the user has specified a non-zero upper bound to the number of particles, impose it here
        if self.rnnmax:
            rnn = minimum(rnn,self.rnnmax)

        # --- Scale appropriately for cylindrical coordinates
        # --- This accounts the difference in area of a grid cell in
        # --- Cartesian (dx*dy) and cylindrical (2 pi*r*dx).
        if self.lcylindrical:
            if solver.xmmin == 0:
                rnn[0,...] *= 0.25*pi*solver.dx/solver.dy
            else:
                rnn[0,...] *= 2.0*pi*solver.xmmin/solver.dy
            rnn[1:,...] *= 2.0*pi*solver.xmesh[1:,newaxis,newaxis]/solver.dy

        # --- Add a random number to the number of particles injected
        # --- so that the average number of particles injected is
        # --- correct.  For example, if rnn < 1., without the
        # --- addition of the random number, no particles would ever
        # --- be injected.  With the random number, particles will be
        # --- injected but the average number will be less than 1.
        #rnn += where(rnn > 0.,random.random(rnn.shape),0.)
        rnn += random.random(rnn.shape)

        relax = self.getrelax()
        if relax != 1.:
            # self.inj_np holds the previous timestep's rnn data
            rnn = relax*rnn + (1. - relax)*self.inj_np

        # --- Save the number for diagnostics
        self.inj_np = rnn.copy()

        # --- Now create the particles. This is easiest to do in fortran.
        # --- For each dual cell, the particles injected in that cell are
        # --- evenly distributed throughout the cell.
        # --- This also gets the E field at the cell center for each
        # --- of the particles.
        nn = sum(aint(rnn))
        if nn == 0: return
        xx,yy,zz = zeros((3,nn),'d')
        ex,ey,ez = zeros((3,nn),'d')
        pp = zeros(nn,'d')
        good = ones(nn,bool)
        nxp = rnn.shape[0] - 1
        nyp = rnn.shape[1] - 1
        nzp = rnn.shape[2] - 1
        createparticlesincells(nxp,nyp,nzp,rnn,Ex,Ey,Ez,self.grid.isinside,
                               self.lcylindrical,
                               dx,dy,dz,nn,xx,yy,zz,ex,ey,ez,pp)
        xx += solver.xmminp
        if not self.l_2d:
            yy += solver.ymminp
        zz += solver.zmminp

        # --- Give particles a thermal velocity.
        # --- This now ignores the fact the roughly half the particles will be
        # --- headed back into the conductor.
        if self.vthermal > 0.:
            vx = random.normal(0.,self.vthermal,nn)
            vy = random.normal(0.,self.vthermal,nn)
            vz = random.normal(0.,self.vthermal,nn)
        else:
            vx = zeros(nn)
            vy = zeros(nn)
            vz = zeros(nn)

        # --- The location of the virtual surface, for each particle
        xv,yv,zv = zeros((3,nn),'d')

        # --- The E field at the cell centers is used to provide a direction
        # --- for the projection onto the surfaces of the conductors.

        # --- Loop over the conductors, handling all of the particles for each
        # --- conductor at once.
        for c in self.conductors:

            # --- Get the particles near the conductor c
            ii = compress(pp == c.condid,arange(nn))
            if len(ii) == 0: continue
            xc = xx[ii]
            yc = yy[ii]
            zc = zz[ii]
            exc = ex[ii]
            eyc = ey[ii]
            ezc = ez[ii]

            # --- Get a velocity from the E fields. Based on the way intercept
            # --- works, the velocity needs to be pointing away from the
            # --- surface. So converting from E fields to velocity depends on
            # --- the charge of the injected particles and whether the starting
            # --- positions are inside or outside. Note that the magnitude
            # --- of the velocity doesn't matter since it would only set the
            # --- time scale of the interception and that is ignored.
            if top.pgroup.sq[self.js] > 0.:
                vxc,vyc,vzc = exc,eyc,ezc
            else:
                vxc,vyc,vzc = -exc,-eyc,-ezc

            isinside = c.isinside(xc,yc,zc).isinside
            vxc = where(isinside,-vxc,vxc)
            vyc = where(isinside,-vyc,vyc)
            vzc = where(isinside,-vzc,vzc)

            # --- Now the intercept can be calculated.
            intercept = c.intercept(xc,yc,zc,vxc,vyc,vzc)
            xi = intercept.xi
            yi = intercept.yi
            zi = intercept.zi
            itheta = intercept.itheta
            iphi = intercept.iphi

            # --- There are particles that could not be projected to the surface.
            # --- Reject them or replace their position with the original.
            lbadparticles = ((xi-xc)**2+(yi-yc)**2+(zi-zc)**2 > dx**2+dy**2+dz**2)
            itheta = where(lbadparticles,0.,itheta)
            iphi = where(lbadparticles,0.,iphi)
            xi = where(lbadparticles,xc,xi)
            yi = where(lbadparticles,yc,yi)
            zi = where(lbadparticles,zc,zi)

            if self.includebadparticles:
                # --- Reject only particles inside of the conductor
                iibad = ii[lbadparticles]
                ddbad = c.distance(xi[lbadparticles],yi[lbadparticles],zi[lbadparticles]).distance
                good[iibad] = (ddbad >= 0.)

            else:
                # --- Reject all bad particles
                # --- They are only removed from the list after the loop over conductors.
                good[lbadparticles] = False

            # --- Get the virtual locations, one grid cell away from the surface.
            # --- The direction of the surface is along the E field lines.
            qsign = sign(top.pgroup.sq[self.js])
            emag = sqrt(exc**2 + eyc**2 + ezc**2)
            de = qsign*self.getinj_d()*dxyz/dvnz(emag)
            xvc = xi + de*exc
            yvc = yi + de*eyc
            zvc = zi + de*ezc

            # --- This alternative uses the intercept surface normal angle. This doesn't give
            # --- nice results near corners.
            #dv = self.getinj_d()*dxyz
            #ct,st = cos(itheta),sin(itheta)
            #cp,sp = cos(iphi),sin(iphi)
            #xvc = xi - cp*st*dv
            #yvc = yi + sp*st*dv
            #zvc = zi + ct*dv

            # --- Now replace the positions with the projected positions
            xx[ii] = xi
            yy[ii] = yi
            zz[ii] = zi
            xv[ii] = xvc
            yv[ii] = yvc
            zv[ii] = zvc

            # --- Set the velocity so that it is only moving away from the
            # --- surface. Use w3d.l_inj_addtempz_abs by default if the
            # --- flag was not set.
            if self.l_inj_addtempz_abs is None:
                addtempz_abs = w3d.l_inj_addtempz_abs
            else:
                addtempz_abs = self.l_inj_addtempz_abs
            if addtempz_abs:
                # --- The velocity is treated as if it is in the frame with z parallel
                # --- to the surface normal. First, set vz to be positive, and
                # -- then transform the velocity to the lab frame.
                vxc = vx[ii]
                vyc = vy[ii]
                vzc = vz[ii]
                vzc = abs(vzc)
                ct,st = cos(itheta),sin(itheta)
                cp,sp = cos(iphi),sin(iphi)
                vxclab = cp*ct*vxc - sp*vyc + cp*st*vzc
                vyclab = sp*ct*vxc + cp*vyc + sp*st*vzc
                vzclab =   -st*vxc +             ct*vzc
                vx[ii] = vxclab
                vy[ii] = vyclab
                vz[ii] = vzclab

        if top.lrelativ:
            gaminv = sqrt(1. - (vx**2 + vy**2 + vz**2)/clight**2)
            gamma = 1./gaminv
            ux = vx*gamma
            uy = vy*gamma
            uz = vz*gamma
        else:
            gaminv = ones(nn,'d')
            ux = vx
            uy = vy
            uz = vz

        # --- Remove the bad particles
        xx = xx[good]
        yy = yy[good]
        zz = zz[good]
        ux = ux[good]
        uy = uy[good]
        uz = uz[good]
        gaminv = gaminv[good]
        xv = xv[good]
        yv = yv[good]
        zv = zv[good]
        nn = len(xx)

        # --- Get the self and applied fields.
        # --- The self fields are gathered at the virtual surface.
        ex,ey,ez = zeros((3,nn),'d')
        bx,by,bz = zeros((3,nn),'d')
        fetche3dfrompositions(top.pgroup.sid[self.js],top.pgroup.ndts[self.js],nn,xv,yv,zv,ex,ey,ez,bx,by,bz)
        fetchb3dfrompositions(top.pgroup.sid[self.js],top.pgroup.ndts[self.js],nn,xv,yv,zv,bx,by,bz)
        exap,eyap,ezap,bxap,byap,bzap = getappliedfields(xx,yy,zz,time=top.time,js=self.js)
        ex += exap
        ey += eyap
        ez += ezap
        bx += bxap
        by += byap
        bz += bzap

        # --- Give the particles a time step size uniformly distributed
        # --- between 0 and top.dt.
        dt = top.dt*random.random(nn)

        # --- Do a half split leap-frog advance.
        # --- Note that this does the advance in place, directly changing the
        # --- input arrays.
        q = top.pgroup.sq[self.js]
        m = top.pgroup.sm[self.js]
        epusht3d(nn,ux,uy,uz,ex,ey,ez,q,m,dt,0.5)
        gammaadv(nn,gaminv,ux,uy,uz,top.gamadv,top.lrelativ)
        bpusht3d(nn,ux,uy,uz,gaminv,bx,by,bz,q,m,dt,0.5,top.ibpush)
        xpusht3d(nn,xx,yy,zz,ux,uy,uz,gaminv,dt)

        # --- Now add the new particles to the simulation.
        nbefore = top.pgroup.nps[self.js]
        addparticles(xx,yy,zz,ux,uy,uz,gi=gaminv,
                     ex=ex,ey=ey,ez=ez,bx=bx,by=by,bz=bz,
                     js=self.js,lmomentum=true,
                     pidpairs=[[self.injectedparticlesid+1,dt]])
        nafter = top.pgroup.nps[self.js]
        top.npinje_s[self.js] = nafter - nbefore

    def finishinjection(self):
        """Complete the advance of the velocity, so that the time is at the same level as existing particles."""
        # --- This is true when the egun model is being used
        if top.inject == 100: return

        if top.pgroup.npmax == 0: return
        q = top.pgroup.sq[self.js]
        m = top.pgroup.sm[self.js]

        i1 = top.pgroup.ins[self.js] - 1
        i2 = i1 + top.pgroup.nps[self.js]
        ii = (top.pgroup.pid[i1:i2,self.injectedparticlesid] > 0.)
        dt = top.pgroup.pid[i1:i2,self.injectedparticlesid][ii]
        top.pgroup.pid[i1:i2,self.injectedparticlesid][ii] = 0.

        xx = top.pgroup.xp[i1:i2][ii]
        yy = top.pgroup.yp[i1:i2][ii]
        zz = top.pgroup.zp[i1:i2][ii]
        ux = top.pgroup.uxp[i1:i2][ii]
        uy = top.pgroup.uyp[i1:i2][ii]
        uz = top.pgroup.uzp[i1:i2][ii]
        gaminv = top.pgroup.gaminv[i1:i2][ii]
        nn = len(xx)
        if nn == 0: return

        ex,ey,ez = zeros((3,nn),'d')
        bx,by,bz = zeros((3,nn),'d')
        fetche3dfrompositions(top.pgroup.sid[self.js],top.pgroup.ndts[self.js],nn,xx,yy,zz,ex,ey,ez,bx,by,bz)
        fetchb3dfrompositions(top.pgroup.sid[self.js],top.pgroup.ndts[self.js],nn,xx,yy,zz,bx,by,bz)
        exap,eyap,ezap,bxap,byap,bzap = getappliedfields(xx,yy,zz,time=top.time,js=self.js)
        ex += exap
        ey += eyap
        ez += ezap
        bx += bxap
        by += byap
        bz += bzap

        # --- Synch the velocities with the positions.
        bpusht3d(nn,ux,uy,uz,gaminv,bx,by,bz,q,m,dt,0.5,top.ibpush)
        epusht3d(nn,ux,uy,uz,ex,ey,ez,q,m,dt,0.5)
        gammaadv(nn,gaminv,ux,uy,uz,top.gamadv,top.lrelativ)

        # --- Go back one half step to match the existing particles.
        fulldt_s = top.dt*top.pgroup.ndts[self.js]
        bpush3d(nn,ux,uy,uz,gaminv,bx,by,bz,q,m,-0.5*fulldt_s,top.ibpush)
        epush3d(nn,ux,uy,uz,ex,ey,ez,q,m,-0.5*fulldt_s)
        gammaadv(nn,gaminv,ux,uy,uz,top.gamadv,top.lrelativ)

        top.pgroup.uxp[i1:i2][ii] = ux
        top.pgroup.uyp[i1:i2][ii] = uy
        top.pgroup.uzp[i1:i2][ii] = uz
        top.pgroup.gaminv[i1:i2][ii] = gaminv

    def registerconductors(self,newconductors):
        if not isinstance(newconductors,list):
            newconductors = [newconductors]
        for c in newconductors:
            assert c.condid != 0,"The conductor id must be nonzero in order for the particle scraping to work."
            self.conductors.append(c)

    def unregisterconductors(self,conductor,nooverlap=0):
        self.conductors.remove(conductor)
        if self.grid is not None:
            if not nooverlap:
                # --- This is horribly inefficient!!!
                self.grid.resetgrid()
                self.updateconductors()
            else:
                self.grid.removeisinside(conductor)

    def updategrid(self,lforce=0):
        """Update the grid to match any changes to the underlying grid,
for example after load balancing.
        """
        if self.grid is None: lforce = 1
        if self.usergrid and not lforce: return
        solvers = self.getsolvers()
        solver = solvers[0]
        if solver is w3d:
            solvertop = top
        else:
            solvertop = solver
        # --- Check if self.grid.decomp is defined. If not, then force
        # --- an update.
        try:
            self.grid.decomp
        except AttributeError:
            lforce = 1
        if not lforce:
            # --- Check if the solver's grid has changed. If not, then
            # --- return since nothing needs to be done.
            gdc = self.grid.decomp
            tdc = solvertop.ppdecomp
            if (self.grid.nxlocal == solver.nxp and
                self.grid.nylocal == solver.nyp and
                self.grid.nzlocal == solver.nzp and
                self.grid.xmmin == solver.xmmin and
                self.grid.xmmax == solver.xmmax and
                self.grid.ymmin == solver.ymmin and
                self.grid.ymmax == solver.ymmax and
                self.grid.zmmin == solver.zmmin and
                self.grid.zmmax == solver.zmmax and
                gdc.ix[gdc.ixproc] == tdc.ix[tdc.ixproc] and
                gdc.iy[gdc.iyproc] == tdc.iy[tdc.iyproc] and
                gdc.iz[gdc.izproc] == tdc.iz[tdc.izproc] and
                gdc.nx[gdc.ixproc] == tdc.nx[tdc.ixproc] and
                gdc.ny[gdc.iyproc] == tdc.ny[tdc.iyproc] and
                gdc.nz[gdc.izproc] == tdc.nz[tdc.izproc]):
                return

        # --- Note that a copy of the decomposition is passed in.
        # --- The decomposition in top may be changed the next time
        # --- loadbalancing is done, but the decomposition in self.grid
        # --- should not be changed. Instead, a whole new grid is created.
        self.grid = Grid(decomp=copy.deepcopy(solvertop.ppdecomp))
        self.updateconductors()

    def updateconductors(self):
        aura = max(self.grid.dx,self.grid.dy,self.grid.dz)
       #if self.lsmooth121:
       #    for c in self.conductors:
       #        self.grid.getisinside(c,aura=-aura)
       #    self.isdeepinside = self.grid.isinside.copy()
        for c in self.conductors:
            self.grid.getisinside(c,aura=aura)

    def smooth121(self,Q):
        Qcopy = Q.copy()
        smooth121nonzero(Qcopy,Q,Q.shape[0]-1,Q.shape[1]-1,Q.shape[2]-1)
        return
        '''
        if Q.shape[1] == 1:
          Q[1:-1,0,1:-1] = (0.0625*(Q[:-2,0,:-2] + Q[2:,0,:-2] +
                                    Q[:-2,0,2:]  + Q[2:,0,2:]) +
                            0.1250*(Q[:-2,0,1:-1] + Q[2:,0,1:-1] +
                                    Q[1:-1,0,:-2] + Q[1:-1,0,2:]) +
                            0.25*Q[1:-1,0,1:-1])
        else:
          Q[1:-1,1:-1,1:-1] = (  0.125*(Q[:-2,:-2,:-2] + Q[2:,:-2,:-2] +
                                        Q[:-2,:-2,2:] + Q[2:,:-2,2:] +
                                        Q[:-2,2:,:-2] + Q[2:,2:,:-2] +
                                        Q[:-2,2:,2:] + Q[2:,2:,2:]) +
                                 0.25*( Q[1:-1,:-2,:-2] + Q[1:-1,2:,:-2] +
                                        Q[1:-1,:-2,2:] + Q[1:-1,2:,2:] +
                                        Q[:-2,1:-1,:-2] + Q[2:,1:-1,:-2] +
                                        Q[:-2,1:-1,2:] + Q[2:,1:-1,2:] +
                                        Q[:-2,:-2,1:-1] + Q[2:,:-2,1:-1] +
                                        Q[:-2,2:,1:-1] + Q[2:,2:,1:-1]) +
                                  0.5*( Q[:-2,1:-1,1:-1] + Q[2:,1:-1,1:-1] +
                                        Q[1:-1,:-2,1:-1] + Q[1:-1,2:,1:-1] +
                                        Q[1:-1,1:-1,:-2] + Q[1:-1,1:-1,2:]) +
                                      ( Q[1:-1,1:-1,1:-1]))*0.125
        '''
