"""General injection algorithms
"""
__all__ = ['InjectionGaussLaw',
           'particleinjection_version','particleinjection_doc']
from warp import *
import generateconductors
import copy

particleinjection_version = "$Id: particleinjection.py,v 1.3 2008/05/12 16:32:01 dave Exp $"
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

After an instance is created, additional conductors can be added by calling
the method registerconductors which takes either a conductor or a list of
conductors are an argument.
    """
    def __init__(self,js=None,conductors=None,vthermal=0.,
                 lcorrectede=None,l_inj_addtempz_abs=None,grid=None):
        self.vthermal = vthermal
        self.lcorrectede = lcorrectede
        self.l_inj_addtempz_abs = l_inj_addtempz_abs

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

        self.enable()

    def enable(self):
        installuserinjection(self.doinjection)

    def disable(self):
        if installeduserinjection(self.doinjection):
            uninstalluserinjection(self.doinjection)

    def getlcorrectede(self):
        if self.lcorrectede is not None:
            return self.lcorrectede
        else:
            return f3d.lcorrectede

    def getEfields(self,solver):
        # --- This routine could do the appropriate sums if mesh refinement
        # --- is being used. It can also call a routine which includes
        # --- the conductor when calculating the field.

        if solver is w3d: phip = solver.phip
        else:             phip = solver.potentialp

        if self.getlcorrectede():

            # --- The shape includes a guard cell in the axis parallel
            # --- to the E field. The calculation of s assumes that there
            # --- is one guard cell on each boundary.
            s = array(phip.shape) - 2
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

            setupconductorfielddata(solver.nxp,solver.nyp,solver.nzp,solver.nz,
                                    solver.dx,solver.dy,solver.dz,
                                    conductorobject,
                                    solvertop.my_index,solvertop.nslaves,
                                    solvertop.izpslave,solvertop.nzpslave)
            sete3dongridwithconductor(conductorobject,phip,
                                      solver.xmmin,solver.ymmin,solver.zmmin,
                                      solver.dx,solver.dy,solver.dz,
                                      solver.nxp,solver.nyp,solver.nzp,
                                      Ex,Ey,Ez,
                                      1,1,1,solverf3d.bounds)

        else:

            # --- Calculate E's directly from grid.
            dx = solver.dx
            dy = solver.dy
            dz = solver.dz

            xslice = slice(1,-1)
            yslice = slice(1,-1)
            zslice = slice(1,-1)

            # --- Make sure that phip is 3D
            if len(phip) == 2:
                tphip = transpose(phip)
                tphip.shape = [tphip.shape[0],0,tphip.shape[1]]
                phip = transpse(tphip)
                yslice = 0

            Ex = ((phip[:-1,yslice,zslice] - phip[1:,yslice,zslice])/dx)
            if yslice != 0:
                Ey = ((phip[xslice,:-1,zslice] - phip[xslice,1:,zslice])/dy)
            else:
                Ey = zeros((phip.shape[0]-2,2,phip.shape[2]-2),'d')
            Ez = ((phip[xslice,yslice,:-1] - phip[xslice,yslice,1:])/dz)

        #self.Ex = Ex
        #self.Ey = Ey
        #self.Ez = Ez
        return Ex,Ey,Ez

    def getintegratedcharge(self,solver):
        dx = solver.dx
        dy = solver.dy
        dz = solver.dz
        if solver is w3d: rhop = solver.rhop
        else:             rhop = solver.sourcep
        Irho = rhop*dx*dy*dz

        return Irho

    def doinjection(self):
        solver = getregisteredsolver()
        if solver is None: solver = w3d
        dx = solver.dx
        dy = solver.dy
        dz = solver.dz
        nxp = solver.nxp
        nyp = solver.nyp
        nzp = solver.nzp

        Ex,Ey,Ez = self.getEfields(solver)
        Qold = self.getintegratedcharge(solver)

        # --- Sum up the integrals of E normal over the sides of the dual cell
        Enorm  = Ex[1:,:,:]*dy*dz
        Enorm -= Ex[:-1,:,:]*dy*dz
        Enorm += Ey[:,1:,:]*dx*dz
        Enorm -= Ey[:,:-1,:]*dx*dz
        Enorm += Ez[:,:,1:]*dx*dy
        Enorm -= Ez[:,:,:-1]*dx*dy
        Enorm *= eps0

        Qnew = Enorm - Qold

        # --- Only inject particle for cells in or near conductors.
        Qnew = where(self.grid.isinside == 0.,0.,Qnew)

        # --- Calculate the number of new particles to add at each grid cell.
        rnn = Qnew/(top.pgroup.sq[self.js]*top.pgroup.sw[self.js])

        # --- Make sure it is positive or zero
        rnn = maximum(rnn,0.)

        # --- Save the number for diagnostics
        self.inj_np = rnn.copy()

        # --- Add a random number to the number of particles injected
        # --- so that the average number of particles injected is
        # --- correct.  For example, if rnn < 1., without the
        # --- addition of the random number, no particles would ever
        # --- be injected.  With the random number, particles will be
        # --- injected and but the average number will be less than 1.
        rnn += where(rnn > 0.,random.random(rnn.shape),0.)

        # --- Now create the particles. This is easiest to do in fortran.
        # --- This also gets the E field at the cell center for each
        # --- of the particles.
        nn = sum(int(rnn))
        if nn == 0: return
        xx,yy,zz = zeros((3,nn),'d')
        ex,ey,ez = zeros((3,nn),'d')
        pp = zeros(nn,'d')
        createparticlesincells(nxp,nyp,nzp,rnn,Ex,Ey,Ez,self.grid.isinside,
                               dx,dy,dz,
                               nn,xx,yy,zz,ex,ey,ez,pp)
        xx += solver.xmmin
        yy += solver.ymmin
        zz += solver.zmminlocal

        # --- Give particles a thermal velocity.
        # --- This now ignores the fact the roughly half the particles will be
        # --- headed back into the conductor.
        vx = random.normal(0.,self.vthermal,nn)
        vy = random.normal(0.,self.vthermal,nn)
        vz = random.normal(0.,self.vthermal,nn)

        # --- Note on the above - a nicer way to get the E field would be to
        # --- get it at each particles position, so that the projection below
        # --- might give smoother results. But, this could be a major problem
        # --- since for many of the particles, the positions chosen above
        # --- would be inside of the conductor so calculating the E field
        # --- would be problematic.

        # --- The E field at the cell centers is used to provide a direction
        # --- for the projection onto the surfaces of the conductors.

        # --- Loop over the conductors, handling all of the particles for each
        # --- conductor at once.
        for c in self.conductors:

            # --- Get the particles near the conductor c
            ii = compress(pp == c.condid,arange(nn))
            if len(ii) == 0: continue
            xc = take(xx,ii)
            yc = take(yy,ii)
            zc = take(zz,ii)
            exc = take(ex,ii)
            eyc = take(ey,ii)
            ezc = take(ez,ii)

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

            # --- For now, as a kludge, in case there are particles that could
            # --- not be projected to the surface, replace the position with the
            # --- original.
            itheta = where(xi == largepos,0.,itheta)
            iphi = where(xi == largepos,0.,iphi)
            xi = where(xi == largepos,xc,xi)
            yi = where(yi == largepos,yc,yi)
            zi = where(zi == largepos,zc,zi)

            # --- Now replace the positions with the projected positions
            put(xx,ii,xi)
            put(yy,ii,yi)
            put(zz,ii,zi)

            # --- Set the velocity so that it is only moving away from the
            # --- surface. Use w3d.l_inj_addtempz_abs by default if the
            # --- flag was not set.
            if self.l_inj_addtempz_abs is None:
                addtempz_abs = w3d.l_inj_addtempz_abs
            else:
                addtempz_abs = self.l_inj_addtempz_abs
            if addtempz_abs:
                # --- The velocity is treated as if it is in the frame relative
                # --- to the surface normal. First, set vz to be positive, and
                # -- then transform the velocity to the lab frame.
                vxc = take(vx,ii)
                vyc = take(vy,ii)
                vzc = take(vz,ii)
                vzc = abs(vzc)
                ct,st = cos(itheta),sin(itheta)
                cp,sp = cos(iphi),sin(iphi)
                vxclab =  ct*cp*vxc + st*vyc - ct*sp*vzc
                vyclab = -st*cp*vxc + ct*vyc + st*sp*vzc
                vzclab =     sp*vxc +             cp*vzc
                put(vx,ii,vxclab)
                put(vy,ii,vyclab)
                put(vz,ii,vzclab)

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

        # --- Get the applied fields.
        exap,eyap,ezap,bx,by,bz = getappliedfields(xx,yy,zz,
                                                   time=top.time,js=self.js)
        ex += exap
        ey += eyap
        ez += ezap

        # --- Give the particles a time step size uniformly distributed
        # --- between 0 and top.dt.
        dt = top.dt*random.random(nn)

        # --- Do a half split leap-frog advance.
        # --- Note that this does the advance in place, directly changing the
        # --- input arrays.
        q = top.pgroup.sq[self.js]
        m = top.pgroup.sm[self.js]
        bpusht3d(nn,ux,uy,uz,gaminv,bx,by,bz,q,m,dt,0.5,top.ibpush)
        epusht3d(nn,ux,uy,uz,ex,ey,ez,q,m,dt,0.5)
        gammaadv(nn,gaminv,ux,uy,uz,top.gamadv,top.lrelativ)
        xpusht3d(nn,xx,yy,zz,ux,uy,uz,gaminv,dt)

        # --- Now add the new particles to the simulation.
        addparticles(xx,yy,zz,ux,uy,uz,gi=gaminv,
                     ex=ex,ey=ey,ez=ez,bx=bx,by=by,bz=bz,
                     js=self.js,lmomentum=true)

    def registerconductors(self,newconductors):
        if not isinstance(newconductors,ListType):
            newconductors = [newconductors]
        for c in newconductors:
            assert c.condid != 0,"The conductor id must be nonzero in order for the particle scraping to work."
            self.conductors.append(c)
        self.updategrid()

    def unregisterconductors(self,conductor,nooverlap=0):
        self.conductors.remove(conductor)
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
        solver = getregisteredsolver()
        if solver is None:
          solver = w3d
          solvertop = top
        else:
          solvertop = solver
        if lparallel: nzlocal = solvertop.nzpslave[me]
        else:         nzlocal = solver.nzlocal
        if (not lforce and (self.grid.nx == solver.nx and
                            self.grid.ny == solver.ny and
                            self.grid.nzlocal == nzlocal and
                            self.grid.xmmin == solver.xmmin and
                            self.grid.xmmax == solver.xmmax and
                            self.grid.ymmin == solver.ymmin and
                            self.grid.ymmax == solver.ymmax and
                            self.grid.zmmin == solver.zmmin and
                            self.grid.zmmax == solver.zmmax and
                            self.grid.izslave[me] == solvertop.izpslave[me])):
          return
        # --- Note that copies of the slave arrays are passed in.
        # --- The arrays in top may be changed the next time loadbalancing is
        # --- done, but the arrays in self.grid should not be changed. Instead,
        # --- a whole new grid is created.
        self.grid = Grid(solver=solver,
                         izslave=solvertop.izpslave.copy(),
                         nzslave=solvertop.nzpslave.copy())
        self.updateconductors()

    def updateconductors(self):
        aura = max(self.grid.dx,self.grid.dy,self.grid.dz)
        for c in self.conductors:
            self.grid.getisinside(c,aura=aura)

