"""LoadBalancer: class wrapping particle load balancing. Sets up automatic
                 periodic load balancing of particles.
"""
__all__ = ['LoadBalancer']
from warp import *
import time

loadbalance_version = "$Id: loadbalance.py,v 1.62 2008/11/19 18:29:59 dave Exp $"

def loadbalancedoc():
    import loadbalance
    print loadbalance.__doc__

#########################################################################
#########################################################################
class LoadBalancer:
    """
Installs load balancer.
Creation arguments:
 - when: dictionary of when to do the load balancing. Keys are time step
         numbers, values are frequency of loadbalancing when top.it is less
         than key. Default is {10:1,100:10,1000000:20}
 - padright,padupperx,paduppery,padupperz=0: 
             Amount of space added to upper end of grid. When not specified,
             it is product of max(v)*top.dt*2 and the number of steps between
             load balances. If not given, the x,y,z values default to padright.
 - padleft,padlowerx,padlowery,padlowerz=0: 
              Amount of space added to lower end of grid. If not given, the
              x,y,z values default to padleft.
 - doloadrho=0: Specifies whether the charge density is recalculated
 - dofs=0: Specifies whether the fields are recalculated
 - verbose=0: Prints output
 - nxguard,nyguard,nzguard=0: Number of extra guard cells to include in the
              field arrays for the particles. Only needed in special cases,
              possibly when the interpolated mover is used since intermediate
              positions in the algorithm may be out of bounds otherwise.
 - spreadx,spready,spreadz=1.: The fraction of processors to spread the work
                               over. Do not use this unless you really know
                               what it means!

Note, if particles only cover a few grid cells, then the distribution is
recalculated on a finer mesh to give better balancing.
    """
    def __init__(self,when=None,padright=None,padleft=None,
                 padupperx=None,paduppery=None,padupperz=None,
                 padlowerx=None,padlowery=None,padlowerz=None,
                 doitnow=0,doloadrho=0,dofs=0,verbose=0,
                 nxguard=0,nyguard=0,nzguard=0,
                 spreadx=1.,spready=1.,spreadz=1.):
        if when is None:
            self.when = {10:1,100:10,1000000:20}
        else:
            self.when = when

        self.padright = padright
        self.padleft = padleft
        if padupperx is None: padupperx = padright
        if paduppery is None: paduppery = padright
        if padupperz is None: padupperz = padright
        self.padupperx = padupperx
        self.paduppery = paduppery
        self.padupperz = padupperz
        if padlowerx is None: padlowerx = padleft
        if padlowery is None: padlowery = padleft
        if padlowerz is None: padlowerz = padleft
        self.padlowerx = padlowerx
        self.padlowery = padlowery
        self.padlowerz = padlowerz

        self.doloadrho = doloadrho
        self.dofs = dofs
        self.verbose = verbose

        self.nxguard = nxguard
        self.nyguard = nyguard
        self.nzguard = nzguard

        self.spreadx = spreadx
        self.spready = spready
        self.spreadz = spreadz

        self.runtime = 0.

        if not lparallel: return

        if doitnow: self.doloadbalance()

        installafterstep(self.doloadbalance)

    def __setstate__(self,dict):
        self.__dict__.update(dict)
        if not isinstalledafterstep(self.doloadbalance):
            installafterstep(self.doloadbalance)

    def doloadbalance(self,lforce=0,doloadrho=None,dofs=None,reorg=None):
        starttime = time.time()

        # --- Set lloadbalanced flag to false. It will be set to true below
        # --- if the load balancing will be done.
        top.lloadbalanced = false

        if not lparallel:
            if self.verbose:
                print "Skipping loadbalance since running in serial"
            endtime = time.time()
            self.runtime += (endtime - starttime)
            return

        # --- Get the current number of live particles.
        if not top.lmoments or top.ifzmmnt == 0 or top.laccumulate_zmoments:
            # --- Sum up the total number of particles calculated.
            nplive = globalsum(top.pgroup.nps)
        else:
            # --- Use the value already calculated.
            nplive = top.pnum[0,-1]

        # --- Check if there are any particles anywhere, and return if not.
        if nplive == 0:
            if self.verbose:
                print "Skipping loadbalance since there are no particles"
            endtime = time.time()
            self.runtime += (endtime - starttime)
            return

        # --- Get the range of particles in each dimension.
        if not top.lmoments or top.ifzmmnt == 0 or top.laccumulate_zmoments:
            # --- If the moments were not calculated, then top.zminp and
            # --- top.zmaxp are not reliable and so need to be calculated.
            # --- In the dimensions where there is no decomposition, skip
            # --- the calculation to not waste time.
            xminp = +largepos
            xmaxp = -largepos
            yminp = +largepos
            ymaxp = -largepos
            zminp = +largepos
            zmaxp = -largepos
            for js in range(top.pgroup.ns):
                if top.pgroup.nps[js] == 0: continue
                i1 = top.pgroup.ins[js] - 1
                i2 = i1 + top.pgroup.nps[js]
                if top.nxprocs > 1:
                    xx = top.pgroup.xp[i1:i2]
                    xminp = min(xminp,min(xx))
                    xmaxp = max(xmaxp,max(xx))
                else:
                    xminp = w3d.xmmin
                    xmaxp = w3d.xmmax
                if top.nyprocs > 1:
                    yy = top.pgroup.yp[i1:i2]
                    yminp = min(yminp,min(yy))
                    ymaxp = max(ymaxp,max(yy))
                else:
                    yminp = w3d.ymmin
                    ymaxp = w3d.ymmax
                if top.nzprocs > 1:
                    zz = top.pgroup.zp[i1:i2]
                    zminp = min(zminp,min(zz))
                    zmaxp = max(zmaxp,max(zz))
                else:
                    zminp = w3d.zmmin
                    zmaxp = w3d.zmmax
            xminp,yminp,zminp = parallelmin([xminp,yminp,zminp])
            xmaxp,ymaxp,zmaxp = parallelmax([xmaxp,ymaxp,zmaxp])
        else:
            # --- Otherwise, use the values already calculated.
            xminp = top.xminp[-1]
            xmaxp = top.xmaxp[-1]
            yminp = top.yminp[-1]
            ymaxp = top.ymaxp[-1]
            zminp = top.zminp[-1]
            zmaxp = top.zmaxp[-1]

        # --- Special check when injection is turned on
        if top.inject:
            # --- Make sure that all of the injection sources are included.
            # --- These are crude estimates of the min and max when xpinject
            # --- and ypinject are nonzero.
            rinj = sqrt(top.ainject**2 + top.binject**2)
            rpinj = sqrt(top.xpinject**2 + top.ypinject**2)
            zinjectmin = max(w3d.zmmin,min(top.zinject - rinj*rpinj))
            zinjectmax = min(w3d.zmmax,max(top.zinject + rinj*rpinj))
            # --- Add in the term accounting for the curvature of the source
            rmax = maximum(top.ainject,top.binject)
            injdepth = max(rmax**2/(top.rinject+sqrt(top.rinject**2-rmax**2)))
            injdepth = max(injdepth,maxnd(w3d.inj_grid))
            if min(top.inj_d) < 0.: zinjectmin = zinjectmin - injdepth
            if max(top.inj_d) > 0.: zinjectmax = zinjectmax + injdepth
            # --- Also make sure that the injection virtual surface is included.
            zinjectmin = zinjectmin + min(0.,min(top.inj_d)*w3d.dz)
            zinjectmax = zinjectmax + max(0.,max(top.inj_d)*w3d.dz)
            zminp = minimum(zinjectmin,zminp)
            zmaxp = maximum(zinjectmax,zmaxp)
            # --- Transverse dimensions
            xminp = minimum(xminp,max(w3d.xmmin,min(top.xinject-rmax)))
            xmaxp = maximum(xmaxp,min(w3d.xmmax,min(top.xinject+rmax)))
            yminp = minimum(yminp,max(w3d.ymmin,min(top.yinject-rmax)))
            ymaxp = maximum(ymaxp,min(w3d.ymmax,min(top.yinject+rmax)))

        if top.tinject:
            zinjectmin = min(top.ztinjmn) - w3d.dz
            zinjectmax = max(top.ztinjmx) + w3d.dz
            if w3d.solvergeom in [w3d.XYZgeom,w3d.XZgeom]:
              maxa = maximum.reduce(top.atinjectz,axis=0)
              maxb = maximum.reduce(top.btinjectz,axis=0)
              xinjectmin = minnd(-maxa+top.xtinject) - 2*w3d.inj_dx
              xinjectmax = maxnd(+maxa+top.xtinject) + 2*w3d.inj_dx
              yinjectmin = minnd(-maxb+top.ytinject) - 2*w3d.inj_dy
              yinjectmax = maxnd(+maxb+top.ytinject) + 2*w3d.inj_dy
            elif w3d.solvergeom == w3d.RZgeom:
              maxa = maximum.reduce(top.atinjectz,axis=0)
              mina = minimum.reduce(top.atinjectz,axis=0)
              xinjectmin = minnd(mina) - 2*w3d.inj_dx
              xinjectmax = maxnd(maxa) + 2*w3d.inj_dx
              yinjectmin = 0.
              yinjectmax = 0.
            xminp = minimum(xminp,max(w3d.xmmin,xinjectmin))
            xmaxp = maximum(xmaxp,min(w3d.xmmax,xinjectmax))
            yminp = minimum(yminp,max(w3d.ymmin,yinjectmin))
            ymaxp = maximum(ymaxp,min(w3d.ymmax,yinjectmax))
            zminp = minimum(zminp,max(w3d.zmmin,zinjectmin))
            zmaxp = maximum(zmaxp,min(w3d.zmmax,zinjectmax))

        # --- Shift into the grid frame
        xminp = xminp - w3d.xmmin
        xmaxp = xmaxp - w3d.xmmin
        yminp = yminp - w3d.ymmin
        ymaxp = ymaxp - w3d.ymmin
        zminp = zminp - w3d.zmmin - top.zbeam
        zmaxp = zmaxp - w3d.zmmin - top.zbeam

        ppdecomp = top.ppdecomp

        # --- Check if uppermost particle is close to edge of last processor
        # --- If so, then force a reloadbalance.
        if ppdecomp.nxprocs > 1 and ppdecomp.xmax[-1] < w3d.xmmax-0.5*w3d.dx:
            if xmaxp > ppdecomp.xmax[-1]-2*w3d.dx:
                lforce = true
                if self.verbose:
                    print "Load balancing since particles near upper end ",
                    print "of mesh in x ",ppdecomp.xmax[-1],w3d.xmmax,xmaxp,
                    print ppdecomp.xmax[-1]-2*w3d.dx

        # --- Check if lowermost particle is close to edge of last processor
        # --- If so, then force a reloadbalance.
        if ppdecomp.nxprocs > 1 and ppdecomp.xmin[0] > w3d.xmmin+0.5*w3d.dx:
            if xminp < ppdecomp.xmin[0]+2*w3d.dx:
                lforce = true
                if self.verbose:
                    print "Load balancing since particles near lower end ",
                    print "of mesh in x ",ppdecomp.xmin[0],w3d.xmmin,xminp,
                    print ppdecomp.xmin[0]+2*w3d.dx

        # --- Check if uppermost particle is close to edge of last processor
        # --- If so, then force a reloadbalance.
        if ppdecomp.nyprocs > 1 and ppdecomp.ymax[-1] < w3d.ymmax-0.5*w3d.dy:
            if ymaxp > ppdecomp.ymax[-1]-2*w3d.dy:
                lforce = true
                if self.verbose:
                    print "Load balancing since particles near upper end ",
                    print "of mesh in y ",ppdecomp.ymax[-1],w3d.ymmax,ymaxp,
                    print ppdecomp.ymax[-1]-2*w3d.dy

        # --- Check if lowermost particle is close to edge of last processor
        # --- If so, then force a reloadbalance.
        if ppdecomp.nyprocs > 1 and ppdecomp.ymin[0] > w3d.ymmin+0.5*w3d.dy:
            if yminp < ppdecomp.ymin[0]+2*w3d.dy:
                lforce = true
                if self.verbose:
                    print "Load balancing since particles near lower end ",
                    print "of mesh in y ",ppdecomp.ymin[0],w3d.ymmin,yminp,
                    print ppdecomp.ymin[0]+2*w3d.dy

        # --- Check if uppermost particle is close to edge of last processor
        # --- If so, then force a reloadbalance.
        if ppdecomp.nzprocs > 1 and ppdecomp.zmax[-1] < w3d.zmmax-0.5*w3d.dz:
            if zmaxp > ppdecomp.zmax[-1]-2*w3d.dz + top.zbeam:
                lforce = true
                if self.verbose:
                    print "Load balancing since particles near upper end ",
                    print "of mesh in z ",ppdecomp.zmax[-1],w3d.zmmax,zmaxp,
                    print ppdecomp.zmax[-1]-2*w3d.dz

        # --- Check if lowermost particle is close to edge of last processor
        # --- If so, then force a reloadbalance.
        if ppdecomp.nzprocs > 1 and ppdecomp.zmin[0] > w3d.zmmin+0.5*w3d.dz:
            if zminp < ppdecomp.zmin[0]+2*w3d.dz + top.zbeam:
                lforce = true
                if self.verbose:
                    print "Load balancing since particles near lower end ",
                    print "of mesh in z ",ppdecomp.zmin[0],w3d.zmmin,zminp,
                    print ppdecomp.zmin[0]+2*w3d.dz

        # --- Find frequency of load balancing
        ii = max(self.when.values())
        for key,value in self.when.items():
            if top.it < key: ii = min(ii,value)

        # --- Just return if load balancing not done now.
        if not lforce and (top.it%ii) != 0:
            if self.verbose:
                print "Skipping loadbalance since it is not time for it"
            endtime = time.time()
            self.runtime += (endtime - starttime)
            return

        if (top.it%ii) == 0 and self.verbose:
            print "Load balancing based on frequency"

        if top.nxprocs > 1:
            self.dodecomposition(0,ii,xminp,xmaxp,self.spreadx,
                                 self.padlowerx,self.padupperx,
                                 w3d.xmmin,w3d.xmmax,w3d.dx,0.,top.nxprocs,
                                 top.pgroup.getpyobject('xp'),
                                 top.pgroup.getpyobject('uxp'),
                                 ppdecomp.nxglobal,self.nxguard,
                                 ppdecomp.xmin,ppdecomp.xmax,
                                 ppdecomp.ix,ppdecomp.nx)
            top.xpminlocal = ppdecomp.xmin[top.ixproc]
            top.xpmaxlocal = ppdecomp.xmax[top.ixproc]
            w3d.xmminp = w3d.xmmin + ppdecomp.ix[top.ixproc]*w3d.dx
            w3d.xmmaxp = w3d.xmmin + (ppdecomp.ix[top.ixproc] +
                                      ppdecomp.nx[top.ixproc])*w3d.dx

        if top.nyprocs > 1:
            self.dodecomposition(1,ii,yminp,ymaxp,self.spready,
                                 self.padlowery,self.paduppery,
                                 w3d.ymmin,w3d.ymmax,w3d.dy,0.,top.nyprocs,
                                 top.pgroup.getpyobject('yp'),
                                 top.pgroup.getpyobject('uyp'),
                                 ppdecomp.nyglobal,self.nyguard,
                                 ppdecomp.ymin,ppdecomp.ymax,
                                 ppdecomp.iy,ppdecomp.ny)
            top.ypminlocal = ppdecomp.ymin[top.iyproc]
            top.ypmaxlocal = ppdecomp.ymax[top.iyproc]
            w3d.ymminp = w3d.ymmin + ppdecomp.iy[top.iyproc]*w3d.dy
            w3d.ymmaxp = w3d.ymmin + (ppdecomp.iy[top.iyproc] +
                                      ppdecomp.ny[top.iyproc])*w3d.dy

        if top.nzprocs > 1:
            self.dodecomposition(2,ii,zminp,zmaxp,self.spreadz,
                                 self.padlowerz,self.padupperz,
                                 w3d.zmmin,w3d.zmmax,w3d.dz,0.,top.nzprocs,
                                 top.pgroup.getpyobject('zp'),
                                 top.pgroup.getpyobject('uzp'),
                                 ppdecomp.nzglobal,self.nzguard,
                                 ppdecomp.zmin,ppdecomp.zmax,
                                 ppdecomp.iz,ppdecomp.nz)
            top.zpminlocal = ppdecomp.zmin[top.izproc]
            top.zpmaxlocal = ppdecomp.zmax[top.izproc]
            w3d.zmminp = w3d.zmmin + ppdecomp.iz[top.izproc]*w3d.dz
            w3d.zmmaxp = w3d.zmmin + (ppdecomp.iz[top.izproc] +
                                      ppdecomp.nz[top.izproc])*w3d.dz

            top.izpslave[:] = ppdecomp.iz
            top.nzpslave[:] = ppdecomp.nz
            top.zpslmin[:] =  ppdecomp.zmin
            top.zpslmax[:] =  ppdecomp.zmax

        # --- Reorganize the particles
        # --- On step zero, a complete reorganization is done so the reorg flag
        # --- is set to true to use the particle sorter which is more efficient
        # --- in that case.
        if reorg is None:
          reorg = (top.it==1)
        if reorg:
          reorgparticles(top.pgroup,w3d.l4symtry,w3d.l2symtry,
                         w3d.solvergeom==w3d.RZgeom)
        else:
          particlegridboundaries3d(top.pgroup,-1)

        # --- Update sizes of grids for particles
        w3d.nxp = ppdecomp.nx[top.ixproc]
        w3d.nyp = ppdecomp.ny[top.iyproc]
        w3d.nzp = ppdecomp.nz[top.izproc]
        if dofs is None: dofs = self.dofs
        solver = getregisteredsolver()
        if solver is not None:
            try:
                solver.resetparticledomains()
            except AttributeError:
                print "Field solver does not have a setparticledomains method"
            # --- Zero out the source that is used for the fieldsolver. This is
            # --- done in case some region of source is no longer covered by
            # --- sourcep.
            try:
                solver.zerosource()
            except AttributeError:
                print "Field solver does not have a zerosource method"
        else:
            if(w3d.solvergeom == w3d.XYZgeom):
                # --- Allocate space with updated nxp, nyp and nzp
                gchange("Fields3dParticles")
            else:
                gchange_rhop_phip_rz()
            # --- Redistribute phi to the particle arrays if a field solve is
            # --- not done.
            if not dofs:
                if getregisteredsolver() is None:
                    for i in range(getnsndtsforsubcycling()):
                        getphipforparticles(i)

        # --- Do some additional work if requested
        if doloadrho is None: doloadrho = self.doloadrho
        if doloadrho: loadrho()
        if dofs: fieldsol(0)

        top.lloadbalanced = true

        endtime = time.time()
        self.runtime += (endtime - starttime)

    def dodecomposition(self,axis,ii,minp,maxp,spread,padlower,padupper,
                   mmin,mmax,dd,beam,nprocs,pp,uu,
                   nnglobal,nguard,ppdecompmin,ppdecompmax,
                   ppdecompii,ppdecompnn):
        if (axis < 2 or (maxp - minp)/dd < 10 or
            not top.lmoments or top.ifzmmnt == 0 or top.laccumulate_zmoments):
            # --- If the particles only extend over a few grid cells,
            # --- recalculate the distribution on a finer grid to get better
            # --- loadbalancing.
            # --- Also, calculate the distribution if the moments were not
            # --- calculated on the most recent step.
            pnum = zeros(1001,'d')
            pmin = max(0.,minp-dd)
            pmax = min(mmax-mmin,maxp+dd)
            for js in range(top.pgroup.ns):
                if top.pgroup.nps[js] == 0: continue
                i1 = top.pgroup.ins[js] - 1
                i2 = i1 + top.pgroup.nps[js]
                setgrid1d(top.pgroup.nps[js],pp[i1:i2],1000,pnum,
                          pmin+beam+mmin,pmax+beam+mmin)
            pnum = parallelsum(pnum)
            pdd = (pmax - pmin)/1000.
        else:
            # --- Otherwise use the already calculated z-moment
            pnum = top.pnumz[:,-1]
            pmin = 0.
            pdd = w3d.dz

        assert max(pnum) > 0.,"No particles found during decomposition"

        # --- Add fictitious data so that actual work is spread only to the
        # --- requested fraction of the processors.
        assert (0. < spread <= 1.),"spread must be between 0 and 1 or 1."
        avepnum = ave(pnum)
        pnum = pnum + avepnum*(1./spread - 1.)

        # --- Convert the number of particles to a decomposition
        domain = decompose(pnum,nprocs,lfullcoverage=0)
        domain = domain*pdd + pmin
        domain[0] = min(domain[0],minp)
        domain[-1] = max(domain[-1],maxp)

        # --- Set domain of each processor.
        ppdecompmin[:] = mmin + domain[:-1]
        ppdecompmax[:] = mmin + domain[1:]

        padlower = self.calcpadlower(axis,ii,padlower,uu,dd)
        padupper = self.calcpadupper(axis,ii,padupper,uu,dd)

        ppdecompmin[0] = max(mmin,ppdecompmin[0] - padlower)
        ppdecompmax[-1] = min(mmax,ppdecompmax[-1] + padupper)

        domaindecomposeparticles(nnglobal,nprocs,nguard,mmin,mmax,dd,
                                 zeros(nprocs,'d'),true,
                                 ppdecompii,ppdecompnn,ppdecompmin,ppdecompmax)

    def calcpadupper(self,axis,ii,padupper,vv,dd):
        # --- Calculate the padding on the upper edge.
        if padupper is None:
            if axis < 2 or not top.lmoments or top.ifzmmnt == 0:
                vmaxp = -largepos
                for js in range(top.pgroup.ns):
                    if top.pgroup.nps[js] == 0: continue
                    i1 = top.pgroup.ins[js] - 1
                    i2 = i1 + top.pgroup.nps[js]
                    vv = vv[i1:i2]*top.pgroup.gaminv[i1:i2]
                    vmaxp = max(vmaxp,max(vv))
                vmaxp = globalmax(vmaxp)
            else:
                vmaxp = max(top.vzmaxp)
            if vmaxp > 0.: padupper = vmaxp*top.dt*ii*2
            else:          padupper = ii*dd
        if self.verbose:
            print "Load balancing padupper%s = "%(['x','y','z'][axis]),padupper
        return padupper

    def calcpadlower(self,axis,ii,padlower,vv,dd):
        # --- Calculate the padding on the lower edge.
        if padlower is None:
            if axis < 2 or not top.lmoments or top.ifzmmnt == 0:
                vminp = +largepos
                for js in range(top.pgroup.ns):
                    if top.pgroup.nps[js] == 0: continue
                    i1 = top.pgroup.ins[js] - 1
                    i2 = i1 + top.pgroup.nps[js]
                    vv = vv[i1:i2]*top.pgroup.gaminv[i1:i2]
                    vminp = min(vminp,min(vv))
                vminp = globalmin(vminp)
            else:
                vminp = min(top.vzminp)
            if vminp < 0.: padlower = -vminp*top.dt*ii*2
            else:          padlower = ii*dd
        if self.verbose:
            print "Load balancing padlower%s = "%(['x','y','z'][axis]),padlower
        return padlower

#########################################################################
#########################################################################
def setparticledomains(zslave,doloadrho=1,dofs=1,padleft=0.,padright=0.,
                       reorg=0,nxguard=0,nyguard=0,nzguard=0):
    """
Sets the particles domains from the input, zslave, in the same way as done
with top.zslave during the generate. This is only meant to be used after
that has already been done.
 - zslave: list of starting locations of the domains in grid cell units
 - doloadrho=1: when true, the charge density is redeposited
 - dofs=1: when true, the fields are recalculated
 - padleft=0, padright=0: extra space added on to leftmost and rightmost
                          domains (up to edge of system) in units of meters
 - reorg=0: when true, call reorg_particles which is fastest when particles
            are to be shifted across multiple processors, otherwise use
            zpartbnd which is fastest when particles are to be shifted only to
            nearest neighbors.
 - nzguard=0: Number of extra guard cells to include in the field arrays for
              the particles.
    """
    if not lparallel: return
    # --- It is assumed that the user supplied decomposition is specified
    # --- in the array zslave.

    # --- All values of zslave must be > 0.
    assert min(zslave[1:]-zslave[:-1]) > 0.,"The length of all particle domains must be > 0."

    ppdecomp = top.ppdecomp

    # --- Set domain of each processor.
    for i in range(npes):
        ppdecomp.zmin[i] = w3d.zmmin + zslave[i]*w3d.dz
        ppdecomp.zmax[i] = w3d.zmmin + zslave[i+1]*w3d.dz

    ppdecomp.zmin[0] = max(w3d.zmmin,ppdecomp.zmin[0] - padleft)
    ppdecomp.zmax[-1] = min(w3d.zmmax,ppdecomp.zmax[-1] + padright)

    """
    domaindecomposeparticles(ppdecomp.nxglobal,ppdecomp.nxprocs,nxguard,
                             w3d.xmmin,w3d.xmmax,w3d.dx,
                             zeros(ppdecomp.nxprocs,'d'),true,
                             ppdecomp.ix,ppdecomp.nx,
                             ppdecomp.xmin,ppdecomp.xmax)

    domaindecomposeparticles(ppdecomp.nyglobal,ppdecomp.nyprocs,nyguard,
                             w3d.ymmin,w3d.ymmax,w3d.dy,
                             zeros(ppdecomp.nyprocs,'d'),true,
                             ppdecomp.iy,ppdecomp.ny,
                             ppdecomp.ymin,ppdecomp.ymax)
    """

    domaindecomposeparticles(ppdecomp.nzglobal,ppdecomp.nzprocs,nzguard,
                             w3d.zmmin,w3d.zmmax,w3d.dz,
                             zeros(ppdecomp.nzprocs,'d'),true,
                             ppdecomp.iz,ppdecomp.nz,
                             ppdecomp.zmin,ppdecomp.zmax)

    top.zpminlocal = ppdecomp.zmin[me]
    top.zpmaxlocal = ppdecomp.zmax[me]

    top.izpslave[:] = ppdecomp.iz
    top.nzpslave[:] = ppdecomp.nz
    top.zpslmin[:] =  ppdecomp.zmin
    top.zpslmax[:] =  ppdecomp.zmax

    # --- Reorganize the particles
    if reorg:
        reorgparticles(top.pgroup,w3d.l4symtry,w3d.l2symtry,
                       w3d.solvergeom==w3d.RZgeom)
    else:
        particlegridboundaries3d(top.pgroup,-1)

    # --- Update sizes of grids for particles
    solver = getregisteredsolver()
    if solver is not None:
        try:
            solver.resetparticledomains()
        except AttributeError:
            print "Field solver does not have a setparticledomains method"
        # --- Zero out the source that is used for the fieldsolver. This is
        # --- done in case some region of source is no longer covered by
        # --- sourcep.
        try:
            solver.zerosource()
        except AttributeError:
            print "Field solver does not have a zerosource method"
    else:
        if(w3d.solvergeom == w3d.XYZgeom):
            w3d.nzp = ppdecomp.nz[me]
            w3d.zmminp = w3d.zmmin + ppdecomp.iz[me]*w3d.dz
            w3d.zmmaxp = w3d.zmminp + w3d.nzp*w3d.dz
            gchange("Fields3dParticles")
        else:
            gchange_rhop_phip_rz()

    # --- Do some additional work if requested
    if doloadrho: loadrho()
    if dofs: fieldsol(0)


#########################################################################
def loadbalanceparticles(doloadrho=1,dofs=1,spread=1.,padleft=0.,padright=0.,
                         reorg=0,pnumz=None,zmin=None,dz=None,
                         zminp=None,zmaxp=None,verbose=0,
                         nxguard=0,nyguard=0,nzguard=0):
    """
Load balances the particles as evenly as possible. The load balancing is
based off of the data in top.pnumz which of course must already have
been calculated. The number density is assumed to vary linearly between
grid points.
 - doloadrho=1: when true, the charge density is redoposited
 - dofs=1: when true, the fields are recalculated
 - spread=1.: fraction of processors to spread the work among
 - padleft=0, padright=0: extra space added on to leftmost and rightmost
                          domains (up to edge of system) in units of meters
 - reorg=0: when true, call reorg_particles which is fastest when particles
            are to be shifted across multiple processors, otherwise use
            zpartbnd which is fastest when particles are to be shifted only to
            nearest neighbors.
 - pnumz=top.pznum: the particle distribution to base the load balancing on
 - zmin=None: optional z-minimum of the grid
 - dz=None: optional grid cell size
 - zminp,zmaxp=None: optional min and max of the region that must be included
                     in the decomposition
 - verbose=0: when true, prints out timing information
 - nzguard=0: Number of extra guard cells to include in the field arrays for
              the particles.
    """

    # --- Convert the number of particles to a decomposition
    zslave = decompose(pnumz,npes,lfullcoverage=0)

    # --- Scale to specified grid if zmin and/or dz input
    if dz is not None: zslave = zslave*dz
    if zmin is not None: zslave = zslave + zmin
    if zminp is not None: zslave[0] = min(zslave[0],zminp)
    if zmaxp is not None: zslave[-1] = max(zslave[-1],zmaxp)
    if dz is not None: zslave = zslave/w3d.dz

    # --- Apply the new domain decomposition.
    setparticledomains(zslave,doloadrho=doloadrho,dofs=dofs,
                       padleft=padleft,padright=padright,reorg=reorg,
                       nxguard=nxguard,nyguard=nyguard,nzguard=nzguard)
    endtime = wtime()
    if verbose: print "Load balance time = ",endtime - starttime

#########################################################################
def loadbalancesor(sgweight=7.0,condweight=2.0):
    """
Load balance the SOR field solver based off of the current timings. This is
needed since some processors may have more conductor points than others.
 - sgweight=7.: weight (in timing) of subgrid points relative to weight of a
                grid cell
 - condweight=2.: weight (in timing) of a conductor points relative to weight
                  of a grid cell
    """
    if not lparallel: return
    # --- Save the old values
    oldizfs = top.izfsslave + 0
    oldnzfs = top.nzfsslave + 0
    oldphi = w3d.phi + 0.
    oldrho = w3d.rho + 0.

    # --- Make sure that the conductor arrays are allocated.
    if f3d.ncondmax == 0: f3d.ncondmax = 1
    if f3d.ncndmax == 0: f3d.ncndmax = 1
    gchange("Conductor3d")
      
    # --- Gather the field solve weights. For each z plane, sum the number of
    # --- grid cells, subgrid points, and conductor points, appropriately
    # --- weighted.
    weight = zeros(top.nzfsslave[me]+1,'d')
    for iz in iota(w3d.izfsmin,w3d.izfsmax):
        nec = len(oldnonzero(logical_not(f3d.iecndz[:f3d.necndbdy]-iz)))
        noc = len(oldnonzero(logical_not(f3d.iocndz[:f3d.nocndbdy]-iz)))
        nc  = len(oldnonzero(logical_not(f3d.izcond[:f3d.ncond]-iz)))
        weight[iz-w3d.izfsmin] = ((w3d.nx+1)*(w3d.ny+1) +
                                  sgweight*(nec + noc) +
                                  condweight*nc)
    weight = gatherallzfsarray(weight)

    # --- Convert to a decomposition
    zslave = decompose(weight,npes,lfullcoverage=1)

    # --- Set domain of each processor.
    # --- This coding ensures that all of the processors have nzfsslave at least
    # --- 2 or greater and that the last processor isn't left with too few cells.
    # --- In cases where all of the processors have values only slightly
    # --- greater than 2, the last processor will likely end up with too few
    # --- cells because of the accumulation of rounding up. Find the processors
    # --- which have the largest amount of roundup and take away one of their
    # --- grid cells until the last processor has enough. Do the same for the
    # --- case where the last processor has too many.
    top.izfsslave[0] = 0
    top.nzfsslave[0] = max(nint(zslave[1]) + 1,2)
    for i in range(1,npes):
        top.izfsslave[i] = top.izfsslave[i-1] + top.nzfsslave[i-1] - 1
        top.nzfsslave[i] = max(nint(zslave[i+1]-zslave[i]) + 1,2)
    top.nzfsslave[-1] = w3d.nz - top.izfsslave[-1]
    while ((zslave[-1]-top.nzfsslave[-1]) > max(zslave[:-1]-top.nzfsslave[:-1])
            or top.nzfsslave[-1] < 2):
        i = argmax(where(greater(top.nzfsslave[:-1],2),
                         top.nzfsslave[:-1]-zslave[:-1],-10000.))
        top.nzfsslave[i] = top.nzfsslave[i] - 1
        top.izfsslave[i+1:] = top.izfsslave[i+1:] - 1
        top.nzfsslave[-1] = w3d.nz - top.izfsslave[-1]
    while (top.nzfsslave[-1]-zslave[-1]) > max(top.nzfsslave[:-1]-zslave[:-1]):
        i = argmax(zslave[:-1]-top.nzfsslave[:-1])
        top.nzfsslave[i] = top.nzfsslave[i] + 1
        top.izfsslave[i+1:] = top.izfsslave[i+1:] + 1
        top.nzfsslave[-1] = w3d.nz - top.izfsslave[-1]

    # --- Adjust the Z data
    _adjustz()

    # --- Shift the existing charge density and phi
    izstart = max(oldizfs[me],top.izfsslave[me])
    izend = min(oldizfs[me]+oldnzfs[me],top.izfsslave[me]+top.nzfsslave[me])
    newiz1 = izstart - top.izfsslave[me]
    newiz2 = izend - top.izfsslave[me] + 1
    oldiz1 = izstart - oldiz[me]
    oldiz2 = izend - oldiz[me] + 1
    w3d.phi[:,:,newiz1+1:newiz2+1] = oldphi[:,:,oldiz1+1:oldiz2+1]
    w3d.rho[:,:,newiz1:newiz2] = oldrho[:,:,oldiz1:oldiz2]

    # --- Correct the locations of conductor points for the field-solver.
    newiz = top.izfsslave
    newnz = top.nzfsslave
    newizfs = top.izfsslave
    newnzfs = top.nzfsslave
    reorgconductors(oldiz,oldnz,oldizfs,oldnzfs,
                    newiz,newnz,newizfs,newnzfs)

    # --- Correct location of injection source.
    if top.inject > 0:
        w3d.inj_grid[:,:,:] = w3d.inj_grid + oldiz[me] - newiz[me]

#########################################################################
def decompose(weight,npes,lfullcoverage=0):
    """
Converts a weight into the size of the domains.
 - weight: array of relative weights of the work done by each processor
 - npes: number of processors
 - lfullcoverage=0: when true, the domains cover the full extent of
                    the system
Returns an array of the same length which is the relative length of each
of the domains.
    """
    assert max(weight) > 0.,"weight must have some positive elements"
    # --- Integrate weight, assuming linear variation between grid points
    nn = len(weight) - 1
    np = 0.5*weight[0] + sum(weight[1:-1]) + 0.5*weight[-1]
    npperpe = 1.*np/npes

    domain = zeros(npes+1,'d')
    ii = 0
    if not lfullcoverage:
        # --- First first non-zero weight, making sure to check first cell too.
        while weight[ii] == 0. and weight[ii+1] == 0.: ii = ii + 1
    delta = 0.
    domain[0] = ii
    for ip in xrange(1,npes):
        fract = 0.
        npint = 0.
        npnext = (weight[ii  ]*((1.-delta)+0.5*(delta**2-1.)) +
                  weight[ii+1]*0.5*(1. - delta**2))
        # --- Get the remaining bit from the previous cell if it is not too much.
        if npnext < npperpe:
            fract = 1. - delta
            ii = ii + 1
            delta = 0.
            npint = npnext
        # --- Keep adding cells until the number per processor is reached.
        while npint + 0.5*(weight[ii]+weight[ii+1]) < npperpe:
            fract = fract + 1.
            delta = 0.
            npint = npint + 0.5*(weight[ii]+weight[ii+1])
            ii = ii + 1
            if ii == nn+1: break
        if ii == nn+1: break
        # --- Add the last little bit to get to exactly npperpe.
        delta1 = delta
        a = 0.5*weight[ii] - 0.5*weight[ii+1]
        b = weight[ii]
        c = (weight[ii]*(delta1 - 0.5*delta1**2) + 0.5*weight[ii+1]*delta1**2 +
             npperpe - npint)
        if b != 0.:
            delta = 2.*c/(sqrt(b**2 - 4.*a*c) + b)
        else:
            delta = sqrt(-c/a)
        fract = fract + delta - delta1
        domain[ip] = domain[ip-1] + fract

    # --- Set the end of the last domain
    if not lfullcoverage:
        # --- Find the last place with non-zero weight, and give the last processor
        # --- everything up to that point.
        for ii in xrange(ii,nn):
          if weight[ii] > 0.: domain[-1] = ii+1
    else:
        domain[-1] = nn
      
    return domain

#########################################################################
def _adjustz():

    #---------------------------------------------------------------------------
    # --- Reset local values
    w3d.nzlocal = top.nzfsslave[me]
    zpmin = w3d.zmmin + top.ppdecomp.iz[me]*w3d.dz
    zpmax = (top.ppdecomp.iz[me]+top.ppdecomp.nz[me])*w3d.dz + w3d.zmmin
    w3d.izfsmin = 0.
    w3d.izfsmax = top.nzfsslave[me]
    w3d.zmminlocal = top.izfsslave[me]*w3d.dz + w3d.zmmin
    w3d.zmmaxlocal = (top.izfsslave[me]+top.nzfsslave[me])*w3d.dz+w3d.zmmin

    # --- Change the alocation of everything effected are reset the meshes.
    gchange("Fields3d")
    gchange("Z_Moments")
    gchange("Hist")
    w3d.zmesh[:] = w3d.zmmin + iota(0,w3d.nz)*w3d.dz
    w3d.zmeshlocal[:] = w3d.zmminlocal + iota(0,w3d.nzlocal)*w3d.dz
    
    # --- Reset the lattice
    setlatt()

#########################################################################
#########################################################################
# --- These are the messy routines for reorganizing the conductor data
def reorgconductors(oldiz,oldnz,oldizfs,oldnzfs,
                    newiz,newnz,newizfs,newnzfs):
    if globalsum(f3d.ncond) > 0:
        # --- Make things easier to deal with by ensuring that all arrays
        # --- are allocated.
        f3d.ncondmax = f3d.ncond + 1
        gchange("Conductor3d")

        # --- Shift the data to be relative to the global system
        f3d.izcond[:] = f3d.izcond[:] + oldiz[me]

        # --- Do the work
        results = _reorgconductorarrays([f3d.ixcond[:f3d.ncond],
                                         f3d.iycond[:f3d.ncond],
                                         f3d.izcond[:f3d.ncond],
                                         f3d.condvolt[:f3d.ncond],
                                         f3d.condnumb[:f3d.ncond]],
                                        f3d.izcond[:f3d.ncond]+0,
                                        oldiz,oldnz,oldizfs,oldnzfs,
                                        newiz,newnz,newizfs,newnzfs)

        # --- Change array sizes and copy the data, localizing it.
        f3d.ncond = len(results[0])
        f3d.ncondmax = f3d.ncond
        gchange("Conductor3d")
        if f3d.ncond > 0:
            f3d.ixcond[:] = results[0]
            f3d.iycond[:] = results[1]
            f3d.izcond[:] = results[2] - newiz[me]
            f3d.condvolt[:] = results[3]
            f3d.condnumb[:] = results[4]

    if globalsum(f3d.necndbdy) > 0:
        # --- Make things easier to deal with by ensuring that all arrays
        # --- are allocated.
        f3d.ncndmax = max(f3d.necndbdy,f3d.nocndbdy) + 1
        gchange("Conductor3d")

        # --- Shift the data to be relative to the global system
        f3d.iecndz[:] = f3d.iecndz[:] + oldiz[me]

        # --- Do the work
        results = _reorgconductorarrays([f3d.iecndx[:f3d.necndbdy],
                                         f3d.iecndy[:f3d.necndbdy],
                                         f3d.iecndz[:f3d.necndbdy],
                                         f3d.ecdelmx[:f3d.necndbdy],
                                         f3d.ecdelmy[:f3d.necndbdy],
                                         f3d.ecdelmz[:f3d.necndbdy],
                                         f3d.ecdelpx[:f3d.necndbdy],
                                         f3d.ecdelpy[:f3d.necndbdy],
                                         f3d.ecdelpz[:f3d.necndbdy],
                                         f3d.ecvolt[:f3d.necndbdy],
                                         f3d.ecnumb[:f3d.necndbdy]],
                                        f3d.iecndz[:f3d.necndbdy]+0,
                                        oldiz,oldnz,oldizfs,oldnzfs,
                                        newiz,newnz,newizfs,newnzfs)

        # --- Change array sizes and copy the data, localizing it.
        f3d.necndbdy = len(results[0])
        if f3d.necndbdy > f3d.ncndmax:
            f3d.ncndmax = max(f3d.necndbdy,f3d.nocndbdy)
            gchange("Conductor3d")
        if f3d.necndbdy > 0:
            f3d.iecndx[:f3d.necndbdy] = results[0]
            f3d.iecndy[:f3d.necndbdy] = results[1]
            f3d.iecndz[:f3d.necndbdy] = results[2] - newiz[me]
            f3d.ecdelmx[:f3d.necndbdy] = results[3]
            f3d.ecdelmy[:f3d.necndbdy] = results[4]
            f3d.ecdelmz[:f3d.necndbdy] = results[5]
            f3d.ecdelpx[:f3d.necndbdy] = results[6]
            f3d.ecdelpy[:f3d.necndbdy] = results[7]
            f3d.ecdelpz[:f3d.necndbdy] = results[8]
            f3d.ecvolt[:f3d.necndbdy] = results[9]
            f3d.ecnumb[:f3d.necndbdy] = results[10]

    if globalsum(f3d.nocndbdy) > 0:
        # --- Make things easier to deal with by ensuring that all arrays
        # --- are allocated.
        f3d.ncndmax = max(f3d.necndbdy,f3d.nocndbdy) + 1
        gchange("Conductor3d")

        # --- Shift the data to be relative to the global system
        f3d.iocndz[:] = f3d.iocndz[:] + oldiz[me]

        # --- Do the work
        results = _reorgconductorarrays([f3d.iocndx[:f3d.nocndbdy],
                                         f3d.iocndy[:f3d.nocndbdy],
                                         f3d.iocndz[:f3d.nocndbdy],
                                         f3d.ocdelmx[:f3d.nocndbdy],
                                         f3d.ocdelmy[:f3d.nocndbdy],
                                         f3d.ocdelmz[:f3d.nocndbdy],
                                         f3d.ocdelpx[:f3d.nocndbdy],
                                         f3d.ocdelpy[:f3d.nocndbdy],
                                         f3d.ocdelpz[:f3d.nocndbdy],
                                         f3d.ocvolt[:f3d.nocndbdy],
                                         f3d.ocnumb[:f3d.nocndbdy]],
                                        f3d.iocndz[:f3d.nocndbdy]+0,
                                        oldiz,oldnz,oldizfs,oldnzfs,
                                        newiz,newnz,newizfs,newnzfs)

        # --- Change array sizes and copy the data, localizing it.
        f3d.nocndbdy = len(results[0])
        if f3d.nocndbdy > f3d.ncndmax:
            f3d.ncndmax = max(f3d.necndbdy,f3d.nocndbdy)
            gchange("Conductor3d")
        if f3d.nocndbdy > 0:
            f3d.iocndx[:f3d.nocndbdy] = results[0]
            f3d.iocndy[:f3d.nocndbdy] = results[1]
            f3d.iocndz[:f3d.nocndbdy] = results[2] - newiz[me]
            f3d.ocdelmx[:f3d.nocndbdy] = results[3]
            f3d.ocdelmy[:f3d.nocndbdy] = results[4]
            f3d.ocdelmz[:f3d.nocndbdy] = results[5]
            f3d.ocdelpx[:f3d.nocndbdy] = results[6]
            f3d.ocdelpy[:f3d.nocndbdy] = results[7]
            f3d.ocdelpz[:f3d.nocndbdy] = results[8]
            f3d.ocvolt[:f3d.nocndbdy] = results[9]
            f3d.ocnumb[:f3d.nocndbdy] = results[10]

#-------------------------------------------------------------------------
def _reorgconductorarrays(arrays,z,oldiz,oldnz,oldizfs,oldnzfs,
                                   newiz,newnz,newizfs,newnzfs):
    # --- Create list to save the incoming data in.
    results = len(arrays)*[[]]

    # --- Loop over global extent of grid, gathering data that is needed
    for iz in range(0,w3d.nz+1):

        # --- If me has the data then get the indices of it.
        if (oldizfs[me] <= iz <= oldizfs[me]+oldnzfs[me]):
            ii = compress(equal(iz,z),arange(len(z)))

            # --- If the data is needed by me, just copy it.
            if (newizfs[me] <= iz <= newizfs[me]+newnzfs[me]):
                for i in range(len(arrays)):
                    results[i] = results[i] + list(take(arrays[i],ii))

        # --- Get the processor which "owns" the data, relative to the old
        # --- grid extents.
        pe = compress(logical_and(less_equal(oldizfs,iz),
                      less_equal(iz,oldizfs+oldnzfs)),arange(npes))[-1]

        if me == pe:
            # --- Loop over processors to check which ones need data
            # --- If the data is needed by others, then send it.
            for ip in range(npes):
                if (not (oldizfs[ip] <= iz <= oldizfs[ip]+oldnzfs[ip]) and
                        (newizfs[ip] <= iz <= newizfs[ip]+newnzfs[ip])):
                    for i in range(len(arrays)):
                        temp = getarray(me,take(arrays[i],ii),ip)

        else:
            # --- Get the data that is needed from other procesors.
            if (not (oldizfs[me] <= iz <= oldizfs[me]+oldnzfs[me]) and
                    (newizfs[me] <= iz <= newizfs[me]+newnzfs[me])):
                for i in range(len(arrays)):
                    results[i] = results[i] + list(getarray(pe,0,me))

    # --- Make sure all processors are done before continuing
    barrier()

    for i in range(len(results)): results[i] = array(results[i])
    return results



