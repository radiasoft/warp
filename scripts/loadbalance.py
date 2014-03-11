"""LoadBalancer: class wrapping particle load balancing. Sets up automatic
                 periodic load balancing of particles.
"""
__all__ = ['LoadBalancer']
from warp import *
import time


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
 - spreadx,spready,spreadz=1.: The fraction of processors to spread the work
                               over. Do not use this unless you really know
                               what it means!
 - laligntogrid=False: When true, the resulting domains will align to the grid.
                       This must be set to True when using the EM solver.
 - mincellsperdomain=2: When aligning domains to the grid, this specifies the
                        minimum number of grid cells per domain.
 - loadbalancefieldsolver=False: When true, include the field solver in the
                                 load balancing (giving it the same
                                 decomposition as the particles).
                                 This must be set to True when using the
                                 EM solver. The field solver domains will
                                 end up the same as the particle domains.
 - fieldtoparticleeffortratio=0.1: Estimate of the amount of effort per grid
                                   cell per field solve compared to the
                                   amount of effort per particle per step.

Note, if particles only cover a few grid cells, then the distribution is
recalculated on a finer mesh to give better balancing.
    """
    def __init__(self,when=None,padright=None,padleft=None,
                 padupperx=None,paduppery=None,padupperz=None,
                 padlowerx=None,padlowery=None,padlowerz=None,
                 doitnow=0,doloadrho=0,dofs=0,verbose=0,
                 spreadx=1.,spready=1.,spreadz=1.,
                 laligntogrid=False,mincellsperdomain=2,
                 loadbalancefieldsolver=False,fieldtoparticleeffortratio=0.1):
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

        self.spreadx = spreadx
        self.spready = spready
        self.spreadz = spreadz

        # --- If load balancing the field solver, then the particle domains
        # --- must be aligned to the grid. This is primarily a requirement
        # --- of the EM solver.
        if loadbalancefieldsolver: laligntogrid = True

        self.laligntogrid = laligntogrid
        self.mincellsperdomain = mincellsperdomain
        self.loadbalancefieldsolver = loadbalancefieldsolver
        self.fieldtoparticleeffortratio = fieldtoparticleeffortratio

        self.runtime = 0.

        self.reorg = None
        self.doparticleboundaries = None

        if not lparallel: return

        if doitnow: self.doloadbalance(lforce=True)

        # --- Just before scraping is probably the best place to do load
        # --- balancing. The scraping will do the shuffling of the particles
        # --- among the processors, so a special call to do that is not
        # --- needed. Right after the scraping, the charge is deposited, so
        # --- reshuffling of the charge density (and potential) is not needed.
        # --- Also, if loadbalancing were done every step, this would
        # --- gaurantee that particles would never be accidently lost
        # --- since the particles would not move in between load balances.
        installbeforescraper(self.doloadbalance)

    def __setstate__(self,dict):
        self.__dict__.update(dict)
        if not isinstalledbeforescraper(self.doloadbalance):
            installbeforescraper(self.doloadbalance)

    def doloadbalance(self,lforce=False,doloadrho=None,dofs=None,
                      doparticleboundaries=0,reorg=None):
        starttime = time.time()
        if self.when is None: return

        # --- Set lloadbalanced flag to false. It will be set to true below
        # --- if the load balancing will be done.
        top.lloadbalanced = false

        if not lparallel:
            if self.verbose:
                print "Skipping loadbalance since running in serial"
            endtime = time.time()
            self.runtime += (endtime - starttime)
            return

        # --- Find frequency of load balancing
        ii = max(self.when.values())
        for key,value in self.when.iteritems():
            if top.it < key: ii = min(ii,value)

        ppdecomp = top.ppdecomp

        # --- Just return if load balancing not done now and the particle decomposition
        # --- fully covers the extent of the field grid. No checks are needed to ensure
        # --- that particles will go outside of the particle domains since those particles
        # --- are scraped anyway.
        if (ppdecomp.ix[0] == 0 and ppdecomp.ix[-1]+ppdecomp.nx[-1] == ppdecomp.nxglobal and
            ppdecomp.iy[0] == 0 and ppdecomp.iy[-1]+ppdecomp.ny[-1] == ppdecomp.nyglobal and
            ppdecomp.iz[0] == 0 and ppdecomp.iz[-1]+ppdecomp.nz[-1] == ppdecomp.nzglobal):
            if not lforce and (top.it%ii) != 0:
                if self.verbose:
                    print "Skipping loadbalance since it is not time for it"
                endtime = time.time()
                self.runtime += (endtime - starttime)
                return

        # --- Older version sometimes used the precalculated moments to
        # --- avoid recalculating the quantities.
        #usemoments = not (not top.lmoments or top.ifzmmnt == 0
        #                  or top.laccumulate_zmoments)
        # --- The new version, with this called before scraping, always
        # --- recalculates the quantities needed since the zmoments will
        # --- not be up to date.
        usemoments = False

        # --- Get the current number of live particles.
        if not usemoments:
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
        if not usemoments:
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
                    xminp = min(xminp,minnd(xx))
                    xmaxp = max(xmaxp,maxnd(xx))
                else:
                    xminp = w3d.xmmin
                    xmaxp = w3d.xmmax
                if top.nyprocs > 1:
                    yy = top.pgroup.yp[i1:i2]
                    yminp = min(yminp,minnd(yy))
                    ymaxp = max(ymaxp,maxnd(yy))
                else:
                    yminp = w3d.ymmin
                    ymaxp = w3d.ymmax
                if top.nzprocs > 1:
                    zz = top.pgroup.zp[i1:i2]
                    zminp = min(zminp,minnd(zz))
                    zmaxp = max(zmaxp,maxnd(zz))
                else:
                    zminp = w3d.zmmin + top.zbeam
                    zmaxp = w3d.zmmax + top.zbeam
            xminp,yminp,zminp = parallelmin([xminp,yminp,zminp])
            xmaxp,ymaxp,zmaxp = parallelmax([xmaxp,ymaxp,zmaxp])
            # --- Make sure that the mins and maxes are within the bounds
            # --- of the grid. This is needed since there may be some
            # --- particles that are out of bounds (since this happens just
            # --- before the particle scraping).
            xminp = max(xminp,w3d.xmmin)
            xmaxp = min(xmaxp,w3d.xmmax)
            yminp = max(yminp,w3d.ymmin)
            ymaxp = min(ymaxp,w3d.ymmax)
            zminp = max(zminp,w3d.zmmin + top.zbeam)
            zmaxp = min(zmaxp,w3d.zmmax + top.zbeam)
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

        # --- Check if uppermost particle is close to edge of last processor
        # --- If so, then force a reloadbalance.
        if ppdecomp.nxprocs > 1 and ppdecomp.xmax[-1] < w3d.xmmax-0.5*w3d.dx:
            if xmaxp > ppdecomp.xmax[-1]-2*w3d.dx:
                lforce = True
                if self.verbose:
                    print "Load balancing since particles near upper end ",
                    print "of mesh in x ",ppdecomp.xmax[-1],w3d.xmmax,xmaxp,
                    print ppdecomp.xmax[-1]-2*w3d.dx

        # --- Check if lowermost particle is close to edge of last processor
        # --- If so, then force a reloadbalance.
        if ppdecomp.nxprocs > 1 and ppdecomp.xmin[0] > w3d.xmmin+0.5*w3d.dx:
            if xminp < ppdecomp.xmin[0]+2*w3d.dx:
                lforce = True
                if self.verbose:
                    print "Load balancing since particles near lower end ",
                    print "of mesh in x ",ppdecomp.xmin[0],w3d.xmmin,xminp,
                    print ppdecomp.xmin[0]+2*w3d.dx

        # --- Check if uppermost particle is close to edge of last processor
        # --- If so, then force a reloadbalance.
        if ppdecomp.nyprocs > 1 and ppdecomp.ymax[-1] < w3d.ymmax-0.5*w3d.dy:
            if ymaxp > ppdecomp.ymax[-1]-2*w3d.dy:
                lforce = True
                if self.verbose:
                    print "Load balancing since particles near upper end ",
                    print "of mesh in y ",ppdecomp.ymax[-1],w3d.ymmax,ymaxp,
                    print ppdecomp.ymax[-1]-2*w3d.dy

        # --- Check if lowermost particle is close to edge of last processor
        # --- If so, then force a reloadbalance.
        if ppdecomp.nyprocs > 1 and ppdecomp.ymin[0] > w3d.ymmin+0.5*w3d.dy:
            if yminp < ppdecomp.ymin[0]+2*w3d.dy:
                lforce = True
                if self.verbose:
                    print "Load balancing since particles near lower end ",
                    print "of mesh in y ",ppdecomp.ymin[0],w3d.ymmin,yminp,
                    print ppdecomp.ymin[0]+2*w3d.dy

        # --- Check if uppermost particle is close to edge of last processor
        # --- If so, then force a reloadbalance.
        if ppdecomp.nzprocs > 1 and ppdecomp.zmax[-1] < w3d.zmmax-0.5*w3d.dz:
            if zmaxp > ppdecomp.zmax[-1]-2*w3d.dz + top.zbeam:
                lforce = True
                if self.verbose:
                    print "Load balancing since particles near upper end ",
                    print "of mesh in z ",ppdecomp.zmax[-1],w3d.zmmax,zmaxp,
                    print ppdecomp.zmax[-1]-2*w3d.dz

        # --- Check if lowermost particle is close to edge of last processor
        # --- If so, then force a reloadbalance.
        if ppdecomp.nzprocs > 1 and ppdecomp.zmin[0] > w3d.zmmin+0.5*w3d.dz:
            if zminp < ppdecomp.zmin[0]+2*w3d.dz + top.zbeam:
                lforce = True
                if self.verbose:
                    print "Load balancing since particles near lower end ",
                    print "of mesh in z ",ppdecomp.zmin[0],w3d.zmmin,zminp,
                    print ppdecomp.zmin[0]+2*w3d.dz

        # --- Shift into the grid frame
        xminp = xminp - w3d.xmmin
        xmaxp = xmaxp - w3d.xmmin
        yminp = yminp - w3d.ymmin
        ymaxp = ymaxp - w3d.ymmin
        zminp = zminp - w3d.zmmin - top.zbeam
        zmaxp = zmaxp - w3d.zmmin - top.zbeam

        # --- Add the padding to the lower and upper end of
        # --- the grid. This is done here now so that if a
        # --- spread is given, it will include the domain within
        # --- the padding.
        padlowerx = self.calcpadlower(0,ii,self.padlowerx,
                                     top.pgroup.getpyobject('uxp'),
                                     w3d.dx,usemoments)
        padupperx = self.calcpadupper(0,ii,self.padupperx,
                                     top.pgroup.getpyobject('uxp'),
                                     w3d.dx,usemoments)
        padlowery = self.calcpadlower(0,ii,self.padlowery,
                                     top.pgroup.getpyobject('uyp'),
                                     w3d.dy,usemoments)
        paduppery = self.calcpadupper(0,ii,self.paduppery,
                                     top.pgroup.getpyobject('uyp'),
                                     w3d.dy,usemoments)
        padlowerz = self.calcpadlower(0,ii,self.padlowerz,
                                     top.pgroup.getpyobject('uzp'),
                                     w3d.dz,usemoments)
        padupperz = self.calcpadupper(0,ii,self.padupperz,
                                     top.pgroup.getpyobject('uzp'),
                                     w3d.dz,usemoments)

        xminp = max(0.,xminp - padlowerx)
        xmaxp = min(w3d.xmmax - w3d.xmmin,xmaxp + padupperx)
        yminp = max(0.,yminp - padlowery)
        ymaxp = min(w3d.ymmax - w3d.ymmin,ymaxp + paduppery)
        zminp = max(0.,zminp - padlowerz)
        zmaxp = min(w3d.zmmax - w3d.zmmin,zmaxp + padupperz)

        # --- Just return if load balancing not done now.
        if not lforce and (top.it%ii) != 0:
            if self.verbose:
                print "Skipping loadbalance since it is not time for it"
            endtime = time.time()
            self.runtime += (endtime - starttime)
            return

        if (top.it%ii) == 0 and self.verbose:
            print "Load balancing based on frequency"

        # --- If including the field solver, setup the spread to include it.
        if self.loadbalancefieldsolver:
            Ng = (1. + w3d.nx)*(1. + w3d.ny)*(1. + w3d.nz)
            Np = globalsum(top.pgroup.nps)
            spread = 1./(Ng*self.fieldtoparticleeffortratio/Np + 1.)
            self.spreadx = spread
            self.spready = spread
            self.spreadz = spread

        if top.nxprocs > 1:
            self.dodecomposition(0,ii,xminp,xmaxp,self.spreadx,
                                 w3d.xmmin,w3d.xmmax,w3d.dx,0.,top.nxprocs,
                                 top.pgroup.getpyobject('xp'),
                                 top.pgroup.getpyobject('uxp'),
                                 ppdecomp.nxglobal,w3d.nxpextra,
                                 ppdecomp.xmin,ppdecomp.xmax,
                                 ppdecomp.ix,ppdecomp.nx,usemoments,self.laligntogrid)
            top.xpminlocal = ppdecomp.xmin[top.ixproc]
            top.xpmaxlocal = ppdecomp.xmax[top.ixproc]
            top.xpmin = ppdecomp.xmin[0]
            top.xpmax = ppdecomp.xmax[-1]
            w3d.xmminp = w3d.xmmin + ppdecomp.ix[top.ixproc]*w3d.dx
            w3d.xmmaxp = w3d.xmmin + (ppdecomp.ix[top.ixproc] +
                                      ppdecomp.nx[top.ixproc])*w3d.dx

        if top.nyprocs > 1:
            self.dodecomposition(1,ii,yminp,ymaxp,self.spready,
                                 w3d.ymmin,w3d.ymmax,w3d.dy,0.,top.nyprocs,
                                 top.pgroup.getpyobject('yp'),
                                 top.pgroup.getpyobject('uyp'),
                                 ppdecomp.nyglobal,w3d.nypextra,
                                 ppdecomp.ymin,ppdecomp.ymax,
                                 ppdecomp.iy,ppdecomp.ny,usemoments,self.laligntogrid)
            top.ypminlocal = ppdecomp.ymin[top.iyproc]
            top.ypmaxlocal = ppdecomp.ymax[top.iyproc]
            top.ypmin = ppdecomp.ymin[0]
            top.ypmax = ppdecomp.ymax[-1]
            w3d.ymminp = w3d.ymmin + ppdecomp.iy[top.iyproc]*w3d.dy
            w3d.ymmaxp = w3d.ymmin + (ppdecomp.iy[top.iyproc] +
                                      ppdecomp.ny[top.iyproc])*w3d.dy

        if top.nzprocs > 1:
            self.dodecomposition(2,ii,zminp,zmaxp,self.spreadz,
                                 w3d.zmmin,w3d.zmmax,w3d.dz,top.zbeam,
                                 top.nzprocs,
                                 top.pgroup.getpyobject('zp'),
                                 top.pgroup.getpyobject('uzp'),
                                 ppdecomp.nzglobal,w3d.nzpextra,
                                 ppdecomp.zmin,ppdecomp.zmax,
                                 ppdecomp.iz,ppdecomp.nz,usemoments,self.laligntogrid)
            top.zpminlocal = ppdecomp.zmin[top.izproc]
            top.zpmaxlocal = ppdecomp.zmax[top.izproc]
            top.zpmin = ppdecomp.zmin[0]
            top.zpmax = ppdecomp.zmax[-1]
            w3d.zmminp = w3d.zmmin + ppdecomp.iz[top.izproc]*w3d.dz
            w3d.zmmaxp = w3d.zmmin + (ppdecomp.iz[top.izproc] +
                                      ppdecomp.nz[top.izproc])*w3d.dz

            top.izpslave[:] = ppdecomp.iz
            top.nzpslave[:] = ppdecomp.nz
            top.zpslmin[:] =  ppdecomp.zmin
            top.zpslmax[:] =  ppdecomp.zmax

        # --- Reorganize the particles only if requested to do so.
        # --- This is otherwise taken care of by the normal particle boundary
        # --- condition calls during the time step.
        # --- On step zero, a complete reorganization is done so the reorg flag
        # --- is set to true to use the particle sorter which is more efficient
        # --- in that case.
        if self.reorg is not None: reorg = self.reorg
        if self.doparticleboundaries is not None: doparticleboundaries = self.doparticleboundaries
        if reorg or (top.it==1 and doparticleboundaries):
            reorgparticles(top.pgroup,w3d.l4symtry,w3d.l2symtry,
                           w3d.solvergeom==w3d.RZgeom)
        elif doparticleboundaries:
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
            ## --- Redistribute phi to the particle arrays if a field solve is
            ## --- not done.
            #if not dofs:
            #    for i in range(getnsndtsforsubcycling()):
            #        getphipforparticles(i)

        if self.loadbalancefieldsolver:
            self.doloadbalancefieldsolver()

        # --- Do some additional work if requested
        if doloadrho is None: doloadrho = self.doloadrho
        if doloadrho: loadrho()
        if dofs: fieldsol(0)

        top.lloadbalanced = true

        endtime = time.time()
        self.runtime += (endtime - starttime)

    def dodecomposition(self,axis,ii,minp,maxp,spread,
                        mmin,mmax,dd,beam,nprocs,pp,uu,
                        nnglobal,npextra,ppdecompmin,ppdecompmax,
                        ppdecompii,ppdecompnn,usemoments,laligntogrid):
        if (axis < 2 or (maxp - minp)/dd < 10 or not usemoments):
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

        #assert max(pnum) > 0.,"No particles found during decomposition"
        if pnum.max() == 0.: return

        # --- Add fictitious data so that actual work is spread only to the
        # --- requested fraction of the processors.
        assert (0. < spread <= 1.),"spread must be between 0 and 1 or 1."
        avepnum = ave(pnum)
        pnum = pnum + avepnum*(1./spread - 1.)

        self.dodecompositionusingpnum(pnum,axis,ii,minp,maxp,
                                      mmin,mmax,dd,pdd,pmin,nprocs,uu,
                                      nnglobal,npextra,ppdecompmin,ppdecompmax,
                                      ppdecompii,ppdecompnn,usemoments,laligntogrid)

    def dodecompositionusingpnum(self,pnum,axis,ii,minp,maxp,
                                 mmin,mmax,dd,pdd,pmin,nprocs,uu,
                                 nnglobal,npextra,ppdecompmin,ppdecompmax,
                                 ppdecompii,ppdecompnn,usemoments,laligntogrid):
        # --- Convert the number of particles to a decomposition
        domain = self.decompose(pnum,nprocs,lfullcoverage=0)
        domain = domain*pdd + pmin
        domain[0] = min(domain[0],minp)
        domain[-1] = max(domain[-1],maxp)

        if laligntogrid:
            # --- Align the domains to the grid.
            # --- It makes sure that all domains are at least mincellsperdomain long.
            idomain = domain/dd
            idomain[0] = nint(idomain[0])
            for i in range(1,npes+1):
                idomain[i] = nint(idomain[i])
                if idomain[i] - idomain[i-1] < self.mincellsperdomain:
                    idomain[i] = idomain[i-1] + self.mincellsperdomain
            if idomain[-1] > nnglobal:
                idomain[-1] = nnglobal
                for i in range(npes-1,-1,-1):
                    if idomain[i+1] - idomain[i] < self.mincellsperdomain:
                        idomain[i] = idomain[i+1] - self.mincellsperdomain

            domain = dd*idomain

        # --- Set domain of each processor.
        ppdecompmin[:] = mmin + domain[:-1]
        ppdecompmax[:] = mmin + domain[1:]

        domaindecomposeparticles(nnglobal,nprocs,npextra,mmin,dd,
                                 zeros(nprocs,'d'),true,
                                 ppdecompii,ppdecompnn,ppdecompmin,ppdecompmax)

    def calcpadupper(self,axis,ii,padupper,uu,dd,usemoments):
        # --- Calculate the padding on the upper edge.
        if padupper is None:
            if axis < 2 or not usemoments:
                vmaxp = -largepos
                for js in range(top.pgroup.ns):
                    if top.pgroup.nps[js] == 0: continue
                    i1 = top.pgroup.ins[js] - 1
                    i2 = i1 + top.pgroup.nps[js]
                    vv = uu[i1:i2]*top.pgroup.gaminv[i1:i2]
                    vmaxp = max(vmaxp,maxnd(vv))
                vmaxp = globalmax(vmaxp)
            else:
                vmaxp = max(top.vzmaxp)
            if vmaxp > 0.: padupper = vmaxp*top.dt*ii*2
            else:          padupper = ii*dd
        if self.verbose:
            print "Load balancing padupper%s = "%(['x','y','z'][axis]),padupper
        return padupper

    def calcpadlower(self,axis,ii,padlower,uu,dd,usemoments):
        # --- Calculate the padding on the lower edge.
        if padlower is None:
            if axis < 2 or not usemoments:
                vminp = +largepos
                for js in range(top.pgroup.ns):
                    if top.pgroup.nps[js] == 0: continue
                    i1 = top.pgroup.ins[js] - 1
                    i2 = i1 + top.pgroup.nps[js]
                    vv = uu[i1:i2]*top.pgroup.gaminv[i1:i2]
                    vminp = min(vminp,minnd(vv))
                vminp = globalmin(vminp)
            else:
                vminp = min(top.vzminp)
            if vminp < 0.: padlower = -vminp*top.dt*ii*2
            else:          padlower = ii*dd
        if self.verbose:
            print "Load balancing padlower%s = "%(['x','y','z'][axis]),padlower
        return padlower

    def decompose(self,weight,npes,lfullcoverage=0):
        """
Converts a weight into the size of the domains.
 - weight: array of relative weights of the work done by each processor
 - npes: number of processors
 - lfullcoverage=0: when true, the domains cover the full extent of
                    the system
Returns an array of the same length which is the relative length of each
of the domains.
        """
        assert maxnd(weight) > 0.,"weight must have some positive elements"
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
        for ip in range(1,npes):
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
                delta = 2.*c/(sqrt(max(0.,b**2 - 4.*a*c)) + b)
            else:
                delta = sqrt(-c/a)
            fract = fract + delta - delta1
            domain[ip] = domain[ip-1] + fract

        # --- Set the end of the last domain
        if not lfullcoverage:
            # --- Find the last place with non-zero weight, and give the last processor
            # --- everything up to that point.
            for ii in range(ii,nn):
              if weight[ii] > 0.: domain[-1] = ii+1
        else:
            domain[-1] = nn

        return domain

    def doloadbalancefieldsolver(self):

        self.domaindecomposefieldstomatchparticles()
        solver = getregisteredsolver()

        if isinstance(solver,EM3D):
            # --- Save a reference to the old data so that it can be copied to
            # --- the new arrays.
            savedblock = solver.block
            savedfsdecomp = solver.fsdecomp

        solver.setupdecompparallel(userfsdecompnx=top.ppdecomp.nx,
                                   userfsdecompny=top.ppdecomp.ny,
                                   userfsdecompnz=top.ppdecomp.nz)

        if isinstance(solver,EM3D):
            solver.allocatefieldarrays()

            # --- This needs to be called since it creates temporaries that
            # --- are dependent on the decomposition
            solver.setuplaser()

            # --- The data needs to be redistributed
            self.redistributeemsolverdata(solver,savedblock,savedfsdecomp)

            # --- Regenerate the data for the particles. This needs to be done
            # --- locally, from the redistributed field grid data.
            solver.setebp()
            if top.efetch[0] != 4:solver.yee2node3d()
            if solver.l_smooth_particle_fields and any(solver.npass_smooth>0):
                solver.smoothfields()

        else:
            # --- For static solvers, the arrays only need to be reallocated.
            solver.allocatedataarrays()

    def domaindecomposefieldstomatchparticles(self):
        top.fsdecomp.nx = top.ppdecomp.nx.copy()
        domaindecomposefields(w3d.nx,top.nxprocs,false,
                              top.fsdecomp.ix,top.fsdecomp.nx,top.grid_overlap)
        top.fsdecomp.ny = top.ppdecomp.ny.copy()
        domaindecomposefields(w3d.ny,top.nyprocs,false,
                              top.fsdecomp.iy,top.fsdecomp.ny,top.grid_overlap)
        top.fsdecomp.nz = top.ppdecomp.nz.copy()
        domaindecomposefields(w3d.nz,top.nzprocs,false,
                              top.fsdecomp.iz,top.fsdecomp.nz,top.grid_overlap)
        top.fsdecomp.xmin = w3d.xmmin + top.fsdecomp.ix*w3d.dx
        top.fsdecomp.ymin = w3d.ymmin + top.fsdecomp.iy*w3d.dy
        top.fsdecomp.zmin = w3d.zmmin + top.fsdecomp.iz*w3d.dz
        top.fsdecomp.xmax = w3d.xmmin + (top.fsdecomp.ix + top.fsdecomp.nx)*w3d.dx
        top.fsdecomp.ymax = w3d.ymmin + (top.fsdecomp.iy + top.fsdecomp.ny)*w3d.dy
        top.fsdecomp.zmax = w3d.zmmin + (top.fsdecomp.iz + top.fsdecomp.nz)*w3d.dz
        w3d.nxlocal = top.fsdecomp.nx[top.ixproc]
        w3d.nylocal = top.fsdecomp.ny[top.iyproc]
        w3d.nzlocal = top.fsdecomp.nz[top.izproc]
        w3d.xmminlocal = top.fsdecomp.xmin[top.ixproc]
        w3d.ymminlocal = top.fsdecomp.ymin[top.iyproc]
        w3d.zmminlocal = top.fsdecomp.zmin[top.izproc]
        w3d.xmmaxlocal = top.fsdecomp.xmax[top.ixproc]
        w3d.ymmaxlocal = top.fsdecomp.ymax[top.iyproc]
        w3d.zmmaxlocal = top.fsdecomp.zmax[top.izproc]
        top.izfsslave = top.fsdecomp.iz
        top.nzfsslave = top.fsdecomp.nz

    def copydecomposition(self,d1,d2):
        # --- The items commented out are not used.
        d2.my_index = d1.my_index
        #d2.nxglobal = d1.nxglobal
        #d2.nyglobal = d1.nyglobal
        #d2.nzglobal = d1.nzglobal
        d2.ixproc = d1.ixproc
        d2.iyproc = d1.iyproc
        d2.izproc = d1.izproc
        d2.nxprocs = d1.nxprocs
        d2.nyprocs = d1.nyprocs
        d2.nzprocs = d1.nzprocs
        d2.mpi_comm = d1.mpi_comm
        #d2.mpi_comm_x = d1.mpi_comm_x
        #d2.mpi_comm_y = d1.mpi_comm_y
        #d2.mpi_comm_z = d1.mpi_comm_z
        d2.gchange()
        #d2.iprocgrid[:] = d1.iprocgrid
        #d2.nprocgrid[:] = d1.nprocgrid
        d2.ix[:] = d1.ix
        d2.nx[:] = d1.nx
        #d2.xmin[:] = d1.xmin
        #d2.xmax[:] = d1.xmax
        d2.iy[:] = d1.iy
        d2.ny[:] = d1.ny
        #d2.ymin[:] = d1.ymin
        #d2.ymax[:] = d1.ymax
        d2.iz[:] = d1.iz
        d2.nz[:] = d1.nz
        #d2.zmin[:] = d1.zmin
        #d2.zmax[:] = d1.zmax
        #d2.mpistatex[:] = d1.mpistatex
        #d2.mpistatey[:] = d1.mpistatey
        #d2.mpistatez[:] = d1.mpistatez

    def redistributeemsolverdata(self,solver,savedblock,savedfsdecomp):
        # --- Create copies of the decomposition objects.
        # --- These are used since various things are changed depending
        # --- on what data is being transferred.
        olddecomp = Decomposition()
        newdecomp = Decomposition()
        self.copydecomposition(savedfsdecomp,olddecomp)
        self.copydecomposition(solver.fsdecomp,newdecomp)

        def transferarray(old,new,name,nonnodeaxis,lsendguards):
            # --- Handy function to transfer the array name from the
            # --- old to the new field type instance.
            # --- When the data is centered between nodes, the number
            # --- of data values is 1 less along that dimension.
            # --- Note that that decomp.ix,iy,iz are still the same.
            # --- Note that only processors that are exchanging data need to
            # --- to call the transfer routine since the transfers are done
            # --- pair-wise, and not as a global alltoall.
            if 'x' in nonnodeaxis and old.nx > 0:
                olddecomp.nx -= 1
                newdecomp.nx -= 1
            if 'y' in nonnodeaxis and old.ny > 0:
                olddecomp.ny -= 1
                newdecomp.ny -= 1
            if 'z' in nonnodeaxis and old.nz > 0:
                olddecomp.nz -= 1
                newdecomp.nz -= 1
            transferarray1toarray23d(old.nx,old.ny,old.nz,getattr(old,name),
                                     new.nx,new.ny,new.nz,getattr(new,name),
                                     old.nxguard,old.nyguard,old.nzguard,lsendguards,
                                     olddecomp,newdecomp)
            if 'x' in nonnodeaxis and old.nx > 0:
                olddecomp.nx += 1
                newdecomp.nx += 1
            if 'y' in nonnodeaxis and old.ny > 0:
                olddecomp.ny += 1
                newdecomp.ny += 1
            if 'z' in nonnodeaxis and old.nz > 0:
                olddecomp.nz += 1
                newdecomp.nz += 1

        # --- Core
        transferarray(savedblock.core.yf,solver.block.core.yf,'Ex',['x'],false)
        transferarray(savedblock.core.yf,solver.block.core.yf,'Ey',['y'],false)
        transferarray(savedblock.core.yf,solver.block.core.yf,'Ez',['z'],false)
        transferarray(savedblock.core.yf,solver.block.core.yf,'Bx',['y','z'],false)
        transferarray(savedblock.core.yf,solver.block.core.yf,'By',['x','z'],false)
        transferarray(savedblock.core.yf,solver.block.core.yf,'Bz',['x','y'],false)

        # --- Note that the p arrays are regenerated afterwards
        # --- This code won't work anyway.
        #transferarray(savedblock.core.yf,solver.block.core.yf,'Exp',['x'],false)
        #transferarray(savedblock.core.yf,solver.block.core.yf,'Eyp',['y'],false)
        #transferarray(savedblock.core.yf,solver.block.core.yf,'Ezp',['z'],false)
        #transferarray(savedblock.core.yf,solver.block.core.yf,'Bxp',['y','z'],false)
        #transferarray(savedblock.core.yf,solver.block.core.yf,'Byp',['x','z'],false)
        #transferarray(savedblock.core.yf,solver.block.core.yf,'Bzp',['x','y'],false)

        # --- If doing damping, then the bar and old need to be exchanged also
        if solver.theta_damp != 0.:
            transferarray(savedblock.core.yf,solver.block.core.yf,'Exold',['x'],false)
            transferarray(savedblock.core.yf,solver.block.core.yf,'Eyold',['y'],false)
            transferarray(savedblock.core.yf,solver.block.core.yf,'Ezold',['z'],false)
            transferarray(savedblock.core.yf,solver.block.core.yf,'Exbar',['x'],false)
            transferarray(savedblock.core.yf,solver.block.core.yf,'Eybar',['y'],false)
            transferarray(savedblock.core.yf,solver.block.core.yf,'Ezbar',['z'],false)

        # --- Exchange the PML data. For each direction, the decomposition objects
        # --- are modified to only include the domains along the appropriate
        # --- boundary.
        def transfersplityee(old,new):
            transferarray(old.syf,new.syf,'exx',['nx'],false)
            transferarray(old.syf,new.syf,'exy',['nx'],false)
            transferarray(old.syf,new.syf,'exz',['nx'],false)
            transferarray(old.syf,new.syf,'eyx',['ny'],false)
            transferarray(old.syf,new.syf,'eyy',['ny'],false)
            transferarray(old.syf,new.syf,'eyz',['ny'],false)
            transferarray(old.syf,new.syf,'ezx',['nz'],false)
            transferarray(old.syf,new.syf,'ezy',['nz'],false)
            transferarray(old.syf,new.syf,'ezz',['nz'],false)
            transferarray(old.syf,new.syf,'bxy',['ny','nz'],false)
            transferarray(old.syf,new.syf,'bxz',['ny','nz'],false)
            transferarray(old.syf,new.syf,'byx',['nx','nz'],false)
            transferarray(old.syf,new.syf,'byz',['nx','nz'],false)
            transferarray(old.syf,new.syf,'bzx',['nx','ny'],false)
            transferarray(old.syf,new.syf,'bzy',['nx','ny'],false)

        ixproc = newdecomp.ixproc
        iyproc = newdecomp.iyproc
        izproc = newdecomp.izproc
        nxprocs = newdecomp.nxprocs
        nyprocs = newdecomp.nyprocs
        nzprocs = newdecomp.nzprocs

        # --- Sides
        def transferside(d1,d2,d3,fieldtype):
            # --- d1 is the normal to the side, d2 and d3 are in the plane
            eold = getattr(savedblock,fieldtype)
            enew = getattr(solver.block,fieldtype)
            # --- The new my_index is calculated relative to the single plane
            # --- taking part in the exchange.
            i2proc = getattr(newdecomp,'i%sproc'%d2)
            i3proc = getattr(newdecomp,'i%sproc'%d3)
            n2procs = getattr(newdecomp,'n%sprocs'%d2)
            olddecomp.my_index = newdecomp.my_index = i2proc + i3proc*n2procs
            # --- The number of processors normal to the plane is set to 1.
            # --- All processors have the iproc=0 and nprocs=1,
            # --- the i=0 and the n=n in the normal direction.
            setattr(olddecomp,'i%sproc'%d1,0)
            setattr(newdecomp,'i%sproc'%d1,0)
            setattr(olddecomp,'n%sprocs'%d1,1)
            setattr(newdecomp,'n%sprocs'%d1,1)
            setattr(olddecomp,'i%s'%d1,zeros(1))
            setattr(newdecomp,'i%s'%d1,zeros(1))
            setattr(olddecomp,'n%s'%d1,zeros(getattr(eold.syf,'n%s'%d1)))
            setattr(newdecomp,'n%s'%d1,zeros(getattr(enew.syf,'n%s'%d1)))
            # --- Exchange all of the arrays
            transfersplityee(eold,enew)
            # --- Reset the decomp objects to the initial values
            self.copydecomposition(savedfsdecomp,olddecomp)
            self.copydecomposition(solver.fsdecomp,newdecomp)

        if solver.bounds[0] == openbc and ixproc == 0 and (nyprocs > 1 or nzprocs > 1):
            transferside('x','y','z','sidexl')
        if solver.bounds[1] == openbc and ixproc == nxprocs-1 and (nyprocs > 1 or nzprocs > 1):
            transferside('x','y','z','sidexr')
        if solver.bounds[2] == openbc and iyproc == 0 and (nxprocs > 1 or nzprocs > 1):
            transferside('y','x','z','sideyl')
        if solver.bounds[3] == openbc and iyproc == nyprocs-1 and (nxprocs > 1 or nzprocs > 1):
            transferside('y','x','z','sideyr')
        if solver.bounds[4] == openbc and izproc == 0 and (nxprocs > 1 or nyprocs > 1):
            transferside('z','x','y','sidezl')
        if solver.bounds[5] == openbc and izproc == nzprocs-1 and (nxprocs > 1 or nyprocs > 1):
            transferside('z','x','y','sidezr')

        # --- Edges
        def transferedge(d1,d2,d3,fieldtype):
            # --- d1 and d2 are the normals to the edge, and d3 is parallel to it
            eold = getattr(savedblock,fieldtype)
            enew = getattr(solver.block,fieldtype)
            olddecomp.my_index = newdecomp.my_index = getattr(newdecomp,'i%sproc'%d3)
            setattr(olddecomp,'i%sproc'%d1,0)
            setattr(newdecomp,'i%sproc'%d1,0)
            setattr(olddecomp,'i%sproc'%d2,0)
            setattr(newdecomp,'i%sproc'%d2,0)
            setattr(olddecomp,'n%sprocs'%d1,1)
            setattr(newdecomp,'n%sprocs'%d1,1)
            setattr(olddecomp,'n%sprocs'%d2,1)
            setattr(newdecomp,'n%sprocs'%d2,1)
            setattr(olddecomp,'i%s'%d1,zeros(1))
            setattr(newdecomp,'i%s'%d1,zeros(1))
            setattr(olddecomp,'n%s'%d1,zeros(getattr(eold.syf,'n%s'%d1)))
            setattr(newdecomp,'n%s'%d1,zeros(getattr(enew.syf,'n%s'%d1)))
            setattr(olddecomp,'i%s'%d2,zeros(1))
            setattr(newdecomp,'i%s'%d2,zeros(1))
            setattr(olddecomp,'n%s'%d2,zeros(getattr(eold.syf,'n%s'%d2)))
            setattr(newdecomp,'n%s'%d2,zeros(getattr(enew.syf,'n%s'%d2)))
            transfersplityee(eold,enew)
            self.copydecomposition(savedfsdecomp,olddecomp)
            self.copydecomposition(solver.fsdecomp,newdecomp)

        if (solver.bounds[0] == openbc and solver.bounds[2] == openbc and ixproc == 0 and iyproc == 0 and nzprocs > 1):
            transferedge('x','y','z','edgexlyl')
        if (solver.bounds[1] == openbc and solver.bounds[2] == openbc and ixproc == nxprocs-1 and iyproc == 0 and nzprocs > 1):
            transferedge('x','y','z','edgexryl')
        if (solver.bounds[0] == openbc and solver.bounds[3] == openbc and ixproc == 0 and iyproc == nyprocs-1 and nzprocs > 1):
            transferedge('x','y','z','edgexlyr')
        if (solver.bounds[1] == openbc and solver.bounds[3] == openbc and ixproc == nxprocs-1 and iyproc == nyprocs-1 and nzprocs > 1):
            transferedge('x','y','z','edgexryr')

        if (solver.bounds[0] == openbc and solver.bounds[4] == openbc and ixproc == 0 and izproc == 0 and nyprocs > 1):
            transferedge('x','z','y','edgexlzl')
        if (solver.bounds[1] == openbc and solver.bounds[4] == openbc and ixproc == nxprocs-1 and izproc == 0 and nyprocs > 1):
            transferedge('x','z','y','edgexrzl')
        if (solver.bounds[0] == openbc and solver.bounds[5] == openbc and ixproc == 0 and izproc == nzprocs-1 and nyprocs > 1):
            transferedge('x','z','y','edgexlzr')
        if (solver.bounds[1] == openbc and solver.bounds[5] == openbc and ixproc == nxprocs-1 and izproc == nzprocs-1 and nyprocs > 1):
            transferedge('x','z','y','edgexrzr')

        if (solver.bounds[2] == openbc and solver.bounds[4] == openbc and iyproc == 0 and izproc == 0 and nxprocs > 1):
            transferedge('y','z','x','edgeylzl')
        if (solver.bounds[3] == openbc and solver.bounds[4] == openbc and iyproc == nyprocs-1 and izproc == 0 and nxprocs > 1):
            transferedge('y','z','x','edgeyrzl')
        if (solver.bounds[2] == openbc and solver.bounds[5] == openbc and iyproc == 0 and izproc == nzprocs-1 and nxprocs > 1):
            transferedge('y','z','x','edgeylzr')
        if (solver.bounds[3] == openbc and solver.bounds[5] == openbc and iyproc == nyprocs-1 and izproc == nzprocs-1 and nxprocs > 1):
            transferedge('y','z','x','edgeyrzr')

        # --- Corners - nothing needs to be done, thank goodness

