"""Class for doing complete multigrid field solve"""
# ToDo:
#  - modify setrho to check if particles are within grid
#  - incorporate instances into the particle mover, so charge is deposited and
#    the E fields gather appropriately.
from warp import *
from generateconductors import installconductors
import MA


##############################################################################
class MultiGrid:
  
  __w3dinputs__ = ['nx','ny','nz',
                   'xmmin','xmmax','ymmin','ymmax','zmmin','zmmax',
                   'bound0','boundnz','boundxy','l2symtry','l4symtry',
                   'solvergeom']
  __f3dinputs__ = ['gridmode','mgparam','downpasses','uppasses',
                   'mgmaxiters','mgtol','mgmaxlevels','mgform',
                   'lcndbndy','icndbndy','laddconductor'] 

  def __init__(self,**kw):

    # --- Kludge - make sure that the multigrid3df routines never sets up
    # --- any conductors.
    f3d.gridmode = 1

    # --- Save input parameters
    for name in MultiGrid.__w3dinputs__:
      if name not in self.__dict__:
        #self.__dict__[name] = kw.pop(name,getattr(w3d,name)) # Python2.3
        self.__dict__[name] = kw.get(name,getattr(w3d,name))
        if kw.has_key(name): del kw[name]
    for name in MultiGrid.__f3dinputs__:
      if name not in self.__dict__:
        #self.__dict__[name] = kw.pop(name,getattr(f3d,name)) # Python2.3
        self.__dict__[name] = kw.get(name,getattr(f3d,name))
        if kw.has_key(name): del kw[name]

    # --- If there are any remaning keyword arguments, raise an error.
    assert len(kw.keys()) == 0,"Bad keyword arguemnts %s"%kw.keys()

    # --- Set set parallel related paramters and calculate mesh sizes
    self.nzfull = self.nz
    self.zmminglobal = self.zmmin
    self.zmmaxglobal = self.zmmax
    self.dx = (self.xmmax - self.xmmin)/self.nx
    self.dy = (self.ymmax - self.ymmin)/self.ny
    self.dz = (self.zmmax - self.zmmin)/self.nz

    self.xmesh = self.xmmin + arange(0,self.nx+1)*self.dx
    self.ymesh = self.ymmin + arange(0,self.ny+1)*self.dy
    self.zmesh = self.zmmin + arange(0,self.nz+1)*self.dz

    # --- Create phi and rho arrays and other arrays. These are created
    # --- with fortran ordering so no transpose and copy is needed when
    # --- they are passed to fortran.
    self.rho = fzeros((1+self.nx,1+self.ny,1+self.nz),'d')
    self.phi = fzeros((1+self.nx,1+self.ny,3+self.nz),'d')
    self.rstar = fzeros(3+self.nz,'d')

    # --- Create a conductor object, which by default is empty.
    self.conductors = ConductorType()

    self.mgiters = 0
    self.mgerror = 0.

    # --- Note that at this time, bends are not supported.
    self.linbend = false

    # --- Turn of build quads option
    self.lbuildquads = false

  def setrho(self,x,y,z,uz,q,w):
    n = len(x)
    setrho3d(self.rho,self.rho,n,x,y,z,top.zgrid,uz,q,w,top.depos,
             self.nx,self.ny,self.nz,self.dx,self.dy,self.dz,
             self.xmmin,self.ymmin,self.zmmin,self.l2symtry,self.l4symtry)

  def gete(self,x,y,z,ex,ey,ez):
    n = len(x)
    sete3d(self.phi,0.,n,x,y,z,top.zgrid,self.xmmin,self.ymmin,self.zmmin,
           self.dx,self.dy,self.dz,self.nx,self.ny,self.nz,1,ex,ey,ez,
           self.l2symtry,self.l4symtry)

  def setphi(self,x,y,z,phi):
    n = len(x)
    getgrid3d(n,x,y,z,phi,self.nx,self.ny,self.nz,self.phi[:,:,1:-1],
              self.xmmin,self.xmmax,self.ymmin,self.ymmax,self.zmmin,self.zmmax,
              self.l2symtry,self.l4symtry)

  def loadrho(self,ins_i=-1,nps_i=-1,is_i=-1,lzero=true):
    if lzero: self.rho[...] = 0.
    for i,n,q,w in zip(top.ins,top.nps,top.sq,top.sw):
      self.setrho(top.xp[i:i+n],top.yp[i:i+n],top.zp[i:i+n],top.uzp[i:i+n],q,w)

  def fetche(self,ipmin,ip,js,ex,ey,ez):
    for i,n in zip(top.ins,top.nps):
      self.gete(top.xp[i:i+n],top.yp[i:i+n],top.zp[i:i+n],ex,ey,ez)

  def installconductor(self,conductor,
                            xmin=None,xmax=None,
                            ymin=None,ymax=None,
                            zmin=None,zmax=None,
                            dfill=top.largepos):
    installconductors(conductor,xmin,xmax,ymin,ymax,zmin,zmax,dfill,
                      top.zbeam,
                      self.nx,self.ny,self.nz,self.nzfull,
                      self.xmmin,self.xmmax,self.ymmin,self.ymmax,
                      self.zmmin,self.zmmax,self.l2symtry,self.l4symtry,
                      solvergeom=self.solvergeom,
                      conductors=self.conductors)

  def solve(self,iwhich=0):
    multigrid3dsolve(iwhich,self.nx,self.ny,self.nz,self.nzfull,
                     self.dx,self.dy,self.dz,self.phi,self.rho,
                     self.rstar,self.linbend,
                     self.bound0,self.boundnz,self.boundxy,
                     self.l2symtry,self.l4symtry,
                     self.xmmin,self.ymmin,self.zmmin,top.zbeam,top.zgrid,
                     self.mgparam,self.mgform,self.mgiters,self.mgmaxiters,
                     self.mgmaxlevels,self.mgerror,self.mgtol,
                     self.downpasses,self.uppasses,
                     self.lcndbndy,self.laddconductor,self.icndbndy,
                     self.lbuildquads,
                     self.gridmode,
                     self.conductors)

  ##########################################################################
  # Define the basic plot commands
  def genericpf(self,kw,pffunc):
    kw['conductors'] = self.conductors
    kw['w3dgrid'] = self
    pffunc(**kw)
  def pfxy(self,**kw): self.genericpf(kw,pfxy)
  def pfzx(self,**kw): self.genericpf(kw,pfzx)
  def pfzy(self,**kw): self.genericpf(kw,pfzy)
  def pfxyg(self,**kw): self.genericpf(kw,pfxyg)
  def pfzxg(self,**kw): self.genericpf(kw,pfzxg)
  def pfzyg(self,**kw): self.genericpf(kw,pfzyg)








  #===========================================================================
  def solve1(self,iwhich=0):
    # --- No initialization needed
    if iwhich == 1: return

    # --- Create temp arrays
    phisave = fzeros(shape(self.phi),'d')
    bendx = fzeros(((self.nx+1)*(self.ny+1)),'d')

    # --- Initialize temporaries
    nxy    = (self.nx+1)*(self.ny+1)
    nxyz   = (self.nx+1)*(self.ny+1)*(self.nz+1)
    dxsqi  = 1./self.dx**2
    dysqi  = 1./self.dy**2
    dzsqi  = 1./self.dz**2
    reps0c = self.mgparam/(eps0*2.*(dxsqi+dysqi+dzsqi))
    rdel   = dzsqi/(dxsqi + dysqi + dzsqi)

    checkconductors(self.nx,self.ny,self.nz,self.nzfull,
                    self.dx,self.dy,self.dz,self.l2symtry,self.l4symtry,
                    self.conductors)

    # --- Preset rho to increase performance (reducing the number of
    # --- multiplies in the main SOR sweep loop).
    if not self.linbend:
      # --- Do the operation in place (to avoid temp arrays)
      multiply(self.rho,reps0c,self.rho)
    else:
      raise "Bends not yet supported"

    # --- If using residual correction form, need to save the original rho.
    # --- Also setup parallel arrays.
    if self.mgform == 2:
      rhosave = self.rho + 0.
      res = fzeros(shape(self.phi),'d')
    #ifdef MPIPARALLEL
    #  lparity = 0
    #  rparity = 0
    #  mggetexchangepes(nslaves,izfsslave,nzfsslave,my_index,
    #                   bound0,boundnz,nzfull,
    #                   lparity,rparity,
    #                   whosendingleft,izsendingleft,
    #                   whosendingright,izsendingright)
    #  localb0 = bound0
    #  localbnz = boundnz
    #  if (izfsslave(my_index) > 0) localb0 = -1
    #  if (izfsslave(my_index)+nz < nzfull) localbnz = -1
    #else
      localb0 = self.bound0
      localbnz = self.boundnz
    #endif

    #   --- Main multigrid v-cycle loop. Calculate error each iteration since
    #   --- very few iterations are done.
    self.mgiters = 0
    self.mgerror = 2.*self.mgtol + 1.
    while (self.mgerror > self.mgtol and self.mgiters < self.mgmaxiters):
      self.mgiters = self.mgiters + 1

      # --- Save current value of phi
      phisave[:,:,:] = self.phi + 0.

      # --- If using residual correction form, calculate the residual and
      # --- copy it into rhosave, zero phisave (the initial error).
      # --- In the calls to cond_potmg and residual, the last argument
      # --- is true, telling the routines to use the actual value of
      # --- voltages rather than zero as is done otherwise for residual
      # --- correction form since it is operating on the error.
      if self.mgform == 2:
        cond_potmg(self.conductors.interior,
                   self.nx,self.ny,self.nz,phisave,0,false,
                   self.mgform,true)
        residual(self.nx,self.ny,self.nz,self.nzfull,dxsqi,dysqi,dzsqi,
                 phisave,rhosave,res,
                 0,localb0,localbnz,self.boundxy,self.l2symtry,self.l4symtry,
                 self.mgparam,self.mgform,true,
                 self.lcndbndy,self.icndbndy,self.conductors)
    #ifdef MPIPARALLEL
    #  mgexchange_phi(nx,ny,nz,nzfull,res,localb0,localbnz,-1,
    #                 my_index,nslaves,izfsslave,nzfsslave,
    #                 whosendingleft,izsendingleft,
    #                 whosendingright,izsendingright)
    #  mgexchange_phiperiodic(nx,ny,nz,nzfull,res,localb0,localbnz,0,
    #                         my_index,nslaves,izfsslave,
    #                         whosendingleft,whosendingright)
    #endif
        self.rho[:,:,:] = res[:,:,1:-1]
        self.phi[:,:,:] = 0.

      # --- Do one vcycle.
      self.vcycle(0,self.nx,self.ny,self.nz,self.nzfull,
                  self.dx,self.dy,self.dz,self.phi,self.rho,
                  self.rstar,self.linbend,self.l2symtry,self.l4symtry,bendx,
                  self.boundxy,self.bound0,self.boundnz,
                  self.mgparam,self.mgform,self.mgmaxlevels,
                  self.downpasses,self.uppasses,self.lcndbndy,
                  self.icndbndy,self.conductors)

      # --- If using residual correction form, add the resulting error to phi.
      if self.mgform == 2: add(self.phi,phisave,self.phi)

      # --- When using residual correction form, the other planes do need
      # --- to be set when using other than Dirichlet boundaries since
      # --- those planes are only set with the error of phi.
      if self.mgform == 2:
        if localb0  == 1: self.phi[:,:,0] = self.phi[:,:,2]
        if localbnz == 1: self.phi[:,:,-1] = self.phi[:,:,-3]
      #ifndef MPIPARALLEL
        if localb0  == 2: self.phi[:,:,0] = self.phi[:,:,-3]
        if localbnz == 2: self.phi[:,:,-1] = self.phi[:,:,2]
      #else
      # mgexchange_phi(nx,ny,nz,nzfull,phi,localb0,localbnz,0,
      #                my_index,nslaves,izfsslave,nzfsslave,
      #                whosendingleft,izsendingleft,
      #                whosendingright,izsendingright)
      # mgexchange_phi(nx,ny,nz,nzfull,phi,localb0,localbnz,-1,
      #                my_index,nslaves,izfsslave,nzfsslave,
      #                whosendingleft,izsendingleft,
      #                whosendingright,izsendingright)
      #endif

      # --- Calculate the change in phi.
      subtract(phisave,self.phi,phisave)
      absolute(phisave,phisave)
      self.mgerror = MA.maximum(phisave)

    #ifdef MPIPARALLEL
    #     --- calculate global sorerror
    #     call parallelmaxrealarray(self.mgerror,1)

    # --- For Dirichlet boundary conditions, copy data into guard planes
    # --- For other boundary conditions, the guard planes are used during
    # --- the solve are so are already set.
    if (self.bound0 == 0): self.phi[:,:,0] = self.phi[:,:,1]
    if (self.boundnz == 0): self.phi[:,:,-1] = self.phi[:,:,-2]

    # --- Make a print out.
    if (self.mgerror > self.mgtol):
      print "Multigrid: Maximum number of iterations reached"
    print ("Multigrid: Error converged to %11.3e in %4d v-cycles"%
           (self.mgerror,self.mgiters))

    # --- If using residual correction form, restore saved rho
    if self.mgform == 2: self.rho[:,:,:] = rhosave

    # --- Restore rho
    if (not self.linbend):
      multiply(self.rho,1./reps0c,self.rho)

  #===========================================================================
  def vcycle(self,mglevel,nx,ny,nz,nzfull,dx,dy,dz,
                  phi,rho,rstar,linbend,l2symtry,l4symtry,
                  bendx,boundxy,globalb0,globalbnz,
                  mgparam,mgform,
                  mgmaxlevels,downpasses,uppasses,
                  lcndbndy,icndbndy,conductors):
   
    res = fzeros(shape(phi),'d')

    dxsqi = 1./dx**2
    dysqi = 1./dy**2
    dzsqi = 1./dz**2

    #ifdef MPIPARALLEL
    #lparityall = 0
    #rparityall = 0
    #mggetexchangepes(nslaves,izfsslave,nzfsslave,my_index,
    #                 globalb0,globalbnz,nzfull,
    #                 lparityall,rparityall,
    #                 whosendingleft,izsendingleft,
    #                 whosendingright,izsendingright)
    #localb0 = globalb0
    #localbnz = globalbnz
    #if (izfsslave(my_index) > 0) localb0 = -1
    #if (izfsslave(my_index)+nz < nzfull) localbnz = -1
    #else
    localb0 = globalb0
    localbnz = globalbnz
    #endif

    for i in range(downpasses):
      self.sorpass3d(mglevel,nx,ny,nz,nzfull,phi,rho,rstar,
                     dxsqi,dysqi,dzsqi,linbend,l2symtry,l4symtry,bendx,
                     localb0,localbnz,boundxy,mgparam,mgform,
                     lcndbndy,icndbndy,conductors)

    # --- Check if this is the finest level. If so, then don't do any further
    # --- coarsening. This is the same check that is done in getmglevels.
    # --- If grid is not at its coarsest level in any of the axis or and
    # --- all dimensions are even, continue the coarsening.
    if ((nx%4) == 0 and (ny%4) == 0 and (nzfull%4) == 0 and
        mglevel < mgmaxlevels):

      # --- Get the residual on the current grid.
      residual(nx,ny,nz,nzfull,dxsqi,dysqi,dzsqi,phi,rho,res,
               mglevel,localb0,localbnz,boundxy,l2symtry,l4symtry,
               mgparam,mgform,false,
               lcndbndy,icndbndy,conductors)
      #ifdef MPIPARALLEL
      #mgexchange_phi(nx,ny,nz,nzfull,res,localb0,localbnz,-1,
      #                  my_index,nslaves,izfsslave,nzfsslave,
      #                  whosendingleft,izsendingleft,
      #                  whosendingright,izsendingright)
      #mgexchange_phiperiodic(nx,ny,nz,nzfull,res,localb0,localbnz,0,
      #                       my_index,nslaves,izfsslave,
      #                       whosendingleft,whosendingright)
      #endif

      # --- If dz > 4/3 dx then only coarsen transversely, otherwise coarsen
      # --- all axis.  This is the same check that is done in getmglevels.
      # --- dz > 4/3 dx <=> (9/16) / dx^2 < 1 / dz^2
      partialcoarsening = (dz > 4./3.*dx)
      #ifdef MPIPARALLEL
      # --- This must be a global operation since, due to roundoff, each
      # --- processor can get a different value if dz == 4./3.*dx.
      #parallellor(partialcoarsening)
      #endif
      if partialcoarsening:

        # --- Allocate new work space
        phi2 = fzeros((1+nx/2,1+ny/2,2+nz+1),'d')
        rho2 = fzeros((1+nx/2,1+ny/2,1+nz),'d')

        # --- Ratio of old to new constant needed to scale the residual for
        # --- the restriction.
        ff = (dxsqi+dysqi+dzsqi)/(dxsqi*0.25 + dysqi*0.25 + dzsqi)
        restrict2d(nx,ny,nz,nzfull,res,rho2,ff,boundxy,localb0,localbnz,
                   l2symtry,l4symtry)

        # --- Continue at the next coarsest level.
        self.vcycle(mglevel+1,nx/2,ny/2,nz,nzfull,
                    dx*2,dy*2,dz,phi2,rho2,
                    rstar,linbend,l2symtry,l4symtry,bendx,
                    boundxy,globalb0,globalbnz,mgparam,mgform,
                    mgmaxlevels,downpasses,uppasses,
                    lcndbndy,icndbndy,conductors)

        # --- Add in resulting error.
        expand2d(nx/2,ny/2,nz,nzfull,phi2,phi,boundxy,localb0,localbnz)

      else:

        #ifdef MPIPARALLEL
        # --- Find domains in coarser grid
        # call mgdividenz(nslaves,izfsslave,nzfsslave,izfsslave2,nzfsslave2,
        #                 nzfull)
        # --- Set new value of nz
        # nznew = nzfsslave2(my_index)
        # --- Difference between starts and ends of coarse and fine grids.
        # --- Should only be in the range 0-2.
        # lparityall = izfsslave - 2*izfsslave2
        # rparityall = 2*(izfsslave2 + nzfsslave2) - (izfsslave + nzfsslave)
        # --- Note that the lparityall and rparityall can only be used in
        # --- MPIPARALLEL sections since they will be unallocated in the
        # --- serial code. So, separate scalars are used in code which is
        # --- used in the serial version.
        # lparity = lparityall(my_index)
        # rparity = rparityall(my_index)
        # --- Get processor with which to exchange data on coarse grid
        # call mggetexchangepes(nslaves,izfsslave2,nzfsslave2,my_index,
        #                       globalb0,globalbnz,nzfull/2,
        #                       lparityall,rparityall,
        #                       whosendingleft2,izsendingleft2,
        #                       whosendingright2,izsendingright2)
        # local2b0 = globalb0
        # local2bnz = globalbnz
        # if (izfsslave2(my_index) > 0) local2b0 = -1
        # if (izfsslave2(my_index) + nznew < nzfull/2) local2bnz = -1
        #else
        nznew = nz/2
        lparity = 0
        rparity = 0
        local2b0 = globalb0
        local2bnz = globalbnz
        #endif

        # --- Alloate new work space
        phi2 = fzeros((1+nx/2,1+ny/2,2+nznew+1),'d')
        rho2 = fzeros((1+nx/2,1+ny/2,1+nznew),'d')

        # --- Restriction - note that scaling factor for residual is always
        # --- 4 for full-coarsening and is compiled into the restriction
        # --- routine.
        restrict3d(nx,ny,nz,nznew,nzfull,res,rho2,boundxy,
                   local2b0,local2bnz,localb0,localbnz,
                   lparity,rparity,l2symtry,l4symtry)
        #ifdef MPIPARALLEL
        # mgexchange_phi(nx/2,ny/2,nznew,nzfull/2,rho2plusmorespace,
        #                local2b0,local2bnz,0,
        #                my_index,nslaves,izfsslave2,nzfsslave2,
        #                whosendingleft2,izsendingleft2,
        #                whosendingright2,izsendingright2)
        #endif

        # --- Continue at the next coarsest level.
        self.vcycle(mglevel+1,nx/2,ny/2,nznew,nzfull/2,
                    dx*2,dy*2,dz*2,phi2,rho2,
                    rstar,linbend,l2symtry,l4symtry,bendx,
                    boundxy,globalb0,globalbnz,mgparam,mgform,
                    mgmaxlevels,downpasses,uppasses,
                    lcndbndy,icndbndy,conductors)

        # --- Add in resulting error.
        expand3d(nx/2,ny/2,nznew,nz,nzfull/2,phi2,phi,
                 boundxy,localb0,localbnz,lparity,rparity)
        #ifdef MPIPARALLEL
        # mgexchange_phiperiodic(nx,ny,nz,nzfull,phi,
        #                        localb0,localbnz,1,
        #                        my_index,nslaves,izfsslave,
        #                        whosendingleft,whosendingright)
        #endif

    # --- Do final SOR passes.
    for i in range(uppasses):
      self.sorpass3d(mglevel,nx,ny,nz,nzfull,phi,rho,rstar,
                     dxsqi,dysqi,dzsqi,linbend,l2symtry,l4symtry,bendx,
                     localb0,localbnz,boundxy,mgparam,mgform,
                     lcndbndy,icndbndy,conductors)

  #===========================================================================
  def sorpass3d(self,mglevel,nx,ny,nz,nzfull,phi,rho,rstar,
                     rdx2,rdy2,rdz2,linbend,l2symtry,l4symtry,bendx,
                     localb0,localbnz,boundxy,mgparam,mgform,
                     lcndbndy,icndbndy,conductors):

    # --- Put desired potential onto conductors in phi array.
    cond_potmg(conductors.interior,nx,ny,nz,phi,mglevel,false,mgform,false)

    # --- Set starting and ending parity.
    #ifdef MPIPARALLEL
    # s_parity = mod(izfsslave(my_index),2)
    # e_parity = mod(s_parity+1,2)
    #else
    s_parity = 0
    e_parity = 1
    #endif

    # --- do loop to cover even and odd points
    for parity in [s_parity,e_parity]:

      sorhalfpass3d(parity,mglevel,nx,ny,nz,nzfull,phi,rho,rstar,
                    rdx2,rdy2,rdz2,linbend,l2symtry,l4symtry,bendx,
                    localb0,localbnz,boundxy,mgparam,mgform,
                    lcndbndy,icndbndy,conductors)

    #ifndef MPIPARALLEL
      if (localb0  == 2): phi[:,:,0] = phi[:,:,-3]
      if (localbnz == 2): phi[:,:,-2:] = phi[:,:,1:3]
    #else
    # mgexchange_phi(nx,ny,nz,nzfull,phi,localb0,localbnz,0,
    #                my_index,nslaves,izfsslave,nzfsslave,
    #                whosendingleft,izsendingleft,
    #                whosendingright,izsendingright)
    # mgexchange_phiperiodic(nx,ny,nz,nzfull,phi,localb0,localbnz,1,
    #                        my_index,nslaves,izfsslave,
    #                        whosendingleft,whosendingright)
    #endif

    #ifdef MPIPARALLEL
    # --- Exchange phi in the z guard planes
    #mgexchange_phi(nx,ny,nz,nzfull,phi,localb0,localbnz,-1,
    #               my_index,nslaves,izfsslave,nzfsslave,
    #               whosendingleft,izsendingleft,
    #               whosendingright,izsendingright)
    #endif

