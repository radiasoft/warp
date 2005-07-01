from warp import *
from multigrid import MultiGrid

class BoltzmanSolver:
  """
 - rhoion: ion number density in the plasma
 - vion: velocity if the ions
         This is not used yet, but will be to control the injected particles
 - te: electron temperature in volts
 - phi0: plasma potential

Method numbers:
 0: Explicitly include the electrons after each full sor iteration.
    This works much the same as methods 2 and 4.
 1: Explicitly includes the electrons after each even and odd
    pass of an sor iteration. This does seem to work OK - in test cases so
    far, it gives similar results to method 3.
 2: Does one iteration of the Newton's method after each full sor iteration.
    Seems to work OK. Convergence is faster the more iterations are done.
    Note that convergence degrades when more
    than one sor cycle is done before applying the Newton's correction.
    The converged value agrees with method 4, but is different than 3. Both 2
    and 4 show a potential greater than phi0 in the plasma region. This does
    not appear in the results with method 1 and 3.
 3: Applies the Newton's method correction after each odd and even pass of the
    sor iteration. This is the method as described by Jin-Soo. It seems to
    work ok. The number of sor iterations can be controlled - taking many at
    the early stage and then decreasing makes convergence vary fast initially,
    and then settles down.
 4: Applies the Newton's method correction after each full sor iteration. Seems
    to work as well as 3, though converged results need to be carefully
    compared. Note that this is the same as method 2, except that here the
    sor passes are called explicitly. This allows the number of iterations
    to be controlled for faster convergence.
 5: Used for test only - Puts a cap on phi of phi0.
 6: FAS multigrid algorithm

Note that in all cases, sormaxit is set to zero so that the functions here
have direct control over the sor iterations.
  """
  def __init__(self,rhoion,vion,te,phi0,n=1,method=3,withegun=1):
    self.rhoion = rhoion
    self.vion = vion
    self.te = te
    self.phi0 = phi0
    self.n = n
    self.method = method
    self.withegun = withegun
    checkconductors(w3d.nx,w3d.ny,w3d.nz,w3d.nzfull,w3d.dx,w3d.dy,w3d.dz,
                    f3d.conductors,0,0,0,0)
    if self.method == 0:
      if not withegun: installafterfs(self.solve0)
      f3d.sormaxit = 0
      f3d.isorerr = 1
      f3d.sorparam = 1.
    elif self.method == 2:
      f3d.sormaxit = 0
      f3d.isorerr = 1
      f3d.sorparam = 1.
      if not withegun: installafterfs(self.solve2)
    elif self.method == 4:
      if not withegun: installafterfs(self.solve4)
      f3d.sormaxit = 0
    elif self.method in [1,3]:
      self.parity = zeros(shape(w3d.rho[:,:,1:-1]))
      self.parity[0::2,0::2,0::2] = 1
      self.parity[1::2,1::2,0::2] = 1
      self.parity[0::2,1::2,1::2] = 1
      self.parity[1::2,0::2,1::2] = 1
      if self.method == 1 and not withegun: installafterfs(self.solve1)
      if self.method == 3 and not withegun: installafterfs(self.solve3)
      f3d.sormaxit = 0
    elif self.method == 5:
      if not withegun: installafterfs(self.solve5)
    elif self.method == 6:
      top.fstype = 13

  def solve0(self):
    # --- These two lines are needed to turn on the field solver, which is
    # --- called from here.
    top.fstype = 3
    f3d.sormaxit = 1
    for i in range(self.n):
      self.phiprev = getphi()[:,:,1:-1].copy()
      fieldsol(-1)
      const = 0.5/(1./w3d.dx**2 + 1./w3d.dy**2 + 1./w3d.dz**2)
      ff = self.phi0 - self.phiprev
      re = self.rhoion*exp(-ff/self.te)
      phi = getphi()[:,:,1:-1]
      phi[:,:,:] = phi - re/eps0*const
      self.phiprev = getphi()[:,:,1:-1].copy()
      #print "PotI,PotE:",maxnd(phi),maxnd(self.rhoelectron/eps0*const)
      #print "DnsI,DnsE:",maxnd(w3d.rho),maxnd(self.rhoelectron)
      cond_potmg(f3d.conductors.interior,w3d.nx,w3d.ny,w3d.nz,w3d.phi,0,1,1)
    # --- Turn the field solver back off.
    f3d.sormaxit = 0
    top.fstype = -1
  
  def solve1(self):
    f3d.sormaxit = 0
    rdx2 = 1./w3d.dx**2
    rdy2 = 1./w3d.dy**2
    rdz2 = 1./w3d.dz**2
    const = 0.5/(1./w3d.dx**2 + 1./w3d.dy**2 + 1./w3d.dz**2)
    reps0c =  1./(eps0*2.*(rdx2+rdy2+rdz2))
    bounds = [w3d.boundxy,w3d.boundxy,w3d.boundxy,w3d.boundxy,
              w3d.bound0,w3d.boundnz]
    for i in range(self.n):
      for iparity in [0,1]:
        self.phiprev = getphi()[:,:,1:-1].copy()
        sorhalfpass3d(iparity,0,w3d.nx,w3d.ny,w3d.nz,w3d.nzfull,
                      w3d.phi,w3d.rho*reps0c,w3d.rstar,
                      rdx2,rdy2,rdz2,false,0.,bounds,
                      1.,1,f3d.lcndbndy,f3d.icndbndy,f3d.conductors)
        ff = self.phi0 - self.phiprev
        re = self.rhoion*exp(-ff/self.te)
        phi = getphi()[:,:,1:-1]
        phi[:,:,:] = where(self.parity==iparity,phi - re/eps0*const,phi)
        cond_potmg(f3d.conductors.interior,w3d.nx,w3d.ny,w3d.nz,w3d.phi,0,1,1)
        #print "PotI,PotE:",maxnd(phi),maxnd(self.rhoelectron/eps0*const)
        #print "DnsI,DnsE:",maxnd(w3d.rho),maxnd(self.rhoelectron)
  
  def solve2(self):
    top.fstype = 3
    rdx2 = 1./w3d.dx**2
    rdy2 = 1./w3d.dy**2
    rdz2 = 1./w3d.dz**2
    const = 0.5/(1./w3d.dx**2 + 1./w3d.dy**2 + 1./w3d.dz**2)
    f3d.sormaxit = 1

    for i in range(self.n):
      fieldsol(-1)
      # --- Treat electron term implicitly by the Newton's Method
      # --- (one iteration is enough)
      Gijk = self.rhoion*(const/eps0)
      phi = getphi()
      numer = Gijk*exp(-(self.phi0 - phi)/self.te)
      denom = 1.0 + numer/self.te
      phinew = phi - numer/denom
      #phi[:,:,:] = where(phi-self.phi0 > -10*self.te,phinew,phi)
      phi[:,:,1:-1] = phinew[:,:,1:-1]
      cond_potmg(f3d.conductors.interior,w3d.nx,w3d.ny,w3d.nz,w3d.phi,0,1,1)

    f3d.sormaxit = 0
    top.fstype = -1

  def solve3(self,n=None):
    if n is None: n = self.n
    f3d.sormaxit = 0
    rdx2 = 1./w3d.dx**2
    rdy2 = 1./w3d.dy**2
    rdz2 = 1./w3d.dz**2
    const = 0.5/(1./w3d.dx**2 + 1./w3d.dy**2 + 1./w3d.dz**2)
    reps0c =  1./(eps0*2.*(rdx2+rdy2+rdz2))
    bounds = [w3d.boundxy,w3d.boundxy,w3d.boundxy,w3d.boundxy,
              w3d.bound0,w3d.boundnz]
    for i in range(n):
      for iparity in [0,1]:
        sorhalfpass3d(iparity,0,w3d.nx,w3d.ny,w3d.nz,w3d.nzfull,
                      w3d.phi,w3d.rho*reps0c,w3d.rstar,
                      rdx2,rdy2,rdz2,false,0.,bounds,
                      1.,1,f3d.lcndbndy,f3d.icndbndy,f3d.conductors)

        # --- Treat electron term implicitly by the Newton's Method
        # --- (one iteration is enough)
        Gijk = self.rhoion*(const/eps0)
        phi = getphi()[:,:,1:-1]
        numer = Gijk*exp(-(self.phi0 - phi)/self.te)
        denom = 1.0 + numer/self.te
        phinew = phi - numer/denom
        phi[:,:,:] = where(self.parity==iparity,phinew,phi)
        cond_potmg(f3d.conductors.interior,w3d.nx,w3d.ny,w3d.nz,w3d.phi,0,1,1)



  def solve4(self,n=None):
    if n is None: n = self.n
    f3d.sormaxit = 0
    rdx2 = 1./w3d.dx**2
    rdy2 = 1./w3d.dy**2
    rdz2 = 1./w3d.dz**2
    const = 0.5/(1./w3d.dx**2 + 1./w3d.dy**2 + 1./w3d.dz**2)
    reps0c =  1./(eps0*2.*(rdx2+rdy2+rdz2))
    bounds = [w3d.boundxy,w3d.boundxy,w3d.boundxy,w3d.boundxy,
              w3d.bound0,w3d.boundnz]
    for i in range(n):
      for iparity in [0,1]:
        sorhalfpass3d(iparity,0,w3d.nx,w3d.ny,w3d.nz,w3d.nzfull,
                      w3d.phi,w3d.rho*reps0c,w3d.rstar,
                      rdx2,rdy2,rdz2,false,0.,bounds,
                      1.,1,f3d.lcndbndy,f3d.icndbndy,f3d.conductors)

      # --- Treat electron term implicitly by the Newton's Method
      # --- (one iteration is enough)
      Gijk = self.rhoion*(const/eps0)
      phi = getphi()
      numer = Gijk*exp(-(self.phi0 - phi)/self.te)
      denom = 1.0 + numer/self.te
      phinew = phi - numer/denom
      phi[:,:,1:-1] = phinew[:,:,1:-1]
      cond_potmg(f3d.conductors.interior,w3d.nx,w3d.ny,w3d.nz,w3d.phi,0,1,1)


  def solve5(self):
    w3d.phi[:,:,:] = minimum(self.phi0,w3d.phi)

  def solve6(self):
    multigridbe3df(-1,w3d.nx,w3d.ny,w3d.nz,w3d.nzfull,
                   w3d.dx,w3d.dy,w3d.dz,w3d.phi,w3d.rho,
                   w3d.rstar,top.linbend,w3d.bound0,w3d.boundnz,w3d.boundxy,
                   w3d.l2symtry,w3d.l4symtry,
                   w3d.xmmin,w3d.ymmin,w3d.zmmin,top.zbeam,top.zgrid,
                   self.rhoion,self.te,self.phi0)



class MultiGridwithBoltzmannElectrons(MultiGrid):
  """
!!!!!!!!!!!!!!!!!!
This is not complete - do not use!!!!!!!!!!
!!!!!!!!!!!!!!!!!!
  """
  def __init__(self,iondensity,electrontemperature,plasmapotential=0.,**kw):
    self.iondensity = iondensity
    self.electrontemperature = electrontemperature
    self.plasmapotential = plasmapotential
    MultiGrid.__init__(self,**kw)


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
                    self.dx,self.dy,self.dz,self.conductors,
                    top.my_index,top.nslaves,top.izfsslave,top.nzfsslave)

    # --- Preset rho to increase performance (reducing the number of
    # --- multiplies in the main SOR sweep loop).
    if not self.linbend:
      # --- Do the operation in place (to avoid temp arrays)
      trho = transpose(self.rho)
      multiply(trho,reps0c,trho)
    else:
      raise "Bends not yet supported"


    #   --- Main multigrid v-cycle loop. Calculate error each iteration since
    #   --- very few iterations are done.
    self.mgiters[0] = 0
    self.mgerror[0] = 2.*self.mgtol + 1.
    while (self.mgerror[0] > self.mgtol and self.mgiters[0] < self.mgmaxiters):
      self.mgiters[0] = self.mgiters[0] + 1

      # --- Save current value of phi
      tphisave = transpose(phisave)
      tphi = transpose(self.phi)
      tphisave[:,:,:] = tphi
      f = fzeros(shape(self.rho),'d')

      # --- Do one vcycle.
      self.vcycle(0,self.nx,self.ny,self.nz,self.nzfull,
                  self.dx,self.dy,self.dz,self.phi,self.rho,f,
                  self.rstar,self.linbend,bendx,self.bounds,
                  self.mgparam,self.mgform,self.mgmaxlevels,
                  self.downpasses,self.uppasses,self.lcndbndy,
                  self.icndbndy,self.conductors)

      # --- Calculate the change in phi.
      subtract(tphisave,tphi,tphisave)
      absolute(tphisave,tphisave)
      self.mgerror[0] = MA.maximum(tphisave)

    # --- For Dirichlet boundary conditions, copy data into guard planes
    # --- For other boundary conditions, the guard planes are used during
    # --- the solve are so are already set.
    if (self.bounds[4] == 0): self.phi[:,:,0] = self.phi[:,:,1]
    if (self.bounds[5] == 0): self.phi[:,:,-1] = self.phi[:,:,-2]

    # --- Make a print out.
    if (self.mgerror[0] > self.mgtol):
      print "Multigrid: Maximum number of iterations reached"
    print ("Multigrid: Error converged to %11.3e in %4d v-cycles"%
           (self.mgerror[0],self.mgiters[0]))

    # --- Restore rho
    if (not self.linbend):
      multiply(self.rho,1./reps0c,self.rho)

  #===========================================================================
  def vcycle(self,mglevel,nx,ny,nz,nzfull,dx,dy,dz,
                  phi,rho,f,rstar,linbend,bendx,bounds,mgparam,mgform,
                  mgmaxlevels,downpasses,uppasses,lcndbndy,icndbndy,conductors):
   
    res = fzeros((3+nx,3+ny,7+nz),'d')

    dxsqi = 1./dx**2
    dysqi = 1./dy**2
    dzsqi = 1./dz**2
    C = 2.*(dxsqi + dysqi + dzsqi)
    phi0 = self.plasmapotential
    te = self.electrontemperature

    localbounds = bounds.copy()

    for i in range(downpasses):
      self.relax3d(mglevel,nx,ny,nz,nzfull,phi,rho,f,res,rstar,
                   dxsqi,dysqi,dzsqi,linbend,bendx,
                   localbounds,mgparam,mgform,lcndbndy,icndbndy,conductors)

    # --- Check if this is the finest level. If so, then don't do any further
    # --- coarsening. This is the same check that is done in getmglevels.
    # --- If grid is not at its coarsest level in any of the axis or and
    # --- all dimensions are even, continue the coarsening.
    if ((nx%4) == 0 and (ny%4) == 0 and (nzfull%4) == 0 and
        mglevel < mgmaxlevels):

      # --- Get the residual on the current grid.
      residual(nx,ny,nz,nzfull,dxsqi,dysqi,dzsqi,phi,rho,res,
               mglevel,localbounds,mgparam,mgform,false,
               lcndbndy,icndbndy,conductors)
      rhoe = self.iondensity*exp((phi-phi0)/te)
      res[1:-1,1:-1,2:-2] = res[1:-1,1:-1,2:-2]*C - rhoe/eps0
      del rhoe

      # --- If dz > 4/3 dx then only coarsen transversely, otherwise coarsen
      # --- all axis.  This is the same check that is done in getmglevels.
      # --- dz > 4/3 dx <=> (9/16) / dx^2 < 1 / dz^2
      partialcoarsening = (dz > 4./3.*dx)
      if partialcoarsening:

        # --- Allocate new work space
        phi2 = fzeros((1+nx/2,1+ny/2,2+nz+1),'d')
        rho2 = fzeros((1+nx/2,1+ny/2,1+nz),'d')

        # --- Ratio of old to new constant needed to scale the residual for
        # --- the restriction.
        ff = (dxsqi+dysqi+dzsqi)/(dxsqi*0.25 + dysqi*0.25 + dzsqi)
        restrict2d(nx,ny,nz,nzfull,res,rho2,ff,localbounds)

        # --- Continue at the next coarsest level.
        self.vcycle(mglevel+1,nx/2,ny/2,nz,nzfull,
                    dx*2,dy*2,dz,phi2,rho2,rstar,linbend,bendx,bounds,
                    mgparam,mgform,mgmaxlevels,downpasses,uppasses,
                    lcndbndy,icndbndy,conductors)

        # --- Add in resulting error.
        expand2d(nx/2,ny/2,nz,nzfull,phi2,phi,localbounds)

      else:

        localboundsH = bounds.copy()

        nznew = nz/2
        lparity = 0
        rparity = 0

        # --- Alloate new work space
        phiH = fzeros((1+nx/2,1+ny/2,2+nznew+1),'d')
        rhoH = fzeros((1+nx/2,1+ny/2,1+nznew),'d')

        rhocopy = fzeros((3+nx,3+ny,7+nz),'d')
        rhocopy[1:-1,1:-1,3:-3] = rho

        restrict3d(nx,ny,nz,nzfull,rhocopy,
                   nx/2,ny/2,nznew,nznew,rhoH,
                   1.,localbounds,localboundsH,0)
        del rhocopy

        res[1:-1,1:-1,3:-3] = res[1:-1,1:-1,3:-3] - f
        dH = fzeros((1+nx/2,1+ny/2,1+nznew),'d')
        restrict3d(nx,ny,nz,nzfull,res,
                   nx/2,ny/2,nznew,nznew,dH,
                   1.,localbounds,localboundsH,0)

        phicopy = fzeros((3+nx,3+ny,7+nz),'d')
        phicopy[1:-1,1:-1,2:-2] = phi
        restrict3d(nx,ny,nz,nzfull,phicopy,
                   nx/2,ny/2,nznew,nznew,phiH[:,:,1:-1],
                   1.,localbounds,localboundsH,0)
        phiHcopy = phiH.copy()
        del phicopy
        resH = fzeros((3+nx/2,3+ny/2,4+nznew+3),'d')
        residual(nx/2,ny/2,nz/2,nzfull/2,0.25*dxsqi,0.25*dysqi,0.25*dzsqi,
                 phiH,rhoH,resH,
                 mglevel+1,localboundsH,mgparam,mgform,false,
                 lcndbndy,icndbndy,conductors)
        rhoe = self.iondensity*exp((phiH[:,:,1:-1]-phi0)/te)
        resH[1:-1,1:-1,3:-3] = resH[1:-1,1:-1,3:-3]*C - rhoe/eps0
        del rhoe

        fH = resH[1:-1,1:-1,3:-3] - dH

        # --- Continue at the next coarsest level.
        self.vcycle(mglevel+1,nx/2,ny/2,nznew,nzfull/2,
                    dx*2,dy*2,dz*2,phiH,rhoH,fH,rstar,linbend,bendx,bounds,
                    mgparam,mgform,mgmaxlevels,downpasses,uppasses,
                    lcndbndy,icndbndy,conductors)

        eH = phiH - phiHcopy

        # --- Add in resulting error.
        expand3d(nx,ny,nz,nz,phi,
                 nx/2,ny/2,nz/2,nz/2,eH,localbounds,0.)

    # --- Do final SOR passes.
    for i in range(uppasses):
      self.relax3d(mglevel,nx,ny,nz,nzfull,phi,rho,f,res,rstar,
                   dxsqi,dysqi,dzsqi,linbend,bendx,localbounds,
                   mgparam,mgform,lcndbndy,icndbndy,conductors)

  #===========================================================================
  def relax3d(self,mglevel,nx,ny,nz,nzfull,phi,rho,f,res,rstar,
                   dxsqi,dysqi,dzsqi,linbend,bendx,bounds,mgparam,mgform,
                   lcndbndy,icndbndy,conductors):

    C = 2.*(dxsqi + dysqi + dzsqi)
    phi0 = self.plasmapotential
    te = self.electrontemperature

    # --- Put desired potential onto conductors in phi array.
    cond_potmg(conductors.interior,nx,ny,nz,phi,mglevel,false,mgform,false)

    # --- Should be done with even-odd ordering to be more correct
    residual(nx,ny,nz,nzfull,dxsqi,dysqi,dzsqi,phi,rho,res,
             mglevel,bounds,mgparam,mgform,false,
             lcndbndy,icndbndy,conductors)
    rhoe = self.iondensity*exp((phi[:,:,1:-1]-phi0)/te)
    res[1:-1,1:-1,3:-3] = res[1:-1,1:-1,3:-3]*C - rhoe/eps0
    phi[:,:,1:-1] = phi[:,:,1:-1] - (res[1:-1,1:-1,3:-3] - f)/(-C - rhoe/(eps0*te))



