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
      if not withegun: installafterfs(self.solve6)

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
    pass

