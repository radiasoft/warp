from warp import *

class BoltzmanSolver:
  """
 - rhoion: ion number density in the plasma
 - vion: velocity if the ions
         This is not used yet, but will be to control the injected particles
 - te: electron temperature in volts
 - phi0: plasma potential

Method numbers:
 0: attempt to explicitly include the electrons after each sor iteration
    Doesn't seem to work.
 1: attempt to explicitly include the electrons after each even and odd
    pass of an sor iteration
    Doesn't seem to work
 2: Does one iteration of the Newton's method after each full sor iteration.
    Seems to work OK, though convergence is slow since there is only one
    application per gun iteration. Note that convergence degrades when more
    than one sor cycle is done before applying the Newton's correction.
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
  """
  def __init__(self,rhoion,vion,te,phi0,n=1,method=3,withegun=1):
    self.rhoion = rhoion
    self.vion = vion
    self.te = te
    self.phi0 = phi0
    self.n = n
    self.method = method
    self.withegun = withegun
    if self.method == 2:
      f3d.sormaxit = 1
      f3d.isorerr = 1
      f3d.sorparam = 1.
      if not withegun: installafterfs(self.solve2)
    elif self.method in [4]:
      if not withegun: installafterfs(self.solve4)
      f3d.sormaxit = 0
    elif self.method in [1,3,4]:
      self.parity = ones(shape(w3d.rho))
      self.parity[0::2,0::2,0::2] = 0
      self.parity[1::2,1::2,0::2] = 0
      self.parity[0::2,1::2,1::2] = 0
      self.parity[1::2,0::2,1::2] = 0
      if self.method == 1: installafterfs(self.solve1)
      if self.method == 3 and not withegun: installafterfs(self.solve3)
      if self.method == 4 and not withegun: installafterfs(self.solve4)
      f3d.sormaxit = 0
    elif self.method == 0:
      installafterfs(self.solve0)

  def solve0(self):
    "This doesn't work"
    for i in range(self.n):
      f3d.sormaxit = 1
      fieldsol(-1)
      phi = w3d.phi[:,:,1:-1]
      self.rhoelectron =1*self.rhoion*exp(-abs(self.phi0-phi)/self.te)
      self.rhoelectron =w3d.rho*exp(-maximum(0,self.phi0-phi)/self.te)
      const = 0.5/(1./w3d.dx**2 + 1./w3d.dy**2 + 1./w3d.dz**2)
      phi[:,:,:] = phi - self.rhoelectron/eps0*const
      print "PotI,PotE:",maxnd(phi),maxnd(self.rhoelectron/eps0*const)
      print "DnsI,DnsE:",maxnd(w3d.rho),maxnd(self.rhoelectron)
      w3d.phi[:,:,0] = w3d.phi[:,:,1]
  
  def solve1(self):
    "This doesn't work"
    f3d.sormaxit = 0
    rdx2 = 1./w3d.dx**2
    rdy2 = 1./w3d.dy**2
    rdz2 = 1./w3d.dz**2
    const = 0.5/(1./w3d.dx**2 + 1./w3d.dy**2 + 1./w3d.dz**2)
    for i in range(self.n):
      for iparity in [0,1]:
        self.phiprev = w3d.phi[:,:,1:-1]
        sorhalfpass3d(iparity,0,w3d.nx,w3d.ny,w3d.nz,w3d.nzfull,
                      w3d.phi,w3d.rho,w3d.rstar,
                      rdx2,rdy2,rdz2,false,0.,[1,1,1,1,0,0],
                      1.,1,f3d.lcndbndy,f3d.icndbndy,f3d.conductors)
        ff = maximum(0,self.phi0 - self.phiprev)
        re = w3d.rho*exp(-ff/self.te)
        self.rhoelectron = where(self.parity==iparity,re,self.rhoelectron)
        phi = w3d.phi[:,:,1:-1]
        phi[:,:,:] = where(self.parity==iparity,phi - re/eps0*const,phi)
        print "PotI,PotE:",maxnd(phi),maxnd(self.rhoelectron/eps0*const)
        print "DnsI,DnsE:",maxnd(w3d.rho),maxnd(self.rhoelectron)
        w3d.phi[:,:,0] = w3d.phi[:,:,1]
  
  def solve2(self):
    rdx2 = 1./w3d.dx**2
    rdy2 = 1./w3d.dy**2
    rdz2 = 1./w3d.dz**2
    const = 0.5/(1./w3d.dx**2 + 1./w3d.dy**2 + 1./w3d.dz**2)

    # --- Treat electron term implicitly by the Newton's Method
    # --- (one iteration is enough)
    Gijk = self.rhoion*(const/eps0)
    phi = getphi()
    numer = Gijk*exp(-(self.phi0 - phi)/self.te)
    denom = 1.0 + numer/self.te
    phinew = phi - numer/denom
    #phi[:,:,:] = where(phi-self.phi0 > -10*self.te,phinew,phi)
    phi[:,:,1:-1] = phinew[:,:,1:-1]

  def solve3(self,n=None):
    if n is None: n = self.n
    f3d.sormaxit = 0
    rdx2 = 1./w3d.dx**2
    rdy2 = 1./w3d.dy**2
    rdz2 = 1./w3d.dz**2
    const = 0.5/(1./w3d.dx**2 + 1./w3d.dy**2 + 1./w3d.dz**2)
    for i in range(n):
      for iparity in [0,1]:
        reps0c =  1./(eps0*2.*(rdx2+rdy2+rdz2))
        sorhalfpass3d(iparity,0,w3d.nx,w3d.ny,w3d.nz,w3d.nzfull,
                      w3d.phi,w3d.rho*reps0c,w3d.rstar,
                      rdx2,rdy2,rdz2,false,0.,[1,1,1,1,0,0],
                      1.,1,f3d.lcndbndy,f3d.icndbndy,f3d.conductors)

        # --- Treat electron term implicitly by the Newton's Method
        # --- (one iteration is enough)
        Gijk = self.rhoion*(const/eps0)
        phi = getphi()
        numer = Gijk*exp(-(self.phi0 - phi)/self.te)
        denom = 1.0 + numer/self.te
        phinew = phi - numer/denom
        phi[:,:,1:-1] = where(self.parity==iparity,phinew,phi)[:,:,1:-1]


  def solve4(self,n=None):
    if n is None: n = self.n
    f3d.sormaxit = 0
    rdx2 = 1./w3d.dx**2
    rdy2 = 1./w3d.dy**2
    rdz2 = 1./w3d.dz**2
    const = 0.5/(1./w3d.dx**2 + 1./w3d.dy**2 + 1./w3d.dz**2)
    for i in range(n):
      for iparity in [0,1]:
        reps0c =  1./(eps0*2.*(rdx2+rdy2+rdz2))
        sorhalfpass3d(iparity,0,w3d.nx,w3d.ny,w3d.nz,w3d.nzfull,
                      w3d.phi,w3d.rho*reps0c,w3d.rstar,
                      rdx2,rdy2,rdz2,false,0.,[1,1,1,1,0,0],
                      1.,1,f3d.lcndbndy,f3d.icndbndy,f3d.conductors)

      # --- Treat electron term implicitly by the Newton's Method
      # --- (one iteration is enough)
      Gijk = self.rhoion*(const/eps0)
      phi = getphi()
      numer = Gijk*exp(-(self.phi0 - phi)/self.te)
      denom = 1.0 + numer/self.te
      phinew = phi - numer/denom
      phi[:,:,1:-1] = phinew[:,:,1:-1]


