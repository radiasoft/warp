from warp import *

class BoltzmanSolver:
  """
 - rhoion: ion number density in the plasma
 - vion: velocity if the ions
         This is not used yet, but will be to control the injected particles
 - te: electron temperature in volts
 - phi0: plasma potential
  """
  def __init__(self,rhoion,vion,te,phi0,n=1,method=2):
    self.rhoion = rhoion
    self.vion = vion
    self.te = te
    self.phi0 = phi0
    self.n = n
    self.method = method
    if self.method == 2:
      installafterfs(self.solve2)
    elif self.method == 1:
      self.rhoelectron = zeros(shape(w3d.rho),'d')
      self.phiprev = zeros(shape(w3d.rho),'d')
      self.parity = ones(shape(w3d.rho))
      self.parity[0::2,0::2,0::2] = 0
      self.parity[1::2,1::2,0::2] = 0
      self.parity[0::2,1::2,1::2] = 0
      self.parity[1::2,0::2,1::2] = 0
      installafterfs(self.solve1)
    elif self.method == 0:
      installafterfs(self.solve0)

  def solve2(self):
    phi = w3d.phi[:,:,2:-1]  # --- This is the updated phi with no electrons

    rdx2 = 1./w3d.dx**2
    rdy2 = 1./w3d.dy**2
    rdz2 = 1./w3d.dz**2
    const = 0.5/(1./w3d.dx**2 + 1./w3d.dy**2 + 1./w3d.dz**2)

    # --- Treat electron term implicitly by the Newton's Method
    # --- (one iteration is enough)
    Gijk = self.rhoion*(const/eps0)
    numer = Gijk*exp((phi - self.phi0)/self.te)
    denom = 1.0 + numer/self.te
    phinew = phi - numer/denom
    phi[:,:,:] = where(phi-self.phi0 > -10*self.te,phinew,phi)

    print "/home/kim/testwarp/boltzmannelectrons_jsk.py"
    print "PotI,PotE:",maxnd(phi),maxnd(self.rhoelectron/eps0*const)
    print "DnsI,DnsE:",maxnd(w3d.rho),maxnd(self.rhoelectron)


  def solve0(self):
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
    f3d.sormaxit = 0
    rdx2 = 1./w3d.dx**2
    rdy2 = 1./w3d.dy**2
    rdz2 = 1./w3d.dz**2
    const = 0.5/(1./w3d.dx**2 + 1./w3d.dy**2 + 1./w3d.dz**2)
    for i in range(self.n):
      for iparity in [0,1]:
        self.phiprev[:,:,:] = w3d.phi[:,:,1:-1]
        sorhalfpass3d(iparity,0,w3d.nx,w3d.ny,w3d.nz,w3d.nzfull,
                      w3d.phi,w3d.rho,w3d.rstar,
                      rdx2,rdy2,rdz2,false,0.,f3d.bounds,
                      1.,1,f3d.lcndbndy,f3d.icndbndy,f3d.conductors)
        ff = maximum(0,self.phi0 - self.phiprev)
        re = w3d.rho*exp(-ff/self.te)
        self.rhoelectron = where(self.parity==iparity,re,self.rhoelectron)
        phi = w3d.phi[:,:,1:-1]
        phi[:,:,:] = where(self.parity==iparity,phi - re/eps0*const,phi)
        print "PotI,PotE:",maxnd(phi),maxnd(self.rhoelectron/eps0*const)
        print "DnsI,DnsE:",maxnd(w3d.rho),maxnd(self.rhoelectron)
        w3d.phi[:,:,0] = w3d.phi[:,:,1]
  

