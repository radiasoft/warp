from lattice import *
generatelattice_version = "$Id: generatelattice.py,v 1.2 2001/09/21 23:22:28 dave Exp $"

# Requires that functions beamduration, beamavesize, and endcondition
# be supplied.

class GenerateLattice:
  def __doc__(s):
    return """
Generates a lattice for final compression.
  - ekin: beam midpulse energy (eV's)
  - totalcharge: beam total charge (Coulombs)
  - aion: mass of beam ions (amu)
  - sigma0: undrepressed phase advance (degrees)
  - occupancy: quadrupole occupancy
  - beamduration: function returning beam duration in seconds
  - beamavesize: function returning average beam transverse size in meters
  - endcondition: function returning true when lattice is complete
The three functions take a single argument, an instance of the
GenerateLattice class.
Function generator returns the lattice in MAD format. The user must call
madtowarp.
"""

  def __init__(s,ekin,totalcharge,aion,sigma0,occupancy,
	       beamduration,beamavesize,endcondition):
    s.ekin = ekin
    s.totalcharge = totalcharge
    s.aion = aion
    s.sigma0 = sigma0*pi/180.
    s.occupancy = occupancy
    s.beamduration = beamduration
    s.beamavesize = beamavesize
    s.endcondition = endcondition
    s.vbeam = s.ekintov(s.ekin,s.aion)

  def ekintov(s,ekin,aion):
    ke = jperev*ekin/(aion*amu*clight**2)
    gammabar = 1. + ke
    return clight*sqrt((2*ke+ke**2)/gammabar**2)

  def pervaence(s):
    return s.totalcharge/(s.beamduration(s)*s.vbeam)/(4.*pi*eps0*s.ekin)

  def generator(s):
    lattice = None
    s.zloc = 0.
    quadsign = 1.
    while not s.endcondition(s):
      hlp = s.beamavesize(s)*sqrt((1.-cos(s.sigma0))/(2.*s.pervaence()))
      quadlen = s.occupancy*hlp
      quadstren = 2.*s.ekin/(hlp**2*s.vbeam)*sqrt(2.*(1.-cos(s.sigma0))/(
                                        s.occupancy**2*(1.-2./3.*s.occupancy)))
      quadaperture = 1.25*s.beamavesize(s) + 0.01
      #print "hlp = ",hlp
      #print "quadlen = ",quadlen
      #print "quadstren = ",quadstren
      #print "quadaperture = ",quadaperture
      #print "pervaence = ",s.pervaence()
      drift = Drft(l=hlp*(1.-s.occupancy)/2.,ap=quadaperture)
      quad = Quad(l=hlp*s.occupancy,db=quadstren*quadsign,ap=quadaperture)
      quadsign = -quadsign
      s.zloc = s.zloc + hlp
      if not lattice:
	lattice = drift + quad + drift
      else:
	lattice = lattice + drift + quad + drift
    return lattice





