from warp import *

def printparameters():

  # --- Exit now if output parameters are not to be printed
  if not top.lprntpara: return

  # --- Formats
  f20 = " %s%11.4e%s\n"
  f40 = " %s%11.4e%s%11.4e%s\n"
  f30 = " %s%8d%s\n"

  textblock = \
      f20%("Atomic number of ion = ",top.aion," ") + \
      f20%("Charge state of ion  = ",top.zion," ") + \
      f40%("Initial X, Y emittances = ",top.emitx,",  ",top.emity," m-rad") + \
      f40%("Initial X,Y envelope radii  = ",top.a0,",  ",top.b0," m") + \
      f40%("Initial X,Y envelope angles = ",top.ap0,",  ",top.bp0," rad")+\
      f20%("Input beam current = ",top.ibeam," amps") + \
      f20%("Current density = ",top.currdens," amps/m**2") + \
      f20%("Charge density = ",top.chrgdens," Coul/m**3") + \
      f20%("Number density = ",top.numbdens," ") + \
      f20%("Plasma frequency     = ",top.omegap," 1/s") + \
      f20%("   times dt          = ",top.omegapdt," ") + \
      f20%("   times quad period = ",top.omegaptq," ") + \
      f20%("Plasma period        = ",top.taup," s") + \
      f40%("X-, Y-Thermal Velocities     = ",top.vthx,",  ",top.vthy," m/s")+\
      f40%("   times dt                  = ",top.vthxdt,",  ",top.vthydt," m")+\
      f40%("   times dt/dx, dt/dy (X, Y) = ",top.vthxdtodx,
           ", ", top.vthydtody," ") + \
      f40%("X-, Y-Debye Wavelengths  = ",top.lamdebx,",  ",top.lamdeby," m")+ \
      f40%("   over dx, dy (X and Y) = ",top.lamdebxodx,",  ",
           top.lamdebyody," ") + \
      f20%("Longitudinal thermal velocity (rms) = ",top.vthz," m/s") + \
      f20%("   times dt                   = ",top.vthzdt," m") + \
      f20%("   times dt/dz                = ",top.vthzdtodz," ") + \
      f20%("Longitudinal Debye wavelength = ",top.lamdebz," m") + \
      f20%("   over dz                    = ",top.lamdebzodz," ")
  plt(textblock,0.12,0.88,justify="LT")
  fma()

  # --- Start a new frame since they all don't fit on one
  textblock = \
      f20%("Beam velocity = ",top.vbeam," m/s") + \
      f20%("   over c     = ",top.vbeamoc," ") + \
      f20%("Kinetic energy = ",top.ekin," eV") + \
      f20%("Weight of simulation particles = ",top.sw[0]," ") + \
      f30%("Number of simulation particles = ",top.np," ") + \
      f20%("Number of real particles = ",top.npreal," ") + \
      f20%("Total mass = ",top.totmass," kg") + \
      f20%("Total charge = ",top.totchrg," Coul") + \
      f20%("Generalized perveance = ",top.genperv," ") + \
      f20%("Characteristic current = ",top.charcurr," amps") + \
      f20%("Budker parameter = ",top.budker," ") + \
      f20%("Timestep size dt = ",top.dt," s") + \
      f20%("Tune length = ",top.tunelen," ") + \
      f40%("Undep. X-, Y-Betatron frequencies  = ",top.omegab0x,",  ",\
           top.omegab0y," 1/s") + \
      f40%("Undep. X-, Y-Betatron periods      = ",top.taub0x,
           ", ",top.taub0y," s") + \
      f40%("Undep. X-, Y-Betatron wavelengths  = ",top.lambdab0x,",  ",\
           top.lambdab0y," m") + \
      f40%("Dep.   X-, Y-Betatron frequencies  = ",top.omegabx,",  ",\
           top.omegaby," 1/s") + \
      f40%("Dep.   X-, Y-Betatron periods      = ",top.taubx,",  ",
           top.tauby," s") + \
      f40%("Dep.   X-, Y-Betatron wavelengths  = ",top.lambdabx,",  ",\
           top.lambdaby," m") + \
      f40%("X-, Y-Tune Depressions (sigma/sigma0) = ",
           top.sigmax/(top.sigma0x+top.smallpos),
           ", ",top.sigmay/(top.sigma0y+top.smallpos)," ") + \
      f20%("Space charge wave velocity = ",top.vwave," m/s") + \
      f20%("Effective wall radius = ",top.rwall," m") + \
      f20%("Geometric factor = ",top.gfactor," ") + \
      f40%("X-, Y-Emittance over Space charge forces = ",top.femtxofscx,
           ", ",top.femtyofscy," ")
  plt(textblock,0.12,0.88,justify="LT")
  fma()

