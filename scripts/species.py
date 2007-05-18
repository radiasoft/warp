"""
This script contains the definition of a number of classes (Particle, Atom, Molecule,
Species), a dictionary of the periodic table of elements, as well as the instanciation of some
usual particles as object (Electron, Positron, Water, atoms from periodic table).
"""
from warp import *
import RandomArray

class Particle(object):
  def __init__(self,mass=None,charge=None,Symbol=None,name=None):
    self.Symbol=Symbol
    self.name = name
    if mass is not None:
      self.mass=mass
      self.M=self.mass
    if charge is not None:
      self.charge=charge
      self.Q=self.charge
    
class Atom(Particle):
  def __init__(self,Symbol,A,Z,Group,Period,name=None):
    Particle.__init__(self,mass=A*amu,Symbol=Symbol,name=name)
    self.A=A
    self.Z=Z
    self.Group=Group
    self.Period=Period

class Molecule(Particle):
  def __init__(self,Symbol,mass,name=None):
    Particle.__init__(self,mass=mass,Symbol=Symbol,name=name)

periodic_table={}
periodic_table['Hydrogen']={'A': 1.0079400000000001, 'Symbol': 'H', 'Z': 1, 'Group': 1, 'Period': 1}
periodic_table['Helium']={'A': 4.0026020000000004, 'Symbol': 'He', 'Z': 2, 'Group': 18, 'Period': 1}
periodic_table['Lithium']={'A': 6.9409999999999998, 'Symbol': 'Li', 'Z': 3, 'Group': 1, 'Period': 2}
periodic_table['Beryllium']={'A': 9.0121819999999992, 'Symbol': 'Be', 'Z': 4, 'Group': 2, 'Period': 2}
periodic_table['Boron']={'A': 10.811, 'Symbol': 'B', 'Z': 5, 'Group': 13, 'Period': 2}
periodic_table['Carbon']={'A': 12.0107, 'Symbol': 'C', 'Z': 6, 'Group': 14, 'Period': 2}
periodic_table['Nitrogen']={'A': 14.0067, 'Symbol': 'N', 'Z': 7, 'Group': 15, 'Period': 2}
periodic_table['Oxygen']={'A': 15.9994, 'Symbol': 'O', 'Z': 8, 'Group': 16, 'Period': 2}
periodic_table['Fluorine']={'A': 18.998403199999998, 'Symbol': 'F', 'Z': 9, 'Group': 17, 'Period': 2}
periodic_table['Neon']={'A': 20.1797, 'Symbol': 'Ne', 'Z': 10, 'Group': 18, 'Period': 2}
periodic_table['Sodium']={'A': 22.98977, 'Symbol': 'Na', 'Z': 11, 'Group': 1, 'Period': 3}
periodic_table['Magnesium']={'A': 24.305, 'Symbol': 'Mg', 'Z': 12, 'Group': 2, 'Period': 3}
periodic_table['Aluminium']={'A': 26.981538, 'Symbol': 'Al', 'Z': 13, 'Group': 13, 'Period': 3}
periodic_table['Silicon']={'A': 28.0855, 'Symbol': 'Si', 'Z': 14, 'Group': 14, 'Period': 3}
periodic_table['Phosphorus']={'A': 30.973761, 'Symbol': 'P', 'Z': 15, 'Group': 15, 'Period': 3}
periodic_table['Sulfur']={'A': 32.064999999999998, 'Symbol': 'S', 'Z': 16, 'Group': 16, 'Period': 3}
periodic_table['Chlorine']={'A': 35.453000000000003, 'Symbol': 'Cl', 'Z': 17, 'Group': 17, 'Period': 3}
periodic_table['Argon']={'A': 39.948, 'Symbol': 'Ar', 'Z': 18, 'Group': 18, 'Period': 3}
periodic_table['Potassium']={'A': 39.098300000000002, 'Symbol': 'K', 'Z': 19, 'Group': 1, 'Period': 4}
periodic_table['Calcium']={'A': 40.078000000000003, 'Symbol': 'Ca', 'Z': 20, 'Group': 2, 'Period': 4}
periodic_table['Scandium']={'A': 44.955910000000003, 'Symbol': 'Sc', 'Z': 21, 'Group': 3, 'Period': 4}
periodic_table['Titanium']={'A': 47.866999999999997, 'Symbol': 'Ti', 'Z': 22, 'Group': 4, 'Period': 4}
periodic_table['Vanadium']={'A': 50.941499999999998, 'Symbol': 'V', 'Z': 23, 'Group': 5, 'Period': 4}
periodic_table['Chromium']={'A': 51.996099999999998, 'Symbol': 'Cr', 'Z': 24, 'Group': 6, 'Period': 4}
periodic_table['Manganese']={'A': 54.938048999999999, 'Symbol': 'Mn', 'Z': 25, 'Group': 7, 'Period': 4}
periodic_table['Iron']={'A': 55.844999999999999, 'Symbol': 'Fe', 'Z': 26, 'Group': 8, 'Period': 4}
periodic_table['Cobalt']={'A': 58.933199999999999, 'Symbol': 'Co', 'Z': 27, 'Group': 9, 'Period': 4}
periodic_table['Nickel']={'A': 58.693399999999997, 'Symbol': 'Ni', 'Z': 28, 'Group': 10, 'Period': 4}
periodic_table['Copper']={'A': 63.545999999999999, 'Symbol': 'Cu', 'Z': 29, 'Group': 11, 'Period': 4}
periodic_table['Zinc']={'A': 65.409000000000006, 'Symbol': 'Zn', 'Z': 30, 'Group': 12, 'Period': 4}
periodic_table['Gallium']={'A': 69.722999999999999, 'Symbol': 'Ga', 'Z': 31, 'Group': 13, 'Period': 4}
periodic_table['Germanium']={'A': 72.640000000000001, 'Symbol': 'Ge', 'Z': 32, 'Group': 14, 'Period': 4}
periodic_table['Arsenic']={'A': 74.921599999999998, 'Symbol': 'As', 'Z': 33, 'Group': 15, 'Period': 4}
periodic_table['Selenium']={'A': 78.959999999999994, 'Symbol': 'Se', 'Z': 34, 'Group': 16, 'Period': 4}
periodic_table['Bromine']={'A': 79.903999999999996, 'Symbol': 'Br', 'Z': 35, 'Group': 17, 'Period': 4}
periodic_table['Krypton']={'A': 83.798000000000002, 'Symbol': 'Kr', 'Z': 36, 'Group': 18, 'Period': 4}
periodic_table['Rubidium']={'A': 85.467799999999997, 'Symbol': 'Rb', 'Z': 37, 'Group': 1, 'Period': 5}
periodic_table['Strontium']={'A': 87.620000000000005, 'Symbol': 'Sr', 'Z': 38, 'Group': 2, 'Period': 5}
periodic_table['Yttrium']={'A': 88.905850000000001, 'Symbol': 'Y', 'Z': 39, 'Group': 3, 'Period': 5}
periodic_table['Zirconium']={'A': 91.224000000000004, 'Symbol': 'Zr', 'Z': 40, 'Group': 4, 'Period': 5}
periodic_table['Niobium']={'A': 92.906379999999999, 'Symbol': 'Nb', 'Z': 41, 'Group': 5, 'Period': 5}
periodic_table['Molybdenum']={'A': 95.939999999999998, 'Symbol': 'Mo', 'Z': 42, 'Group': 6, 'Period': 5}
periodic_table['Technetium']={'A': 98.0, 'Symbol': 'Tc', 'Z': 43, 'Group': 7, 'Period': 5}
periodic_table['Ruthenium']={'A': 101.06999999999999, 'Symbol': 'Ru', 'Z': 44, 'Group': 8, 'Period': 5}
periodic_table['Rhodium']={'A': 102.9055, 'Symbol': 'Rh', 'Z': 45, 'Group': 9, 'Period': 5}
periodic_table['Palladium']={'A': 106.42, 'Symbol': 'Pd', 'Z': 46, 'Group': 10, 'Period': 5}
periodic_table['Silver']={'A': 28.0855, 'Symbol': 'Ag', 'Z': 47, 'Group': 11, 'Period': 5}
periodic_table['Cadmium']={'A': 112.411, 'Symbol': 'Cd', 'Z': 48, 'Group': 12, 'Period': 5}
periodic_table['Indium']={'A': 114.818, 'Symbol': 'In', 'Z': 49, 'Group': 13, 'Period': 5}
periodic_table['Tin']={'A': 118.70999999999999, 'Symbol': 'Sn', 'Z': 50, 'Group': 14, 'Period': 5}
periodic_table['Antimony']={'A': 121.76000000000001, 'Symbol': 'Sb', 'Z': 51, 'Group': 15, 'Period': 5}
periodic_table['Tellurium']={'A': 127.59999999999999, 'Symbol': 'Te', 'Z': 52, 'Group': 16, 'Period': 5}
periodic_table['Iodine']={'A': 126.90447, 'Symbol': 'I', 'Z': 53, 'Group': 17, 'Period': 5}
periodic_table['Xenon']={'A': 131.29300000000001, 'Symbol': 'Xe', 'Z': 54, 'Group': 18, 'Period': 5}
periodic_table['Caesium']={'A': 132.90545, 'Symbol': 'Cs', 'Z': 55, 'Group': 1, 'Period': 6}
periodic_table['Barium']={'A': 137.327, 'Symbol': 'Ba', 'Z': 56, 'Group': 2, 'Period': 6}
periodic_table['Lutetium']={'A': 174.96700000000001, 'Symbol': 'Lu', 'Z': 71, 'Group': 3, 'Period': 6}
periodic_table['Hafnium']={'A': 178.49000000000001, 'Symbol': 'Hf', 'Z': 72, 'Group': 4, 'Period': 6}
periodic_table['Tantalum']={'A': 180.9479, 'Symbol': 'Ta', 'Z': 73, 'Group': 5, 'Period': 6}
periodic_table['Tungsten']={'A': 183.84, 'Symbol': 'W', 'Z': 74, 'Group': 6, 'Period': 6}
periodic_table['Rhenium']={'A': 186.20699999999999, 'Symbol': 'Re', 'Z': 75, 'Group': 7, 'Period': 6}
periodic_table['Osmium']={'A': 190.22999999999999, 'Symbol': 'Os', 'Z': 76, 'Group': 8, 'Period': 6}
periodic_table['Iridium']={'A': 192.21700000000001, 'Symbol': 'Ir', 'Z': 77, 'Group': 9, 'Period': 6}
periodic_table['Platinum']={'A': 195.078, 'Symbol': 'Pt', 'Z': 78, 'Group': 10, 'Period': 6}
periodic_table['Gold']={'A': 196.96655000000001, 'Symbol': 'Au', 'Z': 79, 'Group': 11, 'Period': 6}
periodic_table['Mercury']={'A': 200.59, 'Symbol': 'Hg', 'Z': 80, 'Group': 12, 'Period': 6}
periodic_table['Thallium']={'A': 204.38329999999999, 'Symbol': 'Tl', 'Z': 81, 'Group': 13, 'Period': 6}
periodic_table['Lead']={'A': 207.19999999999999, 'Symbol': 'Pb', 'Z': 82, 'Group': 14, 'Period': 6}
periodic_table['Bismuth']={'A': 208.98038, 'Symbol': 'Bi', 'Z': 83, 'Group': 15, 'Period': 6}
periodic_table['Polonium']={'A': 210.0, 'Symbol': 'Po', 'Z': 84, 'Group': 16, 'Period': 6}
periodic_table['Astatine']={'A': 210.0, 'Symbol': 'At', 'Z': 85, 'Group': 17, 'Period': 6}
periodic_table['Radon']={'A': 220.0, 'Symbol': 'Rn', 'Z': 86, 'Group': 18, 'Period': 6}
periodic_table['Francium']={'A': 223.0, 'Symbol': 'Fr', 'Z': 87, 'Group': 1, 'Period': 7}
periodic_table['Radium']={'A': 226.0, 'Symbol': 'Ra', 'Z': 88, 'Group': 2, 'Period': 7}
periodic_table['Rutherfordium']={'A': 2611.0, 'Symbol': 'Rf', 'Z': 104, 'Group': 4, 'Period': 7}
periodic_table['Lawrencium']={'A': 262.0, 'Symbol': 'Lr', 'Z': 103, 'Group': 3, 'Period': 7}
periodic_table['Dubnium']={'A': 262.0, 'Symbol': 'Db', 'Z': 105, 'Group': 5, 'Period': 7}
periodic_table['Bohrium']={'A': 264.0, 'Symbol': 'Bh', 'Z': 107, 'Group': 7, 'Period': 7}
periodic_table['Seaborgium']={'A': 266.0, 'Symbol': 'Sg', 'Z': 106, 'Group': 6, 'Period': 7}
periodic_table['Meitnerium']={'A': 268.0, 'Symbol': 'Mt', 'Z': 109, 'Group': 9, 'Period': 7}
periodic_table['Darmstadtium']={'A': 271.0, 'Symbol': 'Ds', 'Z': 110, 'Group': 10, 'Period': 7}
periodic_table['Roentgenium']={'A': 272.0, 'Symbol': 'Rg', 'Z': 111, 'Group': 11, 'Period': 7}
periodic_table['Hassium']={'A': 277.0, 'Symbol': 'Hs', 'Z': 108, 'Group': 8, 'Period': 7}
periodic_table['Ununtrium']={'A': 284.0, 'Symbol': 'Uut', 'Z': 113, 'Group': 13, 'Period': 7}
periodic_table['Ununbium']={'A': 285.0, 'Symbol': 'Uub', 'Z': 112, 'Group': 12, 'Period': 7}
periodic_table['Ununpentium']={'A': 288.0, 'Symbol': 'Uup', 'Z': 115, 'Group': 15, 'Period': 7}
periodic_table['Ununquadium']={'A': 289.0, 'Symbol': 'Uuq', 'Z': 114, 'Group': 14, 'Period': 7}
periodic_table['Ununhexium']={'A': 292.0, 'Symbol': 'Uuh', 'Z': 116, 'Group': 16, 'Period': 7}
for k in periodic_table.keys():
  S=periodic_table[k]['Symbol']
  A=periodic_table[k]['A']
  Z=periodic_table[k]['Z']
  G=periodic_table[k]['Group']
  P=periodic_table[k]['Period']
  exec(k+"=Atom(S,A,Z,G,P,name='%s')"%k)
#  exec(k+"=periodic_table['"+k+"']")
del k

Electron=Particle(charge=-echarge,mass=emass,Symbol='e-',name='Electron')
Positron=Particle(charge= echarge,mass=emass,Symbol='e+',name='Positron')

Proton = Particle(mass=  amu, Symbol='P',name='Proton')

Dihydrogen = Molecule(mass= 2.*amu, Symbol='H2',name='Dihydrogen')
Dinitrogen = Molecule(mass=28.*amu, Symbol='N2',name='Dinitrogen')
Dioxygen   = Molecule(mass=32.*amu, Symbol='O2',name='Dioxygen')
Carbon_Monoxide = Molecule(mass=28.*amu, Symbol='CO',name='Carbon_Monoxide')
Carbon_Dioxide  = Molecule(mass=44.*amu, Symbol='CO2',name='Carbon_Dioxide')
Water           = Molecule(mass=18.*amu, Symbol='H2O',name='Water')

class Species:
  """
Creates a new species of particles. All arguments are optional.
  - js=None: global species number (usually this should not be set)
  - type=Electron: one of either an elementary particle, an element, or one of
                   the predefined molecules
  - charge=echarge: charge in Coulombs, will be obtained automatically from
                    type or charge_state
  - mass=emass: mass in kg, will be obtained automatically from type
  - charge_state=0: charge_state of the ion or molecule
  - weight=None: simulation weight, will be obtained automatically from
                 other input
  - name='': species name
  - nautodt=1: number of species to setup for automatic subcycling.
               The slowest will have dt = 2**(nautodt-1)*top.dt.
  - fselfb=0.: z velocity to use when applying relativistic corrections.
               No corrections are done if it is zero.
  - limplicit=false: Flag to turn on the implicit particle advance for this
                     species.
  """
  def __init__(self,js=None,type=Electron,charge=echarge,mass=emass,charge_state=0,weight=None,name='',nautodt=1,
                    fselfb=None,limplicit=None):
    # --- Note that some default arguments are None in case the user had
    # --- set the values in pgroup already, in which case they should not
    # --- be overwritten here unless the inputs are explicitly set.
    self.jslist=[]
    self.type=type
    self.add_group(js,charge=charge,mass=mass,charge_state=charge_state,weight=weight)
    self.charge=top.pgroup.sq[self.jslist[0]]
    self.mass=top.pgroup.sm[self.jslist[0]]
    if type is None or type.__class__ is not Particle:
      self.charge_state=charge_state
    self.name=name
    self.nautodt = nautodt
    if self.nautodt > 1 and top.chdtspid == 0: top.chdtspid = nextpid()
    if fselfb is not None: top.pgroup.fselfb[self.jslist[-1]] = fselfb
    if limplicit is not None: top.pgroup.limplicit[self.jslist[-1]] = limplicit
    for i in range(nautodt-1):
      self.add_group()
      top.pgroup.ndts[self.jslist[-1]]=2*top.pgroup.ndts[self.jslist[-2]]
      if fselfb is not None: top.pgroup.fselfb[self.jslist[-1]] = fselfb
      if limplicit is not None: top.pgroup.limplicit[self.jslist[-1]] = limplicit
      # --- zero out sp_fract for the extra species added with larger ndts
      top.sp_fract[self.jslist[-1]] = 0.
        
  def __setstate__(self,dict):
    self.__dict__.update(dict)
    import species
    # --- Check if the type is already declared in this module.
    # --- If it is, replace it with the one defined here since they
    # --- are supposed to be singletons.
    try:
      self.type = species.__dict__[self.type.name]
    except (KeyError,AttributeError):
      pass

  def add_group(self,js=None,charge=None,mass=None,charge_state=None,weight=None):   
    if js is None:
      if top.pgroup.ns==0:
        condition = 1
      else:
        condition = top.pgroup.sm[0]<>0.
      if condition:
        addspecies()
      js=top.pgroup.ns-1
      if js > 0:
        top.pgroup.ins[js] = top.pgroup.ins[js-1]+top.pgroup.nps[js-1]
      top.pgroup.ipmax_s[js] = top.pgroup.ins[js]+top.pgroup.nps[js]-1
      top.pgroup.sid[js] = js
    self.jslist.append(js)
    type=self.type
    if charge is None:charge=self.charge
    if mass is None:mass=self.mass
    # set charge
    try:
      top.pgroup.sq[js]=type.charge
    except:
      if type is None or type.__class__ is not Particle:
        if charge_state is None:charge_state=self.charge_state
        top.pgroup.sq[js]=echarge*charge_state
      else:
        top.pgroup.sq[js]=charge
    top.zion_s[js]=nint(top.pgroup.sq[js]/echarge)
    # set mass
    try:
      top.pgroup.sm[js]=type.mass
    except:
      try: 
        top.pgroup.sm[js]=type.A*amu
        top.aion_s[js]=type.A
      except:
        top.pgroup.sm[js]=mass
    if weight is not None:
      top.pgroup.sw[js]=weight
    # set atomic number, if any
    try:
      top.nion_s[js]=type.Z
    except:
      pass
     
  def get_density(self,xmin=None,xmax=None,nx=None,ymin=None,ymax=None,ny=None,zmin=None,zmax=None,nz=None,lost=0,charge=0,dens=None):
    if xmin is None:xmin=w3d.xmmin
    if xmax is None:xmax=w3d.xmmax
    if ymin is None:ymin=w3d.ymmin
    if ymax is None:ymax=w3d.ymmax
    if zmin is None:zmin=w3d.zmminglobal
    if zmax is None:zmax=w3d.zmmaxglobal
    if dens is None:
      if nx is None:nx=w3d.nx
      if ny is None:ny=w3d.ny
      if nz is None:nz=w3d.nzfull
      density = fzeros([nx+1,ny+1,nz+1],'d')
      densityc = fzeros([nx+1,ny+1,nz+1],'d')
    else:
      nx = shape(dens)[0]-1
      ny = shape(dens)[1]-1
      nz = shape(dens)[2]-1
      density = dens
      densityc = fzeros([nx+1,ny+1,nz+1],'d')

    for js in self.jslist:
      #print js
      x=getx(js=js,lost=lost,gather=0)
      y=gety(js=js,lost=lost,gather=0)
      z=getz(js=js,lost=lost,gather=0)
      np=len(x)
      if np > 0:
        w=top.pgroup.sw[js]*ones(np,'d')
        if charge:w*=top.pgroup.sq[js]    
        deposgrid3d(1,np,x,y,z,w,nx,ny,nz,density,densityc,xmin,xmax,ymin,ymax,zmin,zmax)
    density[...] *= nx*ny*nz/((xmax-xmin)*(ymax-ymin)*(zmax-zmin))
    density[...] = parallelsum(density)
    if dens is None:return density
  
  def addpart(self,x,y,z,vx,vy,vz,gi=1.,js=None,lmomentum=false,**kw):
    if js is None:
      js=self.jslist[0]
    addparticles(x,y,z,vx,vy,vz,gi=gi,js=js,lmomentum=lmomentum,**kw)
      
  def add_uniform_box(self,np,xmin,xmax,ymin,ymax,zmin,zmax,vthx=0.,
                      vthy=0.,vthz=0.,vxmean=0.,vymean=0.,vzmean=0.,js=None,
                      lmomentum=0,spacing='random',nx=None,ny=None,nz=None,**kw):
    if spacing == 'random':
      x = RandomArray.random(np)
      y = RandomArray.random(np)
      z = RandomArray.random(np)

    elif spacing == 'uniform':
      if ymax > ymin: dims = 3
      else:           dims = 2
      if nx is None: nx = nint(np**(1./dims))
      if ny is None:
        if dims == 3: ny = nint(np**(1./dims))
        else:         ny = 1
      if nz is None: nz = nint(np**(1./dims))
      np = nx*ny*nz
      if dims == 3:
        x,y,z = getmesh3d(0.5/nx,1./nx,nx-1,
                          0.5/ny,1./ny,ny-1,
                          0.5/nz,1./nz,nz-1)
      else:
        x,z = getmesh2d(0.5/nx,1./nx,nx-1,
                        0.5/nz,1./nz,nz-1)
        y = zeros((nx,nz),'d')
      x.shape = (np,)
      y.shape = (np,)
      z.shape = (np,)

    x = xmin + (xmax - xmin)*x
    y = ymin + (ymax - ymin)*y
    z = zmin + (zmax - zmin)*z

    vx=RandomArray.normal(vxmean,vthx,np)
    vy=RandomArray.normal(vymean,vthy,np)
    vz=RandomArray.normal(vzmean,vthz,np)
    if lmomentum:
      gi=1./sqrt(1.+(vx*vx+vy*vy+vz*vz)/clight**2)
    else:
      gi=1.
    self.addpart(x,y,z,vx,vy,vz,js=js,gi=gi,**kw)
    
  def add_uniform_cylinder(self,np,rmax,zmin,zmax,vthx=0.,vthy=0.,vthz=0.,
                           xmean=0.,ymean=0,vxmean=0.,vymean=0.,vzmean=0.,js=None,
                           lmomentum=0,spacing='random',nr=None,nz=None,thetamin=0.,thetamax=2.*pi,
                           **kw):
    """Creates particles, uniformly filling a cylinder.
If top.wpid is nonzero, then the particles are uniformly spaced in radius and the
weights are set appropriately (weight=r/rmax). Otherwise, the particles are spaced uniformly
in radius squared.
 - np: total number of particles to load
 - rmax: radius of cylinder
 - zmin,zmax: z extent of the cylinder
 - vthx, vthy, vthz: thermal velocity, defaults to 0.
 - xmean,ymean: transverse center of the cylinder, defaults to 0.
 - vxmean,vymean,vzmean: directed velocity, defaults to 0.
 - js: particle species number, don't set it unless you mean it
 - lmomentum=false: Set to false when velocities are input as velocities, true
                    when input as massless momentum (as WARP stores them).
                    Only used when top.lrelativ is true.
 - spacing='random': either 'random' or 'uniform' particle spacing. For uniform,
                     r and z are uniform, theta is still random
 - nr,nz: for 'uniform' spacing, number of particles along r and z axis
 - thetamin=0.,thetamax=2.*pi: range of theta
    """

    if spacing == 'random':
      r=RandomArray.random(np)
      z=RandomArray.random(np)
    else:
      if nr is None: nr = nint(np**(1./2.))
      if nz is None: nz = nint(np**(1./2.))
      np = nr*nz
      r,z = getmesh2d(0.5/nr,1./nr,nr-1,
                      0.5/nz,1./nz,nz-1)
      r.shape = (np,)
      z.shape = (np,)
    theta=(thetamax-thetamin)*RandomArray.random(np) + thetamin

    if top.wpid == 0:
      r = sqrt(r)
    else:
      kw['w'] = 2*r
    x=xmean+rmax*r*cos(2.*pi*theta)
    y=ymean+rmax*r*sin(2.*pi*theta)
    z=zmin+(zmax-zmin)*z

    vx=RandomArray.normal(vxmean,vthx,np)
    vy=RandomArray.normal(vymean,vthy,np)
    vz=RandomArray.normal(vzmean,vthz,np)
    if lmomentum:
      gi=1./sqrt(1.+(vx*vx+vy*vy+vz*vz)/clight**2)
    else:
      gi=1
    self.addpart(x,y,z,vx,vy,vz,js=js,gi=gi,**kw)
    
  def add_gaussian_dist(self,np,deltax,deltay,deltaz,vthx=0.,vthy=0.,vthz=0.,
                        xmean=0.,ymean=0.,zmean=0.,vxmean=0.,vymean=0.,vzmean=0.,js=None,**kw):
    x=RandomArray.normal(xmean,deltax,np)
    y=RandomArray.normal(ymean,deltay,np)
    z=RandomArray.normal(zmean,deltaz,np)
    vx=RandomArray.normal(vxmean,vthx,np)
    vy=RandomArray.normal(vymean,vthy,np)
    vz=RandomArray.normal(vzmean,vthz,np)
    self.addpart(x,y,z,vx,vy,vz,js=js,**kw)
    
  def getn(self,**kw):
    return getn(jslist=self.jslist,**kw)
    
  def getx(self,**kw):
    return getx(jslist=self.jslist,**kw)
    
  def gety(self,**kw):
    return gety(jslist=self.jslist,**kw)
    
  def getz(self,**kw):
    return getz(jslist=self.jslist,**kw)
    
  def getr(self,**kw):
    return getr(jslist=self.jslist,**kw)

  def getvx(self,**kw):
    return getvx(jslist=self.jslist,**kw)
    
  def getvy(self,**kw):
    return getvy(jslist=self.jslist,**kw)
    
  def getvz(self,**kw):
    return getvz(jslist=self.jslist,**kw)
    
  def getux(self,**kw):
    return getux(jslist=self.jslist,**kw)
    
  def getuy(self,**kw):
    return getuy(jslist=self.jslist,**kw)
    
  def getuz(self,**kw):
    return getuz(jslist=self.jslist,**kw)
    
  def getxp(self,**kw):
    return getxp(jslist=self.jslist,**kw)
    
  def getyp(self,**kw):
    return getyp(jslist=self.jslist,**kw)
    
  def getrp(self,**kw):
    return getrp(jslist=self.jslist,**kw)
    
  def gettp(self,**kw):
    return gettp(jslist=self.jslist,**kw)
    
  def getgaminv(self,**kw):
    return getgaminv(jslist=self.jslist,**kw)
    
  def getpid(self,**kw):
    return getpid(jslist=self.jslist,**kw)

  def ppxy(self,**kw):
    return ppxy(jslist=self.jslist,**kw)

  def ppxxp(self,**kw):
    return ppxxp(jslist=self.jslist,**kw)

  def ppyyp(self,**kw):
    return ppyyp(jslist=self.jslist,**kw)

  def ppxpyp(self,**kw):
    return ppxpyp(jslist=self.jslist,**kw)

  def ppxvx(self,**kw):
    return ppxvx(jslist=self.jslist,**kw)

  def ppyvy(self,**kw):
    return ppyvy(jslist=self.jslist,**kw)

  def ppxvz(self,**kw):
    return ppxvz(jslist=self.jslist,**kw)

  def ppyvz(self,**kw):
    return ppyvz(jslist=self.jslist,**kw)

  def ppzxy(self,**kw):
    return ppzxy(jslist=self.jslist,**kw)

  def ppzx(self,**kw):
    return ppzx(jslist=self.jslist,**kw)

  def ppzy(self,**kw):
    return ppzy(jslist=self.jslist,**kw)

  def ppzr(self,**kw):
    return ppzr(jslist=self.jslist,**kw)

  def ppzxp(self,**kw):
    return ppzxp(jslist=self.jslist,**kw)

  def ppzvx(self,**kw):
    return ppzvx(jslist=self.jslist,**kw)

  def ppzyp(self,**kw):
    return ppzyp(jslist=self.jslist,**kw)

  def ppzvy(self,**kw):
    return ppzvy(jslist=self.jslist,**kw)

  def ppzvz(self,**kw):
    return ppzvz(jslist=self.jslist,**kw)

  def ppzrp(self,**kw):
    return ppzrp(jslist=self.jslist,**kw)

  def ppzvr(self,**kw):
    return ppzvr(jslist=self.jslist,**kw)

  def ppzvperp(self,**kw):
    return ppzvperp(jslist=self.jslist,**kw)

  def ppvzvperp(self,**kw):
    return ppvzvperp(jslist=self.jslist,**kw)

  def pptrace(self,**kw):
    return pptrace(jslist=self.jslist,**kw)

  def pprrp(self,**kw):
    return pprrp(jslist=self.jslist,**kw)

  def pprtp(self,**kw):
    return pprtp(jslist=self.jslist,**kw)

  def pprvz(self,**kw):
    return pprvz(jslist=self.jslist,**kw)
