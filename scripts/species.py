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
  def __init__(self,Symbol,composition,name=None):
    self.composition=composition
    mass=0.
    self.A=0
    self.Z=0
    for c in composition:
      mass+=c.mass
      self.A+=c.A
      self.Z+=c.Z
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

Dihydrogen = Molecule(composition=[Hydrogen,Hydrogen], Symbol='H2',name='Dihydrogen')
Dinitrogen = Molecule(composition=[Nitrogen,Nitrogen], Symbol='N2',name='Dinitrogen')
Dioxygen   = Molecule(composition=[Oxygen,Oxygen], Symbol='O2',name='Dioxygen')
Carbon_Monoxide = Molecule(composition=[Carbon,Oxygen], Symbol='CO',name='Carbon_Monoxide')
Carbon_Dioxide  = Molecule(composition=[Carbon,Oxygen,Oxygen], Symbol='CO2',name='Carbon_Dioxide')
Water           = Molecule(composition=[Hydrogen,Hydrogen,Oxygen], Symbol='H2O',name='Water')

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
  - efetch=1: Method to use for fetching the self fields. See documentation
              of top.efetch for more info.
  - fselfb=0.: z velocity to use when applying relativistic corrections.
               No corrections are done if it is zero.
  - limplicit=false: Flag to turn on the implicit particle advance for this
                     species.
  """
  def __init__(self,js=None,pgroup=top.pgroup,
                    type=Electron,charge=echarge,mass=emass,charge_state=0,
                    weight=None,name='',nautodt=1,
                    efetch=None,fselfb=None,limplicit=None):
    # --- Note that some default arguments are None in case the user had
    # --- set the values in pgroup already, in which case they should not
    # --- be overwritten here unless the inputs are explicitly set.
    self.jslist=[]
    self.pgroup=pgroup
    self.type=type
    self.add_group(js,charge=charge,mass=mass,charge_state=charge_state,weight=weight)
    self.charge=top.pgroup.sq[self.jslist[0]]
    self.mass=top.pgroup.sm[self.jslist[0]]
    if type is None or type.__class__ is not Particle:
      self.charge_state=charge_state
    self.name=name
    self.nautodt = nautodt
    if self.nautodt > 1 and top.chdtspid == 0: top.chdtspid = nextpid()
    if efetch is not None: top.efetch[self.jslist[-1]] = efetch
    if fselfb is not None: top.pgroup.fselfb[self.jslist[-1]] = fselfb
    if limplicit is not None: top.pgroup.limplicit[self.jslist[-1]] = limplicit
    for i in xrange(nautodt-1):
      self.add_group()
      top.pgroup.ndts[self.jslist[-1]]=2*top.pgroup.ndts[self.jslist[-2]]
      if efetch is not None: top.efetch[self.jslist[-1]] = efetch
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
     
  def get_density(self,xmin=None,xmax=None,nx=None,ymin=None,ymax=None,ny=None,zmin=None,zmax=None,
                  nz=None,lost=0,charge=0,dens=None,l_minmax_grid=true,l_dividebyvolume=1,l4symtry=None,l2symtry=None):
    if l_minmax_grid:
      if xmin is None:xmin=w3d.xmmin
      if xmax is None:xmax=w3d.xmmax
      if ymin is None:ymin=w3d.ymmin
      if ymax is None:ymax=w3d.ymmax
      if zmin is None:zmin=w3d.zmmin
      if zmax is None:zmax=w3d.zmmax
      if l4symtry is None:l4symtry=w3d.l4symtry
      if l2symtry is None:l2symtry=w3d.l2symtry
    else:
      if xmin is None:xmin=min(self.getx())
      if xmax is None:xmax=max(self.getx())
      if ymin is None:ymin=min(self.gety())
      if ymax is None:ymax=max(self.gety())
      if zmin is None:zmin=min(self.getz())
      if zmax is None:zmax=max(self.getz())
      if l4symtry is None:l4symtry=false
      if l2symtry is None:l2symtry=false
    if dens is None:
      if nx is None:nx=w3d.nx
      if ny is None:ny=w3d.ny
      if nz is None:nz=w3d.nz
      density = fzeros([nx+1,ny+1,nz+1],'d')
      densityc = fzeros([nx+1,ny+1,nz+1],'d')
    else:
      nx = shape(dens)[0]-1
      ny = shape(dens)[1]-1
      nz = shape(dens)[2]-1
      density = dens
      densityc = fzeros([nx+1,ny+1,nz+1],'d')

    np=0
    for js in self.jslist:
      np+=getn(js=js)
    if np==0:
      if dens is None:
        return density
      else:
        return
    for js in self.jslist:
      x=getx(js=js,lost=lost,gather=0)
      y=gety(js=js,lost=lost,gather=0)
      z=getz(js=js,lost=lost,gather=0)
      np=shape(x)[0]
      if np>0:
        if top.wpid==0:
          w=top.pgroup.sw[js]*ones(np,'d')
        else:
          w=top.pgroup.sw[js]*getpid(js=js,id=top.wpid-1,gather=0)
        if charge:w*=top.pgroup.sq[js]   
        deposgrid3d(1,np,x,y,z,w,nx,ny,nz,density,densityc,xmin,xmax,ymin,ymax,zmin,zmax)
    if l_dividebyvolume:
      density*=nx*ny*nz/((xmax-xmin)*(ymax-ymin)*(zmax-zmin))
      if l4symtry:
        density[0,:,:] *= 2
        density[:,0,:] *= 2
      if l2symtry:
        density[:,0,:] *= 2
    density[...] = parallelsum(density)
    if dens is None:return density
     
  def addpart(self,x,y,z,vx,vy,vz,gi=1.,js=None,lmomentum=false,**kw):
    if js is None:
      js=self.jslist[0]
    addparticles(x,y,z,vx,vy,vz,gi=gi,js=js,lmomentum=lmomentum,**kw)
      
  def add_uniform_box(self,np,xmin,xmax,ymin,ymax,zmin,zmax,vthx=0.,
                      vthy=0.,vthz=0.,vxmean=0.,vymean=0.,vzmean=0.,js=None,
                      lmomentum=0,spacing='random',nx=None,ny=None,nz=None,
                      lallindomain=0,**kw):
    """
 - lallindomain=0: If true, the code only loads particles within the domain. This
                   only matters when parallel.
    """

    if lallindomain:
      # --- Crop the zmin and zmax to be within the local domain
      zminp = max(zmin,min(zmax,top.zpminlocal+top.zgrid))
      zmaxp = min(zmax,max(zmin,top.zpmaxlocal+top.zgrid))
    else:
      # --- Crop the zmin and zmax to be within the global domain
      zminp = max(zmin,w3d.zmmin+top.zgrid)
      zmaxp = min(zmax,w3d.zmmax+top.zgrid)

    if spacing == 'random':
      # --- Adjust the number of particles to load to based on the
      # --- width of the cropped zmin and max and the original
      if zmin==zmax:
        if lallindomain:
          if zmin>top.zpmaxlocal+top.zgrid or zmax<top.zpminlocal+top.zgrid:
            return
        else:
          if zmin>w3d.zmmax+top.zgrid or zmax<w3d.zmmin+top.zgrid:return
      else:
        np = nint((zmaxp - zminp)/(zmax - zmin)*np)
      if np == 0: return
      x = RandomArray.random(np)
      y = RandomArray.random(np)
      z = RandomArray.random(np)
      z = (zminp + (zmaxp - zminp)*z - zmin)/(zmax - zmin)

    elif spacing == 'uniform':
      if ymax > ymin and zmax > zmin: dims = 3
      else:           dims = 2
      if nx is None: nx = nint(np**(1./dims))
      if ny is None:
        if dims == 3 or ymax>ymin: ny = nint(np**(1./dims))
        else:         ny = 1
      if nz is None:
        if dims == 3 or zmax>zmin: nz = nint(np**(1./dims))
        else:         nz = 1

      # --- Find the range of particle z locations within the cropped
      # --- zmin and max.
      if zmin==zmax:
        if lallindomain:
          if zmin>top.zpmaxlocal+top.zgrid or zmax<top.zpminlocal+top.zgrid:
            return
        else:
          if zmin>w3d.zmmax+top.zgrid or zmax<w3d.zmmin+top.zgrid:return
      else:
        dz = (zmax - zmin)/nz
        izminp = int((zminp - zmin)/dz + 0.5)
        izmaxp = int((zmaxp - zmin)/dz + 0.5)
        nzp = max(0,izmaxp - izminp)
        np = nx*ny*nzp
      if np == 0: return

      if dims == 3:
        x,y,z = getmesh3d(0.5/nx,1./nx,nx-1,
                          0.5/ny,1./ny,ny-1,
                          (izminp + 0.5)/nz,1./nz,nzp-1)
      else:
        if ymin==ymax:
          x,z = getmesh2d(0.5/nx,1./nx,nx-1,
                          (izminp + 0.5)/nz,1./nz,nzp-1)
          y = zeros((nx,nz),'d')
        if zmin==zmax:
          x,y = getmesh2d(0.5/nx,1./nx,nx-1,
                          0.5/ny,1./ny,ny-1)
          z = zeros((nx,ny),'d')

      # --- Perform a transpose so that the data is ordered with increasing z.
      # --- The copy is needed since transposed arrays cannot be reshaped.
      x = transpose(x).copy()
      y = transpose(y).copy()
      z = transpose(z).copy()
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

    kw['lallindomain'] = lallindomain
    self.addpart(x,y,z,vx,vy,vz,js=js,gi=gi,**kw)
    
  def add_uniform_cylinder(self,np,rmax,zmin,zmax,vthx=0.,vthy=0.,vthz=0.,
                           xmean=0.,ymean=0,zmean=0,vxmean=0.,vymean=0.,vzmean=0.,
                           theta=0.,phi=0.,
                           js=None,
                           lmomentum=0,spacing='random',nr=None,nz=None,
                           thetamin=0.,thetamax=2.*pi,
                           lvariableweights=None,lallindomain=0,
                           **kw):
    """Creates particles, uniformly filling a cylinder.
If top.wpid is nonzero, then the particles are uniformly spaced in radius and the
weights are set appropriately (weight=r/rmax). Otherwise, the particles are spaced uniformly
in radius squared.
 - np: total number of particles to load
 - rmax: radius of cylinder
 - zmin,zmax: z extent of the cylinder
 - vthx, vthy, vthz: thermal velocity, defaults to 0.
 - xmean,ymean,zmean: center of the cylinder, defaults to 0.
 - vxmean,vymean,vzmean: directed velocity, defaults to 0.
 - theta,phi: angle of cylinder about the center, 
              The transformation to lab frame is a rotation about the x axis by phi,
              followed by a rotation about the new y axis by theta.
              Note that the position mins and maxs and the velocity averages and means
              are relative to the rotated frame and are automatically transformed
              into the lab frame. The position means are relative to the lab frame.
 - js: particle species number, don't set it unless you mean it
 - lmomentum=false: Set to false when velocities are input as velocities, true
                    when input as massless momentum (as WARP stores them).
                    Only used when top.lrelativ is true.
 - spacing='random': either 'random' or 'uniform' particle spacing. For uniform,
                     r and z are uniform, theta is still random
 - nr,nz: for 'uniform' spacing, number of particles along r and z axis
 - thetamin=0.,thetamax=2.*pi: range of theta around the cylinder
 - lvariableweights: By default, if wpid is set, then the particles will be
                     weighted according to their radius, otherwise all will
                     have the same weight. Use this option to override the
                     default.
 - lallindomain=0: If true, the code only loads particles within the domain. This
                   only matters when parallel.
    """

    if theta == 0. and phi == 0.:
      # --- When no angle is specified, then the clipping to the domain
      # --- can be done directly on the zmin and zmax
      if lallindomain:
        # --- Crop the zmin and zmax to be within the local domain
        zminp = max(zmin,min(zmax,top.zpminlocal))
        zmaxp = min(zmax,max(zmin,top.zpmaxlocal))
      else:
        # --- Crop the zmin and zmax to be within the global domain
        zminp = max(zmin,w3d.zmmin)
        zmaxp = min(zmax,w3d.zmmax)
    else:
      # --- When angles are specified, the clipping must be done on a particle
      # --- by particle basis. Copy the zmin and zmax directly over to start.
      zminp = zmin
      zmaxp = zmax

    if spacing == 'random':
      # --- Adjust the number of particles to load to based on the
      # --- width of the cropped zmin and max and the original
      np = nint((zmaxp - zminp)/(zmax - zmin)*np)
      if np == 0: return
      r = RandomArray.random(np)
      z = RandomArray.random(np)
      z = (zminp + (zmaxp - zminp)*z - zmin)/(zmax - zmin)
    else:
      if nr is None: nr = nint(np**(1./2.))
      if nz is None: nz = nint(np**(1./2.))

      # --- Find the range of particle z locations within the cropped
      # --- zmin and max.
      dz = (zmax - zmin)/nz
      izminp = int((zminp - zmin)/dz + 0.5)
      izmaxp = int((zmaxp - zmin)/dz + 0.5)
      nzp = max(0,izmaxp - izminp)
      np = nr*nzp
      if np == 0: return

      r,z = getmesh2d(0.5/nr,1./nr,nr-1,
                      (izminp + 0.5)/nz,1./nz,nzp-1)

      # --- Perform a transpose so that the data is ordered with increasing z.
      # --- The copy is needed since transposed arrays cannot be reshaped.
      r = transpose(r).copy()
      z = transpose(z).copy()
      r.shape = (np,)
      z.shape = (np,)

    thetap=(thetamax-thetamin)*RandomArray.random(np) + thetamin

    if lvariableweights is None:
      lvariableweights = (top.wpid != 0)
    if lvariableweights and top.wpid == 0:
      top.wpid = nextpid()
    if not lvariableweights:
      r = sqrt(r)
      if top.wpid != 0: kw['w'] = 1.
    else:
      kw['w'] = 2*r

    x = rmax*r*cos(2.*pi*thetap)
    y = rmax*r*sin(2.*pi*thetap)
    z = zmin+(zmax-zmin)*z

    if theta != 0. or phi != 0.:
      # --- Transform positions from rotated frame into the lab frame.
      ct = cos(theta)
      st = sin(theta)
      cp = cos(phi)
      sp = sin(phi)
      x1 = +x*ct - y*st*sp + z*st*cp
      y1 =       + y*cp    + z*sp
      z1 = -x*st - y*ct*sp + z*ct*cp
      x,y,z = x1,y1,z1

    # --- Now add the means, after the rotation transform.
    x += xmean
    y += ymean
    z += zmean

    if theta != 0. or phi != 0.:
      # --- When angles are specified, the clipping must be done on a particle
      # --- by particle basis.
      if lallindomain:
        # --- Crop the z's to be within the local particle domain
        indomain = logical_and(top.zpminlocal<=z,z<top.zpmaxlocal)
      else:
        # --- Crop the z's to be within the global grid domain
        indomain = logical_and(w3d.zmmin<=z,z<w3d.zmmax)
      x = compress(indomain,x)
      y = compress(indomain,y)
      z = compress(indomain,z)
      np = len(z)
      if np == 0: return

    # --- Now the velocities are generated (after clipping the positions).
    vx = RandomArray.normal(vxmean,vthx,np)
    vy = RandomArray.normal(vymean,vthy,np)
    vz = RandomArray.normal(vzmean,vthz,np)

    if theta != 0. or phi != 0.:
      # --- Transform velocities from rotated frame into the lab frame.
      vx1 = +vx*ct - vy*st*sp + vz*st*cp
      vy1 =        + vy*cp    + vz*sp
      vz1 = -vx*st - vy*ct*sp + vz*ct*cp
      vx,vy,vz = vx1,vy1,vz1

    if lmomentum:
      gi=1./sqrt(1.+(vx*vx+vy*vy+vz*vz)/clight**2)
    else:
      gi=1

    kw['lallindomain'] = lallindomain
    self.addpart(x,y,z,vx,vy,vz,js=js,gi=gi,**kw)
    
  def add_gaussian_dist(self,np,deltax,deltay,deltaz,vthx=0.,vthy=0.,vthz=0.,
                        xmean=0.,ymean=0.,zmean=0.,vxmean=0.,vymean=0.,vzmean=0.,
                        zdist='random',nz=1000,fourfold=0,js=None,**kw):
    if fourfold:np=nint(float(np)/4)
    if zdist=='random':
      x=RandomArray.normal(0.,deltax,np)
      y=RandomArray.normal(0.,deltay,np)
      z=RandomArray.normal(0.,deltaz,np)
      vx=RandomArray.normal(0.,vthx,np)
      vy=RandomArray.normal(0.,vthy,np)
      vz=RandomArray.normal(0.,vthz,np)
      self.addpart(x+xmean,y+ymean,z+zmean,vx+vxmean,vy+vymean,vz+vzmean,js=js,**kw)
      if fourfold:
        self.addpart(x+xmean,-y+ymean,z+zmean,vx+vxmean,-vy+vymean,vz+vzmean,lallindomain=true)
        self.addpart(-x+xmean,y+ymean,z+zmean,-vx+vxmean,vy+vymean,vz+vzmean,lallindomain=true)
        self.addpart(-x+xmean,-y+ymean,z+zmean,-vx+vxmean,-vy+vymean,vz+vzmean,lallindomain=true)
    if zdist=='regular': 
      dz=16.*deltaz/nz
      zmin=-(float(nz/2)-0.5)*dz
      for i in xrange(nz):
        zadd=zmin+i*dz
        Nadd =max(0.,np*dz*exp(-0.5*(zadd/deltaz)**2)/(sqrt(2.*pi)*deltaz))
        if ranf()<(Nadd-int(Nadd)):
          Nadd=int(Nadd)+1
        else:
          Nadd=int(Nadd)
        if Nadd>0:
          x=RandomArray.normal(0.,deltax,Nadd)
          y=RandomArray.normal(0.,deltay,Nadd)
          z=zadd+dz*(RandomArray.random(Nadd)-0.5)
          vx=RandomArray.normal(0.,vthx,Nadd)
          vy=RandomArray.normal(0.,vthy,Nadd)
          vz=RandomArray.normal(0.,vthz,Nadd)
          self.addpart(x+xmean,y+ymean,z+zmean,vx+vxmean,vy+vymean,vz+vzmean,lallindomain=true)
    if fourfold:
      self.addpart(x+xmean,-y+ymean,z+zmean,vx+vxmean,-vy+vymean,vz+vzmean,lallindomain=true)
      self.addpart(-x+xmean,y+ymean,z+zmean,-vx+vxmean,vy+vymean,vz+vzmean,lallindomain=true)
      self.addpart(-x+xmean,-y+ymean,z+zmean,-vx+vxmean,-vy+vymean,vz+vzmean,lallindomain=true)
    
  def gather_zmmnts_locs(self):
    get_zmmnts_stations(len(self.jslist),
                              array(self.jslist),
                              self.pgroup,
                              self.nzmmnts_locs,
                              self.zmmnts_locs[0],
                              self.zmmnts_locs[-1],
                              self.zmmnts_locs_vfrm,
                              self.zmmnts_locs_pnum,
                              self.zmmnts_locs_xbar,
                              self.zmmnts_locs_ybar,
                              self.zmmnts_locs_xpbar,
                              self.zmmnts_locs_ypbar,
                              self.zmmnts_locs_x2,
                              self.zmmnts_locs_y2,
                              self.zmmnts_locs_xp2,
                              self.zmmnts_locs_yp2,
                              self.zmmnts_locs_xxp,
                              self.zmmnts_locs_yyp)
    self.zmmnts_gathered=false
    
  def gatherall_zmmnts_locs(self):
    self.zmmnts_locs_pnum = parallelsum(self.zmmnts_locs_pnum)
    self.zmmnts_locs_xbar = parallelsum(self.zmmnts_locs_xbar)
    self.zmmnts_locs_ybar = parallelsum(self.zmmnts_locs_ybar)
    self.zmmnts_locs_xpbar = parallelsum(self.zmmnts_locs_xpbar)
    self.zmmnts_locs_ypbar = parallelsum(self.zmmnts_locs_ypbar)
    self.zmmnts_locs_x2 = parallelsum(self.zmmnts_locs_x2)
    self.zmmnts_locs_y2 = parallelsum(self.zmmnts_locs_y2)
    self.zmmnts_locs_xp2 = parallelsum(self.zmmnts_locs_xp2)
    self.zmmnts_locs_yp2 = parallelsum(self.zmmnts_locs_yp2)
    self.zmmnts_locs_xxp = parallelsum(self.zmmnts_locs_xxp)
    self.zmmnts_locs_yyp = parallelsum(self.zmmnts_locs_yyp)
    if me>0:
      self.zmmnts_locs_pnum[...]=0.
      self.zmmnts_locs_xbar[...]=0.
      self.zmmnts_locs_ybar[...]=0.
      self.zmmnts_locs_xpbar[...]=0.
      self.zmmnts_locs_ypbar[...]=0.
      self.zmmnts_locs_x2[...]=0.
      self.zmmnts_locs_y2[...]=0.
      self.zmmnts_locs_xp2[...]=0.
      self.zmmnts_locs_yp2[...]=0.
      self.zmmnts_locs_xxp[...]=0.
      self.zmmnts_locs_yyp[...]=0.
    self.zmmnts_locs_gathered=true

  def set_zmmnts_locs(self,zmin,zmax,nz,vfrm=0.):
    self.zmmnts_locs  = getmeshcoordinates([zmin],[(zmax-zmin)/(nz-1)],[nz-1])[0]
    self.nzmmnts_locs = nz
    self.zmmnts_locs_vfrm = vfrm
    self.dzmmnts_locs = (zmax-zmin)/(nz-1)
    self.zmmnts_locs_pnum = zeros(self.nzmmnts_locs,'d')
    self.zmmnts_locs_xbar = zeros(self.nzmmnts_locs,'d')
    self.zmmnts_locs_ybar = zeros(self.nzmmnts_locs,'d')
    self.zmmnts_locs_zbar = zeros(self.nzmmnts_locs,'d')
    self.zmmnts_locs_xpbar = zeros(self.nzmmnts_locs,'d')
    self.zmmnts_locs_ypbar = zeros(self.nzmmnts_locs,'d')
    self.zmmnts_locs_x2 = zeros(self.nzmmnts_locs,'d')
    self.zmmnts_locs_y2 = zeros(self.nzmmnts_locs,'d')
    self.zmmnts_locs_z2 = zeros(self.nzmmnts_locs,'d')
    self.zmmnts_locs_xp2 = zeros(self.nzmmnts_locs,'d')
    self.zmmnts_locs_yp2 = zeros(self.nzmmnts_locs,'d')
    self.zmmnts_locs_xxp = zeros(self.nzmmnts_locs,'d')
    self.zmmnts_locs_yyp = zeros(self.nzmmnts_locs,'d')
    if top.xoldpid==0:top.xoldpid=nextpid()
    if top.yoldpid==0:top.yoldpid=nextpid()
    if top.zoldpid==0:top.zoldpid=nextpid()
    if top.uxoldpid==0:top.uxoldpid=nextpid()
    if top.uyoldpid==0:top.uyoldpid=nextpid()
    if top.uzoldpid==0:top.uzoldpid=nextpid()
    self.zmmnts_locs_gathered=true
    installafterstep(self.gather_zmmnts_locs)
   
  def set_zmmnts(self):
    self.zmmnts_pnum = AppendableArray(typecode='d')
    self.zmmnts_xbar =AppendableArray(typecode='d')
    self.zmmnts_ybar =AppendableArray(typecode='d')
    self.zmmnts_zbar =AppendableArray(typecode='d')
    self.zmmnts_xpbar =AppendableArray(typecode='d')
    self.zmmnts_ypbar =AppendableArray(typecode='d')
    self.zmmnts_xpnbar =AppendableArray(typecode='d')
    self.zmmnts_ypnbar =AppendableArray(typecode='d')
    self.zmmnts_x2 =AppendableArray(typecode='d')
    self.zmmnts_y2 =AppendableArray(typecode='d')
    self.zmmnts_z2 =AppendableArray(typecode='d')
    self.zmmnts_xp2 =AppendableArray(typecode='d')
    self.zmmnts_yp2 =AppendableArray(typecode='d')
    self.zmmnts_xxp =AppendableArray(typecode='d')
    self.zmmnts_yyp =AppendableArray(typecode='d')
    self.zmmnts_xpn2 =AppendableArray(typecode='d')
    self.zmmnts_ypn2 =AppendableArray(typecode='d')
    self.zmmnts_xxpn =AppendableArray(typecode='d')
    self.zmmnts_yypn =AppendableArray(typecode='d')
    self.zmmnts_gathered=true
    installafterstep(self.gather_zmmnts)
   
  def gather_zmmnts(self):
    print 'gather_zmmnts'
    self.zmmnts_pnum.append(self.getn(gather=0))
    xbar = 0.
    ybar = 0.
    zbar = 0.
    xpbar = 0.
    ypbar = 0.
    xpnbar = 0.
    ypnbar = 0.
    x2 = 0.
    y2 = 0.
    z2 = 0.
    xp2 = 0.
    yp2 = 0.
    xxpbar = 0.
    yypbar = 0.
    xpn2 = 0.
    ypn2 = 0.
    xxpnbar = 0.
    yypnbar = 0.
    pg=self.pgroup
    nparpgrp = 1024
    xpa = zeros(nparpgrp,'d')
    ypa = zeros(nparpgrp,'d')
    xpna = zeros(nparpgrp,'d')
    ypna = zeros(nparpgrp,'d')
    betaa = zeros(nparpgrp,'d')
    for js in self.jslist:
      if pg.nps[js]==0:continue
      ng = 1+pg.nps[js]/nparpgrp
      for ig in range(ng):
        il = pg.ins[js]-1+nparpgrp*ig
        iu = min(il+nparpgrp,pg.ins[js]-1+pg.nps[js])
        np = iu-il
        x = pg.xp[il:iu]
        y = pg.yp[il:iu]
        z = pg.zp[il:iu]        
        xpa[:np] = pg.uxp[il:iu]/pg.uzp[il:iu]
        ypa[:np] = pg.uyp[il:iu]/pg.uzp[il:iu]
        gaminv =  pg.gaminv[il:iu]
        betaa[:np] = (1.-gaminv)*(1.+gaminv)
        xpna[:np] = xpa[:np]*betaa[:np]/gaminv
        ypna[:np] = ypa[:np]*betaa[:np]*gaminv
        xp=xpa[:np]
        yp=ypa[:np]
        xpn=xpna[:np]
        ypn=ypna[:np]
        beta=betaa[:np]
        xbar+=sum(x)      
        ybar+=sum(y)      
        zbar+=sum(z)      
        xpbar+=sum(xp)      
        ypbar+=sum(yp)      
        xpnbar+=sum(xpn)      
        ypnbar+=sum(ypn)      
        x2+=sum(x*x)      
        y2+=sum(y*y)      
        z2+=sum(z*z)      
        xp2+=sum(xp*xp)      
        yp2+=sum(yp*yp)      
        xpn2+=sum(xpn*xpn)      
        ypn2+=sum(ypn*ypn)      
        xxpbar+=sum(x*xp)      
        yypbar+=sum(y*yp)        
        xxpnbar+=sum(x*xpn)      
        yypnbar+=sum(y*ypn)      
    self.zmmnts_xbar.append(xbar)
    self.zmmnts_ybar.append(ybar)
    self.zmmnts_zbar.append(zbar)
    self.zmmnts_xpbar.append(xpbar)
    self.zmmnts_ypbar.append(ypbar)
    self.zmmnts_x2.append(x2)
    self.zmmnts_y2.append(y2)
    self.zmmnts_z2.append(z2)
    self.zmmnts_xp2.append(xp2)
    self.zmmnts_yp2.append(yp2)
    self.zmmnts_xxp.append(xxpbar)
    self.zmmnts_yyp.append(yypbar)
    self.zmmnts_xpnbar.append(xpnbar)
    self.zmmnts_ypnbar.append(ypnbar)
    self.zmmnts_xpn2.append(xpn2)
    self.zmmnts_ypn2.append(ypn2)
    self.zmmnts_xxpn.append(xxpnbar)
    self.zmmnts_yypn.append(yypnbar)
    self.zmmnts_gathered=false
    print 'gather_zmmnts_done'
    
  def gatherall_zmmnts(self):
    self.zmmnts_pnum.data()[...] = parallelsum(self.zmmnts_pnum.data())
    self.zmmnts_xbar.data()[...] = parallelsum(self.zmmnts_xbar.data())
    self.zmmnts_ybar.data()[...] = parallelsum(self.zmmnts_ybar.data())
    self.zmmnts_xpbar.data()[...] = parallelsum(self.zmmnts_xpbar.data())
    self.zmmnts_ypbar.data()[...] = parallelsum(self.zmmnts_ypbar.data())
    self.zmmnts_xpnbar.data()[...] = parallelsum(self.zmmnts_xpnbar.data())
    self.zmmnts_ypnbar.data()[...] = parallelsum(self.zmmnts_ypnbar.data())
    self.zmmnts_x2.data()[...] = parallelsum(self.zmmnts_x2.data())
    self.zmmnts_y2.data()[...] = parallelsum(self.zmmnts_y2.data())
    self.zmmnts_xp2.data()[...] = parallelsum(self.zmmnts_xp2.data())
    self.zmmnts_yp2.data()[...] = parallelsum(self.zmmnts_yp2.data())
    self.zmmnts_xxp.data()[...] = parallelsum(self.zmmnts_xxp.data())
    self.zmmnts_yyp.data()[...] = parallelsum(self.zmmnts_yyp.data())
    self.zmmnts_xpn2.data()[...] = parallelsum(self.zmmnts_xpn2.data())
    self.zmmnts_ypn2.data()[...] = parallelsum(self.zmmnts_ypn2.data())
    self.zmmnts_xxpn.data()[...] = parallelsum(self.zmmnts_xxpn.data())
    self.zmmnts_yypn.data()[...] = parallelsum(self.zmmnts_yypn.data())
    if me>0:
      self.zmmnts_pnum.data()[...]=0.
      self.zmmnts_xbar.data()[...]=0.
      self.zmmnts_ybar.data()[...]=0.
      self.zmmnts_xpbar.data()[...]=0.
      self.zmmnts_ypbar.data()[...]=0.
      self.zmmnts_xpnbar.data()[...]=0.
      self.zmmnts_ypnbar.data()[...]=0.
      self.zmmnts_x2.data()[...]=0.
      self.zmmnts_y2.data()[...]=0.
      self.zmmnts_xp2.data()[...]=0.
      self.zmmnts_yp2.data()[...]=0.
      self.zmmnts_xxp.data()[...]=0.
      self.zmmnts_yyp.data()[...]=0.
      self.zmmnts_xpn2.data()[...]=0.
      self.zmmnts_ypn2.data()[...]=0.
      self.zmmnts_xxpn.data()[...]=0.
      self.zmmnts_yypn.data()[...]=0.
    self.zmmnts_gathered=true

  def getzmmnts_pnum(self): 
    if not self.zmmnts_gathered:self.gatherall_zmmnts()
    if me==0:
      return self.zmmnts_pnum.data()
    else:
      return None
      
  def getzmmnts_xrms(self): 
    if not self.zmmnts_gathered:self.gatherall_zmmnts()
    if me==0:
      pnum = where(self.zmmnts_pnum.data()>smallpos,self.zmmnts_pnum.data(),1.)
      return sqrt(self.zmmnts_x2.data()/pnum)
    else:
      return None
      
  def getzmmnts_yrms(self): 
    if not self.zmmnts_gathered:self.gatherall_zmmnts()
    if me==0:
      pnum = where(self.zmmnts_pnum.data()>smallpos,self.zmmnts_pnum.data(),1.)
      return sqrt(self.zmmnts_y2.data()/pnum)
    else:
      return None
      
  def getzmmnts_xprms(self): 
    if not self.zmmnts_gathered:self.gatherall_zmmnts()
    if me==0:
      pnum = where(self.zmmnts_pnum.data()>smallpos,self.zmmnts_pnum.data(),1.)
      return sqrt(self.zmmnts_xp2.data()/pnum)
    else:
      return None
      
  def getzmmnts_yprms(self): 
    if not self.zmmnts_gathered:self.gatherall_zmmnts()
    if me==0:
      pnum = where(self.zmmnts_pnum.data()>smallpos,self.zmmnts_pnum.data(),1.)
      return sqrt(self.zmmnts_yp2.data()/pnum)
    else:
      return None
      
  def getzmmnts_xbar(self): 
    if not self.zmmnts_gathered:self.gatherall_zmmnts()
    if me==0:
      pnum = where(self.zmmnts_pnum.data()>smallpos,self.zmmnts_pnum.data(),1.)
      return self.zmmnts_xbar.data()/pnum
    else:
      return None
      
  def getzmmnts_ybar(self): 
    if not self.zmmnts_gathered:self.gatherall_zmmnts()
    if me==0:
      pnum = where(self.zmmnts_pnum.data()>smallpos,self.zmmnts_pnum.data(),1.)
      return self.zmmnts_ybar.data()/pnum
    else:
      return None
      
  def getzmmnts_xpbar(self): 
    if not self.zmmnts_gathered:self.gatherall_zmmnts()
    if me==0:
      pnum = where(self.zmmnts_pnum.data()>smallpos,self.zmmnts_pnum.data(),1.)
      return self.zmmnts_xpbar.data()/pnum
    else:
      return None
      
  def getzmmnts_ypbar(self): 
    if not self.zmmnts_gathered:self.gatherall_zmmnts()
    if me==0:
      pnum = where(self.zmmnts_pnum.data()>smallpos,self.zmmnts_pnum.data(),1.)
      return self.zmmnts_ypbar.data()/pnum
    else:
      return None
      
  def getzmmnts_xxpbar(self): 
    if not self.zmmnts_gathered:self.gatherall_zmmnts()
    if me==0:
      pnum = where(self.zmmnts_pnum.data()>smallpos,self.zmmnts_pnum.data(),1.)
      return self.zmmnts_xxp.data()/pnum
    else:
      return None
      
  def getzmmnts_yypbar(self): 
    if not self.zmmnts_gathered:self.gatherall_zmmnts()
    if me==0:
      pnum = where(self.zmmnts_pnum.data()>smallpos,self.zmmnts_pnum.data(),1.)
      return self.zmmnts_yyp.data()/pnum
    else:
      return None
      
  def getzmmnts_emitxrms(self):
    if not self.zmmnts_gathered:self.gatherall_zmmnts()
    if me>0:return None
    pnum = where(self.zmmnts_pnum.data()>smallpos,self.zmmnts_pnum.data(),1.)
    xbar = self.getzmmnts_xbar()
    xpbar = self.getzmmnts_xpbar()
    xxp = self.zmmnts_xxp.data()/pnum
    x2  = self.zmmnts_x2.data()/pnum
    xp2 = self.zmmnts_xp2.data()/pnum
    return sqrt((x2-xbar*xbar)*(xp2-xpbar*xpbar)-(xxp-xbar*xpbar)**2)

  def getzmmnts_emityrms(self):
    if not self.zmmnts_gathered:self.gatherall_zmmnts()
    if me>0:return None
    pnum = where(self.zmmnts_pnum.data()>smallpos,self.zmmnts_pnum.data(),1.)
    ybar = self.getzmmnts_ybar()
    ypbar = self.getzmmnts_ypbar()
    yyp = self.zmmnts_yyp.data()/pnum
    y2  = self.zmmnts_y2.data()/pnum
    yp2 = self.zmmnts_yp2.data()/pnum
    return sqrt((y2-ybar*ybar)*(yp2-ypbar*ypbar)-(yyp-ybar*ypbar)**2)
      
  def getzmmnts_emitxnrms(self):
    if not self.zmmnts_gathered:self.gatherall_zmmnts()
    if me>0:return None
    pnum = where(self.zmmnts_pnum.data()>smallpos,self.zmmnts_pnum.data(),1.)
    xbar = self.getzmmnts_xbar()
    xpnbar = self.zmmnts_xpnbar.data()/pnum
    xxpn = self.zmmnts_xxpn.data()/pnum
    x2  = self.zmmnts_x2.data()/pnum
    xpn2 = self.zmmnts_xpn2.data()/pnum
    return sqrt((x2-xbar*xbar)*(xpn2-xpnbar*xpnbar)-(xxpn-xbar*xpnbar)**2)

  def getzmmnts_emitynrms(self):
    if not self.zmmnts_gathered:self.gatherall_zmmnts()
    if me>0:return None
    pnum = where(self.zmmnts_pnum.data()>smallpos,self.zmmnts_pnum.data(),1.)
    ybar = self.getzmmnts_ybar()
    ypnbar = self.zmmnts_ypnbar.data()/pnum
    yypn = self.zmmnts_yypn.data()/pnum
    y2  = self.zmmnts_y2.data()/pnum
    ypn2 = self.zmmnts_ypn2.data()/pnum
    return sqrt((y2-ybar*ybar)*(ypn2-ypnbar*ypnbar)-(yypn-ybar*ypnbar)**2)

  def plzmmnts_data(self,x,color=black,width=1,type='solid',xscale=1.,yscale=1.,xoffset=0.):
    if me==0:
      zst = self.zmmnts_locs
      pla(yscale*x,xscale*(zst-xoffset),color=color,width=width)

  def plzmmnts_xbar(self,color=black,width=1,type='solid',xscale=1.,yscale=1.,xoffset=0.):  
    self.plzmmnts_data(self.getzmmnts_xbar(),color,width,type,xscale,yscale,xoffset)
    
  def plzmmnts_ybar(self,color=black,width=1,type='solid',xscale=1.,yscale=1.,xoffset=0.):  
    self.plzmmnts_data(self.getzmmnts_ybar(),color,width,type,xscale,yscale,xoffset)
    
  def plzmmnts_xrms(self,color=black,width=1,type='solid',xscale=1.,yscale=1.,xoffset=0.):  
    self.plzmmnts_data(self.getzmmnts_xrms(),color,width,type,xscale,yscale,xoffset)
    
  def plzmmnts_yrms(self,color=black,width=1,type='solid',xscale=1.,yscale=1.,xoffset=0.):  
    self.plzmmnts_data(self.getzmmnts_yrms(),color,width,type,xscale,yscale,xoffset)
    
  def plzmmnts_xprms(self,color=black,width=1,type='solid',xscale=1.,yscale=1.,xoffset=0.):  
    self.plzmmnts_data(self.getzmmnts_xprms(),color,width,type,xscale,yscale,xoffset)
    
  def plzmmnts_yprms(self,color=black,width=1,type='solid',xscale=1.,yscale=1.,xoffset=0.):  
    self.plzmmnts_data(self.getzmmnts_yprms(),color,width,type,xscale,yscale,xoffset)
    
  def plzmmnts_xpbar(self,color=black,width=1,type='solid',xscale=1.,yscale=1.,xoffset=0.):  
    self.plzmmnts_data(self.getzmmnts_xpbar(),color,width,type,xscale,yscale,xoffset)
    
  def plzmmnts_ypbar(self,color=black,width=1,type='solid',xscale=1.,yscale=1.,xoffset=0.):  
    self.plzmmnts_data(self.getzmmnts_ypbar(),color,width,type,xscale,yscale,xoffset)
    
  def plzmmnts_xxpbar(self,color=black,width=1,type='solid',xscale=1.,yscale=1.,xoffset=0.):  
    self.plzmmnts_data(self.getzmmnts_xxpbar(),color,width,type,xscale,yscale,xoffset)
    
  def plzmmnts_yypbar(self,color=black,width=1,type='solid',xscale=1.,yscale=1.,xoffset=0.):  
    self.plzmmnts_data(self.getzmmnts_yypbar(),color,width,type,xscale,yscale,xoffset)
    
  def plzmmnts_emitx(self,color=black,width=1,type='solid',xscale=1.,yscale=1.,xoffset=0.):  
    self.plzmmnts_data(self.getzmmnts_emitxrms(),color,width,type,xscale,yscale,xoffset)
    
  def plzmmnts_emity(self,color=black,width=1,type='solid',xscale=1.,yscale=1.,xoffset=0.):  
    self.plzmmnts_data(self.getzmmnts_emityrms(),color,width,type,xscale,yscale,xoffset)
    
  def getzmmnts_locs_pnum(self): 
    if not self.zmmnts_locs_gathered:self.gatherall_zmmnts()
    if me==0:
      return self.zmmnts_locs_pnum
    else:
      return None
      
  def getzmmnts_locs_xrms(self): 
    if not self.zmmnts_locs_gathered:self.gatherall_zmmnts()
    if me==0:
      pnum = where(self.zmmnts_locs_pnum>smallpos,self.zmmnts_locs_pnum,1.)
      return sqrt(self.zmmnts_locs_x2/pnum)
    else:
      return None
      
  def getzmmnts_locs_yrms(self): 
    if not self.zmmnts_locs_gathered:self.gatherall_zmmnts()
    if me==0:
      pnum = where(self.zmmnts_locs_pnum>smallpos,self.zmmnts_locs_pnum,1.)
      return sqrt(self.zmmnts_locs_y2/pnum)
    else:
      return None
      
  def getzmmnts_locs_xprms(self): 
    if not self.zmmnts_locs_gathered:self.gatherall_zmmnts()
    if me==0:
      pnum = where(self.zmmnts_locs_pnum>smallpos,self.zmmnts_locs_pnum,1.)
      return sqrt(self.zmmnts_locs_xp2/pnum)
    else:
      return None
      
  def getzmmnts_locs_yprms(self): 
    if not self.zmmnts_locs_gathered:self.gatherall_zmmnts()
    if me==0:
      pnum = where(self.zmmnts_locs_pnum>smallpos,self.zmmnts_locs_pnum,1.)
      return sqrt(self.zmmnts_locs_yp2/pnum)
    else:
      return None
      
  def getzmmnts_locs_xbar(self): 
    if not self.zmmnts_locs_gathered:self.gatherall_zmmnts()
    if me==0:
      pnum = where(self.zmmnts_locs_pnum>smallpos,self.zmmnts_locs_pnum,1.)
      return self.zmmnts_locs_xbar/pnum
    else:
      return None
      
  def getzmmnts_locs_ybar(self): 
    if not self.zmmnts_locs_gathered:self.gatherall_zmmnts()
    if me==0:
      pnum = where(self.zmmnts_locs_pnum>smallpos,self.zmmnts_locs_pnum,1.)
      return self.zmmnts_locs_ybar/pnum
    else:
      return None
      
  def getzmmnts_locs_xpbar(self): 
    if not self.zmmnts_locs_gathered:self.gatherall_zmmnts()
    if me==0:
      pnum = where(self.zmmnts_locs_pnum>smallpos,self.zmmnts_locs_pnum,1.)
      return self.zmmnts_locs_xpbar/pnum
    else:
      return None
      
  def getzmmnts_locs_ypbar(self): 
    if not self.zmmnts_locs_gathered:self.gatherall_zmmnts()
    if me==0:
      pnum = where(self.zmmnts_locs_pnum>smallpos,self.zmmnts_locs_pnum,1.)
      return self.zmmnts_locs_ypbar/pnum
    else:
      return None
      
  def getzmmnts_locs_xxpbar(self): 
    if not self.zmmnts_locs_gathered:self.gatherall_zmmnts()
    if me==0:
      pnum = where(self.zmmnts_locs_pnum>smallpos,self.zmmnts_locs_pnum,1.)
      return self.zmmnts_locs_xxp/pnum
    else:
      return None
      
  def getzmmnts_locs_yypbar(self): 
    if not self.zmmnts_locs_gathered:self.gatherall_zmmnts()
    if me==0:
      pnum = where(self.zmmnts_locs_pnum>smallpos,self.zmmnts_locs_pnum,1.)
      return self.zmmnts_locs_yyp/pnum
    else:
      return None
      
  def getzmmnts_locs_emitxrms(self):
    if not self.zmmnts_locs_gathered:self.gatherall_zmmnts()
    if me>0:return None
    pnum = where(self.zmmnts_locs_pnum>smallpos,self.zmmnts_locs_pnum,1.)
    xbar = self.getzmmnts_locs_xbar()
    xpbar = self.getzmmnts_locs_xpbar()
    xxp = self.zmmnts_locs_xxp/pnum
    x2  = self.zmmnts_locs_x2/pnum
    xp2 = self.zmmnts_locs_xp2/pnum
    return sqrt((x2-xbar*xbar)*(xp2-xpbar*xpbar)-(xxp-xbar*xpbar)**2)

  def getzmmnts_locs_emityrms(self):
    if not self.zmmnts_locs_gathered:self.gatherall_zmmnts()
    if me>0:return None
    pnum = where(self.zmmnts_locs_pnum>smallpos,self.zmmnts_locs_pnum,1.)
    ybar = self.getzmmnts_locs_ybar()
    ypbar = self.getzmmnts_locs_ypbar()
    yyp = self.getzmmnts_locs_yypbar()
    y2  = self.zmmnts_locs_y2/pnum
    yp2 = self.zmmnts_locs_yp2/pnum
    return sqrt((y2-ybar*ybar)*(yp2-ypbar*ypbar)-(yyp-ybar*ypbar)**2)

  def plzmmnts_locs_data(self,x,color=black,width=1,type='solid',xscale=1.,yscale=1.,xoffset=0.):
    if me==0:
      zst = self.zmmnts_locs
      pla(yscale*x,xscale*(zst-xoffset),color=color,width=width)

  def plzmmnts_locs_xbar(self,color=black,width=1,type='solid',xscale=1.,yscale=1.,xoffset=0.):  
    self.plzmmnts_locs_data(self.getzmmnts_locs_xbar(),color,width,type,xscale,yscale,xoffset)
    
  def plzmmnts_locs_ybar(self,color=black,width=1,type='solid',xscale=1.,yscale=1.,xoffset=0.):  
    self.plzmmnts_locs_data(self.getzmmnts_locs_ybar(),color,width,type,xscale,yscale,xoffset)
    
  def plzmmnts_locs_xrms(self,color=black,width=1,type='solid',xscale=1.,yscale=1.,xoffset=0.):  
    self.plzmmnts_locs_data(self.getzmmnts_locs_xrms(),color,width,type,xscale,yscale,xoffset)
    
  def plzmmnts_locs_yrms(self,color=black,width=1,type='solid',xscale=1.,yscale=1.,xoffset=0.):  
    self.plzmmnts_locs_data(self.getzmmnts_locs_yrms(),color,width,type,xscale,yscale,xoffset)
    
  def plzmmnts_locs_xprms(self,color=black,width=1,type='solid',xscale=1.,yscale=1.,xoffset=0.):  
    self.plzmmnts_locs_data(self.getzmmnts_locs_xprms(),color,width,type,xscale,yscale,xoffset)
    
  def plzmmnts_locs_yprms(self,color=black,width=1,type='solid',xscale=1.,yscale=1.,xoffset=0.):  
    self.plzmmnts_locs_data(self.getzmmnts_locs_yprms(),color,width,type,xscale,yscale,xoffset)
    
  def plzmmnts_locs_xpbar(self,color=black,width=1,type='solid',xscale=1.,yscale=1.,xoffset=0.):  
    self.plzmmnts_locs_data(self.getzmmnts_locs_xpbar(),color,width,type,xscale,yscale,xoffset)
    
  def plzmmnts_locs_ypbar(self,color=black,width=1,type='solid',xscale=1.,yscale=1.,xoffset=0.):  
    self.plzmmnts_locs_data(self.getzmmnts_locs_ypbar(),color,width,type,xscale,yscale,xoffset)
    
  def plzmmnts_locs_xxpbar(self,color=black,width=1,type='solid',xscale=1.,yscale=1.,xoffset=0.):  
    self.plzmmnts_locs_data(self.getzmmnts_locs_xxpbar(),color,width,type,xscale,yscale,xoffset)
    
  def plzmmnts_locs_yypbar(self,color=black,width=1,type='solid',xscale=1.,yscale=1.,xoffset=0.):  
    self.plzmmnts_locs_data(self.getzmmnts_locs_yypbar(),color,width,type,xscale,yscale,xoffset)
    
  def plzmmnts_locs_emitx(self,color=black,width=1,type='solid',xscale=1.,yscale=1.,xoffset=0.):  
    self.plzmmnts_locs_data(self.getzmmnts_locs_emitxrms(),color,width,type,xscale,yscale,xoffset)
    
  def plzmmnts_locs_emity(self,color=black,width=1,type='solid',xscale=1.,yscale=1.,xoffset=0.):  
    self.plzmmnts_locs_data(self.getzmmnts_locs_emityrms(),color,width,type,xscale,yscale,xoffset)
    
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
    
  def getex(self,**kw):
    return getex(jslist=self.jslist,**kw)
    
  def getey(self,**kw):
    return getey(jslist=self.jslist,**kw)
    
  def getez(self,**kw):
    return getez(jslist=self.jslist,**kw)
    
  def getbx(self,**kw):
    return getbx(jslist=self.jslist,**kw)
    
  def getby(self,**kw):
    return getby(jslist=self.jslist,**kw)
    
  def getbz(self,**kw):
    return getbz(jslist=self.jslist,**kw)
    
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

  def ppxex(self,**kw):
    return ppxex(jslist=self.jslist,**kw)

  def ppxey(self,**kw):
    return ppxey(jslist=self.jslist,**kw)

  def ppxez(self,**kw):
    return ppxez(jslist=self.jslist,**kw)

  def ppyex(self,**kw):
    return ppyex(jslist=self.jslist,**kw)

  def ppyey(self,**kw):
    return ppyey(jslist=self.jslist,**kw)

  def ppyez(self,**kw):
    return ppyez(jslist=self.jslist,**kw)

  def ppzex(self,**kw):
    return ppzex(jslist=self.jslist,**kw)

  def ppzey(self,**kw):
    return ppzey(jslist=self.jslist,**kw)

  def ppzez(self,**kw):
    return ppzez(jslist=self.jslist,**kw)

  def ppxbx(self,**kw):
    return ppxbx(jslist=self.jslist,**kw)

  def ppxby(self,**kw):
    return ppxby(jslist=self.jslist,**kw)

  def ppxbz(self,**kw):
    return ppxbz(jslist=self.jslist,**kw)

  def ppybx(self,**kw):
    return ppybx(jslist=self.jslist,**kw)

  def ppyby(self,**kw):
    return ppyby(jslist=self.jslist,**kw)

  def ppybz(self,**kw):
    return ppybz(jslist=self.jslist,**kw)

  def ppzbx(self,**kw):
    return ppzbx(jslist=self.jslist,**kw)

  def ppzby(self,**kw):
    return ppzby(jslist=self.jslist,**kw)

  def ppzbz(self,**kw):
    return ppzbz(jslist=self.jslist,**kw)
