"""Ionization: class for generating particles from impact ionization and other inelestic collisions.
"""
__all__ = ['Ionization']
from warp import *
import time
try:
  from txphysics import txionpack
  l_txphysics = True
except:
  l_txphysics = False

def ionizationdoc():
  import Ionization
  print Ionization.__doc__

class Ionization:
  """
Class for generating particles from impact ionization
 - incident_species:   incident species
 - target_species:   target species
 - emitted_species:   created species
 - cross_section: cross sections
 - l_verbose=False: When True, prints information about what collision occur
 - stride=100: The stride in the particle loops, controlling how many particles
               can under go a collision each time step.
 - nx,ny,nz: Dimensions of the grid used when calculating the density of a
             target species. Defaults to the values from w3d.
 - xmin,xmax,ymin,ymax,zmin,zmax: Extent of the grid used when calculating the
                                  density of a target species. Defaults to the
                                  values from w3d.
 - l_timing=False: Turns on timing of the ionization events
  """
  def __init__(self,l_verbose=False,stride=100,nx=None,ny=None,nz=None,xmin=None,xmax=None,ymin=None,ymax=None,zmin=None,zmax=None,l_timing=False):
    self.nx = int(where(nx is None,w3d.nx,nx))
    self.ny = int(where(ny is None,w3d.ny,ny))
    self.nz = int(where(nz is None,w3d.nzlocal,nz))
    self.xmin = float(where(xmin is None,w3d.xmmin,xmin))
    self.xmax = float(where(xmax is None,w3d.xmmax,xmax))
    self.zmin = float(where(zmin is None,w3d.zmminlocal,zmin))
    self.zmax = float(where(zmax is None,w3d.zmmaxlocal,zmax))
    self.dx=(self.xmax-self.xmin)/self.nx
    self.dz=(self.zmax-self.zmin)/self.nz
    if w3d.solvergeom == w3d.RZgeom:
      self.ymin=self.xmin
      self.ymax=self.xmax
      self.dy=self.dx
      self.ny=self.nx
    else:
      self.ymin = float(where(ymin is None,w3d.ymmin,ymin))
      self.ymax = float(where(ymax is None,w3d.ymmax,ymax))
      self.dy=(self.ymax-self.ymin)/self.ny
    self.ndensc=fzeros((self.nx+1,self.ny+1,self.nz+1),'d')
    self.invvol=1./(self.dx*self.dy*self.dz)
    self.l_verbose=l_verbose
    self.stride=stride
    self.inter={}
    self.target_dens={}
    self.npmax=4096
    self.nps={}
    self.x={}
    self.y={}
    self.z={}
    self.ux={}
    self.uy={}
    self.uz={}
    self.gi={}
    self.w={}
    self.pidtag={}
    # --- This is kind of messy, but the injpid must be handled. When injection
    # --- is being done, it needs to be set properly in order for the emitted
    # --- particles to have the correctly calculated E fields when they are created
    # --- near an emitting surface. The injpid of the incident particles are passed to
    # --- the emitted particles.
    self.injpid={}
    self.emitted_id = None
    self.l_timing=l_timing
    self.install()

  def setupcross_section(self,incident_species,emitted_species,cross_section,target_species):
    """If the cross section was specified, this method does nothing. Otherwise
it setups up one of several cross section functions, depending on which
interactions are supported by the libraries.
    """
    if cross_section is not None:
      if isinstance(cross_section,ArrayType):
        # --- Return an instance of a tablized cross section
        return TableCrossSection(incident_species,cross_section)
      else:
        return cross_section

    # --- Check if collision is supported in the Aladdin database
    try:
      cross_section = AladdinCrossSection([incident_species,target_species],emitted_species)
      return cross_section
    except AladdinCrossSectionNotSupported:
      # --- Collision not found in the database
      pass

    if l_txphysics:
      # --- Check if txphysics supports the collision
      try:
        cross_section = TXIonPackCrossSection(incident_species,target_species)
        return cross_section
      except TXIonPackCrossSectionNotSupport:
        # --- Collision is not supported
        pass

    raise CrossSectionNotSupport('Collision not found in the database - it must be supplied')

  def getcross_section(self,cross_section,vi):
    """The cross section can either be a scalar or a functional that returns
a value given a velocity. This returns the Appropriate value.
    """
    try:
      result = cross_section(vi)
    except TypeError:
      result = cross_section
    return result

  def add(self,incident_species,emitted_species,cross_section=None,
          target_species=None,ndens=None,target_fluidvel=None,
          emitted_energy0=None,emitted_energy_sigma=None,
          incident_pgroup=top.pgroup,target_pgroup=top.pgroup,emitted_pgroup=top.pgroup,
          l_remove_incident=None,l_remove_target=None,emitted_tag=None):
    if incident_species not in self.inter:
        self.inter[incident_species]={}
        for key in ['target_species','emitted_species','cross_section','ndens','target_fluidvel',
                    'remove_incident','remove_target',
                    'emitted_energy0','emitted_energy_sigma','emitted_tag',
                    'incident_pgroup','target_pgroup','emitted_pgroup']:
          self.inter[incident_species][key]=[]
    if not iterable(emitted_species): emitted_species=[emitted_species]
    # --- This if block was removed. It is unclear what the reason was for the block.
    # --- It precluded the charge exchange interaction, since in that case the
    # --- incident and emitted species are the same, and the if block removed
    # --- emission of the emitted species.
    #if emitted_species[0] in [incident_species,target_species]:
    #  self.inter[incident_species]['emitted_species'] +=[emitted_species[1:]]
    #else:
    #  self.inter[incident_species]['emitted_species'] +=[emitted_species]

    cross_section = self.setupcross_section(incident_species,emitted_species,cross_section,target_species)

    # --- Only include Species instances in the emitted species. It may have contained
    # --- Particle types that were needed for setting up the cross section calculation.
    e_species = []
    for es in emitted_species:
      if isinstance(es,Species): e_species.append(es)

    # --- Same for target_species - only include it if its an instance of species.
    if isinstance(target_species,Particle): target_species = None

    self.inter[incident_species]['target_species']  +=[target_species]
    self.inter[incident_species]['emitted_species'] +=[e_species]
    self.inter[incident_species]['cross_section']   +=[cross_section]
    self.inter[incident_species]['ndens']           +=[ndens]
    self.inter[incident_species]['target_fluidvel'] +=[target_fluidvel]
    if l_remove_incident is None:
      if incident_species.type is e_species[0].type:
        self.inter[incident_species]['remove_incident']+=[1]
      else:
        self.inter[incident_species]['remove_incident']+=[0]
    else:
      self.inter[incident_species]['remove_incident']+=[l_remove_incident]
    if l_remove_target is None and target_species is not None:
      if target_species.type is e_species[0].type:
        self.inter[incident_species]['remove_target']+=[1]
      else:
        self.inter[incident_species]['remove_target']+=[0]
    else:
      self.inter[incident_species]['remove_target']+=[l_remove_target]
    if emitted_energy0 is None and not self.inter[incident_species]['remove_incident'][-1]:
      # --- If the incident species is not being removed, then the emitted
      # --- particles are drawn from a random distribution. If not specified,
      # --- the default energy of the emitted particles is zero.
      emitted_energy0 = 0.
      emitted_energy_sigma = 0.
    try:
      len(emitted_energy0)
    except TypeError:
      # --- Make sure the emitted_energy0 has the same length as e_species
      emitted_energy0 = len(e_species)*[emitted_energy0]
    try:
      len(emitted_energy_sigma)
    except TypeError:
      # --- Make sure the emitted_energy_sigma has the same length as e_species
      emitted_energy_sigma = len(e_species)*[emitted_energy_sigma]
    self.inter[incident_species]['emitted_energy0']   +=[emitted_energy0]
    self.inter[incident_species]['emitted_energy_sigma']   +=[emitted_energy_sigma]
    self.inter[incident_species]['emitted_tag']   +=[emitted_tag]
    if emitted_tag is not None and self.emitted_id is None:
      self.emitted_id = nextpid()
    self.inter[incident_species]['incident_pgroup']=incident_pgroup
    self.inter[incident_species]['target_pgroup']  =target_pgroup
    self.inter[incident_species]['emitted_pgroup'] =emitted_pgroup
    if target_species is not None:
      if target_species not in self.target_dens:
        self.target_dens[target_species]={}
        for key in ['ndens','ndens_updated']:
          self.target_dens[target_species][key]=[]
        self.target_dens[target_species]['ndens']           =fzeros((self.nx+1,self.ny+1,self.nz+1),'d')
        self.target_dens[target_species]['target_fluidvel'] =fzeros((self.nx+1,self.ny+1,self.nz+1,3),'d')
        self.target_dens[target_species]['ndens_updated']   =0

    for e in e_species:
      js=e.jslist[0]
      if emitted_pgroup not in self.x:
        self.nps[emitted_pgroup]={}
        self.x[emitted_pgroup]={}
        self.y[emitted_pgroup]={}
        self.z[emitted_pgroup]={}
        self.ux[emitted_pgroup]={}
        self.uy[emitted_pgroup]={}
        self.uz[emitted_pgroup]={}
        self.gi[emitted_pgroup]={}
        self.w[emitted_pgroup]={}
        self.pidtag[emitted_pgroup]={}
        self.injpid[emitted_pgroup]={}
      if js not in self.x[emitted_pgroup]:
        self.nps[emitted_pgroup][js]=0
        self.x[emitted_pgroup][js]=fzeros(self.npmax,'d')
        self.y[emitted_pgroup][js]=fzeros(self.npmax,'d')
        self.z[emitted_pgroup][js]=fzeros(self.npmax,'d')
        self.ux[emitted_pgroup][js]=fzeros(self.npmax,'d')
        self.uy[emitted_pgroup][js]=fzeros(self.npmax,'d')
        self.uz[emitted_pgroup][js]=fzeros(self.npmax,'d')
        self.gi[emitted_pgroup][js]=fzeros(self.npmax,'d')
        if top.wpid > 0:
          self.w[emitted_pgroup][js]=fzeros(self.npmax,'d')
        if emitted_tag is not None:
          self.pidtag[emitted_pgroup][js]=fzeros(self.npmax,'d')
        if top.injpid > 0:
          self.injpid[emitted_pgroup][js]=fzeros(self.npmax,'d')

  def add_ionization(self,incident_species,emitted_species,target_species=None,cross_section=None,**kw):
    """The incident species will knock an electron off of the target species.
If the cross section is not given, it will be obtained from the databse.
The reduced charge state of the target must be included in the emitted species
for the proper cross section to be determined.
    """
    # --- Generate a temporary list of emitted species, making sure that it includes all of the
    # --- appropriate species so that the proper cross section can be found. This includes
    # --- the incident species, the reduced charge species and an electron. The user must have
    # --- provided the reduced charge species.
    e_species = copy.copy(emitted_species)
    if not iterable(e_species): e_species = [e_species]
    l_need_electron = True
    for es in e_species:
        if es is Electron: l_need_electron = False
        if isinstance(es,Species) and es.type is Electron: l_need_electron = False
    if l_need_electron: e_species.append(Electron)
    if incident_species not in e_species: e_species.append(incident_species)

    cross_section = self.setupcross_section(incident_species,e_species,cross_section,target_species)

    self.add(incident_species,emitted_species,cross_section,target_species,
             l_remove_target=False,**kw)

  def add_chargeexchange(self,incident_species,target_species=None,cross_section=None,**kw):
    """An electron will be exchanged between the incident and target species. A particle of the incident species will be emitted at the energy emitted_energy_sigma. The target species should be a particle type, or if it is not given, the cross section must be.
    """
    cross_section = self.setupcross_section(incident_species,[target_species,incident_species],
                                            cross_section,target_species)
    self.add(incident_species,incident_species,cross_section,target_species,
             l_remove_incident=True,**kw)

  def add_detachment(self,incident_species,target_species=None,emitted_species=None,cross_section=None,**kw):
    """An electron will be detached from the incident species when it
collides with the target species, resulting in a reduced charge state
particle and an electron.  If the cross section is not given, it will be
obtained from the database. The reduced charge state species or its particle
type must be listed in the emitted species so that the proper cross section
can be found. If specified, the emitted species will be emitted with the
velocity of the incident particle.
    """
    # --- Generate a temporary list of emitted species, making sure that it includes all of the
    # --- appropriate species so that the proper cross section can be found. This includes
    # --- the target species, the reduced charge species and an electron. The user must have
    # --- provided the reduced charge species.
    e_species = copy.copy(emitted_species)
    if not iterable(e_species): e_species = [e_species]
    l_need_electron = True
    for es in e_species:
        if es is Electron: l_need_electron = False
        if isinstance(es,Species) and es.type is Electron: l_need_electron = False
    if l_need_electron: e_species.append(Electron)
    if target_species not in e_species: e_species.append(target_species)

    cross_section = self.setupcross_section(incident_species,e_species,cross_section,target_species)

    self.add(incident_species,emitted_species,cross_section,target_species,
             l_remove_incident=True,**kw)

  add_stripping = add_detachment

  def install(self):
    if not isinstalleduserinjection(self.generate):
      installuserinjection(self.generate)

  def addpart(self,nn,x,y,z,ux,uy,uz,gi,pg,js,tag,injpid,w=1.):
    ilf=0
    if injpid is not None:
      # --- This is needed in case injection is setup after the interactions
      # --- are setup. If that happens, the array would not have been allocated.
      try:
        self.injpid[pg][js]
      except KeyError:
        self.injpid[pg][js]=fzeros(self.npmax,'d')
    while self.nps[pg][js]+nn>self.npmax:
      il=self.nps[pg][js]
      iu=min(il+nn,self.npmax)
      nf=iu-il
      self.x [pg][js][il:iu]= x[ilf:ilf+nf]
      self.y [pg][js][il:iu]= y[ilf:ilf+nf]
      self.z [pg][js][il:iu]= z[ilf:ilf+nf]
      self.ux[pg][js][il:iu]=ux[ilf:ilf+nf]
      self.uy[pg][js][il:iu]=uy[ilf:ilf+nf]
      self.uz[pg][js][il:iu]=uz[ilf:ilf+nf]
      self.gi[pg][js][il:iu]=gi[ilf:ilf+nf]
      if top.wpid>0:self.w[pg][js][il:iu]=w[ilf:ilf+nf]
      if tag is not None:
        self.pidtag[pg][js][il:iu]=tag
      if injpid is not None:
        self.injpid[pg][js][il:iu]=injpid[ilf:ilf+nf]
      self.nps[pg][js]+=nf
      self.flushpart(pg,js)
      ilf+=nf
      nn-=nf
    il=self.nps[pg][js]
    iu=il+nn
    self.x [pg][js][il:iu]=x [ilf:]
    self.y [pg][js][il:iu]=y [ilf:]
    self.z [pg][js][il:iu]=z [ilf:]
    self.ux[pg][js][il:iu]=ux[ilf:]
    self.uy[pg][js][il:iu]=uy[ilf:]
    self.uz[pg][js][il:iu]=uz[ilf:]
    self.gi[pg][js][il:iu]=gi[ilf:]
    if top.wpid>0:self.w[pg][js][il:iu]=w[ilf:]
    if tag is not None:
      self.pidtag[pg][js][il:iu]=tag
    if injpid is not None:
      self.injpid[pg][js][il:iu]=injpid[ilf:]
    self.nps[pg][js]+=nn

  def flushpart(self,pg,js):
      if self.nps[pg][js]>0:
         nn=self.nps[pg][js]
#         condition = (self.x[js][:nn]>w3d.xmmin) & (self.x[js][:nn]<w3d.xmmax) & \
#                     (self.y[js][:nn]>w3d.ymmin) & (self.y[js][:nn]<w3d.ymmax) & \
#                     (self.z[js][:nn]>w3d.zmminlocal) & (self.z[js][:nn]<w3d.zmmaxlocal)
#         if sum(condition)<nn:
#           ic = compress(condition==0,arange(nn))
#           print 'ioniz: out of bound: ',js,self.x[js][ic],self.y[js][ic],self.z[js][ic]
#           f=PW.PW('pos.pdb')
#           f.x=self.x[js][:nn]
#           f.y=self.y[js][:nn]
#           f.z=self.z[js][:nn]
#           f.close()
#           raise Exception('')
#         window(5);ppg(self.y[pg][js][:nn],self.x[pg][js][:nn]);limits(w3d.xmmin,w3d.xmmax,w3d.ymmin,w3d.ymmax);refresh()
         pidpairs = []
         try:
           pidpairs.append([self.emitted_id,self.pidtag[pg][js][:nn]])
         except KeyError:
           pass
         if top.injpid > 0:
           pidpairs.append([top.injpid,self.injpid[pg][js][:nn]])
         if top.wpid==0:
           w=1.
         else:
           w=self.w[pg][js][:nn]
         addparticles(x=self.x[pg][js][:nn],
                      y=self.y[pg][js][:nn],
                      z=self.z[pg][js][:nn],
                      vx=self.ux[pg][js][:nn],
                      vy=self.uy[pg][js][:nn],
                      vz=self.uz[pg][js][:nn],
                      gi=self.gi[pg][js][:nn],
                      w=w,
                      pidpairs=pidpairs,
                      lmomentum=True,
                      pgroup=pg,
                      lallindomain=True,
                      js=js)
         self.nps[pg][js]=0

  def printall(self,l_cgm=0):
    swidth=0
    twidth=0
    ewidth={}
    title='        *** Particle/particle interactions ***\n'
    textblock=''
    for incident_species in self.inter:
      swidth=max(swidth,len(incident_species.name))
      for it,target_species in enumerate(self.inter[incident_species]['target_species']):
        if target_species is None:
          tname='background gas'
        else:
          tname=target_species.name
        twidth=max(twidth,len(tname))
        firstelement = []
        #if not self.inter[incident_species]['remove_incident'][it]:
        #  firstelement = [incident_species]
        #if not self.inter[incident_species]['remove_target'][it]:
        #  firstelement = [target_species]
        for ie,emitted_species in enumerate(firstelement+self.inter[incident_species]['emitted_species'][it]):
          ename=emitted_species.name
          if ie not in ewidth:
            ewidth[ie]=len(ename)
          else:
            ewidth[ie]=max(ewidth[ie],len(ename))
    fs='%%-%gs'%swidth
    ft='%%-%gs'%twidth
    fe={}
    for ie in ewidth:
      fe[ie]='%%-%gs'%ewidth[ie]
    for incident_species in self.inter:
      sname=fs%incident_species.name[:swidth]
      textblock+='\n'
      for it,target_species in enumerate(self.inter[incident_species]['target_species']):
        if target_species is None:
          tname=ft%'background gas'[:twidth]
        else:
          tname=ft%target_species.name[:twidth]
        if callable(self.inter[incident_species]['cross_section'][it]):
          try:
            cname = self.inter[incident_species]['cross_section'][it].name
          except AttributeError:
            cname = ''
        else:
          cname='(CS=%g)'%(self.inter[incident_species]['cross_section'][it])
        firstelement = []
        #if not self.inter[incident_species]['remove_incident'][it]:
        #  firstelement = [incident_species]
        #if not self.inter[incident_species]['remove_target'][it]:
        #  firstelement = [target_species]
        for ie,emitted_species in enumerate(firstelement+self.inter[incident_species]['emitted_species'][it]):
          if ie==0:
            ename=fe[ie]%emitted_species.name[:ewidth[ie]]
          else:
            ename+=' + '+fe[ie]%emitted_species.name[:ewidth[ie]]
        textblock += sname+' + '+tname+' => '+ename+' '+cname+'\n'
    if l_cgm:
      plt(title,0.20,0.905,justify="LT",height=14)
      plt(textblock,0.13,0.88,justify="LT",height=9,font='courier')
      fma()
    else:
      print title
      print textblock


  def cleanuppgroups(self):
    """Resets all of the pgroups to top.pgroup"""
    for inter in self.inter.values():
      inter['incident_pgroup']= top.pgroup
      inter['target_pgroup']  = top.pgroup
      inter['emitted_pgroup'] = top.pgroup

    def _cleanit(x):
      #if len(x) > 1: return
      # --- This assumes that there is only one key
      kk = x.keys()
      for pg in kk:
        x[top.pgroup] = x[pg]
        del x[pg]

    _cleanit(self.nps)
    _cleanit(self.x)
    _cleanit(self.y)
    _cleanit(self.z)
    _cleanit(self.ux)
    _cleanit(self.uy)
    _cleanit(self.uz)
    _cleanit(self.gi)
    _cleanit(self.w)
    _cleanit(self.pidtag)
    _cleanit(self.injpid)

#printall(io,l_cgm=1)

  def generate(self,dt=None):
    if dt is None:dt=top.dt
    if self.l_timing:t1 = time.clock()
    for target_species in self.target_dens:
      self.target_dens[target_species]['ndens_updated']=0
    for incident_species in self.inter:
      npinc = 0
      ispushed=0
      ipg=self.inter[incident_species]['incident_pgroup']
      tpg=self.inter[incident_species]['target_pgroup']
      epg=self.inter[incident_species]['emitted_pgroup']
      for js in incident_species.jslist:
        npinc+=ipg.nps[js]
        if ipg.ldts[js]:ispushed=1
      if npinc==0 or not ispushed:continue
      for it,target_species in enumerate(self.inter[incident_species]['target_species']):
        ndens=self.inter[incident_species]['ndens'][it]
        target_fluidvel=self.inter[incident_species]['target_fluidvel'][it]
        if ndens is not None:
          continue
        else:
          if self.target_dens[target_species]['ndens_updated']:
            continue
          else:
            self.target_dens[target_species]['ndens_updated']=1
          ndens = self.target_dens[target_species]['ndens']
          target_fluidvel=self.target_dens[target_species]['target_fluidvel']
          nptarget=0
          for jstarget in target_species.jslist:
            nptarget+=tpg.nps[jstarget]
          if nptarget==0:continue
          self.ndensc[...]=0.
          ndens[...]=0.
          for jstarget in target_species.jslist:
            i1 = tpg.ins[jstarget] - 1
            i2 = tpg.ins[jstarget] + tpg.nps[jstarget] - 1
            xt=tpg.xp[i1:i2]
            yt=tpg.yp[i1:i2]
            zt=tpg.zp[i1:i2]
            git=tpg.gaminv[i1:i2]
            vxt=tpg.uxp[i1:i2]*git
            vyt=tpg.uyp[i1:i2]*git
            vzt=tpg.uzp[i1:i2]*git
            fact=1.
            if w3d.l4symtry:
              xt=abs(xt)
              fact=0.25
            elif w3d.l2symtry:
              fact=0.5
            if w3d.l2symtry or w3d.l4symtry:yt=abs(yt)
            if top.wpid==0:
              weights=ones(tpg.nps[jstarget],'d')
            else:
              weights=tpg.pid[i1:i2,top.wpid-1]
            # --- deposit density
            deposgrid3d(1,tpg.nps[jstarget],xt,yt,zt,
                        tpg.sw[jstarget]*self.invvol*fact*weights,
                        self.nx,self.ny,self.nz,ndens,self.ndensc,
                        self.xmin,self.xmax,self.ymin,self.ymax,
                        self.zmin,self.zmax)
            # --- computes target fluid velocity
            deposgrid3dvect(0,tpg.nps[jstarget],xt,yt,zt,vxt,vyt,vzt,
                        tpg.sw[jstarget]*self.invvol*fact*weights,
                        self.nx,self.ny,self.nz,target_fluidvel,self.ndensc,
                        self.xmin,self.xmax,self.ymin,self.ymax,
                        self.zmin,self.zmax)

#          if w3d.l2symtry or w3d.l4symtry:self.ndens[:,0,:]*=2.
#          if w3d.l4symtry:self.ndens[0,:,:]*=2.

    for incident_species in self.inter:
      npinc = 0
      ispushed=0
      ipg=self.inter[incident_species]['incident_pgroup']
      tpg=self.inter[incident_species]['target_pgroup']
      epg=self.inter[incident_species]['emitted_pgroup']
      for js in incident_species.jslist:
        npinc+=ipg.nps[js]
        if ipg.ldts[js]:ispushed=1
      if npinc==0 or not ispushed:continue
      for it,target_species in enumerate(self.inter[incident_species]['target_species']):
        ndens=self.inter[incident_species]['ndens'][it]
        target_fluidvel=self.inter[incident_species]['target_fluidvel'][it]
#        ndens=1.e23
        for js in incident_species.jslist:
          i1 = ipg.ins[js] - 1 + top.it%self.stride
          i2 = ipg.ins[js] + ipg.nps[js] - 1
          xi=ipg.xp[i1:i2:self.stride]#.copy()
          yi=ipg.yp[i1:i2:self.stride]#.copy()
          zi=ipg.zp[i1:i2:self.stride]#.copy()
          ni = shape(xi)[0]
          gaminvi=ipg.gaminv[i1:i2:self.stride]#.copy()
          uxi=ipg.uxp[i1:i2:self.stride]#.copy()
          uyi=ipg.uyp[i1:i2:self.stride]#.copy()
          uzi=ipg.uzp[i1:i2:self.stride]#.copy()
          if top.wpid > 0:
            # --- Save the wpid of the incident particles so that it can be
            # --- passed to the emitted particles.
            wi = ipg.pid[i1:i2:self.stride,top.wpid-1]
          else:
            wi = 1.
          if top.injpid > 0:
            # --- Save the injpid of the incident particles so that it can be
            # --- passed to the emitted particles.
            injpid = ipg.pid[i1:i2:self.stride,top.injpid-1]
          else:
            injpid = None
          # --- get velocity in lab frame if using a boosted frame of reference
          if top.boost_gamma>1.:
            uzboost = clight*sqrt(top.boost_gamma**2-1.)
            setu_in_uzboosted_frame3d(ni,uxi,uyi,uzi,gaminvi,
                                      -uzboost,
                                      top.boost_gamma)
          vxi=uxi*gaminvi
          vyi=uyi*gaminvi
          vzi=uzi*gaminvi
          # --- get local target density
          if ndens is None:
            ndens = self.target_dens[target_species]['ndens']
          if type(ndens) in [type(0.),type(0)]:
            dp=ones(ni,'d')*ndens
            if target_fluidvel is None:
              xmin=self.xmin
              xmax=self.xmax
              ymin=self.ymin
              ymax=self.ymax
              zmin=self.zmin
              zmax=self.zmax
            else:
              vxtf = target_fluidvel[0]
              vytf = target_fluidvel[1]
              vztf = target_fluidvel[2]
              xmin=self.xmin+vxtf*top.time
              xmax=self.xmax+vxtf*top.time
              ymin=self.ymin+vytf*top.time
              ymax=self.ymax+vytf*top.time
              zmin=self.zmin+vztf*top.time
              zmax=self.zmax+vztf*top.time
            if w3d.solvergeom==w3d.RZgeom:
              ri=sqrt(xi*xi+yi*yi)
              dp=where((ri>=xmin) & (ri<=xmax) & \
                       (zi>=zmin) & (zi<=zmax),dp,0.)
            else:
              dp=where((xi>=xmin) & (xi<=xmax) & \
                       (yi>=ymin) & (yi<=ymax) & \
                       (zi>=zmin) & (zi<=zmax),dp,0.)
          else:
            dp=zeros(ni,'d')
            getgrid3d(ni,xi,yi,zi,dp,
                        self.nx,self.ny,self.nz,ndens, \
                        self.xmin,self.xmax,self.ymin,self.ymax,self.zmin,self.zmax,
                        w3d.l2symtry,w3d.l4symtry)
          # --- get local target fluid velocity
          if target_fluidvel is None:
            if target_species is None:
              target_fluidvel = [0.,0.,0.]
            else:
              target_fluidvel = self.target_dens[target_species]['target_fluidvel']
          if type(target_fluidvel) in [type([])]:
            vxtf = target_fluidvel[0]
            vytf = target_fluidvel[1]
            vztf = target_fluidvel[2]
          else:
            vxtf=zeros(ni,'d')
            vytf=zeros(ni,'d')
            vztf=zeros(ni,'d')
            getgrid3d(ni,xi,yi,zi,vxtf,
                        self.nx,self.ny,self.nz,target_fluidvel[...,0], \
                        self.xmin,self.xmax,self.ymin,self.ymax,self.zmin,self.zmax,
                        w3d.l2symtry,w3d.l4symtry)
            getgrid3d(ni,xi,yi,zi,vytf,
                        self.nx,self.ny,self.nz,target_fluidvel[...,1], \
                        self.xmin,self.xmax,self.ymin,self.ymax,self.zmin,self.zmax,
                        w3d.l2symtry,w3d.l4symtry)
            getgrid3d(ni,xi,yi,zi,vztf,
                        self.nx,self.ny,self.nz,target_fluidvel[...,2], \
                        self.xmin,self.xmax,self.ymin,self.ymax,self.zmin,self.zmax,
                        w3d.l2symtry,w3d.l4symtry)

          # compute the relative velocity
          # NOTE that at this point, the target species is assumed to have a negligible velocity.
          # this needs to be modified if this approximation is not valid.
          vxr = vxi-vxtf
          vyr = vyi-vytf
          vzr = vzi-vztf
          vi=sqrt(vxr*vxr+vyr*vyr+vzr*vzr)
          cross_section = self.getcross_section(self.inter[incident_species]['cross_section'][it],vi)

          # probability
          ncol = dp*cross_section*vi*dt*ipg.ndts[js]*self.stride
          if top.boost_gamma>1.:ncol*=top.gammabar_lab/top.gammabar

          # --- If the incident species is being removed, then only one collision event can happen.
          # --- Otherwise, a single particle may collide multiple times in a time step if ncol > 1.
          if self.inter[incident_species]['remove_incident'][it]:
            # --- Note that ncol is being set to slightly less than one. In the code below, adding ranf
            # --- to it will guarantee that it will end up as one. (Is this needed?)
            ncol = where(ncol>=1.,1.-1.e-10,ncol)

          # --- Get a count of the number of collisions for each particle. A random number is added to
          # --- ncol so that a fractional value has chance to result in a collision.
          ncoli=int(ncol + ranf(ncol))

          # --- Select the particles that will collide
          io=compress(ncoli>0,arange(ni))
          nnew = len(io)

          if None in self.inter[incident_species]['emitted_energy0'][it]:
            # --- When emitted_energy0 is not specified, use the velocity of
            # --- the incident particles for the emitted particles.
            uxnewsave = uxi
            uynewsave = uyi
            uznewsave = uzi

          if self.inter[incident_species]['remove_incident'][it]:
            # --- if projectile is modified, then need to delete it
            put(ipg.gaminv,array(io)*self.stride+i1,0.)

          # --- The position of the incident particle is at or near the incident particle
          xnew = xi
          ynew = yi
          znew = zi

          # --- Loop until there are no more collision events that need handling
          while(nnew>0):

            # --- The emitted particles positions, in some cases, are slightly
            # --- offset from the incident
            xnewp = xnew[io]
            ynewp = ynew[io]
            znewp = znew[io]
            xnew = xnewp+(ranf(xnewp)-0.5)*1.e-10*self.dx
            ynew = ynewp+(ranf(ynewp)-0.5)*1.e-10*self.dy
            znew = znewp+(ranf(znewp)-0.5)*1.e-10*self.dz
            if top.wpid == 0: 
              w = 1.
            else:
              w = wi[io]

            # --- The injpid value needs to be copied to the emitted particles
            # --- so that they are handled properly in the region near the source.
            if top.injpid > 0: injpid = injpid[io]

            # --- If the emitted energy was not specified, the emitted particle will be
            # --- given the same velocity of the incident particle.
            if None in self.inter[incident_species]['emitted_energy0'][it]:
              uxnewsave = uxnewsave[io]
              uynewsave = uynewsave[io]
              uznewsave = uznewsave[io]

            for ie,emitted_species in enumerate(self.inter[incident_species]['emitted_species'][it]):

              if self.inter[incident_species]['emitted_energy0'][it][ie] is not None:
                # --- Create new velocities for the emitted particles.
                ek0ionel = self.inter[incident_species]['emitted_energy0'][it][ie]
                esigionel = self.inter[incident_species]['emitted_energy_sigma'][it][ie]
                if esigionel==0.:
                  ek = zeros(nnew)
                else:
                  ek = SpRandom(0.,esigionel,nnew) #kinetic energy
                ek=abs(ek+ek0ionel) #kinetic energy
                fact = jperev/(emass*clight**2)
                gamma=ek*fact+1.
                u=clight*sqrt(ek*fact*(gamma+1.))
                # velocity direction: random in (x-y) plane plus small longitudinal component:
                phi=2.*pi*ranf(u)
                vx=cos(phi); vy=sin(phi); vz=0.01*ranf(u)
                # convert into a unit vector:
                vu=sqrt(vx**2+vy**2+vz**2)
                # renormalize:
                vx/=vu; vy/=vu; vz/=vu
                # find components of v*gamma:
                uxnew=u*vx
                uynew=u*vy
                uznew=u*vz
              else:
                uxnew = uxnewsave
                uynew = uynewsave
                uznew = uznewsave

              ginew = 1./sqrt(1.+(uxnew**2+uynew**2+uznew**2)/clight**2)
              # --- get velocity in boosted frame if using a boosted frame of reference
              if top.boost_gamma>1.:
                setu_in_uzboosted_frame3d(shape(ginew)[0],uxnew,uynew,uznew,ginew,
                                          uzboost,
                                          top.boost_gamma)

              if self.l_verbose:print 'add ',nnew, emitted_species.name,' from by impact ionization:',incident_species.name,'+',((target_species is None and 'background gas') or target_species.name)
              if self.inter[incident_species]['remove_incident'][it] and (emitted_species.type is incident_species.type):
                self.addpart(nnew,xnewp,ynewp,znewp,uxnew,uynew,uznew,ginew,epg,emitted_species.jslist[0],
                             self.inter[incident_species]['emitted_tag'][it],injpid,w)
              else:
                self.addpart(nnew,xnew,ynew,znew,uxnew,uynew,uznew,ginew,epg,emitted_species.jslist[0],
                             self.inter[incident_species]['emitted_tag'][it],injpid,w)
            ncoli = ncoli[io] - 1
            io = arange(nnew)[ncoli>0]
            nnew = len(io)

    # make sure that all particles are added and cleared
    for pg in self.x:
      for js in self.x[pg]:
        self.flushpart(pg,js)
 
    for incident_species in self.inter:
      for js in incident_species.jslist:
        processlostpart(top.pgroup,js+1,top.clearlostpart,top.time,top.zbeam)

    if self.l_timing:print 'time ionization = ',time.clock()-t1,'s'

class CrossSectionNotSupport(Exception):
    pass

class TXIonPackCrossSectionNotSupport(Exception):
    pass

class TXIonPackCrossSection(object):
    """
For supported collisions, retrieves the cross section from the txionpack library.
    """
    # Initialize the kinds of gases
    #
    # 0 is H (protons only)
    # 1 is H2
    # 2 is He
    # 3 is CO2
    # 4 is CO
    # 5 is O2
    # 6 is N2
    # 7-8 for electrons only
    # 7 is Ar
    # 8 is Ne
    txphysics_targets={Hydrogen:0,
                       Dihydrogen:1,
                       Helium:2,
                       Carbon_Dioxide:3,
                       Carbon_Monoxide:4,
                       Dioxygen:5,
                       Dinitrogen:6,
                       Argon:7,
                       Neon:8}

    def __init__(self,incident_species,target_species):
        self.incident_species = incident_species
        self.target_species = target_species
        if isinstance(incident_species,Species):
            itype = incident_species.type
        elif isinstance(incident_species,Particle):
            itype = incident_species
        if isinstance(target_species,Species):
            ttype = target_species.type
        elif isinstance(target_species,Particle):
            ttype = target_species
        else:
            raise TXIonPackCrossSectionNotSupport('target_species must be given')
        try:
            self.gastype = TXIonPackCrossSection.txphysics_targets[ttype]
        except KeyError:
            raise TXIonPackCrossSectionNotSupport('target species type not supported in txphysics package')
        # This is an integer flag to specify electron (0) or ion (1)
        if incident_species.type in (Electron,Positron):
            self.incident = 0
        else:
            self.incident = 1
        self.name = 'TxPhysics ' + itype.Symbol + ' ' + ttype.Symbol

    def __call__(self,vi):
        return txionpack.get_sigma_impact_array(vi, self.gastype, self.incident)

class AladdinCrossSectionNotSupported(Exception):
    pass

class AladdinCrossSection(object):
    """
This accesses the Aladdin cross section database to find an approrpriate
cross section for the given collision event. Currently, it only supports
events in section 8 (Hyrogen and Helium collions), and only those that use a
Chebyshev fitting.
    """
    def __init__(self,incident_species,emitted_species):
        for i_species in incident_species:
            if (i_species is Electron or
                (isinstance(i_species,Species) and i_species.type is Electron)):
                raise AladdinCrossSectionNotSupported('incident electron not supported in the database')
        self.incident_species = incident_species
        self.emitted_species = emitted_species
        # --- This is not needed, but save it for later maybe.
        #self.header = re.compile('\$ (?P<itype>[A-Z]+)  (?P<inc1>[a-zA-z]+)(?P<inc1n>\{[0-9]\})* (?P<inc1chg>\[[+-][0-9]\])  (?P<inc2>[a-zA-z]+)(?P<inc2n>\{[0-9]\})* (?P<inc2chg>\[[+-][0-9]\])   (?P<emit1>[a-zA-z]+)(?P<emit1n>\{[0-9]\})* (?P<emit1chg>\[[+-][0-9]\])  (?P<emit2>[a-zA-z]+)(?P<emit2n>\{[0-9]\})* (?P<emit2chg>\[[+-][0-9]\])(?P<emite>  [0-9]e)*')
        self.data = self.findcollision()
        if len(self.data) == 0:
            raise AladdinCrossSectionNotSupported('collision not found in database')
        self.Emin = self.data[-2]
        self.Emax = self.data[-1]
        self.logEmin = log(self.data[-2])
        self.logEmax = log(self.data[-1])
        self.coeffs = self.data[:-2]
        self.mass = self.incident_species[0].mass

        # --- Calculate the cross sections at Emin*2 and Emax/2. These are
        # --- used to provide an approximate scaling for E < Emax and E > Emax.
        # --- The cross section is assumed to scale linearly in log-log space
        # --- beyond the Emin and Emax, except that if the cross section is
        # --- increasing, it is assumed fixed instead.
        self.logcsEmin = log(self.chebyshev(self.Emin))
        self.logcsEmin2 = log(max(self.chebyshev(self.Emin*2.),self.chebyshev(self.Emin)))
        self.logcsslopemin = (self.logcsEmin2 - self.logcsEmin)/(log(self.Emin*2.) - self.logEmin)
        self.logcsEmax = log(self.chebyshev(self.Emax))
        self.logcsEmax2 = log(max(self.chebyshev(self.Emax/2.),self.chebyshev(self.Emax)))
        self.logcsslopemax = (self.logcsEmax2 - self.logcsEmax)/(log(self.Emax/2.) - self.logEmax)

    def __call__(self,vi):
        gamma = sqrt(1./(1. - vi**2/clight**2))
        E = self.mass*clight**2*(gamma - 1.)/jperev
        result = self.chebyshev(E)
        result[E==0] = 0.
        ii = nonzero(logical_and(0. < E,E < self.Emin))[0]
        if len(ii) > 0:
            result[ii] = exp(self.logcsslopemin*(log(E[ii]) - self.logEmin) + self.logcsEmin)
        ii = nonzero(E > self.Emax)[0]
        if len(ii) > 0:
            result[ii] = exp(self.logcsslopemax*(log(E[ii]) - self.logEmax) + self.logcsEmax)
        return result

    def constructheadertext(self,species):
        if isinstance(species,Species):
            stype = species.type
            charge_state = species.charge_state
        elif isinstance(species,Particle):
            stype = species
            charge_state = 0
        else:
            return ''
        t = stype.Symbol
        if t[-1] == '2': t = t[:-1] + '\{2\}'
        t += ' \[\%+d\]'%charge_state
        return t

    def constuctheaders(self):
        h_inc0 = self.constructheadertext(self.incident_species[0])
        if len(self.incident_species) == 2:
            h_inc1 = self.constructheadertext(self.incident_species[1])
        else:
            h_inc1 = ''
        h_emit = []
        h_electron = ''
        for s in self.emitted_species:
            if isinstance(s,Species):
                etype = s.type
            elif isinstance(s,Particle):
                etype = s
            if etype == Electron:
                if h_electron == '':
                    h_electron = '  e'
                else:
                    if h_electron == '  e': ie = 1
                    else:                   ie = int(h_electron[2])
                    ie += 1
                    h_electron = '  %de'%ie
            else:
                h_emit.append(self.constructheadertext(s))
        h_emit0 = h_emit[0]
        if len(h_emit) == 2:
            h_emit1 = h_emit[1]
        else:
            h_emit1 = ''
        #print h_inc0 + '  ' + h_inc1 + '   ' + h_emit0 + '  ' + h_emit1 + h_electron
        #print h_inc1 + '  ' + h_inc0 + '   ' + h_emit0 + '  ' + h_emit1 + h_electron
        #print h_inc0 + '  ' + h_inc1 + '   ' + h_emit1 + '  ' + h_emit0 + h_electron
        #print h_inc1 + '  ' + h_inc0 + '   ' + h_emit1 + '  ' + h_emit0 + h_electron
        h1 = re.compile(h_inc0 + '  ' + h_inc1 + '   ' + h_emit0 + '  ' + h_emit1 + h_electron)
        h2 = re.compile(h_inc1 + '  ' + h_inc0 + '   ' + h_emit0 + '  ' + h_emit1 + h_electron)
        h3 = re.compile(h_inc0 + '  ' + h_inc1 + '   ' + h_emit1 + '  ' + h_emit0 + h_electron)
        h4 = re.compile(h_inc1 + '  ' + h_inc0 + '   ' + h_emit1 + '  ' + h_emit0 + h_electron)
        return h1,h2,h3,h4

    def datafilename(self):
        import ionization
        return os.path.join(os.path.dirname(ionization.__file__),'aladdin_8.txt')

    def findcollision(self):
        h1,h2,h3,h4 = self.constuctheaders()
        found = False
        data = []
        with open(self.datafilename(),'r') as f:
            for line in f:
                if found:
                    if line[0] == '&':
                        eval_func = line.split()[-1]
                        if eval_func != '#CHEB':
                            found = False
                        continue
                    if len(line) < 5: break
                    data.extend([float(x) for x in line.split()])
                if h1.search(line):
                    self.name = 'Aladdin ' + h1.pattern.replace('\\','')
                    found = True
                if h2.search(line):
                    self.name = 'Aladdin ' + h2.pattern.replace('\\','')
                    found = True
                if h3.search(line):
                    self.name = 'Aladdin ' + h3.pattern.replace('\\','')
                    found = True
                if h4.search(line):
                    self.name = 'Aladdin ' + h4.pattern.replace('\\','')
                    found = True
        return data

    def chebyshev(self,E):
        try:
            Eclipped = E.clip(self.Emin,self.Emax)
        except AttributeError:
            # --- clip is not defined - hopefully it means that E is a scalar
            Eclipped = max(self.Emin,min(self.Emax,E))
        logE = log(Eclipped)
        x = ((logE - self.logEmin) - (self.logEmax - logE))/(self.logEmax - self.logEmin)
        result = self.coeffs[-1]
        prev2 = 0.
        for c in self.coeffs[-2:0:-1]:
            prev = result
            result = c + 2*x*prev - prev2
            prev2 = prev
        prev = result
        result = 0.5*self.coeffs[0] + x*prev - prev2
        return exp(result)*cm**2

class TableCrossSection(object):
    """Will generate a cross section given a table of data. The data must be
a 2-D array with the following format.

data = array([[E1,CS1],[E2,CS2],[E3,CS3],...])

The energies, E1, E2 etc, must be in increasing order.

At least for now, linear interpolation in log-log space is used for simplicity.
    """
    def __init__(self,incident_species,data):
        self.data = data
        self.logdata = log(where(data>0.,data,smallpos))
        self.name = 'Table'
        self.incident_species = incident_species
        self.mass = self.incident_species.mass

        # --- For E outside of the range of data, the cross section is assumed to scale
        # --- linearly in log-log space, except that if the cross section is
        # --- increasing, it is assumed fixed instead.
        self.logcsEmin = self.logdata[0,1]
        self.logcsEmin2 = max(self.logdata[1,1],self.logdata[0,1])
        self.logcsslopemin = (self.logcsEmin2 - self.logcsEmin)/(self.logdata[1,0] - self.logdata[0,0])
        self.logcsEmax = self.logdata[-1,1]
        self.logcsEmax2 = max(self.logdata[-2,1],self.logdata[-1,1])
        self.logcsslopemax = (self.logcsEmax2 - self.logcsEmax)/(self.logdata[-2,0] - self.logdata[-1,0])

    def __call__(self,vi):
        gamma = sqrt(1./(1. - vi**2/clight**2))
        E = self.mass*clight**2*(gamma - 1.)/jperev
        Earray = E
        if not iterable(Earray): Earray = array([Earray])
        Eclipped = Earray.clip(self.data[0,0],self.data[-1,0])

        ii = searchsorted(self.data[1:,0],Eclipped)
        w = (log(Eclipped) - self.logdata[ii,0])/(self.logdata[ii+1,0] - self.logdata[ii,0])
        result = exp(self.logdata[ii,1]*(1. - w) + self.logdata[ii+1,1]*w)

        result[Earray==0] = 0.
        ii = nonzero(logical_and(0. < Earray,Earray < self.data[0,0]))[0]
        if len(ii) > 0:
            result[ii] = exp(self.logcsslopemin*(log(Earray[ii]) - self.logdata[0,0]) + self.logcsEmin)
        ii = nonzero(Earray > self.data[-1,0])[0]
        if len(ii) > 0:
            result[ii] = exp(self.logcsslopemax*(log(Earray[ii]) - self.logdata[-1,0]) + self.logcsEmax)

        if not iterable(E): result = result[0]

        return result


