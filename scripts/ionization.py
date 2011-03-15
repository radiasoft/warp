"""
Ionization: class for generating particles from impact ionization.
"""
__all__ = ['Ionization']
from warp import *
import time
try:
  from txphysics import txionpack
  l_txphysics=1
except:
  l_txphysics=0


ionization_version = "$Id: ionization.py,v 1.17 2011/03/15 23:12:31 grote Exp $"
def ionizationdoc():
  import Ionization
  print Ionization.__doc__

class Ionization:
  """
Class for generating particles from impact ionization. 
 - incident_species:   incident species
 - target_species:   target species
 - emitted_species:   created species
 - cross_section: cross sections
  """
  def __init__(self,l_verbose=0,stride=100,nx=None,ny=None,nz=None,xmin=None,xmax=None,ymin=None,ymax=None,zmin=None,zmax=None,l_timing=false):
    self.nx = int(where(nx is None,w3d.nx,nx))
    self.ny = int(where(ny is None,w3d.ny,ny))
    self.nz = int(where(nz is None,w3d.nzlocal,nz))
    self.xmin = float(where(xmin is None,w3d.xmmin,xmin))
    self.xmax = float(where(xmax is None,w3d.xmmax,xmax))
    self.ymin = float(where(ymin is None,w3d.ymmin,ymin))
    self.ymax = float(where(ymax is None,w3d.ymmax,ymax))
    self.zmin = float(where(zmin is None,w3d.zmminlocal,zmin))
    self.zmax = float(where(zmax is None,w3d.zmmaxlocal,zmax))   
    self.dx=(self.xmax-self.xmin)/self.nx
    self.dy=(self.ymax-self.ymin)/self.ny
    self.dz=(self.zmax-self.zmin)/self.nz
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
    self.l_timing=l_timing
    if l_txphysics:
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
      self.txphysics_targets={Hydrogen:0, \
                              Dihydrogen:1, \
                              Helium:2, \
                              Carbon_Dioxide:3, \
                              Carbon_Monoxide:4, \
                              Dioxygen:5, \
                              Dinitrogen:6, \
                              Argon:7, \
                              Neon:8}
    self.install()

  def add(self,incident_species,emitted_species,cross_section=None,
          target_species=None,ndens=None,emitted_energy0=0.,emitted_energy_sigma=0.,
          incident_pgroup=top.pgroup,target_pgroup=top.pgroup,emitted_pgroup=top.pgroup):
    if not self.inter.has_key(incident_species):
        self.inter[incident_species]={}
        for key in ['target_species','emitted_species','cross_section','ndens',
                    'remove_incident','remove_target','emitted_energy0','emitted_energy_sigma',
                    'incident_pgroup','target_pgroup','emitted_pgroup']:
          self.inter[incident_species][key]=[]
    if type(emitted_species)<>type([]):emitted_species=[emitted_species]
    self.inter[incident_species]['target_species']  +=[target_species]
    if emitted_species[0] in [incident_species,target_species]:
      self.inter[incident_species]['emitted_species'] +=[emitted_species[1:]]
    else:
      self.inter[incident_species]['emitted_species'] +=[emitted_species]
    self.inter[incident_species]['cross_section']   +=[cross_section]
    self.inter[incident_species]['ndens']           +=[ndens]
    if incident_species is emitted_species[0]:
      self.inter[incident_species]['remove_incident']+=[1]
    else:
      self.inter[incident_species]['remove_incident']+=[0]
    if target_species is emitted_species[0]:
      self.inter[incident_species]['remove_target']+=[1]
    else:
      self.inter[incident_species]['remove_target']+=[0]
    self.inter[incident_species]['emitted_energy0']   +=[emitted_energy0]
    self.inter[incident_species]['emitted_energy_sigma']   +=[emitted_energy_sigma]
    self.inter[incident_species]['incident_pgroup']=incident_pgroup
    self.inter[incident_species]['target_pgroup']  =target_pgroup
    self.inter[incident_species]['emitted_pgroup'] =emitted_pgroup
    if target_species is not None:
      if not self.target_dens.has_key(target_species):
        self.target_dens[target_species]={}
        for key in ['ndens','ndens_updated']:
          self.target_dens[target_species][key]=[]
        self.target_dens[target_species]['ndens']           =fzeros((self.nx+1,self.ny+1,self.nz+1),'d')
        self.target_dens[target_species]['ndens_updated']   =0
      
    for e in emitted_species:
      js=e.jslist[0]
      if not self.x.has_key(emitted_pgroup):
        self.nps[emitted_pgroup]={}
        self.x[emitted_pgroup]={}
        self.y[emitted_pgroup]={}
        self.z[emitted_pgroup]={}
        self.ux[emitted_pgroup]={}
        self.uy[emitted_pgroup]={}
        self.uz[emitted_pgroup]={}
        self.gi[emitted_pgroup]={}
      if not self.x[emitted_pgroup].has_key(js):
        self.nps[emitted_pgroup][js]=0
        self.x[emitted_pgroup][js]=fzeros(self.npmax,'d')
        self.y[emitted_pgroup][js]=fzeros(self.npmax,'d')
        self.z[emitted_pgroup][js]=fzeros(self.npmax,'d')
        self.ux[emitted_pgroup][js]=fzeros(self.npmax,'d')
        self.uy[emitted_pgroup][js]=fzeros(self.npmax,'d')
        self.uz[emitted_pgroup][js]=fzeros(self.npmax,'d')
        self.gi[emitted_pgroup][js]=fzeros(self.npmax,'d')

  def install(self):
    if not isinstalleduserinjection(self.generate):
      installuserinjection(self.generate)

  def addpart(self,nn,x,y,z,ux,uy,uz,gi,pg,js):
    ilf=0
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
#           raise('')
         addparticles(x=self.x[pg][js][:nn],
                      y=self.y[pg][js][:nn],
                      z=self.z[pg][js][:nn],
                      vx=self.ux[pg][js][:nn],
                      vy=self.uy[pg][js][:nn],
                      vz=self.uz[pg][js][:nn],
                      gi=self.gi[pg][js][:nn],
                      lmomentum=True,
                      pgroup=pg,
                      js=js)
         self.nps[pg][js]=0

  def printall(self,l_cgm=0):
    swidth=0
    twidth=0
    ewidth={}
    title='        *** Particle/particle interactions ***\n'
    textblock=''
    for incident_species in self.inter.keys():
      swidth=max(swidth,len(incident_species.name))
      for it,target_species in enumerate(self.inter[incident_species]['target_species']):
        if target_species is None:
          tname='background gas'
        else:
          tname=target_species.name
        twidth=max(twidth,len(tname))
        firstelement = []
        if not self.inter[incident_species]['remove_incident'][it]:
          firstelement = [incident_species]
        if not self.inter[incident_species]['remove_target'][it]:
          firstelement = [target_species]
        for ie,emitted_species in enumerate(firstelement+self.inter[incident_species]['emitted_species'][it]):
          ename=emitted_species.name
          if not ewidth.has_key(ie):
            ewidth[ie]=len(ename)
          else:
            ewidth[ie]=max(ewidth[ie],len(ename))
    fs='%%-%gs'%swidth
    ft='%%-%gs'%twidth
    fe={}
    for ie in ewidth.keys():
      fe[ie]='%%-%gs'%ewidth[ie]
    for incident_species in self.inter.keys():
      sname=fs%incident_species.name[:swidth]
      textblock+='\n'
      for it,target_species in enumerate(self.inter[incident_species]['target_species']):
        if target_species is None:
          tname=ft%'background gas'[:twidth]
        else:
          tname=ft%target_species.name[:twidth]
        if self.inter[incident_species]['cross_section'][it] is None:
          cname = ''
        else:
          cname='(CS=%g)'%(self.inter[incident_species]['cross_section'][it])
        firstelement = []
        if not self.inter[incident_species]['remove_incident'][it]:
          firstelement = [incident_species]
        if not self.inter[incident_species]['remove_target'][it]:
          firstelement = [target_species]
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

#printall(io,l_cgm=1)

  def generate(self,dt=None):
    if dt is None:dt=top.dt
    if self.l_timing:t1 = time.clock()
    for target_species in self.target_dens.keys():
      self.target_dens[target_species]['ndens_updated']=0    
    for incident_species in self.inter.keys():
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
        if ndens is not None:
          continue
        else:
          if self.target_dens[target_species]['ndens_updated']:
            continue
          else:
            self.target_dens[target_species]['ndens_updated']=1
          ndens = self.target_dens[target_species]['ndens']
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
            fact=1.
            if w3d.l4symtry:
              xt=abs(xt)
              fact=0.25
            elif w3d.l2symtry:
              fact=0.5
            if w3d.l2symtry or w3d.l4symtry:yt=abs(yt)
            deposgrid3d(1,tpg.nps[jstarget],xt,yt,zt,
                        tpg.sw[jstarget]*self.invvol*fact*ones(tpg.nps[jstarget],'d'),
                        self.nx,self.ny,self.nz,ndens,self.ndensc,
                        self.xmin,self.xmax,self.ymin,self.ymax,
                        self.zmin,self.zmax)    

#          if w3d.l2symtry or w3d.l4symtry:self.ndens[:,0,:]*=2.
#          if w3d.l4symtry:self.ndens[0,:,:]*=2.
      
    for incident_species in self.inter.keys():
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
        for js in incident_species.jslist:
          i1 = ipg.ins[js] - 1 + top.it%self.stride
          i2 = ipg.ins[js] + ipg.nps[js] - 1
          xi=ipg.xp[i1:i2:self.stride]#.copy()
          yi=ipg.yp[i1:i2:self.stride]#.copy()
          zi=ipg.zp[i1:i2:self.stride]#.copy()
          ni = shape(xi)[0]
          gaminvi=ipg.gaminv[i1:i2:self.stride].copy()
          uxi=ipg.uxp[i1:i2:self.stride].copy()
          uyi=ipg.uyp[i1:i2:self.stride].copy()
          uzi=ipg.uzp[i1:i2:self.stride].copy()
          # --- get velocity in lab frame if using a boosted frame of reference
          if top.boost_gamma>1.:
            uzboost = clight*sqrt(top.boost_gamma**2-1.)
            setu_in_uzboosted_frame3d(ni,uxi,uyi,uzi,gaminvi,
                                      -uzboost,
                                      top.boost_gamma)
          vxi=uxi*gaminvi
          vyi=uyi*gaminvi
          vzi=uzi*gaminvi
          # compute the relative velocity
          # NOTE that at this point, the target species is assumed to have a negligible velocity.
          # this needs to be modified if this approximation is not valid.
          vi=sqrt(vxi*vxi+vyi*vyi+vzi*vzi)
          if ndens is None:
            ndens = self.target_dens[target_species]['ndens']
          if type(ndens) in [type(0.),type(0)]:
            dp=ndens
          else:
            dp=zeros(ni,'d')
            getgrid3d(ni,xi,yi,zi,dp,
                        self.nx,self.ny,self.nz,ndens, \
                        self.xmin,self.xmax,self.ymin,self.ymax,self.zmin,self.zmax, 
                        w3d.l2symtry,w3d.l4symtry)    
          # if no cross-section is provided, then attempt to get it from txphysics
          cross_section = self.inter[incident_species]['cross_section'][it] # cross-section
          if cross_section is None:
            try:
              gastype = self.txphysics_targets[target_species.type]
              if self.l_verbose:print incident_species.name+' on '+target_species.name+':', gastype
            except:                    
              raise('Error in ionization: cross section of '+incident_species.name+' on '+target_species.name+' is not available. Please provide.')
            # This is an integer flag to specify electron (0) or ion (1)
            if incident_species.type in (Electron,Positron):
              incident = 0
            else:
              incident = 1
            cross_section = txionpack.get_sigma_impact_array(vi, gastype, incident)

          # probability
          ncol = dp*cross_section*vi*dt*ipg.ndts[js]*self.stride
          if top.boost_gamma>1.:ncol*=top.gammabar_lab/top.gammabar
          l_ionization_projectile=0
          if self.inter[incident_species]['remove_incident'][it]:
            ncol = where(ncol>=1.,1.-1.e-10,ncol)
          ncoli=int(ncol)
          ncol=ncol-ncoli
          r=ranf(ncol)
          ncoli=where(r<ncol,ncoli+1,ncoli)
          io=compress(ncoli>0,arange(ni))
          nnew = len(io)
          if self.inter[incident_species]['remove_incident'][it]:
            uxnew = uxi
            uynew = uyi
            uznew = uzi
            # if projectile is modified, then need to delete it
            put(ipg.gaminv,array(io)*self.stride+i1,0.)
          xnew = xi
          ynew = yi
          znew = zi
          ifg = 0
          while(nnew>0):
            xnewp = take(xnew,io)
            ynewp = take(ynew,io)
            znewp = take(znew,io)
            ifg+=1
            xnew = xnewp+(ranf(xnewp)-0.5)*1.e-10*self.dx
            ynew = ynewp+(ranf(ynewp)-0.5)*1.e-10*self.dy
            znew = znewp+(ranf(znewp)-0.5)*1.e-10*self.dz
            for emitted_species in self.inter[incident_species]['emitted_species'][it]:
              if self.inter[incident_species]['remove_incident'][it]:
                uxnew = take(uxnew,io)
                uynew = take(uynew,io)
                uznew = take(uznew,io)
              else:
#                uxnew = zeros(nnew,float64)+1.e-10
#                uynew = zeros(nnew,float64)+1.e-10
#                uznew = zeros(nnew,float64)+1.e-10
                 ek0ionel = self.inter[incident_species]['emitted_energy0'][it]
                 esigionel = self.inter[incident_species]['emitted_energy_sigma'][it]
                 if esigionel==0.:
                   ek = zeros(nnew)
                 else:
                   ek =SpRandom(0.,esigionel,nnew)	#kinetic energy
                 ek=abs(ek+ek0ionel)	#kinetic energy
                 fact = echarge/(emass*clight**2)
                 gamma=ek*fact+1.		
                 u=clight*sqrt(ek*fact*(gamma+1.))	
                 # velocity direction: random in (x-y) plane plus small longitudianl component:
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

              ginew = 1./sqrt(1.+(uxnew**2+uynew**2+uznew**2)/clight**2)
              # --- get velocity in boosted frame if using a boosted frame of reference
              if top.boost_gamma>1.:
                setu_in_uzboosted_frame3d(shape(ginew)[0],uxnew,uynew,uznew,ginew,
                                          uzboost,
                                          top.boost_gamma)

              if self.l_verbose:print 'add ',nnew, emitted_species.name,' from by impact ionization:',incident_species.name,'+',target_species.name 
              if self.inter[incident_species]['remove_incident'][it] and (emitted_species.type is incident_species.type):
                self.addpart(nnew,xnewp,ynewp,znewp,uxnew,uynew,uznew,ginew,epg,emitted_species.jslist[0])
              else:
                self.addpart(nnew,xnew,ynew,znew,uxnew,uynew,uznew,ginew,epg,emitted_species.jslist[0])
            ncoli=take(ncoli,io)-1
            io=compress(ncoli>0,arange(nnew))
            nnew = len(io)
                                           
    # make sure that all particles are added and cleared 
    for pg in self.x.keys():
      for js in self.x[pg].keys():
        self.flushpart(pg,js)
        processlostpart(pg,js+1,top.clearlostpart,top.time,top.zbeam)
                       
    if self.l_timing:print 'time ionization = ',time.clock()-t1,'s'

