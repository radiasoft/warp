"""
Ionization: class for generating particles from impact ionization.
"""
from warp import *
from species import *
import time
try:
  from txphysics import txionpack
  l_txphysics=1
except:
  l_txphysics=0

ionization_version = "$Id: ionization.py,v 1.14 2009/02/05 18:18:40 dave Exp $"
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
    self.vx={}
    self.vy={}
    self.vz={}
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
    
  def add(self,incident_species,emitted_species,cross_section=None,target_species=None,ndens=None):
    if not self.inter.has_key(incident_species):
        self.inter[incident_species]={}
        for key in ['target_species','emitted_species','cross_section','ndens','remove_incident','remove_target']:
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
      self.inter[incident_species]['remove_incident']+=[0]
    else:
      self.inter[incident_species]['remove_incident']+=[1]
    if target_species is emitted_species[0]:
      self.inter[incident_species]['remove_target']+=[0]
    else:
      self.inter[incident_species]['remove_target']+=[1]
    if target_species is not None:
      if not self.target_dens.has_key(target_species):
        self.target_dens[target_species]={}
        for key in ['ndens','ndens_updated']:
          self.target_dens[target_species][key]=[]
        self.target_dens[target_species]['ndens']           =fzeros((self.nx+1,self.ny+1,self.nz+1),'d')
        self.target_dens[target_species]['ndens_updated']   =0
      
    for e in emitted_species:
      js=e.jslist[0]
      if not self.x.has_key(js):
        self.nps[js]=0
        self.x[js]=fzeros(self.npmax,'d')
        self.y[js]=fzeros(self.npmax,'d')
        self.z[js]=fzeros(self.npmax,'d')
        self.vx[js]=fzeros(self.npmax,'d')
        self.vy[js]=fzeros(self.npmax,'d')
        self.vz[js]=fzeros(self.npmax,'d')

  def install(self):
    if not isinstalleduserinjection(self.generate):
      installuserinjection(self.generate)

  def addpart(self,nn,x,y,z,vx,vy,vz,js):
    ilf=0
    while self.nps[js]+nn>self.npmax:
      il=self.nps[js]
      iu=min(il+nn,self.npmax)
      nf=iu-il
      self.x [js][il:iu]= x[ilf:ilf+nf]
      self.y [js][il:iu]= y[ilf:ilf+nf]
      self.z [js][il:iu]= z[ilf:ilf+nf]
      self.vx[js][il:iu]=vx[ilf:ilf+nf]
      self.vy[js][il:iu]=vy[ilf:ilf+nf]
      self.vz[js][il:iu]=vz[ilf:ilf+nf]
      self.nps[js]+=nf
      self.flushpart(js)
      ilf+=nf
      nn-=nf
    il=self.nps[js]
    iu=il+nn
    self.x [js][il:iu]=x [ilf:]
    self.y [js][il:iu]=y [ilf:]
    self.z [js][il:iu]=z [ilf:]
    self.vx[js][il:iu]=vx[ilf:]
    self.vy[js][il:iu]=vy[ilf:]
    self.vz[js][il:iu]=vz[ilf:]
    self.nps[js]+=nn
      
  def flushpart(self,js):
      if self.nps[js]>0:
         nn=self.nps[js]
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
         addparticles(x=self.x[js][:nn],
                      y=self.y[js][:nn],
                      z=self.z[js][:nn],
                      vx=self.vx[js][:nn],
                      vy=self.vy[js][:nn],
                      vz=self.vz[js][:nn],
                      js=js)
         self.nps[js]=0

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

  def generate(self):
    if self.l_timing:t1 = time.clock()
    for target_species in self.target_dens.keys():
      self.target_dens[target_species]['ndens_updated']=0    
    for incident_species in self.inter.keys():
      npinc = 0
      ispushed=0
      for js in incident_species.jslist:
        npinc+=top.pgroup.nps[js]
        if top.pgroup.ldts[js]:ispushed=1
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
            nptarget+=top.pgroup.nps[jstarget]
          if nptarget==0:continue
          self.ndensc[...]=0.
          ndens[...]=0.
          for jstarget in target_species.jslist:
            i1 = top.pgroup.ins[jstarget] - 1
            i2 = top.pgroup.ins[jstarget] + top.pgroup.nps[jstarget] - 1
            xt=top.pgroup.xp[i1:i2]
            yt=top.pgroup.yp[i1:i2]
            zt=top.pgroup.zp[i1:i2]
            fact=1.
            if w3d.l4symtry:
              xt=abs(xt)
              fact=0.25
            elif w3d.l2symtry:
              fact=0.5
            if w3d.l2symtry or w3d.l4symtry:yt=abs(yt)
            deposgrid3d(1,top.pgroup.nps[jstarget],xt,yt,zt,
                        top.pgroup.sw[jstarget]*self.invvol*fact*ones(top.pgroup.nps[jstarget],'d'),
                        self.nx,self.ny,self.nz,ndens,self.ndensc,
                        self.xmin,self.xmax,self.ymin,self.ymax,
                        self.zmin,self.zmax)    

#          if w3d.l2symtry or w3d.l4symtry:self.ndens[:,0,:]*=2.
#          if w3d.l4symtry:self.ndens[0,:,:]*=2.
      
    for incident_species in self.inter.keys():
      npinc = 0
      ispushed=0
      for js in incident_species.jslist:
        npinc+=top.pgroup.nps[js]
        if top.pgroup.ldts[js]:ispushed=1
      if npinc==0 or not ispushed:continue
      for it,target_species in enumerate(self.inter[incident_species]['target_species']):
        ndens=self.inter[incident_species]['ndens'][it]
        for js in incident_species.jslist:
          i1 = top.pgroup.ins[js] - 1 + top.it%self.stride
          i2 = top.pgroup.ins[js] + top.pgroup.nps[js] - 1
          xi=top.pgroup.xp[i1:i2:self.stride]#.copy()
          yi=top.pgroup.yp[i1:i2:self.stride]#.copy()
          zi=top.pgroup.zp[i1:i2:self.stride]#.copy()
          ni = shape(xi)[0]
          gaminvi=top.pgroup.gaminv[i1:i2:self.stride]#.copy()
          vxi=top.pgroup.uxp[i1:i2:self.stride]*gaminvi#.copy()
          vyi=top.pgroup.uyp[i1:i2:self.stride]*gaminvi#.copy()
          vzi=top.pgroup.uzp[i1:i2:self.stride]*gaminvi#.copy()
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
          ncol = dp*cross_section*vi*top.dt*top.pgroup.ndts[js]*self.stride
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
            vxnew = vxi
            vynew = vyi
            vznew = vzi
            # if projectile is modified, then need to delete it
            put(top.pgroup.gaminv,array(io)*self.stride+i1,0.)
          xnew = xi
          ynew = yi
          znew = zi
          while(nnew>0):
            #print nnew
            xnewp = take(xnew,io)
            ynewp = take(ynew,io)
            znewp = take(znew,io)
            vnew = zeros(nnew,float64)+1.e-10
            xnew = xnewp+(ranf(vnew)-0.5)*1.e-10*self.dx
            ynew = ynewp+(ranf(vnew)-0.5)*1.e-10*self.dy
            znew = znewp+(ranf(vnew)-0.5)*1.e-10*self.dz
            if self.inter[incident_species]['remove_incident'][it]:
              vxnew = take(vxnew,io)
              vynew = take(vynew,io)
              vznew = take(vznew,io)
            for emitted_species in self.inter[incident_species]['emitted_species'][it]:
              if self.l_verbose:print 'add ',nnew, emitted_species.name,' from by impact ionization:',incident_species.name,'+',target_species.name 
              if self.inter[incident_species]['remove_incident'][it] and (emitted_species.type is incident_species.type):
                self.addpart(nnew,xnewp,ynewp,znewp,vxnew,vynew,vznew,emitted_species.jslist[0])
              else:
                self.addpart(nnew,xnew,ynew,znew,vnew,vnew,vnew,emitted_species.jslist[0])
            ncoli=take(ncoli,io)-1
            io=compress(ncoli>0,arange(nnew))
            nnew = len(io)
                                           
    # make sure that all particles are added and cleared 
    for js in self.x.keys():
      self.flushpart(js)
      processlostpart(top.pgroup,js+1,top.clearlostpart,top.time,top.zbeam)
                       
    if self.l_timing:print 'time ionization = ',time.clock()-t1,'s'

