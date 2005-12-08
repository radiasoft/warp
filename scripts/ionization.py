"""
Ionization: class for generating particles from impact ionization.
"""
from warp import *
import time

ionization_version = "$Id: ionization.py,v 1.2 2005/12/08 23:36:45 jlvay Exp $"
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
    self.nz = int(where(nz is None,w3d.nz,nz))
    self.xmin = float(where(xmin is None,w3d.xmmin,xmin))
    self.xmax = float(where(xmax is None,w3d.xmmax,xmax))
    self.ymin = float(where(ymin is None,w3d.ymmin,ymin))
    self.ymax = float(where(ymax is None,w3d.ymmax,ymax))
    self.zmin = float(where(zmin is None,w3d.zmmin,zmin))
    self.zmax = float(where(zmax is None,w3d.zmmax,zmax))   
    self.dx=(self.xmax-self.xmin)/w3d.nx
    self.dy=(self.ymax-self.ymin)/w3d.ny
    self.dz=(self.zmax-self.zmin)/w3d.nz
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
    self.install()
    
  def add(self,incident_species,emitted_species,cross_section,target_species=None,ndens=None):
    if not self.inter.has_key(incident_species):
        self.inter[incident_species]={}
        for key in ['target_species','emitted_species','cross_section','ndens']:
          self.inter[incident_species][key]=[]
    if type(emitted_species)<>type([]):emitted_species=[emitted_species]
    self.inter[incident_species]['target_species']  +=[target_species]
    self.inter[incident_species]['emitted_species'] +=[emitted_species]
    self.inter[incident_species]['cross_section']   +=[cross_section]
    self.inter[incident_species]['ndens']           +=[ndens]
    if target_species is not None:
      if not self.target_dens.has_key(target_species):
        self.target_dens[target_species]={}
        for key in ['ndens','ndens_updated']:
          self.target_dens[target_species][key]=[]
        self.target_dens[target_species]['ndens']           =fzeros((self.nx+1,self.ny+1,self.nz+1),'d')
        self.target_dens[target_species]['ndens_updated']   =0
      
    self.inter[incident_species]['cross_section']   +=[cross_section]
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
    if not isinstalledafterfs(self.generate):
      installafterfs(self.generate)

  def addpart(self,nn,x,y,z,vx,vy,vz,js):
    if self.nps[js]+nn>self.npmax:self.flushpart(js)
    il=self.nps[js]
    iu=il+nn
    xx=ones(nn,'d')
    self.x[js][il:iu]=x*xx
    self.y[js][il:iu]=y*xx
    self.z[js][il:iu]=z*xx
    self.vx[js][il:iu]=vx
    self.vy[js][il:iu]=vy
    self.vz[js][il:iu]=vz
    self.nps[js]+=nn
      
  def flushpart(self,js):
      if self.nps[js]>0:
         nn=self.nps[js]
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
        for ie,emitted_species in enumerate(self.inter[incident_species]['emitted_species'][it]):
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
        cname='CS=%g'%(self.inter[incident_species]['cross_section'][it])
        for ie,emitted_species in enumerate(self.inter[incident_species]['emitted_species'][it]):
          if ie==0:
            ename=fe[ie]%emitted_species.name[:ewidth[ie]]
          else:
            ename+=' + '+fe[ie]%emitted_species.name[:ewidth[ie]]
        textblock += sname+' + '+tname+' => '+ename+' ('+cname+')\n'
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
        npinc+=top.nps[js]
        if top.it%top.ndts[js]==0:ispushed=1
      if npinc==0 or not ispushed:continue
      for it,target_species in enumerate(self.inter[incident_species]['target_species']):
        cross_section = self.inter[incident_species]['cross_section'][it] # cross-section
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
            nptarget+=top.nps[jstarget]
          if nptarget==0:continue
          self.ndensc[...]=0.
          ndens[...]=0.
          for jstarget in target_species.jslist:
            i1 = top.ins[jstarget] - 1
            i2 = top.ins[jstarget] + top.nps[jstarget] - 1
            xt=top.xp[i1:i2]
            yt=top.yp[i1:i2]
            zt=top.zp[i1:i2]
            fact=1.
            if w3d.l4symtry:
              xt=abs(xt)
              fact=0.25
            elif w3d.l2symtry:
              fact=0.5
            if w3d.l2symtry or w3d.l4symtry:yt=abs(yt)
            deposgrid3d(1,top.nps[jstarget],xt,yt,zt,top.sw[jstarget]*self.invvol*fact*ones(top.nps[jstarget],'d'), \
                        self.nx,self.ny,self.nz,ndens,self.ndensc, \
                        self.xmin,self.xmax,self.ymin,self.ymax,self.zmin,self.zmax)    

#          if w3d.l2symtry or w3d.l4symtry:self.ndens[:,0,:]*=2.
#          if w3d.l4symtry:self.ndens[0,:,:]*=2.
      
    for incident_species in self.inter.keys():
      npinc = 0
      ispushed=0
      for js in incident_species.jslist:
        npinc+=top.nps[js]
        if top.it%top.ndts[js]==0:ispushed=1
      if npinc==0 or not ispushed:continue
      for it,target_species in enumerate(self.inter[incident_species]['target_species']):
        cross_section = self.inter[incident_species]['cross_section'][it] # cross-section
        ndens=self.inter[incident_species]['ndens'][it]
        for js in incident_species.jslist:
          i1 = top.ins[js] - 1 + top.it%self.stride
          i2 = top.ins[js] + top.nps[js] - 1
          xi=top.xp[i1:i2:self.stride]#.copy()
          yi=top.yp[i1:i2:self.stride]#.copy()
          zi=top.zp[i1:i2:self.stride]#.copy()
          ni = shape(xi)[0]
          gaminvi=top.gaminv[i1:i2:self.stride]#.copy()
          vxi=top.uxp[i1:i2:self.stride]*gaminvi#.copy()
          vyi=top.uyp[i1:i2:self.stride]*gaminvi#.copy()
          vzi=top.uzp[i1:i2:self.stride]*gaminvi#.copy()
          vi=sqrt(vxi*vxi+vyi*vyi+vzi*vzi)
          if ndens is not None:
            dp=ndens
          else:
            ndens = self.target_dens[target_species]['ndens']
            dp=zeros(ni,'d')
            getgrid3d(ni,xi,yi,zi,dp,
                        self.nx,self.ny,self.nz,ndens, \
                        self.xmin,self.xmax,self.ymin,self.ymax,self.zmin,self.zmax, 
                        w3d.l2symtry,w3d.l4symtry)    
          # probability
          ncol = dp*cross_section*vi*top.dt*top.ndts[js]*self.stride
          l_ionization_projectile=0
          for emitted_species in self.inter[incident_species]['emitted_species'][it]:
            if emitted_species.type is incident_species.type:l_ionization_projectile=1           
          if l_ionization_projectile:
            ncol = where(ncol>=1.,1.-1.e-10,ncol)
          ncoli=int(ncol)
          ncol=ncol-ncoli
          ncoli+=1
          r=ranf(ncol)
          io=compress(r<ncol,arange(ni))
          nnew = len(io)
          if l_ionization_projectile:
            vxnew = take(vxi,io)
            vynew = take(vyi,io)
            vznew = take(vzi,io)
            put(top.gaminv[i1:i2],io,0.)
          xnew = xi
          ynew = yi
          znew = zi
          while(nnew>0):
            xnewp = take(xnew,io)
            ynewp = take(ynew,io)
            znewp = take(znew,io)
            vnew = zeros(nnew,Float)+1.e-10
            xnew = xnewp+ranf(vnew)*0.25*self.dx
            ynew = ynewp+ranf(vnew)*0.25*self.dy
            znew = znewp+ranf(vnew)*0.25*self.dz
#          print 'add ',nnew, 'particles by impact ionization'
            for emitted_species in self.inter[incident_species]['emitted_species'][it]:
              if emitted_species.type is incident_species.type:
                self.addpart(nnew,xnewp,ynewp,znewp,vxnew,vynew,vznew,emitted_species.jslist[0])
              else:
                self.addpart(nnew,xnew,ynew,znew,vnew,vnew,vnew,emitted_species.jslist[0])
            ncoli=take(ncoli,io)-1
            io=compress(ncoli>0,arange(nnew))
            nnew = len(io)
                                           
    # make sure that all particles are added and cleared 
    for js in self.x.keys():
      self.flushpart(js)
      processlostpart(js+1,top.clearlostpart,top.time,top.zbeam)
                       
    if self.l_timing:print 'time ionization = ',time.clock()-t1,'s'

