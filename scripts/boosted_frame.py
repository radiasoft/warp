from warp import *

class Boosted_Frame(object):
  """
Class transforming particle distribution from lab frame to boosted frame. 
Boosted particles can optionally be injected through a plane. In this case, the 
boosted particles are moved from the main top.pgroup to a separate particle group, 
drift at velocity top.vbeam until they reach the injection plane though which 
they are injected using Warp's injection routines.
In the current implementation, the distribution needs to fit entirely in the 
simulation zone, which is too restrictive for some applications and will need 
to be lifted in the future.
  """
  def __init__(self,gammaframe,direction=1.,l_setselfb=1):
    top.boost_gamma=gammaframe
    self.gammaframe=gammaframe
    self.l_setselfb=l_setselfb
    self.betaframe  = direction*sqrt(1.-1./self.gammaframe**2)
    self.betabeam_lab=top.vbeam/clight
    self.betabeamfrm_lab=top.vbeamfrm/clight
    top.vbeam_lab = top.vbeam
#    top.gammabar_lab = top.gammabar
    top.gammabar_lab=1./sqrt(1.-(top.vbeam_lab/clight)**2)
    top.vbeam=clight*(self.betabeam_lab-self.betaframe)/(1.-self.betabeam_lab*self.betaframe)
    top.vbeamfrm=clight*(self.betabeamfrm_lab-self.betaframe)/(1.-self.betabeamfrm_lab*self.betaframe)
    top.gammabar=1./sqrt(1.-(top.vbeam/clight)**2)
    if self.l_setselfb:
      fselfbcopy = top.fselfb.copy()
      top.fselfb[...]=(top.fselfb[...]-self.betaframe*clight)/(1.-top.fselfb[...]*self.betaframe/clight)
      for js in range(shape(top.pgroup.fselfb)[0]):
        for jss in range(top.nsselfb):
          if top.pgroup.fselfb[js] == fselfbcopy[jss]:
            top.pgroup.fselfb[js]=top.fselfb[jss]

  def boost(self,species,zinject=0.,tinit=0.,l_inject_plane=1,lallindomain=0,
            l_rmzmean=1.,l_deprho=1,l_savezinit=0,l_focus=1,l_project=1):
   print 'enter boost',top.pgroup.nps
   if l_savezinit:
     if top.zbirthlabpid==0:
       top.zbirthlabpid = nextpid()
       top.pgroup.npid = top.npid
       top.pgroup.gchange()
     iupr=-1
     for js in species.jslist:
       ilpr=iupr+1
       iupr=ilpr+getn(js=js,bcast=0,gather=0)
       if getn(js=js,bcast=0,gather=0):
         top.pgroup.pid[ilpr:iupr,top.zbirthlabpid-1] = top.pgroup.zp[ilpr:iupr].copy()
   if l_inject_plane:
    pg = top.pgroup
    self.species=species
    self.zinject=zinject
    self.tinit=tinit
    self.l_rmzmean=l_rmzmean
    self.pgroup = ParticleGroup()
#    self.pgroup.ns = pg.ns#len(species.jslist)
    self.pgroup.ns = len(species.jslist)
    self.pgroup.npmax = species.getn(bcast=0,gather=0)
    self.pgroup.npid = pg.npid
    self.pgroup.lebcancel_pusher = top.pgroup.lebcancel_pusher
    self.pgroup.gchange()
    iupr=-1
    for jspr,js in enumerate(species.jslist):
      ilpr=iupr+1
      iupr=ilpr+getn(pgroup=pg,js=js,bcast=0,gather=0)
      self.pgroup.sq[jspr] = pg.sq[js]
      self.pgroup.sm[jspr] = pg.sm[js]
      self.pgroup.sw[jspr] = pg.sw[js]
      self.pgroup.sid[jspr] = pg.sid[js]
      self.pgroup.ndts[jspr] = pg.ndts[js]
      self.pgroup.ldts[jspr] = pg.ldts[js]
      self.pgroup.lvdts[jspr] = pg.lvdts[js]
      self.pgroup.iselfb[jspr] = pg.iselfb[js]
      self.pgroup.dtscale[jspr] = pg.dtscale[js]
      self.pgroup.limplicit[jspr] = pg.limplicit[js]
      self.pgroup.iimplicit[jspr] = pg.iimplicit[js]
      self.pgroup.zshift[jspr] = pg.zshift[js]
      self.pgroup.ins[jspr]=ilpr+1
      self.pgroup.nps[jspr]=getn(pgroup=pg,js=js,bcast=0,gather=0)
      if l_inject_plane:
       if getn(pgroup=pg,js=js,bcast=0,gather=0)>0: 
        z=getz(pgroup=pg,js=js,bcast=0,gather=0)
       else:
        z=array([])
       if self.l_rmzmean:
         zmean=globalave(z)
       else:
         zmean=0.
       if getn(pgroup=pg,js=js,bcast=0,gather=0)>0: 
        gaminvbeam_lab = getgaminv(pgroup=pg,js=js,bcast=0,gather=0)
        betabeam_lab  = sqrt(1.-gaminvbeam_lab*gaminvbeam_lab)
        betabeam_frame = (betabeam_lab-self.betaframe)/(1.-betabeam_lab*self.betaframe)
        gammabeam_frame  = 1./sqrt(1.-betabeam_frame*betabeam_frame)
        zcopy = z.copy()
        z=z-zmean
        # --- get data at z=0
        vx = getvx(pgroup=pg,js=js,bcast=0,gather=0)
        vy = getvy(pgroup=pg,js=js,bcast=0,gather=0)
        vz = getvz(pgroup=pg,js=js,bcast=0,gather=0)
        t = z/vz
        x = getx(pgroup=pg,js=js,bcast=0,gather=0)
        y = gety(pgroup=pg,js=js,bcast=0,gather=0)
        if not l_project:
          x-=t*vx
          y-=t*vy
        # --- correct for focusing effect from shift from z=0 to zinject
        if l_focus:
          tfoc = -zinject*self.gammaframe/vz#pr
          x = x-tfoc*vx#pr
          y = y-tfoc*vy#pr
        # --- get data in boosted frame
        tpr = -self.gammaframe*t
        zpr = self.gammaframe*self.betaframe*clight*t
       else:
        zpr=array([])
       if top.boost_z0==0.:
         top.boost_z0 = -globalmax(zpr)
       if getn(pgroup=pg,js=js,bcast=0,gather=0)>0: 
        fact = 1./(1.-self.betaframe*vz/clight)
        vxpr = vx*fact/self.gammaframe
        vypr = vy*fact/self.gammaframe
        vzpr = (vz-self.betaframe*clight)*fact
       # --- get data at t=0 in boosted frame
#        zpr = zpr - vzpr*tpr
        zpr = zpr - top.vbeam*tpr
#        zpr = zcopy*top.gammabar_lab/top.gammabar
        # --- make sure that z<=0
#        zpr += top.boost_z0 
        # --- sets location of beam center at t=0 in boosted frame
        gammapr = 1./sqrt(1.-(vxpr*vxpr+vypr*vypr+vzpr*vzpr)/clight**2)
        self.pgroup.uxp[ilpr:iupr]=vxpr*gammapr
        self.pgroup.uyp[ilpr:iupr]=vypr*gammapr
        self.pgroup.uzp[ilpr:iupr]=vzpr*gammapr
        self.pgroup.gaminv[ilpr:iupr]=1./gammapr
        self.pgroup.xp[ilpr:iupr] = x
        self.pgroup.yp[ilpr:iupr] = y
        self.pgroup.zp[ilpr:iupr] = zpr
        if pg.npid>0:self.pgroup.pid[ilpr:iupr,:] = getpid(pgroup=pg,js=js,bcast=0,gather=0,id=-1)
        if top.uxoldpid>0:self.pgroup.pid[ilpr:iupr,top.uxoldpid-1]=self.pgroup.uxp[ilpr:iupr]
        if top.uyoldpid>0:self.pgroup.pid[ilpr:iupr,top.uyoldpid-1]=self.pgroup.uyp[ilpr:iupr]
        if top.uzoldpid>0:self.pgroup.pid[ilpr:iupr,top.uzoldpid-1]=self.pgroup.uzp[ilpr:iupr]
      pg.nps[js]=0
#      if pg.fselfb[js] != 0.:
#        pg.fselfb[js]=(pg.fselfb[js]-self.betaframe*clight)/(1.-pg.fselfb[js]*self.betaframe/clight)
      self.pgroup.fselfb[jspr] = pg.fselfb[js]
    # --- check for particle out of bounds and exchange particles among processors if needed
    top.ns=self.pgroup.ns
#    zpartbnd(self.pgroup,w3d.zmmax,w3d.zmmin,w3d.dz)
    particlegridboundaries3d(top.pgroup,-1)
    top.ns=top.pgroup.ns
    # --- Specify injection of the particles
    top.inject   = 1 #3
    top.injctspc = 1
    top.npinject = 0
    top.zinject  = zinject
    top.ainject  = w3d.xmmax
    top.binject  = w3d.ymmax
    top.apinject = 0.e0
    top.bpinject = 0.e0
    top.lvinject = false  # if false, source conductor input by user
    top.inj_d    = 2.0
    top.inj_f    = 1.0
    top.finject[0][1:]=0.
    top.linj_efromgrid=True
    vbeamfrmtmp = top.vbeamfrm
    injctint(pg)
    top.vbeamfrm = vbeamfrmtmp
    w3d.l_inj_user_particles = true
    w3d.l_inj_user_particles_v = true
    w3d.l_inj_user_particles_dt = true
    w3d.l_inj_zmminmmaxglobal = true
#    installuserparticlesinjection(self.add_boosted_species)
    installbeforestep(self.add_boosted_species)
    if l_deprho:
      self.depos=top.depos.copy()
      top.depos='none'
      installbeforefs(self.add_boosted_rho)
    self.hn = []
    self.hinj = []
    self.hbf = []
   else:
    pg=top.pgroup
    for jspr,js in enumerate(species.jslist):
#      if pg.fselfb[js] != 0.:
#        pg.fselfb[js]=(pg.fselfb[js]-self.betaframe*clight)/(1.-pg.fselfb[js]*self.betaframe/clight)
      il=top.pgroup.ins[js]-1
      iu=il+top.pgroup.nps[js]
      if getn(pgroup=pg,js=js,bcast=0,gather=0)>0: 
        z=getz(pgroup=pg,js=js,bcast=0,gather=0)
      else:
        z=0.
      zmean=globalave(z)
      if getn(pgroup=pg,js=js,bcast=0,gather=0)>0: 
        uzfrm = self.gammaframe*self.betaframe*clight
        tpr =  self.gammaframe*top.time-uzfrm*top.pgroup.zp[il:iu]/clight**2
        top.pgroup.zp[il:iu] = self.gammaframe*top.pgroup.zp[il:iu]-uzfrm*top.time
#        top.pgroup.zp[il:iu]=zmean+(top.pgroup.zp[il:iu]-zmean)/(self.gammaframe*(1.-self.betaframe*self.betabeam_lab))
        vx = getvx(pgroup=pg,js=js,bcast=0,gather=0)
        vy = getvy(pgroup=pg,js=js,bcast=0,gather=0)
        vz = getvz(pgroup=pg,js=js,bcast=0,gather=0)
        fact = 1./(1.-self.betaframe*vz/clight)
        vxpr = vx*fact/self.gammaframe
        vypr = vy*fact/self.gammaframe
        vzpr = (vz-self.betaframe*clight)*fact
        top.pgroup.xp[il:iu] = top.pgroup.xp[il:iu] - tpr*vxpr
        top.pgroup.yp[il:iu] = top.pgroup.yp[il:iu] - tpr*vypr
        top.pgroup.zp[il:iu] = top.pgroup.zp[il:iu] - tpr*vzpr
        gammapr = 1./sqrt(1.-(vxpr*vxpr+vypr*vypr+vzpr*vzpr)/clight**2)
        top.pgroup.uxp[il:iu]=vxpr*gammapr
        top.pgroup.uyp[il:iu]=vypr*gammapr
        top.pgroup.uzp[il:iu]=vzpr*gammapr
        top.pgroup.gaminv[il:iu]=1./gammapr
        if top.uxoldpid>0:top.pgroup.pid[il:iu,top.uxoldpid-1]=top.pgroup.uxp[il:iu]
        if top.uyoldpid>0:top.pgroup.pid[il:iu,top.uyoldpid-1]=top.pgroup.uyp[il:iu]
        if top.uzoldpid>0:top.pgroup.pid[il:iu,top.uzoldpid-1]=top.pgroup.uzp[il:iu]
   if not lallindomain:particleboundaries3d(top.pgroup,-1,False)
   print 'exit boost',top.pgroup.nps
   
  def add_boosted_species(self):
    if self.pgroup.npid != top.pgroup.npid:
      self.pgroup.npid = top.pgroup.npid
      self.pgroup.gchange()
    for js in range(self.pgroup.ns):
     if self.pgroup.nps[js]>0:
      il=self.pgroup.ins[js]-1
      iu=il+self.pgroup.nps[js]
#      self.pgroup.xp[il:iu]+=top.dt*getvx(pgroup=self.pgroup,js=js,bcast=0,gather=0)
#      self.pgroup.yp[il:iu]+=top.dt*getvy(pgroup=self.pgroup,js=js,bcast=0,gather=0)
#      self.pgroup.zp[il:iu]+=top.dt*self.pgroup.uzp[il:iu]*self.pgroup.gaminv[il:iu] # WARNING: this can cause particles to get out of bounds
      if top.zoldpid>0:self.pgroup.pid[il:iu,top.zoldpid-1]=self.pgroup.zp[il:iu].copy()
      self.pgroup.zp[il:iu]+=top.dt*top.vbeam
    if all(self.pgroup.nps==0):
      w3d.npgrp = 0
      gchange("Setpwork3d")
    nps = parallelsum(self.pgroup.nps)
    top.finject[0][1:]=0.
    if all(nps==0):
      top.inject=0
    for js in range(self.pgroup.ns):
     if self.pgroup.nps[js]>0:
      il=self.pgroup.ins[js]-1
      iu=il+self.pgroup.nps[js]
      top.zinject=self.zinject-top.time*self.betaframe*clight
      ii=compress(self.pgroup.zp[il:iu]>top.zinject,il+arange(getn(pgroup=self.pgroup,js=js,bcast=0,gather=0)))
      w3d.npgrp = len(ii)
      w3d.npidgrp = top.npid
      gchange("Setpwork3d")
      if len(ii)>0:
        gi=take(self.pgroup.gaminv,ii)
        vz = take(self.pgroup.uzp,ii)*gi
        w3d.xt = take(self.pgroup.xp,ii).copy()
        w3d.yt = take(self.pgroup.yp,ii).copy()
        w3d.uxt = take(self.pgroup.uxp,ii)*gi
        w3d.uyt = take(self.pgroup.uyp,ii)*gi
        w3d.uzt = vz
#        w3d.bpt = (take(self.pgroup.zp,ii)-top.zinject)/vz
        w3d.bpt = (take(self.pgroup.zp,ii)-top.zinject)/top.vbeam
        for ipid in range(top.npid):
          w3d.pidt[:,ipid] = take(self.pgroup.pid[:,ipid],ii)
#        gi=getgaminv(pgroup=self.pgroup,js=js,bcast=0,gather=0)
        put(self.pgroup.gaminv,ii,0.)     
        npo = self.pgroup.nps[0]
        processlostpart(self.pgroup,js+1,top.clearlostpart,top.time+top.dt*self.pgroup.ndts[js],top.zbeam)
    self.hn.append(getn())
    self.hinj.append(globalsum(w3d.npgrp))
    self.hbf.append(globalsum(self.pgroup.nps[0]))
        
  def pln(self):
    pla(self.hn)
    pla(self.hinj,color=red)
    pla(self.hbf,color=blue)
    pla(cumsum(self.hinj),color=green)

  def add_boosted_speciesold(self):
    for js in range(self.pgroup.ns):
     if self.pgroup.nps[js]>0:
      il=self.pgroup.ins[js]-1
      iu=il+self.pgroup.nps[js]
#      self.pgroup.xp[il:iu]+=top.dt*getvx(pgroup=self.pgroup,js=js,bcast=0,gather=0)
#      self.pgroup.yp[il:iu]+=top.dt*getvy(pgroup=self.pgroup,js=js,bcast=0,gather=0)
#      self.pgroup.zp[il:iu]+=top.dt*getvz(pgroup=self.pgroup,js=js,bcast=0,gather=0) # WARNING: this can cause particles to get out of bounds
      self.pgroup.zp[il:iu]+=top.dt*top.vbeam
    for js in range(self.pgroup.ns):
     if self.pgroup.nps[js]>0:
      il=self.pgroup.ins[js]-1
      iu=il+self.pgroup.nps[js]
      ii=compress(self.pgroup.zp[il:iu]>self.zinject-top.time*self.betaframe*clight,il+arange(getn(pgroup=self.pgroup,js=js,bcast=0,gather=0)))
      if len(ii)>0:
        if self.pgroup.npid>0:
          pid=take(self.pgroup.pid,ii,0)
        else:
          pid=0.
        self.species.addpart(x=take(self.pgroup.xp,ii),
                           y=take(self.pgroup.yp,ii),
                           z=take(self.pgroup.zp,ii),
                           vx=take(self.pgroup.uxp,ii),
                           vy=take(self.pgroup.uyp,ii),
                           vz=take(self.pgroup.uzp,ii),
                           gi=take(self.pgroup.gaminv,ii),
                           pid=pid,
                           lmomentum=true,lallindomain=true)
        gi=getgaminv(pgroup=self.pgroup,js=js,bcast=0,gather=0)
        put(self.pgroup.gaminv,ii,0.)
        processlostpart(self.pgroup,js+1,top.clearlostpart,top.time+top.dt*self.pgroup.ndts[js],top.zbeam)

  def add_boosted_rho(self):
    if self.pgroup.npid != top.pgroup.npid:
      self.pgroup.npid = top.pgroup.npid
      self.pgroup.gchange()
#    if rstrip(top.depos.tostring())=='none': return
#    w3d.lbeforelr=0
    fs=getregisteredsolver()
    pg=top.pgroup
    top.depos=self.depos
    particleboundaries3d(self.pgroup,-1,False)
    fs.loadrho(pgroups=[self.pgroup,top.pgroup])
    top.depos='none'
#    w3d.lbeforelr=1

  def add_boosted_rho_old(self):
    if rstrip(top.depos.tostring())=='none': return
    w3d.lbeforelr=0
    if 1:#getn(pgroup=self.pgroup)>0:
      fs=getregisteredsolver()
      pg=top.pgroup
#      fs.zerosourcep()
#      top.laccumulate_rho=true
      top.depos=self.depos
      fs.loadrho(pgroups=[top.pgroup,self.pgroup])
#      top.laccumulate_rho=false
      self.depos=top.depos.copy()
      top.depos='none'
      fs.aftersetsourcep()
#      fs.aftersetsourcep(lzero=1)
    w3d.lbeforelr=1

  def getn(self):
    n1 = self.species.getn()
    n2 = getn(pgroup=self.pgroup)
    return n1+n2
    
  def getx(self,**kw):
    n1 = self.species.getn(**kw)
    if n1>0:
      z1 = self.species.getx(**kw)
    else:
      z1 = array([])
    n2 = getn(pgroup=self.pgroup,**kw)
    if n2>0:
      z2 = getx(pgroup=self.pgroup,**kw)
    else:
      z2 = array([])
    return concatenate([z1,z2])
    
  def gety(self,**kw):
    n1 = self.species.getn(**kw)
    if n1>0:
      z1 = self.species.gety(**kw)
    else:
      z1 = array([])
    n2 = getn(pgroup=self.pgroup,**kw)
    if n2>0:
      z2 = gety(pgroup=self.pgroup,**kw)
    else:
      z2 = array([])
    return concatenate([z1,z2])
    
  def getz(self,**kw):
    n1 = self.species.getn(**kw)
    if n1>0:
      z1 = self.species.getz(**kw)
    else:
      z1 = array([])
    n2 = getn(pgroup=self.pgroup,**kw)
    if n2>0:
      z2 = getz(pgroup=self.pgroup,**kw)
    else:
      z2 = array([])
    return concatenate([z1,z2])
    
  def getgaminv(self,**kw):
    n1 = self.species.getn(**kw)
    if n1>0:
      z1 = self.species.getgaminv(**kw)
    else:
      z1 = array([])
    n2 = getn(pgroup=self.pgroup,**kw)
    if n2>0:
      z2 = getgaminv(pgroup=self.pgroup,**kw)
    else:
      z2 = array([])
    return concatenate([z1,z2])
        
  def getux(self,**kw):
    n1 = self.species.getn(**kw)
    if n1>0:
      z1 = self.species.getux(**kw)
    else:
      z1 = array([])
    n2 = getn(pgroup=self.pgroup,**kw)
    if n2>0:
      z2 = getux(pgroup=self.pgroup,**kw)
    else:
      z2 = array([])
    return concatenate([z1,z2])
    
  def getuy(self,**kw):
    n1 = self.species.getn(**kw)
    if n1>0:
      z1 = self.species.getuy(**kw)
    else:
      z1 = array([])
    n2 = getn(pgroup=self.pgroup,**kw)
    if n2>0:
      z2 = getuy(pgroup=self.pgroup,**kw)
    else:
      z2 = array([])
    return concatenate([z1,z2])
    
  def getuz(self,**kw):
    n1 = self.species.getn(**kw)
    if n1>0:
      z1 = self.species.getuz(**kw)
    else:
      z1 = array([])
    n2 = getn(pgroup=self.pgroup,**kw)
    if n2>0:
      z2 = getuz(pgroup=self.pgroup,**kw)
    else:
      z2 = array([])
    return concatenate([z1,z2])
    
  def getvx(self,**kw):
    n1 = self.species.getn(**kw)
    if n1>0:
      z1 = self.species.getvx(**kw)
    else:
      z1 = array([])
    n2 = getn(pgroup=self.pgroup,**kw)
    if n2>0:
      z2 = getvx(pgroup=self.pgroup,**kw)
    else:
      z2 = array([])
    return concatenate([z1,z2])
    
  def getvy(self,**kw):
    n1 = self.species.getn(**kw)
    if n1>0:
      z1 = self.species.getvy(**kw)
    else:
      z1 = array([])
    n2 = getn(pgroup=self.pgroup,**kw)
    if n2>0:
      z2 = getvy(pgroup=self.pgroup,**kw)
    else:
      z2 = array([])
    return concatenate([z1,z2])
    
  def getvz(self,**kw):
    n1 = self.species.getn(**kw)
    if n1>0:
      z1 = self.species.getvz(**kw)
    else:
      z1 = array([])
    n2 = getn(pgroup=self.pgroup,**kw)
    if n2>0:
      z2 = getvz(pgroup=self.pgroup,**kw)
    else:
      z2 = array([])
    return concatenate([z1,z2])

  def get_density(self,xmin=None,xmax=None,nx=None,ymin=None,ymax=None,ny=None,zmin=None,zmax=None,
                    nz=None,lost=0,charge=0,dens=None,l_minmax_grid=true,l_dividebyvolume=1,l4symtry=None,l2symtry=None):
        if l_minmax_grid:
            if xmin is None:xmin=w3d.xmmin
            if xmax is None:xmax=w3d.xmmax
            if ymin is None:ymin=w3d.ymmin
            if ymax is None:ymax=w3d.ymmax
            if zmin is None:zmin=w3d.zmmin+top.zgrid
            if zmax is None:zmax=w3d.zmmax+top.zgrid
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
            if w3d.solvergeom is w3d.XYgeom:
                density = fzeros([nx+1,ny+1],'d')
                densityc = fzeros([nx+1,ny+1],'d')
            else:
                if w3d.solvergeom in [w3d.Zgeom]:
                    density = fzeros([nz+1],'d')
                    densityc = fzeros([nz+1],'d')
                elif w3d.solvergeom in [w3d.XZgeom,w3d.RZgeom]:
                    density = fzeros([nx+1,nz+1],'d')
                    densityc = fzeros([nx+1,nz+1],'d')
                else:
                    density = fzeros([nx+1,ny+1,nz+1],'d')
                    densityc = fzeros([nx+1,ny+1,nz+1],'d')
        else:
            if w3d.solvergeom is w3d.XYgeom:
                nx = shape(dens)[0]-1
                ny = shape(dens)[1]-1
            else:
                if w3d.solvergeom in [w3d.Zgeom]:
                    nz = shape(dens)[0]-1
                if w3d.solvergeom in [w3d.XZgeom,w3d.RZgeom]:
                    nx = shape(dens)[0]-1
                    nz = shape(dens)[1]-1
                else:
                    nx = shape(dens)[0]-1
                    ny = shape(dens)[1]-1
                    nz = shape(dens)[2]-1
            density = dens
            densityc = 0.*dens

        np=0
        for pgroup in [top.pgroup,self.pgroup]:
          for js in self.species.jslist:
            np+=getn(pgroup=pgroup,js=js)
        if np == 0:
            if dens is None:
                return density
            else:
                return
        for pgroup in [top.pgroup,self.pgroup]:
          for js in self.species.jslist:
            x=getx(pgroup=pgroup,js=js,lost=lost,gather=0)
            y=gety(pgroup=pgroup,js=js,lost=lost,gather=0)
            z=getz(pgroup=pgroup,js=js,lost=lost,gather=0)
            if w3d.solvergeom==w3d.RZgeom:x=sqrt(x*x+y*y)
            np=shape(x)[0]
            if np > 0:
                if top.wpid == 0:
                    w=self.pgroup.sw[js]*ones(np,'d')
                else:
                    w=self.pgroup.sw[js]*getpid(pgroup=pgroup,js=js,id=top.wpid-1,gather=0)
                if charge:w*=self.pgroup.sq[js]
                if w3d.solvergeom is w3d.Zgeom:
                    deposgrid1d(1,np,z,w,nz,density,densityc,zmin,zmax)
                elif w3d.solvergeom is w3d.XYgeom:
                    deposgrid2d(1,np,x,y,w,nx,ny,density,densityc,xmin,xmax,ymin,ymax)
                elif w3d.solvergeom in [w3d.XZgeom,w3d.RZgeom]:
                    deposgrid2d(1,np,x,z,w,nx,nz,density,densityc,xmin,xmax,zmin,zmax)
                else:
                    deposgrid3d(1,np,x,y,z,w,nx,ny,nz,density,densityc,xmin,xmax,ymin,ymax,zmin,zmax)
        if w3d.solvergeom is w3d.Zgeom:
            if l_dividebyvolume:
                density*=nz/(zmax-zmin)
        elif w3d.solvergeom is w3d.XYgeom:
            if l_dividebyvolume:
                density*=nx*ny/((xmax-xmin)*(ymax-ymin))
                if l4symtry:
                    density[0,:] *= 2
                    density[:,0] *= 2
                if l2symtry:
                    density[:,0] *= 2
                if w3d.boundxy is periodic:
                    density[0,:] += density[-1,:]; density[-1,:]=density[0,:]
                    density[:,0] += density[:,-1]; density[:,-1]=density[:,0]
        else:
            if w3d.solvergeom in [w3d.XZgeom,w3d.RZgeom]:
                if l_dividebyvolume:
                    density*=nx*nz/((xmax-xmin)*(zmax-zmin))
                    if l4symtry:
                        density[0,:] *= 2
                    if w3d.boundxy is periodic:
                        density[0,:] += density[-1,:]; density[-1,:]=density[0,:]
                    if w3d.bound0 is periodic:
                        density[:,0] += density[:,-1]; density[:,-1]=density[:,0]
                    if w3d.solvergeom==w3d.RZgeom:
                        dr = (xmax-xmin)/nx
                        r = arange(nx+1)*dr
                        for j in range(1,nx+1):
                            density[j,:] /= 2.*pi*r[j]
                        density[0,:] /= pi*dr/2
            else:
                if l_dividebyvolume:
                    density*=nx*ny*nz/((xmax-xmin)*(ymax-ymin)*(zmax-zmin))
                    if l4symtry:
                        density[0,:,:] *= 2
                        density[:,0,:] *= 2
                    if l2symtry:
                        density[:,0,:] *= 2
                    if w3d.boundxy is periodic:
                        density[0,:,:] += density[-1,:,:]; density[-1,:,:]=density[0,:,:]
                        density[:,0,:] += density[:,-1,:]; density[:,-1,:]=density[:,0,:]
                    if w3d.bound0 is periodic:
                        density[:,:,0] += density[:,:,-1]; density[:,:,-1]=density[:,:,0]
        density[...] = parallelsum(density)
        if dens is None: return density
      
