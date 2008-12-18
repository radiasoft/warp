from warp import *

class Boosted_Frame(object):
  def __init__(self,gammaframe,direction=1.):
    top.boost_gamma=gammaframe
    self.gammaframe=gammaframe
    self.betaframe  = direction*sqrt(1.-1./self.gammaframe**2)
    self.betabeam_lab=top.vbeam/clight
    self.betabeamfrm_lab=top.vbeamfrm/clight
    top.vbeam_lab = top.vbeam
    top.gammabar_lab = top.gammabar
    top.vbeam=clight*(self.betabeam_lab-self.betaframe)/(1.-self.betabeam_lab*self.betaframe)
    top.vbeamfrm=clight*(self.betabeamfrm_lab-self.betaframe)/(1.-self.betabeamfrm_lab*self.betaframe)
    top.gammabar=1./sqrt(1.-(top.vbeam/clight)**2)
    top.fselfb[...]=where(top.fselfb==0.,0.,(top.fselfb[...]-self.betaframe*clight)/(1.-top.fselfb[...]*self.betaframe/clight))

  def boost(self,species,zinject=0.,tinit=0.,l_inject_plane=1):
   print 'enter boost',top.pgroup.nps
   if l_inject_plane:
    pg = top.pgroup
    self.species=species
    self.zinject=zinject
    self.tinit=tinit
    self.pgroup = ParticleGroup()
#    self.pgroup.ns = pg.ns#len(species.jslist)
    self.pgroup.ns = len(species.jslist)
    self.pgroup.npmax = species.getn(bcast=0,gather=0)
    self.pgroup.npid = pg.npid
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
       zmean=globalave(z)
       if getn(pgroup=pg,js=js,bcast=0,gather=0)>0: 
        gaminvbeam_lab = getgaminv(pgroup=pg,js=js,bcast=0,gather=0)
        betabeam_lab  = sqrt(1.-gaminvbeam_lab**2)
        betabeam_frame = (betabeam_lab-self.betaframe)/(1.-betabeam_lab*self.betaframe)
        gammabeam_frame  = 1./sqrt(1.-betabeam_frame**2)
        z=z-zmean
        # --- get data at z=0
        vx = getvx(pgroup=pg,js=js,bcast=0,gather=0)
        vy = getvy(pgroup=pg,js=js,bcast=0,gather=0)
        vz = getvz(pgroup=pg,js=js,bcast=0,gather=0)
        t = z/vz
        x = getx(pgroup=pg,js=js,bcast=0,gather=0)#-t*vx
        y = gety(pgroup=pg,js=js,bcast=0,gather=0)#-t*vy
        # --- get data in boosted frame
        tpr = -self.gammaframe*t
        zpr = self.gammaframe*self.betaframe*clight*t
        if top.boost_z0==0.:
          top.boost_z0 = -globalmax(zpr)
        fact = 1./(1.-self.betaframe*vz/clight)
        vxpr = vx*fact/self.gammaframe
        vypr = vy*fact/self.gammaframe
        vzpr = (vz-self.betaframe*clight)*fact
        # --- get data at t=0 in boosted frame
        zpr = zpr - vzpr*tpr
        # --- make sure that z<=0
        zpr += top.boost_z0 
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
      if pg.fselfb[js]<>0.:
        pg.fselfb[js]=(pg.fselfb[js]-self.betaframe*clight)/(1.-pg.fselfb[js]*self.betaframe/clight)
      self.pgroup.fselfb[jspr] = pg.fselfb[js]
    # --- check for particle out of bounds and exchange particles among processors if needed
    top.ns=self.pgroup.ns
#    zpartbnd(self.pgroup,w3d.zmmax,w3d.zmmin,w3d.dz)
    particlegridboundaries3d(top.pgroup,-1)
    top.ns=top.pgroup.ns
    self.depos=top.depos.copy()
    top.depos='none'
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
    vbeamfrmtmp = top.vbeamfrm
    injctint(pg)
    top.vbeamfrm = vbeamfrmtmp
    w3d.l_inj_user_particles = true
    w3d.l_inj_user_particles_v = true
    w3d.l_inj_user_particles_dt = true
    w3d.l_inj_zmminmmaxglobal = true
#    installuserparticlesinjection(self.add_boosted_species)
    installbeforestep(self.add_boosted_species)
    installbeforefs(self.add_boosted_rho)
   else:
    pg=top.pgroup
    for jspr,js in enumerate(species.jslist):
      if pg.fselfb[js]<>0.:
        pg.fselfb[js]=(pg.fselfb[js]-self.betaframe*clight)/(1.-pg.fselfb[js]*self.betaframe/clight)
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
   particleboundaries3d(top.pgroup,-1,False)
   print 'exit boost',top.pgroup.nps
   
  def add_boosted_species(self):
    for js in range(self.pgroup.ns):
     if self.pgroup.nps[js]>0:
      il=self.pgroup.ins[js]-1
      iu=il+self.pgroup.nps[js]
#      self.pgroup.xp[il:iu]+=top.dt*getvx(pgroup=self.pgroup,js=js,bcast=0,gather=0)
#      self.pgroup.yp[il:iu]+=top.dt*getvy(pgroup=self.pgroup,js=js,bcast=0,gather=0)
#      self.pgroup.zp[il:iu]+=top.dt*getvz(pgroup=self.pgroup,js=js,bcast=0,gather=0) # WARNING: this can cause particles to get out of bounds
      self.pgroup.zp[il:iu]+=top.dt*top.vbeam
    if all(self.pgroup.nps==0):
      w3d.npgrp = 0
      gchange("Setpwork3d")
    nps = parallelsum(self.pgroup.nps)
    if all(nps==0):
      top.inject=0
    for js in range(self.pgroup.ns):
     if self.pgroup.nps[js]>0:
      il=self.pgroup.ins[js]-1
      iu=il+self.pgroup.nps[js]
      ii=compress(self.pgroup.zp[il:iu]>self.zinject-top.time*self.betaframe*clight,il+arange(getn(pgroup=self.pgroup,js=js,bcast=0,gather=0)))
      w3d.npgrp = len(ii)
      gchange("Setpwork3d")
      top.zinject=self.zinject-top.time*self.betaframe*clight
      if len(ii)>0:
        gi=take(self.pgroup.gaminv,ii)
        vz = take(self.pgroup.uzp,ii)*gi
        w3d.xt = take(self.pgroup.xp,ii)
        w3d.yt = take(self.pgroup.yp,ii)
        w3d.uxt = take(self.pgroup.uxp,ii)*gi
        w3d.uyt = take(self.pgroup.uyp,ii)*gi
        w3d.uzt = vz
        w3d.bpt = (take(self.pgroup.zp,ii)-top.zinject)/vz
        gi=getgaminv(pgroup=self.pgroup,js=js,bcast=0,gather=0)
        put(self.pgroup.gaminv,ii,0.)     
        npo = self.pgroup.nps[0]
        processlostpart(self.pgroup,js+1,top.clearlostpart,top.time+top.dt*self.pgroup.ndts[js],top.zbeam)

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
#    if rstrip(top.depos.tostring())=='none': return
#    w3d.lbeforelr=0
    fs=getregisteredsolver()
    pg=top.pgroup
    top.depos=self.depos
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
