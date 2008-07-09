from warp import *

class Boosted_Frame(object):
  def __init__(self,gammaframe,direction=1.):
    self.gammaframe=gammaframe
    self.betaframe  = direction*sqrt(1.-1./self.gammaframe**2)
    self.betabeam_lab=top.vbeamfrm/clight
    top.vbeamfrm=clight*(self.betabeam_lab-self.betaframe)/(1.-self.betabeam_lab*self.betaframe)
    top.vbeam=top.vbeamfrm
    top.gammabar=1./sqrt(1.-(top.vbeam/clight)**2)
    top.fselfb[...]=(top.fselfb[...]-self.betaframe*clight)/(1.-top.fselfb[...]*self.betaframe/clight)

  def boost(self,species,zadd=0.,tinit=0.,l_inject_plane=1):
   if l_inject_plane:
    pg = top.pgroup
    self.species=species
    self.zadd=zadd
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
        z=0.
       zmean=globalave(z)
       if getn(pgroup=pg,js=js,bcast=0,gather=0)>0: 
        gaminvbeam_lab = getgaminv(pgroup=pg,js=js,bcast=0,gather=0)
        betabeam_lab  = sqrt(1.-gaminvbeam_lab**2)
        betabeam_frame = (betabeam_lab-self.betaframe)/(1.-betabeam_lab*self.betaframe)
        gammabeam_frame  = 1./sqrt(1.-betabeam_frame**2)
        z=z-zmean
        tpr = self.gammaframe*self.betaframe*z/clight-tinit
        vx = getvx(pgroup=pg,js=js,bcast=0,gather=0)
        vy = getvy(pgroup=pg,js=js,bcast=0,gather=0)
        vz = getvz(pgroup=pg,js=js,bcast=0,gather=0)
        print 'vz',min(vz),max(vz)
        fact = 1./(1.-self.betaframe*vz/clight)
        vxpr = vx*fact/self.gammaframe
        vypr = vy*fact/self.gammaframe
        vzpr = (vz-self.betaframe*clight)*fact
        print 'vzpr',min(vzpr),max(vzpr)
        gammapr = 1./sqrt(1.-(vxpr*vxpr+vypr*vypr+vzpr*vzpr)/clight**2)
        self.pgroup.uxp[ilpr:iupr]=vxpr*gammapr
        self.pgroup.uyp[ilpr:iupr]=vypr*gammapr
        self.pgroup.uzp[ilpr:iupr]=vzpr*gammapr
        self.pgroup.gaminv[ilpr:iupr]=1./gammapr
        self.pgroup.xp[ilpr:iupr] = getx(pgroup=pg,js=js,bcast=0,gather=0)#+tpr*vxpr
        self.pgroup.yp[ilpr:iupr] = gety(pgroup=pg,js=js,bcast=0,gather=0)#+tpr*vypr
        self.pgroup.zp[ilpr:iupr] = self.gammaframe*z    +tpr*vzpr
        if pg.npid>0:self.pgroup.pid[ilpr:iupr,:] = getpid(pgroup=pg,js=js,bcast=0,gather=0,id=-1)
        if top.uxoldpid>0:self.pgroup.pid[ilpr:iupr,top.uxoldpid-1]=self.pgroup.uxp[ilpr:iupr]
        if top.uyoldpid>0:self.pgroup.pid[ilpr:iupr,top.uyoldpid-1]=self.pgroup.uyp[ilpr:iupr]
        if top.uzoldpid>0:self.pgroup.pid[ilpr:iupr,top.uzoldpid-1]=self.pgroup.uzp[ilpr:iupr]
      pg.nps[js]=0
      pg.fselfb[js]=(pg.fselfb[js]-self.betaframe*clight)/(1.-pg.fselfb[js]*self.betaframe/clight)
      self.pgroup.fselfb[jspr] = pg.fselfb[js]
    # --- check for particle out of bounds and exchange particles among processors if needed
    top.ns=self.pgroup.ns
    zpartbnd(self.pgroup,w3d.zmmax,w3d.zmmin,w3d.dz)
    top.ns=top.pgroup.ns
    self.depos=top.depos.copy()
    installbeforeloadrho(self.add_boosted_species)
    installbeforefs(self.add_boosted_rho)
   else:
    pg=top.pgroup
    for jspr,js in enumerate(species.jslist):
      pg.fselfb[js]=(pg.fselfb[js]-self.betaframe*clight)/(1.-pg.fselfb[js]*self.betaframe/clight)
      il=top.pgroup.ins[js]-1
      iu=il+top.pgroup.nps[js]
      if getn(pgroup=pg,js=js,bcast=0,gather=0)>0: 
        z=getz(pgroup=pg,js=js,bcast=0,gather=0)
      else:
        z=0.
      zmean=globalave(z)
      if getn(pgroup=pg,js=js,bcast=0,gather=0)>0: 
        top.pgroup.zp[il:iu]=zmean+(top.pgroup.zp[il:iu]-zmean)/(self.gammaframe*(1.-self.betaframe*self.betabeam_lab))
        vx = getvx(pgroup=pg,js=js,bcast=0,gather=0)
        vy = getvy(pgroup=pg,js=js,bcast=0,gather=0)
        vz = getvz(pgroup=pg,js=js,bcast=0,gather=0)
        fact = 1./(1.-self.betaframe*vz/clight)
        vxpr = vx*fact/self.gammaframe
        vypr = vy*fact/self.gammaframe
        vzpr = (vz-self.betaframe*clight)*fact
        gammapr = 1./sqrt(1.-(vxpr*vxpr+vypr*vypr+vzpr*vzpr)/clight**2)
        top.pgroup.uxp[il:iu]=vxpr*gammapr
        top.pgroup.uyp[il:iu]=vypr*gammapr
        top.pgroup.uzp[il:iu]=vzpr*gammapr
        top.pgroup.gaminv[il:iu]=1./gammapr
        if top.uxoldpid>0:top.pgroup.pid[il:iu,top.uxoldpid-1]=top.pgroup.uxp[il:iu]
        if top.uyoldpid>0:top.pgroup.pid[il:iu,top.uyoldpid-1]=top.pgroup.uyp[il:iu]
        if top.uzoldpid>0:top.pgroup.pid[il:iu,top.uzoldpid-1]=top.pgroup.uzp[il:iu]
   
  def add_boosted_species(self):
    for js in range(self.pgroup.ns):
     if self.pgroup.nps[js]>0:
      il=self.pgroup.ins[js]-1
      iu=il+self.pgroup.nps[js]
#      self.pgroup.xp[il:iu]+=top.dt*getvx(pgroup=self.pgroup,js=js,bcast=0,gather=0)
#      self.pgroup.yp[il:iu]+=top.dt*getvy(pgroup=self.pgroup,js=js,bcast=0,gather=0)
      self.pgroup.zp[il:iu]+=top.dt*getvz(pgroup=self.pgroup,js=js,bcast=0,gather=0)
    for js in range(self.pgroup.ns):
     if self.pgroup.nps[js]>0:
      il=self.pgroup.ins[js]-1
      iu=il+self.pgroup.nps[js]
      ii=compress(self.pgroup.zp[il:iu]>self.zadd-top.time*self.betaframe*clight,il+arange(getn(pgroup=self.pgroup,js=js,bcast=0,gather=0)))
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
