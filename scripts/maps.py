from warp import *

class Maps:
  def __init__(self,nux,nuy,C,nstations=1,nparpgrp=top.nparpgrp,l_mode=1,l_verbose=0,betax=None,betay=None):
     self.nux = nux
     self.nuy = nuy
     self.C   = C
     self.nstations = nstations
     # --- computes number of steps
     self.sigmax = 2.*pi*self.nux/self.nstations
     self.sigmay = 2.*pi*self.nuy/self.nstations
     if betax is None:
       self.betax = self.C/(2.*pi*self.nux)
     else:
       self.betax = betax
     if betay is None:
       self.betay = self.C/(2.*pi*self.nuy)
     else:
       self.betay = betay
     self.Mtx = array([[ cos(self.sigmax)           , self.betax*sin(self.sigmax)], \
                       [-sin(self.sigmax)/self.betax,            cos(self.sigmax)]])
     self.Mty = array([[ cos(self.sigmay)           , self.betay*sin(self.sigmay)], \
                       [-sin(self.sigmay)/self.betay,            cos(self.sigmay)]])
     self.l_verbose=l_verbose
     self.l_mode=l_mode
     self.nparpgrp = nparpgrp

  def apply_transfer_map(self,sp):
    js=sp.jslist[0]
    pg=top.pgroup
    ng = 1+pg.nps[js]/self.nparpgrp
    if pg.nps[js]==0:return
    for ig in range(ng):
      il = pg.ins[js]-1+self.nparpgrp*ig
      iu = min(il+self.nparpgrp,pg.ins[js]-1+pg.nps[js])
      np = iu-il
      apply_simple_map(np,
                       pg.xp[il:iu],
                       pg.yp[il:iu],
                       pg.uxp[il:iu],
                       pg.uyp[il:iu],
                       pg.gaminv[il:iu],
                       self.Mtx,
                       self.Mty,
                       top.vbeam)
#      xp = pg.xp[il:iu].copy()
#      yp = pg.yp[il:iu].copy()
#      scf = pg.gaminv[il:iu]/top.vbeam
#      pg.xp [il:iu] = self.Mtx[0,0]*xp     + self.Mtx[0,1]*pg.uxp[il:iu]*scf
#      pg.uxp[il:iu] = self.Mtx[1,0]*xp/scf + self.Mtx[1,1]*pg.uxp[il:iu]
#      pg.yp [il:iu] = self.Mty[0,0]*yp     + self.Mty[0,1]*pg.uyp[il:iu]*scf
#      pg.uyp[il:iu] = self.Mty[1,0]*yp/scf + self.Mty[1,1]*pg.uyp[il:iu]

  def apply_space_charge_kick(self,sp):
    js=sp.jslist[0]
    pg=top.pgroup
    ng = 1+pg.nps[js]/self.nparpgrp
    if pg.nps[js]==0:return
    for ig in range(ng):
      il = pg.ins[js]-1+self.nparpgrp*ig
      iu = min(il+self.nparpgrp,pg.ins[js]-1+pg.nps[js])
      np = iu-il
      if self.l_verbose:print 'fetche3d'
      fselfb=top.fselfb.copy()
      top.fselfb[...]=0.
      top.pgroup.fselfb[...]=0.
      fsolver=getregisteredsolver()
      if fsolver is not None:
        efetchsave = top.efetch[js]+0
        top.efetch[js] = 1
      fetche3d(top.pgroup,il+1,np,js+1)
      if fsolver is not None:
        top.efetch[js] = efetchsave
      top.fselfb[...]=fselfb
      top.pgroup.fselfb[...]=fselfb
      if self.l_mode==2:
        lzeros = where((pg.zp[il:iu]<0.5*w3d.zmmin) | (pg.zp[il:iu]>0.5*w3d.zmmax),1,0)
        pg.ex[il:iu] = where(lzeros,0.,pg.ex[il:iu])
        pg.ey[il:iu] = where(lzeros,0.,pg.ey[il:iu])
        pg.ez[il:iu] = where(lzeros,0.,pg.ez[il:iu])
      if self.l_verbose:print 'epush beam'
      epush3d(np,pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],
              pg.ex[il:iu],pg.ey[il:iu],pg.ez[il:iu],pg.sq[js],pg.sm[js],top.dt)
      gammaadv(np,pg.gaminv[il:iu],pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],
               top.gamadv,top.lrelativ)

  def apply_bnd_conditions(self,sp):
    js=sp.jslist[0]
    pg=top.pgroup
    ng = 1+pg.nps[js]/self.nparpgrp
    if pg.nps[js]==0:return
    for ig in range(ng):
      il = pg.ins[js]-1+self.nparpgrp*ig
      iu = min(il+self.nparpgrp,pg.ins[js]-1+pg.nps[js])
      np = iu-il
      if self.l_verbose:print 'stckxy3d beam'
      zpartbndwithdata(np,pg.zp[il:iu],pg.uzp[il:iu],pg.gaminv[il:iu],
                       w3d.zmmaxlocal,w3d.zmminlocal,w3d.dz,top.zgrid)
      stckxy3d(np,pg.xp[il:iu],w3d.xmmax,w3d.xmmin,w3d.dx,
                  pg.yp[il:iu],w3d.ymmax,w3d.ymmin,w3d.dy,
                  pg.zp[il:iu],w3d.zmminlocal,w3d.dz,
                  pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],pg.gaminv[il:iu],
                  top.zgrid,top.zbeam,w3d.l2symtry,w3d.l4symtry,top.pboundxy,true)
    processlostpart(top.pgroup,js+1,top.clearlostpart,top.time+top.dt*top.pgroup.ndts[js],top.zbeam)

