from warp import *

class Maps:
  def __init__(self,station1,station2,nparpgrp=top.nparpgrp,l_mode=1,l_verbose=0, xtune=0., 
               ytune=0., eta=0., harm_num=0., ring_circum=0., sync_tune=0., cross_zero = 0):
     #define required lattice parameters 
     self.bx1  = station1["betax"]
     self.bx2  = station2["betax"]
     self.ax1  = station1["alphax"]
     self.ax2  = station2["alphax"]
     self.dx1  = station1["dispx"]
     self.dx2  = station2["dispx"]
     self.dpx1 = station1["disppx"]
     self.dpx2 = station2["disppx"]
     self.by1  = station1["betay"]
     self.by2  = station2["betay"]
     self.ay1  = station1["alphay"]
     self.ay2  = station2["alphay"]
     self.dy1  = station1["dispy"]
     self.dy2  = station2["dispy"]
     self.dpy1 = station1["disppy"]
     self.dpy2 = station2["disppy"]     
     self.Qx   = xtune 
     self.Qy   = ytune 
     self.eta  = eta
     self.h    = harm_num
     self.C    = ring_circum 
     self.mus  = sync_tune 
     self.l_verbose = l_verbose
     self.l_mode = l_mode
     self.nparpgrp = nparpgrp
     self.ex = zeros(self.nparpgrp,'d')
     self.ey = zeros(self.nparpgrp,'d')
     self.ez = zeros(self.nparpgrp,'d')
     self.RF1 =  station1["RFkick"]
     self.RF2 =  station2["RFkick"]
     self.ec_flag = station1["ec_flag"] 
     if cross_zero == 0 :
        self.L = station2["s_dist"] - station1["s_dist"]    
        self.phx  = (station2["phasex"] - station1["phasex"])*2*pi
        self.phy  = (station2["phasey"] - station1["phasey"])*2*pi
     else:
        self.L = station2["s_dist"] - station1["s_dist"] + ring_circum
        self.phx  = (station2["phasex"] - station1["phasex"] + xtune)*2*pi
        self.phy  = (station2["phasey"] - station1["phasey"] + ytune)*2*pi
     if self.RF1 == 0 or self.RF2 ==  0:               
        Mx11 = sqrt(self.bx2/self.bx1)*(cos(self.phx)+self.ax1*sin(self.phx))
        Mx12 = sqrt(self.bx2*self.bx1)*sin(self.phx)
        Mx21 = -(1/sqrt(self.bx1*self.bx2)*((self.ax2-self.ax1)*cos(self.phx) + \
                (1+self.ax1*self.ax2)*sin(self.phx)))
        Mx22 = sqrt(self.bx1/self.bx2)*(cos(self.phx) - self.ax2*sin(self.phx)) 
        My11 = sqrt(self.by2/self.by1)*(cos(self.phy)+self.ay1*sin(self.phy))
        My12 = sqrt(self.by2*self.by1)*sin(self.phy)
        My21 = -(1/sqrt(self.by1*self.by2)*((self.ay2-self.ay1)*cos(self.phy) + \
                (1+self.ay1*self.ay2)*sin(self.phy)))
        My22 = sqrt(self.by1/self.by2)*(cos(self.phy) - self.ay2*sin(self.phy))
        Mz11 = 1. 
        Mz12 = -self.eta*self.L
        #Mz12 = 0   #testthis
        Mz21 = 0.
        Mz22 = 1.
        Tx13 = self.dx2  - Mx11*self.dx1 - Mx12*self.dpx1 
        Tx23 = self.dpx2 - Mx21*self.dx1 - Mx22*self.dpx1
        Ty13 = self.dy2  - My11*self.dy1 - My12*self.dpy1 
        Ty23 = self.dpy2 - My21*self.dy1 - My22*self.dpy1
         
        Map = array([[Mx11,Mx12,0.,0.,0.,Tx13],\
                     [Mx21,Mx22,0.,0.,0.,Tx23],\
                     [0.,0.,My11,My12,0.,Ty13],\
                     [0.,0.,My21,My22,0.,Ty23],\
                     [0.,0.,0.,0.,Mz11, Mz12 ],\
                     [0.,0.,0.,0.,Mz21, Mz22 ]])         
     else:
    
         Be   = top.vbeam/top.clight
         Mz21    = (2*pi*self.mus)**2/(self.eta*Be*self.C)
         #Mz21 = 0   #testthis  
  
         Map = array([[1,0.,0.,0.,0.,0.], \
                      [0.,1.,0.,0.,0.,0.],\
                      [0.,0.,1.,0.,0.,0.],\
                      [0.,0.,0.,1.,0.,0.],\
                      [0.,0.,0.,0.,1.,0.],\
                      [0.,0.,0.,0.,Mz21,1.]])   
     
     self.Map = Map  

  def apply_transfer_map(self,sp):
    uzb =  top.gammabar*top.vbeam
    js=sp.jslist[0]
    pg=top.pgroup
    ng = 1+pg.nps[js]/self.nparpgrp
    for ig in range(ng):
      il  =  pg.ins[js]-1+self.nparpgrp*ig
      iu  =  min(il+self.nparpgrp,pg.ins[js]-1+pg.nps[js])
      np  =  iu-il
      xp  =  pg.xp[il:iu].copy()
      yp  =  pg.yp[il:iu].copy()
      zp  =  pg.zp[il:iu].copy()
      scf =  pg.gaminv[il:iu]/top.vbeam
      uxp =  pg.uxp[il:iu]*scf
      uyp =  pg.uyp[il:iu]*scf
      uzp =  (pg.uzp[il:iu] - uzb)/uzb
               
      pg.xp[il:iu]  = self.Map[0,0]*xp   + self.Map[0,1]*uxp \
                    + self.Map[0,2]*yp   + self.Map[0,3]*uyp  \
                    + self.Map[0,4]*zp   + self.Map[0,5]*uzp  
                          
      pg.uxp[il:iu] = (self.Map[1,0]*xp   + self.Map[1,1]*uxp \
                      + self.Map[1,2]*yp   + self.Map[1,3]*uyp \
                      + self.Map[1,4]*zp   + self.Map[1,5]*uzp)/scf 
                      
      pg.yp[il:iu]  = self.Map[2,0]*xp   + self.Map[2,1]*uxp \
                    + self.Map[2,2]*yp   + self.Map[2,3]*uyp  \
                    + self.Map[2,4]*zp   + self.Map[2,5]*uzp        
                       
      pg.uyp[il:iu] = (self.Map[3,0]*xp   + self.Map[3,1]*uxp  \
                      + self.Map[3,2]*yp   + self.Map[3,3]*uyp  \
                      + self.Map[3,4]*zp   + self.Map[3,5]*uzp)/scf    
                      
      pg.zp[il:iu]  = self.Map[4,0]*xp   + self.Map[4,1]*uxp \
                    + self.Map[4,2]*yp   + self.Map[4,3]*uyp \
                    + self.Map[4,4]*zp   + self.Map[4,5]*uzp     
                      
      pg.uzp[il:iu] = (self.Map[5,0]*xp   + self.Map[5,1]*uxp  \
                      + self.Map[5,2]*yp   + self.Map[5,3]*uyp \
                      + self.Map[5,4]*zp   + self.Map[5,5]*uzp)*uzb + uzb     
                     
      gammaadv(np,pg.gaminv[il:iu],pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],
               top.gamadv,top.lrelativ)                      



  def apply_space_charge_kick(self,sp):
    if self.RF1 == 0 or self.RF2 ==  0:               
     top.dt = self.L/top.vbeam
     js=sp.jslist[0]
     pg=top.pgroup
     ng = 1+pg.nps[js]/self.nparpgrp
     for ig in range(ng):
       il = pg.ins[js]-1+self.nparpgrp*ig
       iu = min(il+self.nparpgrp,pg.ins[js]-1+pg.nps[js])
       np = iu-il
       if self.l_verbose:print 'fetche3d'
       fetche3dfrompositions(js+1,pg.ndts,np,
                            pg.xp[il:iu],pg.yp[il:iu],pg.zp[il:iu],
                            self.ex[:np],self.ey[:np],self.ez[:np])
       if self.l_mode==2:
         lzeros = where( (pg.zp[il:iu]<0.5*w3d.zmmin) | (pg.zp[il:iu]>0.5*w3d.zmmax) ,1,0)
         self.ex[:np] = where(lzeros,0.,self.ex[:np])
         self.ey[:np] = where(lzeros,0.,self.ey[:np])
         self.ez[:np] = where(lzeros,0.,self.ez[:np])
       if self.l_verbose:print 'epush beam'
       epush3d(np,pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],
              self.ex[:np],self.ey[:np],self.ez[:np],pg.sq[js],pg.sm[js],top.dt)
       gammaadv(np,pg.gaminv[il:iu],pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],
               top.gamadv,top.lrelativ)
     print 'no RF, so epush'  
    else:
      print 'RF, so no epush' 
  
  def apply_bnd_conditions(self,sp):
    js=sp.jslist[0]
    pg=top.pgroup
    ng = 1+pg.nps[js]/self.nparpgrp
    for ig in range(ng):
      il = pg.ins[js]-1+self.nparpgrp*ig
      iu = min(il+self.nparpgrp,pg.ins[js]-1+pg.nps[js])
      np = iu-il
      if self.l_verbose:print 'stckxy3d beam'
      zpartbndwithdata(np,pg.zp[il:iu],pg.uzp[il:iu],pg.gaminv[il:iu],
                       w3d.zmmax,w3d.zmmin,w3d.dz,top.zgrid)
      stckxy3d(np,pg.xp[il:iu],w3d.xmmax,w3d.xmmin,w3d.dx,
                  pg.yp[il:iu],w3d.ymmax,w3d.ymmin,w3d.dy,
                  pg.zp[il:iu],w3d.zmmin,w3d.dz,
                  pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],pg.gaminv[il:iu],
                  top.zgrid,top.zbeam,w3d.l2symtry,w3d.l4symtry,top.pboundxy,true)
    processlostpart(top.pgroup,js+1,top.clearlostpart,top.time+top.dt*top.pgroup.ndts[js],top.zbeam)

class Maps_simple:
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
                       pg.uzp[il:iu],
                       self.Mtx,
                       self.Mty)
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

