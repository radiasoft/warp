from warp import *

class Maps:
  def __init__(self,nux,nuy,nuz,C,nstations=1,z_rms=0.,eta=0.,xchrom=0.,ychrom=0.,
                    alphax=0.,alphay=0.,dispx=0.,dispy=0.,disppx=0.,disppy=0.,
                    nparpgrp=top.nparpgrp,l_mode=1,l_verbose=0,betax=None,betay=None):
     self.nux = nux
     self.nuy = nuy
     self.nuz = nuz
     self.xchrom=xchrom
     self.ychrom=ychrom
     self.C   = C
     self.z_rms = z_rms
     self.nstations = nstations
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
     self.alphax=alphax
     self.dispx=dispx
     self.disppx=disppx
     self.alphay=alphay
     self.dispy=dispy
     self.disppy=disppy
     self.omega0 = 2.*pi*top.vbeam/self.C
     self.omegax = nux*self.omega0
     self.omegay = nuy*self.omega0
     self.omegaz = nuz*self.omega0
     self.eta = eta
     self.l_verbose=l_verbose
     self.l_mode=l_mode
     self.nparpgrp = nparpgrp

  def apply_transfer_map(self,pg,il,iu,dt,l_push_z=true):
    top.pgroup=pg
    np = iu - il   
    ax1=ax2=self.alphax
    bx1=bx2=self.betax
    dx1=dx2=self.dispx
    dpx1=dpx2=self.disppx
    Qx=self.nux
    xchrom=self.xchrom
    phasex=self.omegax*dt
    eta=self.eta
    ay1=ay2=self.alphay
    by1=by2=self.betay
    dy1=dy2=self.dispy
    dpy1=dpy2=self.disppy
    Qy=self.nuy
    ychrom=self.ychrom
    phasey=self.omegay*dt
    eta=self.eta
    omegaz=self.omegaz
    if l_push_z:
      phz=self.omegaz*dt
    else:
      phz=0.
    apply_linear_map(np,pg.xp[il:iu],pg.yp[il:iu],pg.zp[il:iu],
                        pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],pg.gaminv[il:iu],
                        top.vbeam,top.gammabar,
                        ax1,ax2,bx1,bx2,dx1,dx2,
                        dpx1,dpx2,Qx,xchrom,phasex,
                        ay1,ay2,by1,by2,dy1,dy2,
                        dpy1,dpy2,Qy,ychrom,phasey,eta,omegaz,phz)

class Maps_simple:
  def __init__(self,nux,nuy,nuz,C,nstations=1,z_rms=0.,eta=0.,xchrom=0.,ychrom=0.,
                    nparpgrp=top.nparpgrp,l_mode=1,l_verbose=0,betax=None,betay=None):
     self.nux = nux
     self.nuy = nuy
     self.nuz = nuz
     self.xchrom=xchrom
     self.ychrom=ychrom
     self.C   = C
     self.z_rms = z_rms
     self.nstations = nstations
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
     self.omega0 = 2.*pi*top.vbeam/self.C
     self.omegax = nux*self.omega0
     self.omegay = nuy*self.omega0
     self.omegaz = nuz*self.omega0
     self.eta = eta
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
                       w3d.zmmax,w3d.zmmin,w3d.dz,top.zgrid)
      stckxy3d(np,pg.xp[il:iu],w3d.xmmax,w3d.xmmin,w3d.dx,
                  pg.yp[il:iu],w3d.ymmax,w3d.ymmin,w3d.dy,
                  pg.zp[il:iu],w3d.zmminlocal,w3d.dz,
                  pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],pg.gaminv[il:iu],
                  top.zgrid,top.zbeam,w3d.l2symtry,w3d.l4symtry,top.pboundxy,true)
    processlostpart(top.pgroup,js+1,top.clearlostpart,top.time+top.dt*top.pgroup.ndts[js],top.zbeam)

  def apply_synchrotron_motion(self,pg,il,iu,dt):
   if self.omegaz<>0.:
    uzb =  top.gammabar*top.vbeam
    Beta_z = -self.eta/(self.omegaz*top.gammabar)
    Cn = cos(self.omegaz*dt)
    Sn = sin(self.omegaz*dt)
    zp  =  pg.zp[il:iu].copy()
    uzp =  (pg.uzp[il:iu] - uzb)  

    pg.zp[il:iu]  =  Cn*zp        + Beta_z*Sn*uzp
    pg.uzp[il:iu] = -Sn*zp/Beta_z + Cn*uzp + uzb
    np = iu - il               
    gammaadv(np,pg.gaminv[il:iu],pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],
             top.gamadv,top.lrelativ)
   else:
    pg.zp[il:iu]+=(pg.uzp[il:iu]*pg.gaminv[il:iu]-top.vbeam)*dt

# Headtail
# 	     sz[kmain] = csa * szz0                - eta * C/omegas * ssa * dpp0;
#	     dp[kmain] = omegas/eta/C * ssa * szz0 + csa * dpp0;
# if omegas==0:
#        sz[kmain] = szz0 - eta*circ*dpp0;

  def apply_synchrotron_motion_old(self,pg,il,iu,fmult):
    uzb =  top.gammabar*top.vbeam
    Beta_z = self.z_rms/(uzb*self.dPbyP) 
    Cn = cos(2*top.pi*self.nuz*fmult)
    Sn = sin(2*top.pi*self.nuz*fmult)
    zp  =  pg.zp[il:iu].copy()
    uzp =  (pg.uzp[il:iu] - uzb)   
    np = iu - il               
    zp1=pg.zp[il:iu].copy()
    pg.zp[il:iu] =  Cn*zp + Beta_z*Sn*uzp
    pg.uzp[il:iu] = -Sn*zp/Beta_z + Cn*uzp + uzb
    zp2=pg.zp[il:iu].copy()
    if(max(abs(zp2-zp1))>w3d.dz):
      print 'error',max(abs(zp2-zp1)),w3d.dz
    gammaadv(np,pg.gaminv[il:iu],pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],
             top.gamadv,top.lrelativ)
 
class Maps_lattice:
  def __init__(self, Map, length, s_posn, ecflag = 1, name = "", alphax=0., alphay=0.,
                     dispx=0.,dispy=0.,disppx=0.,disppy=0.):
     self.Map = Map
     self.L   = length
     self.s_posn = s_posn
     self.ecflag = ecflag
     self.name = name
     self.alphax = alphax
     self.alphay = alphay
     self.dispx = dispx
     self.disppx = disppx
     self.dispy = dispy
     self.disppy = disppy

  def apply_transfer_map(self,pg,il,iu):
      apply_map(iu-il,
                       pg.xp[il:iu],
                       pg.yp[il:iu],
                       pg.zp[il:iu],
                       pg.uxp[il:iu],
                       pg.uyp[il:iu],
                       pg.uzp[il:iu],
                       pg.gaminv[il:iu],
                       self.Map,
                       top.vbeam,
                       top.gammabar)

class Maps_Parameters:
  def __init__(self,gamma_t = 0,harmn_num = 0,sync_freq = 0,sync_tune = 0,
                    eta = 0,gamma_b = 0,Ergy = 0,alpha = 0, mo = 938280000., 
                    circum=0):
   if gamma_b == 0:
    self.gamma = Ergy/mo 
   else:
    self.gamma = gamma_b 
   if gamma_t == 0:
    self.eta = alpha**2 - 1/self.gamma**2
   else:
    self.eta = 1/gamma_t**2 - 1/self.gamma**2 #slippage factor   
   if sync_tune == 0:
    self.mus = sync_freq*circum/top.clight
   else:
    self_mus = sync_tune 
   self.h = harmn_num
   self.C= circum  
   self.Beta = sqrt(1 - 1/(self.gamma**2))

  def Thin_quad(self,kl):
    M11 = 1.
    M12 = 0. 
    M13 = 0.
    M14 = 0.
    M15 = 0.
    M16 = 0.

    M21 = -kl
    M22 = 1.
    M23 = 0.
    M24 = 0. 
    M25 = 0.
    M26 = 0.

    M31 = 0.
    M32 = 0.
    M33 = 1.
    M34 = 0.
    M35 = 0.
    M36 = 0.

    M41 = 0.
    M42 = 0.
    M43 = kl 
    M44 = 1.
    M45 = 0.
    M46 = 0.

    M51 = 0.
    M52 = 0.
    M53 = 0.
    M54 = 0.
    M55 = 1.
    M56 = 0.

    M61 = 0.
    M62 = 0.
    M63 = 0.
    M64 = 0.
    M65 = 0.
    M66 = 1.

    Map = array([[M11,M12,M13,M14,M15,M16],\
               [M21,M22,M23,M24,M25,M26],\
               [M31,M32,M33,M34,M35,M36],\
               [M41,M42,M43,M44,M45,M46],\
               [M51,M52,M53,M54,M55,M56],\
               [M61,M62,M63,M64,M65,M66]])

    return Map

  def Thick_quad(self,k,L):
    Be   = self.Beta
    ga = self.gamma
    if (k>0. and L>0.) or (k<0. and L<0.):
      k = abs(k)
      L = abs(L)
      # --- focusing quad
      cx = cos(k*L)
      sx = sin(k*L)/k
      cy = cos(k*L) 
      sy = -sin(k*L)/k 
#      cy = cosh(k*L) 
#      sy = -sinh(k*L)/k 
    if (k>0. and L<0.) or (k<0. and L>0.):
      k = abs(k)
      L = abs(L)
      # --- defocusing quad
#      cx = cosh(k*L)
#      sx = -sinh(k*L)/k
      cx = cos(k*L)
      sx = -sin(k*L)/k
      cy = cos(k*L) 
      sy = sin(k*L)/k 

    M11 = cx
    M12 = sx
    M13 = 0.
    M14 = 0.
    M15 = 0.
    M16 = 0. 

    M21 = -k**2*sx
    M22 = cx
    M23 = 0.
    M24 = 0.
    M25 = 0.
    M26 = 0. 

    M31 = 0.
    M32 = 0.
    M33 = cy
    M34 = sy
    M35 = 0.
    M36 = 0.

    M41 = 0.
    M42 = 0.
    M43 = -k**2*sy
    M44 = cy
    M45 = 0.
    M46 = 0.

    M51 = 0.
    M52 = 0.
    M53 = 0.
    M54 = 0.
    M55 = 1.
    M56 = L/(ga*ga*Be*Be)

    M61 = 0.
    M62 = 0.
    M63 = 0.
    M64 = 0.
    M65 = 0.
    M66 = 1.

    Map = array([[M11,M12,M13,M14,M15,M16],\
               [M21,M22,M23,M24,M25,M26],\
               [M31,M32,M33,M34,M35,M36],\
               [M41,M42,M43,M44,M45,M46],\
               [M51,M52,M53,M54,M55,M56],\
               [M61,M62,M63,M64,M65,M66]])

    return Map

  def Drft(self,L):
    Be   = self.Beta 
    ga = self.gamma 
    M11 = 1.
    M12 = L 
    M13 = 0.
    M14 = 0.
    M15 = 0.
    M16 = 0.

    M21 = 0.
    M22 = 1.
    M23 = 0.
    M24 = 0.
    M25 = 0.
    M26 = 0.

    M31 = 0.
    M32 = 0.
    M33 = 1.
    M34 = L 
    M35 = 0.
    M36 = 0.

    M41 = 0.
    M42 = 0.
    M43 = 0.
    M44 = 1.
    M45 = 0. 
    M46 = 0. 

    M51 = 0.
    M52 = 0.
    M53 = 0.
    M54 = 0.
    M55 = 1.
    M56 = L/(ga*ga*Be*Be) 

    M61 = 0.
    M62 = 0.
    M63 = 0.
    M64 = 0.
    M65 = 0.
    M66 = 1.
         
    Map = array([[M11,M12,M13,M14,M15,M16],\
                [M21,M22,M23,M24,M25,M26],\
                [M31,M32,M33,M34,M35,M36],\
                [M41,M42,M43,M44,M45,M46],\
                [M51,M52,M53,M54,M55,M56],\
                [M61,M62,M63,M64,M65,M66]])

    return Map

  def Bend(self,L,angle):
    Be   = self.Beta 
    ga = self.gamma 
    h = angle/L;
    K = 0.0;
    kx = sqrt(h**2+K);
    cx = cos(kx*L);
    sx = sin(kx*L)/kx;
    dx = (1-cx)/kx**2;
    J1 = (L - sx)/kx**2;
    sy = L;
    cy = 1.0;
    ky = 0;

    M11 = cx 
    M12 = sx 
    M13 = 0.
    M14 = 0.
    M15 = 0.
    M16 = (h/Be)*dx

    M21 = -kx**2*sx 
    M22 = cx 
    M23 = 0.
    M24 = 0.
    M25 = 0.
    M26 = (h/Be)*sx

    M31 = 0.
    M32 = 0.
    M33 = cy 
    M34 = sy 
    M35 = 0.
    M36 = 0.

    M41 = 0.
    M42 = 0.
    M43 = -ky**2*sy
    M44 = cy 
    M45 = 0.
    M46 = 0.

    M51 = -(h/Be)*sx 
    M52 = -(h/Be)*dx
    M53 =  0.
    M54 =  0. 
    M55 =  1.
    M56 =  -(h/Be)**2*J1+L/(Be**2*ga**2)

    M61 = 0.
    M62 = 0.
    M63 = 0. 
    M64 = 0.
    M65 = 0.
    M66 = 1.
 
    Map = array([[M11,M12,M13,M14,M15,M16],\
               [M21,M22,M23,M24,M25,M26],\
               [M31,M32,M33,M34,M35,M36],\
               [M41,M42,M43,M44,M45,M46],\
               [M51,M52,M53,M54,M55,M56],\
               [M61,M62,M63,M64,M65,M66]])

    return Map

  def RFkick(self):
    mus = self.mus
    eta = self.eta
    C = self.C

    M11 = 1.
    M12 = 0.
    M13 = 0.
    M14 = 0.
    M15 = 0.
    M16 = 0.

    M21 = 0. 
    M22 = 1.
    M23 = 0.
    M24 = 0.
    M25 = 0.
    M26 = 0.

    M31 = 0.
    M32 = 0.
    M33 = 1.
    M34 = 0.
    M35 = 0.
    M36 = 0.

    M41 = 0.
    M42 = 0.
    M43 = 0. 
    M44 = 1.
    M45 = 0.
    M46 = 0.

    M51 = 0.
    M52 = 0.
    M53 = 0.
    M54 = 0.
    M55 = 1.
    M56 = 0.

    M61 = 0.
    M62 = 0.
    M63 = 0.
    M64 = 0.
    M65 = (2*pi*mus)**2/(eta*C) 
    M66 = 1.

    Map = array([[M11,M12,M13,M14,M15,M16],\
               [M21,M22,M23,M24,M25,M26],\
               [M31,M32,M33,M34,M35,M36],\
               [M41,M42,M43,M44,M45,M46],\
               [M51,M52,M53,M54,M55,M56],\
               [M61,M62,M63,M64,M65,M66]])

    return Map

class Maps_twiss:
  """
  not complete
  """
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

