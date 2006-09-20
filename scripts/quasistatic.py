"""
This script contains the class for performing runs in quasistatic approximation.
"""
from warp import *

class Quasistatic:
  def __init__(self,slist=None,MRroot=None,l_verbose=0,l_findmgparam=0,
               nparpgrp=top.nparpgrp,Ninit=1000,l_mode=1,l_selfe=1,pboundxy=None):
    frz.init_base(w3d.nx,w3d.ny,w3d.dx,w3d.dy,w3d.xmmin,w3d.ymmin)
    if MRroot is None:
      self.l_MR = 0
      frz.mgridrz_accuracy = f3d.mgtol
    else:
      self.l_MR = 1
      frz.mgridrz_accuracy = MRroot.mgtol
      self.MRroot = MRroot
      b=MRroot
      brz=frz.basegrid
      l_addblock=true
      while l_addblock:
        try:
          b=b.children[0]
        except:
          l_addblock=false
        if l_addblock:
          add_subgrid(brz.gid[0],
                      b.nx,#-2*b.nguard*b.refinement[0],
                      b.ny,#-2*b.nguard*b.refinement[1],
                      b.dx,b.dy,
                      b.xmmin,#+b.nguard*b.dx*b.refinement[0],
                      b.ymmin,#+b.nguard*b.dy*b.refinement[1],
                      b.nguard*b.refinement[0],
                      b.nguard*b.refinement[0],
                      b.nguard*b.refinement[1],
                      b.nguard*b.refinement[1])
          brz=brz.down
    self.dt = w3d.dz/top.vbeam # time step for pushing electrons
    self.nparpgrp = nparpgrp
    self.ex = zeros(self.nparpgrp,'d')
    self.ey = zeros(self.nparpgrp,'d')
    self.ez = zeros(self.nparpgrp,'d')
    self.l_verbose = l_verbose
    self.slist = slist
    self.l_findmgparam = l_findmgparam
    self.Ninit = Ninit
    self.l_mode = l_mode
    self.l_selfe = l_selfe
    if pboundxy is None:
      self.pboundxy = top.pboundxy
    else:
      self.pboundxy = pboundxy
      
  def push(self,l_return_dist=false):
    # --- defines shortcuts
    pg = top.pgroup
    bg = frz.basegrid
    # --- generate electrons
    self.slist[0].add_uniform_box(self.Ninit,w3d.xmmin,w3d.xmmax,w3d.ymmin,w3d.ymmax,0.,0.,1.e-10,1.e-10,1.e-10)
    if l_return_dist:
      dist=[]
      sp=self.slist[0]
      dist.append([sp.getx(),sp.gety(),sp.getvx(),sp.getvy()])
    # --- switch to 2-D solver
    solvergeomcp = w3d.solvergeom
    fstypecp = top.fstype
    w3d.solvergeom = w3d.XYgeom
    # --- loop over 3-D grid mesh in z
    if self.l_mode==2:
      izmin = w3d.nz/4
      nzmax = 2*izmin
    else:
      izmin = 0
      nzmax = w3d.nz
    for izs in range(nzmax,-1,-1):
      iz = izmin+izs
      if self.l_verbose:print 'iz = ',iz
      # --- reset 2-D rho arrays
      reset_rzmgrid_rho()
      # --- loop over species list
      for sp in self.slist:
        # --- loop over species index
        for js in sp.jslist:
          np = pg.nps[js]
          il = pg.ins[js]-1
          iu = il+pg.nps[js]
          if self.l_verbose:print 'deposit rho' # MR OK
          # --- deposit rho
          rhoweightrz(pg.xp[il:iu],pg.yp[il:iu],pg.yp[il:iu],np,pg.sq[js]*pg.sw[js],bg.nr,bg.nz,
                      bg.dr,bg.dz,bg.rmin,top.zgrid)
      # --- distribute rho among patches
      distribute_rho_rz()
      if self.l_verbose:print 'feldsolve'   # MR OK
      if self.l_findmgparam and top.it==0 and iz==w3d.nz/2:frz.find_mgparam_rz(false)
      solve_mgridrz(bg,frz.mgridrz_accuracy,true)
      # --- stack 2-D phi from electrons into 3-D array
      # --- and set 2-D phi from precalculated 3-D potential
      if self.l_MR:
        b3d = self.MRroot
        brz = bg
        for ig in range(frz.ngrids):
          if ig>0:
            b3d = b3d.children[0]
            brz = brz.down
          phitmp = brz.phi[1:-1,1:-1].copy()
          if self.l_selfe:
            brz.phi[1:-1,1:-1] += b3d.phi[:,:,iz+1] 
          else:
            brz.phi[1:-1,1:-1] = b3d.phi[:,:,iz+1].copy()
          b3d.phi[:,:,iz+1] = phitmp.copy()
      else:
        phitmp = bg.phi[1:-1,1:-1].copy()
        if self.l_selfe:
          bg.phi[1:-1,1:-1] += w3d.phi[:,:,iz+1] 
        else:
          bg.phi[1:-1,1:-1] = w3d.phi[:,:,iz+1].copy() 
        w3d.phi[:,:,iz+1] = phitmp.copy()
#      if top.it==0:window(0);self.slist[0].ppxy();pfxy(iz=iz);fma()
      # --- update guard cells of 2-D solver
      updateguardcells2d()
      # --- loop over species list
      for sp in self.slist:
        # --- loop over species index
        for js in sp.jslist:
          ng = 1+pg.nps[js]/self.nparpgrp
          for ig in range(ng):
            il = pg.ins[js]-1+self.nparpgrp*ig
            iu = min(il+self.nparpgrp,pg.ins[js]-1+pg.nps[js])
            np = iu-il
            # --- gather self-forces
            if self.l_verbose:print 'fieldweight' # MR OK
            fieldweightxz(pg.xp[il:iu],pg.yp[il:iu],self.ex[:np],self.ey[:np],np,top.zgrid)
            if self.l_verbose:print 'epush'
            epush3d(np,pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],
                    self.ex[:np],self.ey[:np],self.ez[:np],pg.sq[js],pg.sm[js],0.5*self.dt)
            gammaadv(np,pg.gaminv[il:iu],pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],
                     top.gamadv,top.lrelativ)
            if self.l_verbose:print 'xpush'
            xpush3d(np,pg.xp[il:iu],pg.yp[il:iu],pg.zp[il:iu],
                    pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],pg.gaminv[il:iu],self.dt)
            if self.l_verbose:print 'stckxy3d'
            stckxy3d(np,pg.xp[il:iu],w3d.xmmax,w3d.xmmin,w3d.dx,
                     pg.yp[il:iu],w3d.ymmax,w3d.ymmin,w3d.dy,
                     pg.zp[il:iu],w3d.zmmin,w3d.dz,
                     pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],pg.gaminv[il:iu],
                     top.zgrid,top.zbeam,w3d.l2symtry,w3d.l4symtry,self.pboundxy,true)
          processlostpart(top.pgroup,js+1,top.clearlostpart,top.time+self.dt*top.pgroup.ndts[js],top.zbeam)
          ng = 1+pg.nps[js]/self.nparpgrp
          for ig in range(ng):
            il = pg.ins[js]-1+self.nparpgrp*ig
            iu = min(il+self.nparpgrp,pg.ins[js]-1+pg.nps[js])
            np = iu-il
            if self.l_verbose:print 'fieldweight' # MR OK
            fieldweightxz(pg.xp[il:iu],pg.yp[il:iu],self.ex[:np],self.ey[:np],np,top.zgrid)
            if self.l_verbose:print 'epush'
            epush3d(np,pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],
                    self.ex[:np],self.ey[:np],self.ez[:np],pg.sq[js],pg.sm[js],0.5*self.dt)
            gammaadv(np,pg.gaminv[il:iu],pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],
                     top.gamadv,top.lrelativ)
      if l_return_dist:
        sp=self.slist[0]
        dist.append([sp.getx(),sp.gety(),sp.getvx(),sp.getvy()])
    if self.l_mode==2:
      print izmin,nzmax,self.l_mode,self.l_MR
      if self.l_MR:
        b3d = self.MRroot
        for ig in range(frz.ngrids):
          if ig>0:
            b3d = b3d.children[0]
          for iz in range(izmin+1):
            b3d.phi[:,:,iz] = b3d.phi[:,:,izmin+1]
            b3d.phi[:,:,-iz-1] = b3d.phi[:,:,izmin+nzmax+1]
      else:
        for iz in range(izmin+1):
          w3d.phi[:,:,iz] = w3d.phi[:,:,izmin+1]
          w3d.phi[:,:,-iz-1] = w3d.phi[:,:,izmin+nzmax+1]
    # --- switch back to 3-D solver
    w3d.solvergeom = solvergeomcp
    top.fstype = fstypecp
    if l_return_dist: 
      return dist

class Quasistatic2:
  def __init__(self,slist=None,MRroot=None,l_verbose=0,l_findmgparam=0,nparpgrp=top.nparpgrp,Ninit=1000):
#    wxy.lthick=1
    self.pg = ParticleGroup()
    # --- switch to 2-D solver
    solvergeomcp = w3d.solvergeom
    fstypecp = top.fstype
    w3d.solvergeom = w3d.XYgeom
    package('wxy');generate()
    # --- switch back to 3-D solver
    w3d.solvergeom = solvergeomcp
    top.fstype = fstypecp
    package('w3d')
    if MRroot is None:
      self.l_MR = 0
      frz.mgridrz_accuracy = f3d.mgtol
    else:
      self.l_MR = 1
      frz.mgridrz_accuracy = MRroot.mgtol
      self.MRroot = MRroot
      b=MRroot
      brz=frz.basegrid
      l_addblock=true
      while l_addblock:
        try:
          b=b.children[0]
        except:
          l_addblock=false
        if l_addblock:
          add_subgrid(brz.gid[0],
                      b.nx,#-2*b.nguard*b.refinement[0],
                      b.ny,#-2*b.nguard*b.refinement[1],
                      b.dx,b.dy,
                      b.xmmin,#+b.nguard*b.dx*b.refinement[0],
                      b.ymmin,#+b.nguard*b.dy*b.refinement[1],
                      b.nguard*b.refinement[0],
                      b.nguard*b.refinement[0],
                      b.nguard*b.refinement[1],
                      b.nguard*b.refinement[1])
          brz=brz.down
    self.dt = w3d.dz/top.vbeam # time step for pushing electrons
    self.ex = zeros(top.nparpgrp,'d')
    self.ey = zeros(top.nparpgrp,'d')
    self.ez = zeros(top.nparpgrp,'d')
    self.l_verbose = l_verbose
    self.slist = slist
    self.l_findmgparam = l_findmgparam
    self.nparpgrp = nparpgrp
    self.Ninit = Ninit
    
  def push(self):
    # --- defines shortcuts
    pg = top.pgroup
    bg = frz.basegrid
    # --- generate electrons
    self.slist[0].add_uniform_box(self.Ninit,w3d.xmmin,w3d.xmmax,w3d.ymmin,w3d.ymmax,0.,0.,1.e-10,1.e-10,1.e-10)
    # --- switch to 2-D solver
    solvergeomcp = w3d.solvergeom
    fstypecp = top.fstype
    w3d.solvergeom = w3d.XYgeom
    # --- loop over 3-D grid mesh in z
    for iz in range(w3d.nz,-1,-1):
      if self.l_verbose:print 'iz = ',iz
      # --- reset 2-D rho arrays
      reset_rzmgrid_rho()
      # --- loop over species list
      for sp in self.slist:
        # --- loop over species index
        for js in sp.jslist:
          np = pg.nps[js]
          il = pg.ins[js]-1
          iu = il+pg.nps[js]
          if self.l_verbose:print 'deposit rho' # MR OK
          # --- deposit rho
          rhoweightrz(pg.xp[il:iu],pg.yp[il:iu],pg.yp[il:iu],np,pg.sq[js]*pg.sw[js],bg.nr,bg.nz,
                      bg.dr,bg.dz,bg.rmin,top.zgrid)
      # --- distribute rho among patches
      distribute_rho_rz()
      if self.l_verbose:print 'feldsolve'   # MR OK
      if self.l_findmgparam and top.it==0 and iz==w3d.nz/2:frz.find_mgparam_rz(false)
      solve_mgridrz(bg,frz.mgridrz_accuracy,true)
      # --- stack 2-D phi from electrons into 3-D array
      # --- and set 2-D phi from precalculated 3-D potential
      if self.l_MR:
        b3d = self.MRroot
        brz = bg
        for ig in range(frz.ngrids):
          if ig>0:
            b3d = b3d.children[0]
            brz = brz.down
          phitmp = brz.phi[1:-1,1:-1].copy()
          brz.phi[1:-1,1:-1] += b3d.phi[:,:,iz+1] 
          b3d.phi[:,:,iz+1] = phitmp.copy()
      else:
        phitmp = bg.phi[1:-1,1:-1].copy()
        bg.phi[1:-1,1:-1] += w3d.phi[:,:,iz+1] 
        w3d.phi[:,:,iz+1] = phitmp.copy()
      # --- loop over species list
      for sp in self.slist:
        # --- loop over species index
        for js in sp.jslist:
          ng = 1+pg.nps[js]/self.nparpgrp
          for ig in range(ng):
            il = pg.ins[js]-1+self.nparpgrp*ig
            iu = min(il+self.nparpgrp,pg.ins[js]-1+pg.nps[js])
            np = iu-il
            # --- gather self-forces
            if self.l_verbose:print 'fieldweight' # MR OK
            fieldweightxz(pg.xp[il:iu],pg.yp[il:iu],self.ex[:np],self.ey[:np],np,top.zgrid)
            if self.l_verbose:print 'epush'
            epush3d(np,pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],
                    self.ex[:np],self.ey[:np],self.ez[:np],pg.sq[js],pg.sm[js],0.5*self.dt)
            gammaadv(np,pg.gaminv[il:iu],pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],
                     top.gamadv,top.lrelativ)
            if self.l_verbose:print 'xpush'
            xpush3d(np,pg.xp[il:iu],pg.yp[il:iu],pg.zp[il:iu],
                    pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],pg.gaminv[il:iu],self.dt)
            if self.l_verbose:print 'stckxy3d'
            stckxy3d(np,pg.xp[il:iu],w3d.xmmax,w3d.xmmin,w3d.dx,
                     pg.yp[il:iu],w3d.ymmax,w3d.ymmin,w3d.dy,
                     pg.zp[il:iu],w3d.zmmin,w3d.dz,
                     pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],pg.gaminv[il:iu],
                     top.zgrid,top.zbeam,w3d.l2symtry,w3d.l4symtry,top.pboundxy,true)
          processlostpart(top.pgroup,js+1,top.clearlostpart,top.time+self.dt*top.pgroup.ndts[js],top.zbeam)
          ng = 1+pg.nps[js]/self.nparpgrp
          for ig in range(ng):
            il = pg.ins[js]-1+self.nparpgrp*ig
            iu = min(il+self.nparpgrp,pg.ins[js]-1+pg.nps[js])
            np = iu-il
            if self.l_verbose:print 'fieldweight' # MR OK
            fieldweightxz(pg.xp[il:iu],pg.yp[il:iu],self.ex[:np],self.ey[:np],np,top.zgrid)
            if self.l_verbose:print 'epush'
            epush3d(np,pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],
                    self.ex[:np],self.ey[:np],self.ez[:np],pg.sq[js],pg.sm[js],0.5*self.dt)
            gammaadv(np,pg.gaminv[il:iu],pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],
                     top.gamadv,top.lrelativ)
    # --- switch back to 3-D solver
    w3d.solvergeom = solvergeomcp
    top.fstype = fstypecp

