"""
This script contains the class for performing runs in quasistatic approximation.
"""
from warp import *

class Quasistatic:
  def __init__(self,slist=None,MRroot=None,l_verbose=0,l_findmgparam=0,
               nparpgrp=top.nparpgrp,Ninit=1000,l_mode=1,l_selfe=1,l_maps=0,pboundxy=None,
               conductors=[],l_elecuniform=0,scraper=None):
    if w3d.solvergeom==w3d.RZgeom:
      self.l_rz=1
      self.gridr=GRIDtype()
      gridcp=frz.basegrid
      frz.basegrid=self.gridr
      w3d.solvergeom=w3d.Rgeom
      nguardz=frz.nguardz+0
      init_base(w3d.nx,0,w3d.dx,0.,0.,0.,false)
      frz.basegrid=gridcp
      w3d.solvergeom=w3d.RZgeom
      frz.nguardz=nguardz
      frz.init_base(w3d.nx,w3d.nz,w3d.dx,w3d.dz,w3d.xmmin,w3d.zmmin,false)
      self.l_MR=0
    else:
      self.l_rz=0
      solvergeomcp=w3d.solvergeom
      w3d.solvergeom=w3d.XYgeom
      frz.init_base(w3d.nx,w3d.ny,w3d.dx,w3d.dy,w3d.xmmin,w3d.ymmin,false)
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
          except:
            l_addblock=false
        g=frz.basegrid
        for conductor in conductors:
          if l_addblock:raise('conductors need to be installed in quasistatic MR')
          try:
            cond = conductor.cond
          except AttributeError:
            cond = conductor
          self.conductors=ConductorType()
          installconductors(cond,conductors=self.conductors,gridrz=g,nz=0)
      w3d.solvergeom=solvergeomcp
    self.dt = w3d.dz/top.vbeam # time step for pushing electrons
    self.nparpgrp = nparpgrp
    self.l_verbose = l_verbose
    self.slist = slist
    self.l_findmgparam = l_findmgparam and not self.l_rz
    self.Ninit = Ninit
    self.l_mode = l_mode
    if npes>1:
      self.l_mode=1
    else:
      self.l_mode = l_mode
    self.l_selfe = l_selfe
    self.l_elecuniform=l_elecuniform
    self.scraper=scraper
    self.l_maps=l_maps
    if pboundxy is None:
      self.pboundxy = top.pboundxy
    else:
      self.pboundxy = pboundxy
      
  def push(self,l_return_dist=false,l_return_rho=false,l_plotelec=0):
    if npes>1:self.l_mode=1
    # --- defines shortcuts
    pg = top.pgroup
    bg = frz.basegrid
    # --- generate electrons (only on last processor)
    if me==max(0,npes-1) and (top.it==0  or not l_plotelec):
     pg.nps[1]=0
     if self.l_mode==2:
      zmax = w3d.zmmaxglobal/2+top.zgrid
     else:
      zmax = w3d.zmmaxglobal-1.e-10*w3d.dz+top.zgrid
     if top.prwall<sqrt(w3d.xmmax**2+w3d.ymmax**2) or self.l_rz:
      self.slist[0].add_uniform_cylinder(self.Ninit,min(top.prwall,w3d.xmmax),zmax,zmax,1.e-10,1.e-10,1.e-10,lallindomain=true)
     else:
      if self.l_elecuniform:
        spacing='uniform'
      else:
        spacing='random'
      self.slist[0].add_uniform_box(self.Ninit,w3d.xmmin,w3d.xmmax,
                                    w3d.ymmin,w3d.ymmax,zmax,zmax,
                                    1.e-10,1.e-10,1.e-10,spacing=spacing,lallindomain=true)
    if self.scraper is not None:self.scraper.scrapeall(local=1,clear=1)
    if l_return_dist:
      dist=[]
      sp=self.slist[0]
      dist.append([sp.getx(),sp.gety(),sp.getvx(),sp.getvy()])
    if l_return_rho:
      rho = zeros([w3d.nx+1,w3d.ny+1,w3d.nz+1],'d')
    # --- switch to 2-D solver
    if self.l_rz:
      gridcp=frz.basegrid
      frz.basegrid=self.gridr
      bg=self.gridr
      solvergeomcp = w3d.solvergeom
      fstypecp = top.fstype
      w3d.solvergeom = w3d.Rgeom
    else:
      solvergeomcp = w3d.solvergeom
      fstypecp = top.fstype
      w3d.solvergeom = w3d.XYgeom
#    if npes>1:self.MRroot.phi[...]=0.
    # --- loop over 3-D/R-Z grid mesh in z
    if npes==0 and self.l_mode==2:
      izmin = w3d.nz/4
      nzmax = 2*izmin
    else:
      izmin = 0
#      nzmax = nint((getz(js=self.slist[0])[0]-w3d.zmmin)/w3d.dz)#w3d.nz
      nzmax = w3d.nzp
    for izs in range(nzmax,-1,-1):
      iz = izmin+izs
      if self.l_verbose:print 'me = ',me,';iz = ',iz
      # --- reset 2-D rho arrays
      if self.l_rz:
        self.gridr.rho[...]=0.
      else:
        reset_rzmgrid_rho()
      # --- loop over species list
      for sp in self.slist:
        # --- loop over species index
        for js in sp.jslist:
          np = pg.nps[js]
          if np==0:continue
          il = pg.ins[js]-1
          iu = il+pg.nps[js]
          if self.l_verbose:print 'me = ',me,';deposit rho' # MR OK
          # --- deposit rho
          if self.l_rz:
            rhoweightr(pg.xp[il:iu],pg.yp[il:iu],np,pg.sq[js]*pg.sw[js],
                       bg.nr,bg.dr,bg.rmin)
          else:
            rhoweightrz(pg.xp[il:iu],pg.yp[il:iu],pg.yp[il:iu],np,
                        pg.sq[js]*pg.sw[js],bg.nr,bg.nz,
                        bg.dr,bg.dz,bg.rmin,0.)
      # --- distribute rho among patches
      distribute_rho_rz()
      if l_return_rho:rho[:,:,iz]=frz.basegrid.rho
      if self.l_verbose:print 'me = ',me,';feldsolve'   # MR OK
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
            enzl=b3d.extradimslower[2]
            enzu=b3d.extradimsupper[2]
          else:
            enzl=0
            enzu=0
          if self.l_maps:
            phitmp = brz.phi[1:-1,1:-1].copy()
            if self.l_selfe:
#              brz.phi[1:-1,1:-1] += b3d.phi[...,iz+1+enzl] 
              brz.phi[1:-1,1:-1] += b3d.potentialparray[:,:,iz+1+enzl,0,0]
            else:
              brz.phi[1:-1,1:-1] = b3d.phi[...,iz+1+enzl].copy()
            b3d.potentialparray[:,:,iz+1+enzl,0,0] = phitmp.copy()
#            b3d.phi[:,:,iz+1+enzl] = phitmp.copy()
            del phitmp
          else:
            if self.l_selfe:
#              brz.phi[1:-1,1:-1] += b3d.phi[...,iz+1+enzl] 
#              b3d.phi[:,:,iz+1+enzl] = brz.phi[1:-1,1:-1]
              b3d.potentialparray[:,:,iz+1+enzl,0,1] = brz.phi[1:-1,1:-1]
              brz.phi[1:-1,1:-1] += b3d.potentialparray[:,:,iz+1+enzl,0,0]
              b3d.potentialparray[:,:,iz+1+enzl,0,0]=0. # would be added again by getallpotentialpforparticles() otherwise
            else:
              phitmp = brz.phi[1:-1,1:-1].copy()
              brz.phi[1:-1,1:-1] = b3d.phi[...,iz+1+enzl].copy()
              b3d.phi[:,:,iz+1+enzl] += phitmp.copy()
      else:
       if self.l_rz:
        phitmp = bg.phi.copy()
        if self.l_selfe:
          bg.phi[:,0] += gridcp.phi[...,iz+1] 
        else:
          bg.phi[:,0] = gridcp.phi[...,iz+1].copy() 
        gridcp.phi[...,iz+1] = phitmp[:,0].copy()
       else:
        if self.l_maps:
          phitmp = bg.phi[1:-1,1:-1].copy()
          if self.l_selfe:
            bg.phi[1:-1,1:-1] += w3d.phi[:,:,iz+1] 
          else:
            bg.phi[1:-1,1:-1] = w3d.phi[:,:,iz+1].copy() 
          w3d.phi[:,:,iz+1] = phitmp.copy()
        else:
          raise('this part needs to be programmed: Quasistatic+Boris-l_MR')
          
#     if l_plotelec and iz<w3d.nzfull/max(1,npes):
#        ppgeneric(getvx(js=1),getx(js=1))
#        window(0);self.slist[0].ppxy();        limits(w3d.xmmin,w3d.xmmax,w3d.ymmin,w3d.ymmax)
#        ppxex(js=1)
#        pfxy(iz=iz);
#        fma()
      # --- update guard cells of 2-D solver
      updateguardcells2d()
      if(iz>0): # push particles but on last step
       # --- loop over species list
       for sp in self.slist:
        # --- loop over species index
        for js in sp.jslist:
          ng = 1+pg.nps[js]/self.nparpgrp
          for ig in range(ng):
            il = pg.ins[js]-1+self.nparpgrp*ig
            iu = min(il+self.nparpgrp,pg.ins[js]-1+pg.nps[js])
            np = iu-il
            if np==0:continue
            # --- gather self-forces
            if self.l_verbose:print 'me = ',me,';fieldweight' # MR OK
            if self.l_rz:
              fieldweightr(pg.xp[il:iu],pg.yp[il:iu],pg.ex[il:iu],pg.ey[il:iu],np)
            else:
              fieldweightxz(pg.xp[il:iu],pg.yp[il:iu],pg.ex[il:iu],pg.ey[il:iu],np,0.,top.efetch[js])
#            fetche3d(pg,il+1,np,js+1)
            if self.l_verbose:print 'me = ',me,';epush'
            epush3d(np,pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],
                    pg.ex[il:iu],pg.ey[il:iu],pg.ez[il:iu],pg.sq[js],pg.sm[js],0.5*self.dt)
            gammaadv(np,pg.gaminv[il:iu],pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],
                     top.gamadv,top.lrelativ)
            if self.l_verbose:print 'me = ',me,';xpush'
            xpush3d(np,pg.xp[il:iu],pg.yp[il:iu],pg.zp[il:iu],
                    pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],pg.gaminv[il:iu],self.dt)
            pg.zp[il:iu]-=w3d.dz
#            pli(frz.basegrid.phi);refresh()
#            ppzx(js=1,color=red,msize=2);refresh()
#            window(3);ppzvx(js=1,msize=2);window(0)
            if self.l_verbose:print 'me = ',me,';stckxy3d'
            stckxy3d(np,pg.xp[il:iu],w3d.xmmax,w3d.xmmin,w3d.dx,
                     pg.yp[il:iu],w3d.ymmax,w3d.ymmin,w3d.dy,
                     pg.zp[il:iu],w3d.zmmin,w3d.dz,
                     pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],pg.gaminv[il:iu],
                     top.zgrid,top.zbeam,w3d.l2symtry,w3d.l4symtry,self.pboundxy,true)
          if self.scraper is not None:
            self.scraper.scrapeall(local=1,clear=1)
          else:
            processlostpart(top.pgroup,js+1,top.clearlostpart,top.time+self.dt*top.pgroup.ndts[js],top.zbeam)
          top.npslost=0
          ng = 1+pg.nps[js]/self.nparpgrp
          for ig in range(ng):
            il = pg.ins[js]-1+self.nparpgrp*ig
            iu = min(il+self.nparpgrp,pg.ins[js]-1+pg.nps[js])
            np = iu-il
            if np==0:continue
            if self.l_verbose:print 'me = ',me,';fieldweight' # MR OK
            if self.l_rz:
              fieldweightr(pg.xp[il:iu],pg.yp[il:iu],pg.ex[il:iu],pg.ey[il:iu],np)
            else:
              fieldweightxz(pg.xp[il:iu],pg.yp[il:iu],pg.ex[il:iu],pg.ey[il:iu],np,top.zgrid,top.efetch[js])
#            fetche3d(pg,il+1,np,js+1)
            if self.l_verbose:print 'me = ',me,';epush'
            epush3d(np,pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],
                    pg.ex[il:iu],pg.ey[il:iu],pg.ez[il:iu],pg.sq[js],pg.sm[js],0.5*self.dt)
            gammaadv(np,pg.gaminv[il:iu],pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],
                     top.gamadv,top.lrelativ)
      if l_return_dist:
        sp=self.slist[0]
        dist.append([sp.getx(),sp.gety(),sp.getvx(),sp.getvy()])
#      if pg.nps[1]>0:
#        zel=ave(getz(js=1,gather=0,bcast=0))
#        print 'me,iz,zel',me,iz,zel,nint((zel-w3d.zmminglobal)/w3d.dz)
    # --- switch back to 3-D solver
    w3d.solvergeom = solvergeomcp
    top.fstype = fstypecp
    if self.l_rz:
      frz.basegrid=gridcp
      bg=frz.basegrid
    if self.l_mode==1:
      self.fillphizends(1)
    else:
      self.fillphizends(izmin)
    if self.MRroot is not None:
      self.setidosolve(self.MRroot,0)
#      self.MRroot.getallpotentialpforparticles()    
#      fieldsol()
      self.MRroot.solve()
      self.setidosolve(self.MRroot,1)      
    if me==0:pg.nps[1]=0      
    if npes>1:
      js = 1
      if pg.nps[js]>0:
        il=top.pgroup.ins[js]-1
        iu=il+top.pgroup.nps[js]
        top.pgroup.zp[il:iu]=w3d.zmmin-1.e-10*w3d.dz
      zpartbnd(top.pgroup,w3d.zmmax,w3d.zmmin,w3d.dz)
    if l_return_dist: 
      return dist
    if l_return_rho:
      return rho

  def setefetch(self,fsolver,efetch):
        fsolver.efetch=efetch
        try:
          self.setefetch(fsolver.children[0],efetch)
        except:
          return
 
  def setidosolve(self,fsolver,idosolve):
        fsolver.l_internal_dosolve=idosolve
        try:
          self.setidosolve(fsolver.children[0],idosolve)
        except:
          return
 
  def fillphizends(self,nz):
      if self.l_MR:
        b3d = self.MRroot
        for ig in range(frz.ngrids):
          if ig>0:
            b3d = b3d.children[0]
            enzl=b3d.extradimslower[2]
            enzu=b3d.extradimsupper[2]
          else:
            enzl=0
            enzu=0
          for iz in range(nz+enzl):
            b3d.phi[:,:,iz] = b3d.phi[:,:,nz+enzl]
          for iz in range(nz+enzu):
            b3d.phi[:,:,-iz-1] = b3d.phi[:,:,-nz-1-enzu]
      else:
        if self.l_rz:
          for iz in range(nz):
            bg.phi[:,iz] = bg.phi[:,nz]
            bg.phi[:,-iz-1] = bg.phi[:,-nz-1]
        else:
          for iz in range(nz):
            w3d.phi[:,:,iz] = w3d.phi[:,:,nz]
            w3d.phi[:,:,-iz-1] = w3d.phi[:,:,-nz-1]

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
            fieldweightxz(pg.xp[il:iu],pg.yp[il:iu],self.ex[:np],self.ey[:np],np,top.zgrid,top.efetch[js])
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
            fieldweightxz(pg.xp[il:iu],pg.yp[il:iu],self.ex[:np],self.ey[:np],np,top.zgrid,top.efetch[js])
            if self.l_verbose:print 'epush'
            epush3d(np,pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],
                    self.ex[:np],self.ey[:np],self.ez[:np],pg.sq[js],pg.sm[js],0.5*self.dt)
            gammaadv(np,pg.gaminv[il:iu],pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],
                     top.gamadv,top.lrelativ)
    # --- switch back to 3-D solver
    w3d.solvergeom = solvergeomcp
    top.fstype = fstypecp

