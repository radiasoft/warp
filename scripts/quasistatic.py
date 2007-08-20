"""
This script contains the class for performing runs in quasistatic approximation.
In parallel, there is a caveat that each processor carries particles at different time steps. 
"""
from warp import *
from getzmom import *
from AMR import *
from appendablearray import *
import __main__

class Quasistatic:
  # this class differs from the previous one that it does not use the 3-D solver for the ions, but the 2-D one.
  def __init__(self,ions,MRroot=None,l_verbose=0,l_findmgparam=0,
               nparpgrp=top.nparpgrp,Ninit=1000,l_mode=1,l_selfe=1,maps=None,pboundxy=None,
               conductors=[],l_elecuniform=0,scraper=None,l_weakstrong=0,nelecperiod=1,
               backgroundtype=Electron):
    w3d.solvergeom=w3d.XYgeom
    self.gridelecs=[]
    self.gridions=[]
    self.gridionscp=[]
    for gxy in [self.gridions,self.gridelecs,self.gridionscp]:
     for i in range(2):
      gxy.append(GRIDtype())
      frz.basegrid=gxy[-1]
      frz.init_base(w3d.nx,w3d.ny,w3d.dx,w3d.dy,w3d.xmmin,w3d.ymmin,false)
      if MRroot is None:
        self.l_MR = 0
        frz.mgridrz_accuracy = f3d.mgtol
      else:
        self.l_MR = 1
        frz.mgridrz_accuracy = MRroot.mgtol
        # note that the 3-D MR grid must have been setup by the user (this should be get rid off)
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
    self.MRroot = MRroot
    self.dt = w3d.dz/top.vbeam # time step for pushing electrons
    self.nparpgrp = nparpgrp
    self.l_verbose = l_verbose
    self.l_findmgparam = l_findmgparam
    self.Ninit = Ninit
    self.l_mode = l_mode
    if npes>1:
      self.l_mode=1
    else:
      self.l_mode = l_mode
    self.l_selfe = l_selfe
    self.l_elecuniform=l_elecuniform
    self.l_weakstrong=l_weakstrong
    self.nelecperiod=nelecperiod
    self.scraper=scraper
    if maps is not None:
      self.l_maps=true
      self.maps=maps
    if pboundxy is None:
      self.pboundxy = top.pboundxy
    else:
      self.pboundxy = pboundxy
    # --- sets z range
    if self.l_mode==2:
      self.izmin = w3d.nz/4
      self.izmax = w3d.nz-izmin
    else:
      self.izmin = 0
      self.izmax = w3d.nzp
    self.wz0=self.wz1=None
    self.rhoe = zeros([w3d.nx+1,w3d.ny+1,w3d.nzp+1],'d')
    self.phie = zeros([w3d.nx+3,w3d.ny+3,w3d.nzp+1],'d')
    self.rhoi = zeros([w3d.nx+1,w3d.ny+1,w3d.nzp+1],'d')
    self.phii = zeros([w3d.nx+3,w3d.ny+3,w3d.nzp+1],'d')
    self.pgions = top.pgroup
    self.ions = [ions]
    for iz in range(globalmax(w3d.nzp)-1):
      self.ions.append(Species(type=ions.type,fselfb=top.pgroup.fselfb[0]))
      top.pgroup.sq[-1] = top.pgroup.sq[0]
    self.pgelec = ParticleGroup()
    self.pgelec.gchange()
    top.pgroup = self.pgelec
    self.electrons = Species(type=backgroundtype)
    top.pgroup = self.pgions
    self.ionstoprev = [0]
    self.fullsortdone = 0
    self.pnum = AppendableArray(typecode='d')
    self.xbar = AppendableArray(typecode='d')
    self.ybar = AppendableArray(typecode='d')
    self.zbar = AppendableArray(typecode='d')
    self.xpbar = AppendableArray(typecode='d')
    self.ypbar = AppendableArray(typecode='d')
    self.x2 = AppendableArray(typecode='d')
    self.y2 = AppendableArray(typecode='d')
    self.z2 = AppendableArray(typecode='d')
    self.xp2 = AppendableArray(typecode='d')
    self.yp2 = AppendableArray(typecode='d')
    self.xxp = AppendableArray(typecode='d')
    self.yyp = AppendableArray(typecode='d')
    self.timemmnts = AppendableArray(typecode='d')
    self.iz=-1

  def step(self,nt=1,l_return_dist=false,l_plotelec=0,l_savedist=0):
   # --- if running in parallel, enforcing l_mode=1
   if npes>1:self.l_mode=1
   for it in range(nt):
     self.sort_ions_along_z()
     if top.it==0 or me>=(npes-top.it):
       self.getmmnts()
     self.compute_sw()
     if l_savedist:
       f=PW.PW('dist%g.pdb'%top.it)
       for iz in range(w3d.nzp):
         exec('f.x%i=getx(js=iz)'%iz)
         exec('f.y%i=gety(js=iz)'%iz)
         exec('f.z%i=getz(js=iz)'%iz)
         exec('f.vx%i=getvx(js=iz)'%iz)
         exec('f.vy%i=getvy(js=iz)'%iz)
         exec('f.vz%i=getvz(js=iz)'%iz)
         exec('f.ex%i=getex(js=iz)'%iz)
         exec('f.ey%i=getey(js=iz)'%iz)
         exec('f.ez%i=getez(js=iz)'%iz)
       f.close()
     if (not self.l_weakstrong or top.it==0) and (top.it%self.nelecperiod==0):
       # --- generate electrons (on last processor only)
       if me==max(0,npes-1) and (top.it==0  or not l_plotelec):self.create_electrons()
       # --- push 2-D electrons slice along bunch

       # --- (diagnostic) stacking of electron distribution
       if l_return_dist:
         self.dist=[]
         sp=self.slist[0]
         self.dist.append([sp.getx(),sp.gety(),sp.getvx(),sp.getvy()])

       # --- loop over "3-D" grid mesh in z
       for iz in range(self.izmax-1,self.izmin-1,-1):
         self.iz=iz
         
         if iz==(self.izmax-3):
           self.sendrecv_storedions()
           
         # --- deposit ions charge
         self.deposit_ions()

         # --- call 2-D field solver for ions
         frz.basegrid=self.gridions[1];mk_grids_ptr()
         # --- optimize mgparam
         if self.l_findmgparam and top.it==0 and self.iz in [w3d.nz/2,w3d.nz/2+1]:find_mgparam_rz(false)
         solve_mgridrz(self.gridions[1],frz.mgridrz_accuracy,true)
         if self.iz==0:solve_mgridrz(self.gridions[0],frz.mgridrz_accuracy,true)

         # --- deposit electrons charge
         self.deposit_electrons(1)

         # --- call 2-D field solver for electrons
         frz.basegrid=self.gridelecs[1];mk_grids_ptr()
         # --- optimize mgparam
         if self.l_findmgparam and top.it==0 and self.iz in [w3d.nz/2,w3d.nz/2+1]:find_mgparam_rz(false)
         solve_mgridrz(self.gridelecs[1],frz.mgridrz_accuracy,true)

         # --- add electron and ion fields
         self.add_ei_fields(1)

         self.rhoe[:,:,self.iz+1] = self.gridelecs[1].rho
         self.phie[:,:,self.iz+1] = self.gridelecs[1].phi
         self.rhoi[:,:,self.iz+1] = self.gridions[1].rho
         self.phii[:,:,self.iz+1] = self.gridions[1].phi

         # --- push electrons 
         self.push_electrons()
         if iz==0:
           self.deposit_electrons(0)
           frz.basegrid=self.gridelecs[0];mk_grids_ptr()
           solve_mgridrz(self.gridelecs[0],frz.mgridrz_accuracy,true)
           self.add_ei_fields(0)

         # --- plot electrons
#         if l_plotelec :self.plot_electrons()
         if l_plotelec and iz<w3d.nz/max(1,npes)-1:self.plot_electrons()

         # --- push ions
         if me>=(npes-1-top.it):
           self.push_ions()

         # --- (diagnostic) stacking of electron distribution
         if l_return_dist:
           sp=self.slist[0]
           self.dist.append([sp.getx(),sp.gety(),sp.getvx(),sp.getvy()])

         # --- switch 2-D solvers
         gtmp=self.gridelecs.pop(0)
         self.gridelecs.append(gtmp)
         gtmp=self.gridions.pop(0)
         self.gridions.append(gtmp)
#       for iz in range(globalmax(w3d.nzp)-w3d.nzp):
#         if l_plotelec :self.plot_electrons()
         
       self.rhoe[:,:,0] = self.gridelecs[1].rho
       self.phie[:,:,0] = self.gridelecs[1].phi
       self.rhoi[:,:,0] = self.gridions[1].rho
       self.phii[:,:,0] = self.gridions[1].phi

       self.reset_rho1()
#       if lparallel:mpi.barrier()

       self.store_ionstoprev()
       
     # --- clear electrons on processor 0, shift them to the previous ones for me>0
     self.clear_electrons()
     top.pgroup = self.pgions
     
     # --- update time, time counter
     top.time+=top.dt
     top.it+=1
     print me,top.it,self.iz,'it = %i, time = %gs.'%(top.it,top.time)

  def reset_rho1(self):
    if self.l_verbose:print me,top.it,self.iz,'enter reset_rho1'
    # --- sends rho1 to previous processor
    # --- reset rho1 on proc 0
    if not lparallel:
      bg=frz.basegrid=self.gridions[1]
      mk_grids_ptr()
      reset_rzmgrid_rho()
    else:
      if me>0:
        tosend = []
        g = self.gridions[1]
        for ig in range(frz.ngrids):
          if ig>0:
            g = g.down
          tosend.append(g.rho)
        mpi.send(tosend,me-1)
      if me<npes-1:
        recved = mpirecv(me+1)
        g = self.gridions[1]
        for ig in range(frz.ngrids):
          if ig>0:
            g = g.down
          g.rho[...] = recved[ig]
      if me == npes-1:
        bg=frz.basegrid=self.gridions[1]
        mk_grids_ptr()
        reset_rzmgrid_rho()
    if self.l_verbose:print me,top.it,self.iz,'exir reset_rho1'

  def store_ionstoprev(self):
    if self.l_verbose:print me,top.it,self.iz,'enter store_ionstoprev'
    js = 0
    pg = self.pgions
    il = pg.ins[js]-1
    iu = il+pg.nps[js]
    ii = il+compress(pg.zp[il:iu]<w3d.zmminp,arange(pg.nps[js]))    
    if len(ii)==0:
      self.ionstoprev = [0]
    else:
      self.ionstoprev = [len(ii)]
      self.ionstoprev.append(pg.xp[il:iu].copy)
      self.ionstoprev.append(pg.yp[il:iu].copy)
      self.ionstoprev.append(pg.zp[il:iu].copy)
      self.ionstoprev.append(pg.uxp[il:iu].copy)
      self.ionstoprev.append(pg.uyp[il:iu].copy)
      self.ionstoprev.append(pg.uzp[il:iu].copy)
      self.ionstoprev.append(pg.gaminv[il:iu].copy)
      if pg.npid>0:
        self.ionstoprev.append(pg.pid[il:iu,:].copy)
    put(pg.gaminv,ii,0.)
    if self.l_verbose:print me,top.it,self.iz,'exit store_ionstoprev'
    
  def sendrecv_storedions(self):
    if not lparallel:return
    if self.l_verbose:print me,top.it,self.iz,'enter sendrecv_storedions'
    # --- sends stored ions
    if me>0:
      mpi.send(self.ionstoprev,me-1)
      if self.ionstoprev[0]>0:
        print me, 'sends ',self.ionstoprev[0],' ions to ',me-1
    # --- receives stored ions
    if me<npes-1:
      js = w3d.nzp-1
      pg = self.pgions
      recved = mpirecv(me+1)
      np = recved[0]
      if np>0:
        print me, 'recvs ',np,' ions from ',me+1
        if pg.npid>0:
          pid=self.recved[-1]
        else:
          pid=0.
        addparticles(js = js,
                     x  = recved[1],
                     y  = recved[2],
                     z  = recved[3],
                     vx = recved[4],
                     vy = recved[5],
                     vz = recved[6],
                     gi = recved[7],
                     pid=pid,
                     lmomentum=1,
                     lallindomain=1,
                     pgroup=self.pgions)
        self.set_sw(js)
    if self.l_verbose:print me,top.it,self.iz,'exit sendrecv_storedions'
             
  def sendparticlestonext(self,ii):
    if self.l_verbose:print me,top.it,self.iz,'enter sendparticlestonext'
    tosend = []
    tosend.append(len(ii))
    if len(ii)>0:
      tosend.append(take(pg.xp[ii]))
      tosend.append(take(pg.yp[ii]))
      tosend.append(take(pg.zp[ii]))
      tosend.append(take(pg.uxp[ii]))
      tosend.append(take(pg.uyp[ii]))
      tosend.append(take(pg.uzp[ii]))
      tosend.append(take(pg.gaminv[ii]))
      if top.npid>0:tosend.append(take(pg.pid[ii,:]))
    mpi.send(tosend,me+1)
    if self.ionstoprev[0]>0:
      print me, 'sends ',self.ionstoprev[0],' ions to ',me+1
    if len(ii)>0:
      put(pg.gaminv,ii,0.)
    if self.l_verbose:print me,top.it,self.iz,'exit sendparticlestonext'
    
  def recvparticlesfromprevious(self):
    if self.l_verbose:print me,top.it,self.iz,'enter recvparticlesfromprevious'
    pg = self.pgions
    top.pgroup = pg
    recved = mpirecv(me-1)
    np = recved[0]
    if np>0:
      print me, 'recvs ',np,' ions from ',me-1
      if top.npid>0:
        pid=self.recved[8]
      else:
        pid=0.
      addparticles(js = 0,
                   x  = recved[1],
                   y  = recved[2],
                   z  = recved[3],
                   vx = recved[4],
                   vy = recved[5],
                   vz = recved[6],
                   gi = recved[7],
                   pid=pid,
                   lmomentum=1,
                   lallindomain=1,
                   pgroup=self.pgions)
      self.set_sw(0)
    if self.l_verbose:print me,top.it,self.iz,'exit recvparticlesfromprevious'
             
  def deposit_ions(self):
    if self.l_verbose:print me,top.it,self.iz,'enter deposit_ions'
    pg = self.pgions
    top.pgroup = pg
    js = self.iz
    il = pg.ins[js]-1
    iu = il+pg.nps[js]
    bg=frz.basegrid=self.gridions[0]
    mk_grids_ptr()
    # --- reset 2-D rho arrays
    reset_rzmgrid_rho()
    if iu-il>0:
      xp = pg.xp[il:iu]
      yp = pg.yp[il:iu]
      q = pg.sq[js]*pg.sw[js]/w3d.dz
      # --- deposit ion charge in first grid
      rhoweightrz_weights(xp,yp,yp,self.wz0[js],iu-il,q,bg.nr,bg.nz,bg.dr,bg.dz,bg.rmin,0.)
    if self.iz<self.izmax-1:
      bg=frz.basegrid=self.gridions[1]
      mk_grids_ptr()
      if iu-il>0:
        # --- deposit ion charge in second grid
        rhoweightrz_weights(xp,yp,yp,self.wz1[js],iu-il,q,bg.nr,bg.nz,bg.dr,bg.dz,bg.rmin,0.)
      if self.iz==0:self.deposit_ions_last_step()
      # --- distribute rho among patches (WARNING: must be called only once after everything has been deposited.)
      bg=frz.basegrid=self.gridions[1]
      mk_grids_ptr()
      distribute_rho_rz()
      if self.iz==0:
        bg=frz.basegrid=self.gridions[0]
        mk_grids_ptr()
        distribute_rho_rz()
    if self.l_verbose:print me,top.it,self.iz,'exit deposit_ions'
                      
  def deposit_ions_last_step(self):
  # deposits ions on last cell at last step, using ions newly advanced, from last two cell.
    if self.l_verbose:print me,top.it,self.iz,'enter deposit_ions_last_step'
    pg = self.pgions
    top.pgroup = pg
    iz = w3d.nzp-1
    js = iz
    il = pg.ins[js]-1
    iu = il+pg.nps[js]
    np=iu-il
    # --- reset 2-D rho arrays
    ii=[]
    if np>0:
      q = pg.sq[js]*pg.sw[js]/w3d.dz
      zmin=w3d.zmminp+iz*w3d.dz
      wz1i = (pg.zp[il:iu]-zmin)/w3d.dz
      # --- select particles within zmminp <--> zmmaxp
      ii = compress( (wz1i>=0.) and (wz1i<=1.), arange(np))
      np =len(ii)
      if np>0:
        wz1 = take(wz1i,ii)
        ii+=il
        xp = take(pg.xp,ii)
        yp = take(pg.yp,ii)
        # --- deposit ion charge in first grid
        bg=frz.basegrid=self.gridionscp[0]
        mk_grids_ptr()
        reset_rzmgrid_rho()
        rhoweightrz_weights(xp,yp,yp,wz1,np,q,bg.nr,bg.nz,bg.dr,bg.dz,bg.rmin,0.)
      # --- select particles outside bounds (> zmmaxp)
      ii = compress( wz1i>1., arange(np))
      np =len(ii)
      if np>0:
        wz0 = 1.-take(wz1i,ii)
        ii+=il
        xp = take(pg.xp,ii)
        yp = take(pg.yp,ii)
        # --- deposit ion charge in first grid
        bg=frz.basegrid=self.gridionscp[1]
        mk_grids_ptr()
        reset_rzmgrid_rho()
        rhoweightrz_weights(xp,yp,yp,wz0,np,q,bg.nr,bg.nz,bg.dr,bg.dz,bg.rmin,0.)
    # --- stacks rho to be sent to next processor
    if lparallel and me<npes-1:
      tosend=[]
      g0 = self.gridionscp[0]
      g1 = self.gridionscp[1]
      for ig in range(frz.ngrids):
        if ig>0:
          g0 = g0.down
          g1 = g1.down
        tosend.append(g0.rho)
        tosend.append(g1.rho)
      mpi.send(tosend,me+1)
      del tosend
      self.sendparticlestonext(ii)
    if lparallel and me>0:
      g0 = self.gridions[0]
      g1 = self.gridions[1]
      rhorecved = mpirecv(me-1)
      for ig in range(frz.ngrids):
        if ig>0:
          g0 = g0.down
          g1 = g1.down
        g0.rho += rhorecved[ig*2]
        g1.rho += rhorecved[ig*2+1]
      del rhorecved
      self.recvparticlesfromprevious()
    if self.l_verbose:print me,top.it,self.iz,'exit deposit_ions_last_step'
#    if lparallel:mpi.barrier()
                      
  def push_electrons(self):
       if self.l_verbose:print me,top.it,self.iz,'enter push_electrons'
       pg = self.pgelec
       top.pgroup = self.pgelec
       frz.basegrid = self.gridions[1]
       mk_grids_ptr()
       # --- loop over species list
       for sp in [self.electrons]:
        # --- loop over species index
        for js in sp.jslist:
          ng = 1+pg.nps[js]/self.nparpgrp
          for ig in range(ng):
            il = pg.ins[js]-1+self.nparpgrp*ig
            iu = min(il+self.nparpgrp,pg.ins[js]-1+pg.nps[js])
            np = iu-il
            if np==0:continue
            # --- gather self-forces
            if self.l_verbose:print me,top.it,self.iz,'me = ',me,';fieldweight' # MR OK
            pg.ex[il:iu]=0.
            pg.ey[il:iu]=0.
            pg.ez[il:iu]=0.
            fieldweightxz(pg.xp[il:iu],pg.yp[il:iu],pg.ex[il:iu],pg.ey[il:iu],np,0.,top.efetch[js])
            if self.l_verbose:print me,top.it,self.iz,'me = ',me,';epush'
            epush3d(np,pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],
                    pg.ex[il:iu],pg.ey[il:iu],pg.ez[il:iu],pg.sq[js],pg.sm[js],0.5*self.dt)
            gammaadv(np,pg.gaminv[il:iu],pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],
                     top.gamadv,top.lrelativ)
            if self.l_verbose:print me,top.it,self.iz,'me = ',me,';xpush'
            xpush3d(np,pg.xp[il:iu],pg.yp[il:iu],pg.zp[il:iu],
                    pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],pg.gaminv[il:iu],self.dt)
            pg.zp[il:iu]-=w3d.dz
#            pli(frz.basegrid.phi);refresh()
#            ppzx(js=1,color=red,msize=2);refresh()
#            window(3);ppzvx(js=1,msize=2);window(0)
            if self.l_verbose:print me,top.it,self.iz,'me = ',me,';stckxy3d'
            stckxy3d(np,pg.xp[il:iu],w3d.xmmax,w3d.xmmin,w3d.dx,
                     pg.yp[il:iu],w3d.ymmax,w3d.ymmin,w3d.dy,
                     pg.zp[il:iu],w3d.zmminlocal,w3d.dz,
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
            if self.l_verbose:print me,top.it,self.iz,'me = ',me,';fieldweight' # MR OK
            pg.ex[il:iu]=0.
            pg.ey[il:iu]=0.
            pg.ez[il:iu]=0.
            fieldweightxz(pg.xp[il:iu],pg.yp[il:iu],pg.ex[il:iu],pg.ey[il:iu],np,top.zgrid,top.efetch[js])
#            fetche3d(pg,il+1,np,js+1)
            if self.l_verbose:print me,top.it,self.iz,'me = ',me,';epush'
            epush3d(np,pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],
                    pg.ex[il:iu],pg.ey[il:iu],pg.ez[il:iu],pg.sq[js],pg.sm[js],0.5*self.dt)
            gammaadv(np,pg.gaminv[il:iu],pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],
                     top.gamadv,top.lrelativ)
       if self.l_verbose:print me,top.it,self.iz,'exit push_electrons'

  def sort_ions_along_z(self):
    if self.l_verbose:print me,top.it,self.iz,'enter sort_ions_along_z'
    if not self.fullsortdone:
      self.full_sort_ions_along_z()
      return
    pg = self.pgions
    if sum(pg.nps)==0:return
    for js in range(w3d.nzp):
      il = pg.ins[js]-1
      iu = il+pg.nps[js]
      zmin = js*w3d.dz+w3d.zmminp
      iz = int((pg.zp[il:iu]-zmin)/w3d.dz) 
      izleft = compress(iz<0,arange(pg.nps[js]))
      izright = compress(iz>0,arange(pg.nps[js]))
      if js==0 and len(izleft)>0:raise('Error in sort_ions_along_z:js==0 and len(izleft)>0')   
      if js==w3d.nzp-1 and len(izright)>0:raise('Error in sort_ions_along_z:js==w3d.nzp-1 and len(izright)>0')   
      if len(izleft)>0:
        toaddleft = [take(pg.xp,izleft).copy(),
                     take(pg.yp,izleft).copy(),
                     take(pg.zp,izleft).copy(),
                     take(pg.uxp,izleft).copy(),
                     take(pg.uyp,izleft).copy(),
                     take(pg.uzp,izleft).copy(),
                     take(pg.gaminv,izleft).copy()]
        if pg.npid>0:toaddleft.append(take(pg.pid,izleft).copy())
      if len(izright)>0:
        toaddright = [take(pg.xp,izright).copy(),
                     take(pg.yp,izright).copy(),
                     take(pg.zp,izright).copy(),
                     take(pg.uxp,izright).copy(),
                     take(pg.uyp,izright).copy(),
                     take(pg.uzp,izright).copy(),
                     take(pg.gaminv,izright).copy()]
        if pg.npid>0:toaddright.append(take(pg.pid,izright).copy())
      if len(izleft)>0:
        if pg.npid==0:
          pid = 0.
        else:
          pid = toaddleft[-1]
        addparticles(js=js-1,
                     x = toaddleft[0],
                     y = toaddleft[1],
                     z = toaddleft[2],
                     vx = toaddleft[3],
                     vy = toaddleft[4],
                     vz = toaddleft[5],
                     gi = toaddleft[6],
                     pid = pid,
                     lmomentum=1,
                     lallindomain=1)
      if len(izright)>0:
        if pg.npid==0:
          pid = 0.
        else:
          pid = toaddright[-1]
        addparticles(js=js+1,
                     x = toaddright[0],
                     y = toaddright[1],
                     z = toaddright[2],
                     vx = toaddright[3],
                     vy = toaddright[4],
                     vz = toaddright[5],
                     gi = toaddright[6],
                     pid = pid,
                     lmomentum=1,
                     lallindomain=1)
    if self.l_verbose:print me,top.it,self.iz,'exit sort_ions_along_z'

  def full_sort_ions_along_z(self):
    if self.l_verbose:print me,top.it,self.iz,'enter full_sort_ions_along_z'
    pg = self.pgions
    top.pgroup = self.pgions
    if sum(pg.nps)==0:return
    js = 0
    il = pg.ins[js]-1
    iu = il+pg.nps[js]
    x=pg.xp[il:iu].copy()
    y=pg.yp[il:iu].copy()
    z=pg.zp[il:iu].copy()
    ux=pg.uxp[il:iu].copy()
    uy=pg.uyp[il:iu].copy()
    uz=pg.uzp[il:iu].copy()
    gi=pg.gaminv[il:iu].copy()
    if pg.npid>0:
      pid = pg.pid[il:iu,:].copy()
    pg.nps[0]=0
    for iz in range(w3d.nzp):
      zmin = iz*w3d.dz+w3d.zmminp
      ii = compress((z>=zmin) & (z<zmin+w3d.dz),arange(shape(x)[0]))
      if len(ii)>0:
        if pg.npid>0:
          pida = take(pid,ii)
        else:
          pida = 0.
        addparticles(x=take(x,ii),
                     y=take(y,ii),
                     z=take(z,ii),
                     vx=take(ux,ii),
                     vy=take(uy,ii),
                     vz=take(uz,ii),
                     gi=take(gi,ii),
                     js=iz,
                     pid = pida,
                     lmomentum=1,
                     lallindomain=1)
    self.fullsortdone = 1
    if self.l_verbose:print me,top.it,self.iz,'enter full_sort_ions_along_z'

  def compute_sw(self):
    if self.l_verbose:print me,top.it,self.iz,'enter compute_sw'
    del self.wz0,self.wz1
    self.wz0=[]
    self.wz1=[]
    for iz in range(w3d.nzp):
      wz0,wz1 = self.get_sw(iz)
      self.wz0.append(wz0)
      self.wz1.append(wz1)
    if self.l_verbose:print me,top.it,self.iz,'exit compute_sw'

  def get_sw(self,js):
      if self.l_verbose:print me,top.it,self.iz,'enter get_sw'
      pg = self.pgions
      il = pg.ins[js]-1
      iu = il+pg.nps[js]
      if pg.nps[js]==0:
        return None,None
      else:
        zmin=w3d.zmminp+js*w3d.dz
        wz1 = (pg.zp[il:iu]-zmin)/w3d.dz
        return 1.-wz1,wz1
      if self.l_verbose:print me,top.it,self.iz,'exit get_sw'

  def set_sw(self,js):
      if self.l_verbose:print me,top.it,self.iz,'enter set_sw'
      wz1 = self.get_sw(js)
      self.wz1[js] = wz1
      self.wz0[js] = 1.-wz1
      if self.l_verbose:print me,top.it,self.iz,'exit set_sw'
        
  def push_ions(self):
    if self.l_verbose:print me,top.it,self.iz,'enter push_ions'
    pg = self.pgions
    top.pgroup = pg
    if me>=(npes-1-top.it):
#    if me==(npes-1-top.it):
      js = self.iz
      il = pg.ins[js]-1
      iu = il+pg.nps[js]
      np = pg.nps[js]
      if np>0:
        pg.ex[il:iu]=0.
        pg.ey[il:iu]=0.
        pg.ez[il:iu]=0.
        frz.basegrid = self.gridelecs[1]
        mk_grids_ptr()
        fieldweightxz(pg.xp[il:iu],pg.yp[il:iu],pg.ex[il:iu],pg.ey[il:iu],np,top.zgrid,top.efetch[0])
        pg.ex[il:iu]*=self.wz1[js]
        pg.ey[il:iu]*=self.wz1[js]
        if self.iz==0:
          ex = zeros(np,'d')
          ey = zeros(np,'d')
          frz.basegrid = self.gridelecs[0]
          mk_grids_ptr()
          fieldweightxz(pg.xp[il:iu],pg.yp[il:iu],ex,ey,np,top.zgrid,top.efetch[0])
          ex*=self.wz0[js]
          ey*=self.wz0[js]
          pg.ex[il:iu]+=ex
          pg.ey[il:iu]+=ey
          self.advance_ions(js)
      if self.iz<w3d.nzp-1:        
        js = self.iz+1
        il = pg.ins[js]-1
        iu = il+pg.nps[js]
        np = pg.nps[js]
        if np>0:
          ex = zeros(np,'d')
          ey = zeros(np,'d')
          frz.basegrid = self.gridelecs[1]
          mk_grids_ptr()
          fieldweightxz(pg.xp[il:iu],pg.yp[il:iu],ex,ey,np,top.zgrid,top.efetch[0])
          ex*=self.wz0[js]
          ey*=self.wz0[js]
          pg.ex[il:iu]+=ex
          pg.ey[il:iu]+=ey
          self.advance_ions(js)
    if self.l_verbose:print me,top.it,self.iz,'exit push_ions'
             
  def advance_ions(self,js):
    if self.l_verbose:print me,top.it,self.iz,'enter advance_ions'
    pg = self.pgions
    np = pg.nps[js]
    if np==0:return
    il = pg.ins[js]-1
    iu = il+pg.nps[js]
    if self.l_maps:
      # --- apply space_charge kick
      epush3d(np,pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],
                 pg.ex[il:iu], pg.ey[il:iu], pg.ez[il:iu], 
                 pg.sq[js],pg.sm[js],top.dt)
      gammaadv(np,pg.gaminv[il:iu],pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],
               top.gamadv,top.lrelativ)
            # --- apply transfer map 
#      pg.xp[il:iu]*=0.99
#      pg.yp[il:iu]*=0.99
#      return
      apply_simple_map(np,
                       pg.xp[il:iu],
                       pg.yp[il:iu],
                       pg.uxp[il:iu],
                       pg.uyp[il:iu],
                       pg.gaminv[il:iu],
                       self.maps.Mtx,
                       self.maps.Mty,
                       top.vbeam)
    self.apply_bnd_conditions(js)
    if self.l_verbose:print me,top.it,self.iz,'exit advance_ions'
    
  def apply_bnd_conditions(self,js):
    if self.l_verbose:print me,top.it,self.iz,'enter apply_bnd_conditions'
    pg=self.pgions
    il = pg.ins[js]-1
    iu = il+pg.nps[js]
    stckxy3d(pg.nps[js],pg.xp[il:iu],w3d.xmmax,w3d.xmmin,w3d.dx,
                  pg.yp[il:iu],w3d.ymmax,w3d.ymmin,w3d.dy,
                  pg.zp[il:iu],w3d.zmminlocal,w3d.dz,
                  pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],pg.gaminv[il:iu],
                  top.zgrid,top.zbeam,w3d.l2symtry,w3d.l4symtry,top.pboundxy,true)
    if js==0 or js==w3d.nzp-1:
      if js==0:top.pboundnz=-1
      if js==w3d.nzp-1:top.pbound0=-1
      zpartbndwithdata(pg.nps[js],pg.zp[il:iu],pg.uzp[il:iu],pg.gaminv[il:iu],
                       w3d.zmmaxlocal,w3d.zmminlocal,w3d.dz,top.zgrid)
      if js==0:top.pboundnz=0
      if js==w3d.nzp-1:top.pbound0=0
    if self.scraper is not None:self.scraper.scrape(js)
    processlostpart(pg,js+1,top.clearlostpart,top.time+top.dt*pg.ndts[js],top.zbeam)
    if self.l_verbose:print me,top.it,self.iz,'enter apply_bnd_conditions'

  def clear_electrons(self):
    if self.l_verbose:print me,top.it,self.iz,'enter clear_electrons'
    pg = self.pgelec
    top.pgroup = self.pgelec
    # --- clear electrons on processor 0, shift them to the previous ones for me>0
    if lparallel and me>0:
      js = 0
      tosend=[pg.nps[js]]
      if pg.nps[js]>0:
        il=top.pgroup.ins[js]-1
        iu=il+top.pgroup.nps[js]
        tosend.append(pg.xp[il:iu])
        tosend.append(pg.yp[il:iu])
        tosend.append(pg.zp[il:iu])
        tosend.append(pg.uxp[il:iu])
        tosend.append(pg.uyp[il:iu])
        tosend.append(pg.uzp[il:iu])
        tosend.append(pg.gaminv[il:iu])
        if top.npid>0:tosend.append(pg.pid[:,il:iu])
      mpi.send(tosend,me-1)
    pg.nps[0]=0      
    if lparallel and me<npes-1:
      self.recved = mpirecv(me+1)
      np=self.recved[0]
      if np>0:
        if top.npid>0:
          pid=self.recved[8]
        else:
          pid=0.
        self.electrons.addpart(self.recved[1],
                              self.recved[2],
                              w3d.zmmaxp-1.e-10*w3d.dz,#
                              #self.recved[3],
                              self.recved[4],
                              self.recved[5],
                              self.recved[6],
                              self.recved[7],
                              pid=pid,
                              lmomentum=1,
                              lallindomain=1)
    if self.l_verbose:print me,top.it,self.iz,'exit clear_electrons'
#    if lparallel:mpi.barrier()
    
  def plot_electrons(self):
        if self.l_verbose:print me,top.it,self.iz,'enter plot_electrons'
#        ppgeneric(getvx(js=1),getx(js=1))
        window(0);fma();
#        frz.basegrid=self.gridions[1]
#        mk_grids_ptr()
#        plphirz();refresh()
#        window(0);fma();
#        n=getn(js=1)
#        if n>0:
#          x = self.slist[0].getx()
#          y = self.slist[0].gety()
#          ex = self.slist[0].getex()
#          if me==0:ppco(x,y,ex,msize=3,ncolor=100);refresh()
#        if iz==0:
        self.electrons.ppxex(msize=2);refresh()
#        else:
#          self.slist[0].ppxex(msize=2);
#        limits(w3d.xmmin,w3d.xmmax,w3d.ymmin,w3d.ymmax)
#        ppxex(js=1)
#        pfxy(iz=iz);
#        fma()
        if self.l_verbose:print me,top.it,self.iz,'enter plot_electrons'

  def create_electrons(self):
     if self.l_verbose:print me,top.it,self.iz,'enter create_electrons'
     pg = self.pgelec
     top.pgroup = self.pgelec
     # --- generate electrons (on last processor only)
     pg.nps[0]=0
     if self.l_mode==2:
      zmax = w3d.zmmax/2+top.zgrid
     else:
      zmax = w3d.zmmax-1.e-10*w3d.dz+top.zgrid
     if top.prwall<sqrt(w3d.xmmax**2+w3d.ymmax**2):
      self.electrons.add_uniform_cylinder(self.Ninit,min(top.prwall,w3d.xmmax),zmax,zmax,1.e-10,1.e-10,1.e-10,lallindomain=true)
     else:
      if self.l_elecuniform:
        spacing='uniform'
      else:
        spacing='random'
      self.electrons.add_uniform_box(self.Ninit,w3d.xmmin,w3d.xmmax,
                                    w3d.ymmin,w3d.ymmax,zmax,zmax,
                                    1.e-10,1.e-10,1.e-10,spacing=spacing,lallindomain=true)
     # --- scrape electrons 
     if self.scraper is not None:self.scraper.scrapeall(local=1,clear=1)
#     shrinkpart(pg)
     if self.l_verbose:print me,top.it,self.iz,'exit create_electrons'

  def deposit_electrons(self,i=1):
      if self.l_verbose:print me,top.it,self.iz,'enter deposit_electrons'
      pg = self.pgelec
      top.pgroup = self.pgelec
      bg = frz.basegrid = self.gridelecs[i]
      mk_grids_ptr()
      # --- reset 2-D rho arrays
      reset_rzmgrid_rho()
      # --- loop over species list
      for sp in [self.electrons]:
        # --- loop over species index
        for js in sp.jslist:
          np = pg.nps[js]
          if np==0:continue
          il = pg.ins[js]-1
          iu = il+pg.nps[js]
          if self.l_verbose:print me,top.it,self.iz,'me = ',me,';deposit rho' # MR OK
          # --- deposit rho
          rhoweightrz(pg.xp[il:iu],pg.yp[il:iu],pg.yp[il:iu],np,
                      pg.sq[js]*pg.sw[js],bg.nr,bg.nz,
                      bg.dr,bg.dz,bg.rmin,0.)
      # --- distribute rho among patches
      distribute_rho_rz()
      if self.l_verbose:print me,top.it,self.iz,'exit deposit_electrons'

  def add_ei_fields(self,i=1):
    if self.l_verbose:print me,top.it,self.iz,'enter add_ei_fields'
    if self.l_selfe:
        ge = self.gridelecs[i]
        gi = self.gridions[i]
        for ig in range(frz.ngrids):
          if ig>0:
            ge = ge.down
            gi = gi.down
          gi.phi[...] += ge.phi[...] 
    if frz.l_get_fields_on_grid:
        ge = self.gridelecs[i]
        gi = self.gridions[i]
        for g in [ge,gi]:
          frz.basegrid = g
          mk_grids_ptr()
          getallfieldsfromphip()
    if self.l_verbose:print me,top.it,self.iz,'exit add_ei_fields'

  def getmmnts(self):
    if self.l_verbose:print me,top.it,self.iz,'enter getmmnts'
    pg = self.pgions
    self.pnum.append(sum(self.pgions.nps))
    xbar = 0.
    ybar = 0.
    zbar = 0.
    xpbar = 0.
    ypbar = 0.
    x2 = 0.
    y2 = 0.
    z2 = 0.
    xp2 = 0.
    yp2 = 0.
    xxp = 0.
    yyp = 0.
    for js in range(self.pgions.ns):
      x = getx(js=js,gather=0)
      y = gety(js=js,gather=0)
      z = getz(js=js,gather=0)
      xp = getvx(js=js,gather=0)/getvz(js=js,gather=0)
      yp = getvy(js=js,gather=0)/getvz(js=js,gather=0)
      xbar+=sum(x)      
      ybar+=sum(y)      
      zbar+=sum(z)      
      xpbar+=sum(xp)      
      ypbar+=sum(yp)      
      x2+=sum(x*x)      
      y2+=sum(y*y)      
      z2+=sum(z*z)      
      xp2+=sum(xp*xp)      
      yp2+=sum(yp*yp)      
      xxp+=sum(x*xp)      
      yyp+=sum(y*yp)      
    self.xbar.append(xbar)
    self.ybar.append(ybar)
    self.zbar.append(zbar)
    self.xpbar.append(xpbar)
    self.ypbar.append(ypbar)
    self.x2.append(x2)
    self.y2.append(y2)
    self.z2.append(z2)
    self.xp2.append(xp2)
    self.yp2.append(yp2)
    self.xxp.append(xxp)
    self.yyp.append(yyp)
    if top.it==0 or not lparallel:
      self.timemmnts.append(top.time)
    else:
      self.timemmnts.append(top.time+(me-npes+1)*top.dt)
    if self.l_verbose:print me,top.it,self.iz,'exit getmmnts'
      
  def getpnum(self):
    if not lparallel:
      return self.pnum.data()
    else:
      if me==0:
        return self.pnum.data()   
      else:
        return self.pnum.data()[:-me]  

  def getxrms(self):
    if not lparallel:
      return sqrt(self.x2.data()/self.pnum.data())
    else:
      if me==0:
        return sqrt(parallelsum(self.x2.data()[:])/parallelsum(self.pnum.data()[:]))   
      else:
        return sqrt(parallelsum(self.x2.data()[:-me])/parallelsum(self.pnum.data()[:-me]))   

  def getyrms(self):
    if not lparallel:
      return sqrt(self.y2.data()/self.pnum.data())
    else:
      if me==0:
        return sqrt(parallelsum(self.y2.data()[:])/parallelsum(self.pnum.data()[:]))   
      else:
        return sqrt(parallelsum(self.y2.data()[:-me])/parallelsum(self.pnum.data()[:-me]))   

  def getzrms(self):
    if not lparallel:
      return sqrt(self.z2.data()/self.pnum.data())
    else:
      if me==0:
        return sqrt(parallelsum(self.z2.data()[:])/parallelsum(self.pnum.data()[:]))   
      else:
        return sqrt(parallelsum(self.z2.data()[:-me])/parallelsum(self.pnum.data()[:-me]))   

  def getxbar(self):
    if not lparallel:
      return self.xbar.data()/self.pnum.data()
    else:
      if me==0:
        return (parallelsum(self.xbar.data()[:])/parallelsum(self.pnum.data()[:]))   
      else:
        return (parallelsum(self.xbar.data()[:-me])/parallelsum(self.pnum.data()[:-me]))   

  def getybar(self):
    if not lparallel:
      return self.ybar.data()/self.pnum.data()
    else:
      if me==0:
        return (parallelsum(self.ybar.data()[:])/parallelsum(self.pnum.data()[:]))   
      else:
        return (parallelsum(self.ybar.data()[:-me])/parallelsum(self.pnum.data()[:-me]))   

  def getzbar(self):
    if not lparallel:
      return self.zbar.data()/self.pnum.data()
    else:
      if me==0:
        return (parallelsum(self.zbar.data()[:])/parallelsum(self.pnum.data()[:]))   
      else:
        return (parallelsum(self.zbar.data()[:-me])/parallelsum(self.pnum.data()[:-me]))   

  def getxpbar(self):
    if not lparallel:
      return self.xpbar.data()/self.pnum.data()
    else:
      if me==0:
        return (parallelsum(self.xpbar.data()[:])/parallelsum(self.pnum.data()[:]))   
      else:
        return (parallelsum(self.xpbar.data()[:-me])/parallelsum(self.pnum.data()[:-me]))   

  def getypbar(self):
    if not lparallel:
      return self.ypbar.data()/self.pnum.data()
    else:
      if me==0:
        return (parallelsum(self.ypbar.data()[:])/parallelsum(self.pnum.data()[:]))   
      else:
        return (parallelsum(self.ypbar.data()[:-me])/parallelsum(self.pnum.data()[:-me]))   

  def getemitxrms(self):
    xbar = self.getxbar()
    xpbar = self.getxpbar()
    if not lparallel:
      pnum = self.pnum.data()
      x2  = self.x2.data()/pnum
      xp2 = self.xp2.data()/pnum
      xxp = self.xxp.data()/pnum
    else:
      if me==0:
        pnum = parallelsum(self.pnum.data())
        x2  = parallelsum(self.x2.data())/pnum
        xp2 = parallelsum(self.xp2.data())/pnum
        xxp = parallelsum(self.xxp.data())/pnum
      else:
        pnum = parallelsum(self.pnum.data()[:-me])
        x2  = parallelsum(self.x2.data()[:-me])/pnum
        xp2 = parallelsum(self.xp2.data()[:-me])/pnum
        xxp = parallelsum(self.xxp.data()[:-me])/pnum
    return sqrt((x2-xbar*xbar)*(xp2-xpbar*xpbar)-(xxp-xbar*xpbar)**2)
        
  def getemityrms(self):
    ybar = self.getybar()
    ypbar = self.getypbar()
    if not lparallel:
      pnum = self.pnum.data()
      y2  = self.y2.data()/pnum
      yp2 = self.yp2.data()/pnum
      yyp = self.yyp.data()/pnum
    else:
      if me==0:
        pnum = parallelsum(self.pnum.data())
        y2  = parallelsum(self.y2.data())/pnum
        yp2 = parallelsum(self.yp2.data())/pnum
        yyp = parallelsum(self.yyp.data())/pnum
      else:
        pnum = parallelsum(self.pnum.data()[:-me])
        y2  = parallelsum(self.y2.data()[:-me])/pnum
        yp2 = parallelsum(self.yp2.data()[:-me])/pnum
        yyp = parallelsum(self.yyp.data()[:-me])/pnum
    return sqrt((y2-ybar*ybar)*(yp2-ypbar*ypbar)-(yyp-ybar*ypbar)**2)
        
  def plfrac(self,color=black,width=1,type='solid',xscale=1.,yscale=1.):  
    qspnum = self.getpnum()
    if me==0:
      qst = self.timemmnts.data()
      pla(yscale*qspnum/qspnum[0],xscale*qst,color=color,width=width)

  def plxrms(self,color=black,width=1,type='solid',xscale=1.,yscale=1.):  
    qsxrms = self.getxrms()
    if me==0:
      qst = self.timemmnts.data()
      pla(yscale*qsxrms,xscale*qst,color=color,width=width)

  def plyrms(self,color=black,width=1,type='solid',xscale=1.,yscale=1.):  
    qsyrms = self.getyrms()
    if me==0:
      qst = self.timemmnts.data()
      pla(yscale*qsyrms,xscale*qst,color=color,width=width)

  def plzrms(self,color=black,width=1,type='solid',xscale=1.,yscale=1.):  
    qszrms = self.getzrms()
    if me==0:
      qst = self.timemmnts.data()
      pla(yscale*qszrms,xscale*qst,color=color,width=width)

  def plxbar(self,color=black,width=1,type='solid',xscale=1.,yscale=1.):  
    qsxbar = self.getxbar()
    if me==0:
      qst = self.timemmnts.data()
      pla(yscale*qsxbar,xscale*qst,color=color,width=width)

  def plybar(self,color=black,width=1,type='solid',xscale=1.,yscale=1.):  
    qsybar = self.getybar()
    if me==0:
      qst = self.timemmnts.data()
      pla(yscale*qsybar,xscale*qst,color=color,width=width)

  def plzbar(self,color=black,width=1,type='solid',xscale=1.,yscale=1.):  
    qszbar = self.getzbar()
    if me==0:
      qst = self.timemmnts.data()
      pla(yscale*qszbar,xscale*qst,color=color,width=width)

  def plemitx(self,color=black,width=1,type='solid',xscale=1.,yscale=1.):  
    qsexrms = self.getemitxrms()
    if me==0:
      qst = self.timemmnts.data()
      pla(yscale*qsexrms,xscale*qst,color=color,width=width)

  def plemity(self,color=black,width=1,type='solid',xscale=1.,yscale=1.):  
    qseyrms = self.getemityrms()
    if me==0:
      qst = self.timemmnts.data()
      pla(yscale*qseyrms,xscale*qst,color=color,width=width)

class Quasistaticold:
  def __init__(self,slist=None,MRroot=None,beam=None,l_verbose=0,l_findmgparam=0,
               nparpgrp=top.nparpgrp,Ninit=1000,l_mode=1,l_selfe=1,maps=None,pboundxy=None,
               conductors=[],l_elecuniform=0,scraper=None,l_weakstrong=0,nelecperiod=1):
    assert beam is not None, 'beam must be passed as argument of Quasistatic'
    self.beam=beam
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
      frz.init_base(w3d.nx,w3d.nz,w3d.dx,w3d.dz,w3d.xmmin,w3d.zmminlocal,false)
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
    self.MRroot = MRroot
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
    self.l_weakstrong=l_weakstrong
    self.nelecperiod=nelecperiod
    self.scraper=scraper
    if maps is not None:
      self.l_maps=true
      self.maps=maps
    if pboundxy is None:
      self.pboundxy = top.pboundxy
    else:
      self.pboundxy = pboundxy
      
  def step(self,nt=1,l_plotelec=0):
   # --- if running in parallel, enforcing l_mode=1
   if npes>1:self.l_mode=1
   for it in range(nt):
     if (not self.l_weakstrong or top.it==0) and (top.it%self.nelecperiod==0):
       loadrho()
       fieldsol()
       # --- switch to 2-D solver
       self.solver_3d_to_2d()
       # --- generate electrons (on last processor only)
       if me==max(0,npes-1) and (top.it==0  or not l_plotelec):self.create_electrons()
       # --- push 2-D electrons slice along bunch
       self.push_electrons(l_plotelec=l_plotelec)
       # --- switch back to 3-D solver
       self.solver_2d_to_3d()
     # --- push ions
     self.push_ions()
     # --- update time
     top.time+=top.dt
     # --- compute moments
     if self.l_verbose: print me,top.it,self.iz,'me = ',me,';compute zmmnt'
     zmmnt()
     minidiag(top.it,top.time,top.lspecial)
     # --- update time counter
     top.it+=1

  def push_ions(self):
    if me>=(npes-1-top.it):
#    if me==(npes-1-top.it):
     if self.l_maps:
      if self.l_verbose: print me,top.it,self.iz,'me = ',me,';apply space_charge kick'
      self.maps.apply_space_charge_kick(self.beam)
      if self.l_verbose: print me,top.it,self.iz,'me = ',me,';apply transfer map'
      self.maps.apply_transfer_map(self.beam)
      if self.l_verbose: print me,top.it,self.iz,'me = ',me,';apply bnd conditions'
      self.maps.apply_bnd_conditions(self.beam)

  def save_ions(self,zmax):
    # -- save data for ions with z<zmax
    # --- set shortcuts
    pg = top.pgroup
    js = 0
    il = pg.ins[js]-1
    iu = il+pg.nps[js]
    self.isaved = compress(pg.zp[il:iu]<zmax,il+arange(pg.nps[js]))
    self.xsaved = take(pg.xp[il:iu],self.isaved)
    self.ysaved = take(pg.yp[il:iu],self.isaved)
    self.zsaved = take(pg.zp[il:iu],self.isaved)
    self.uxsaved = take(pg.uxp[il:iu],self.isaved)
    self.uysaved = take(pg.uyp[il:iu],self.isaved)
    self.uzsaved = take(pg.uzp[il:iu],self.isaved)
    self.gisaved = take(pg.gaminv[il:iu],self.isaved)
    if pg.npid>0:self.pidsaved = take(pg.pid[il:iu,:],self.isaved)
    
  def push_electrons(self,l_return_dist=false,l_plotelec=0):
    # --- set shortcuts
    pg = top.pgroup
    bg = frz.basegrid

    # --- (diagnostic) stacking of electron distribution
    if l_return_dist:
      self.dist=[]
      sp=self.slist[0]
      self.dist.append([sp.getx(),sp.gety(),sp.getvx(),sp.getvy()])

    # --- sets z range
    if self.l_mode==2:
      izmin = w3d.nz/4
      izmax = w3d.nz-izmin
    else:
      izmin = 0
      izmax = w3d.nzp

    # --- loop over 3-D grid mesh in z
    for iz in range(izmax,izmin-1,-1):
      self.iz=iz
      if self.l_verbose:print me,top.it,self.iz,'me = ',me,';iz = ',iz
      # --- solve 2d field
      self.solve2dfield()
      # --- switch 2d and 3d slice of potentials
      self.switch2d3d(iz)          
      if l_plotelec and iz<w3d.nz/max(1,npes):self.plot_electrons()
      # --- push particles but on last step
      if(iz>0):self.push_electrons_2d()
      # --- (diagnostic) stacking of electron distribution
      if l_return_dist:
        sp=self.slist[0]
        self.dist.append([sp.getx(),sp.gety(),sp.getvx(),sp.getvy()])

    # --- clear electrons on processor 0, shift them to the previous ones for me>0
    self.clear_electrons()

    # --- fills end of 3-D phi arrays
    if self.l_mode==1:
      self.fillphizends(1)
    else:
      self.fillphizends(izmin)

#    if self.MRroot is not None:
#      self.setidosolve(self.MRroot,0)
##      self.MRroot.getallpotentialpforparticles()    
##      fieldsol()
#      self.MRroot.solve()
#      self.setidosolve(self.MRroot,1)      
      
  def clear_electrons(self):
    # --- clear electrons on processor 0, shift them to the previous ones for me>0
    if me==0:
      pg.nps[1]=0      
    if npes>1:
      js = 1
      if pg.nps[js]>0:
        il=top.pgroup.ins[js]-1
        iu=il+top.pgroup.nps[js]
        top.pgroup.zp[il:iu]=w3d.zmminlocal-1.e-10*w3d.dz
      zpartbnd(top.pgroup,w3d.zmmaxlocal,w3d.zmminlocal,w3d.dz)

  def plot_electrons(self):
#        ppgeneric(getvx(js=1),getx(js=1))
        window(0);fma();
#        ppco(getx(js=1),gety(js=1),getex(js=1),msize=2,ncolor=100);refresh()
#        if iz==0:
        self.slist[0].ppxex(msize=2);refresh()
#        else:
#          self.slist[0].ppxex(msize=2);
#        limits(w3d.xmmin,w3d.xmmax,w3d.ymmin,w3d.ymmax)
#        ppxex(js=1)
#        pfxy(iz=iz);
#        fma()

  def create_electrons(self):
     # --- define shortcuts
     pg = top.pgroup
     bg = frz.basegrid
     # --- generate electrons (on last processor only)
     pg.nps[1]=0
     if self.l_mode==2:
      zmax = w3d.zmmax/2+top.zgrid
     else:
      zmax = w3d.zmmax-1.e-10*w3d.dz+top.zgrid
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
     # --- scrape electrons 
     if self.scraper is not None:self.scraper.scrapeall(local=1,clear=1)
#     shrinkpart(pg)

  def solve2dfield(self,l_findmgparam=false):
      if self.l_verbose:print me,top.it,self.iz,'me = ',me,'; enter fieldsolve'   # MR OK
      # --- define shortcuts
      pg = top.pgroup
      bg = frz.basegrid
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
          if self.l_verbose:print me,top.it,self.iz,'me = ',me,';deposit rho' # MR OK
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
      
      # --- optimize mgparam
      if self.l_findmgparam and top.it==0 and self.iz==w3d.nz/2:find_mgparam_rz(false)

      # --- call 2-D field solver
      solve_mgridrz(bg,frz.mgridrz_accuracy,true)
      if self.l_verbose:print me,top.it,self.iz,'me = ',me,'; exit fieldsolve'   # MR OK

  def solver_3d_to_2d(self):
    # --- switch to 2-D solver
    if self.l_rz:
      self.gridcp=frz.basegrid
      frz.basegrid=self.gridr
      bg=self.gridr
      self.solvergeomcp = w3d.solvergeom
      self.fstypecp = top.fstype
      w3d.solvergeom = w3d.Rgeom
    else:
      self.solvergeomcp = w3d.solvergeom
      self.fstypecp = top.fstype
      w3d.solvergeom = w3d.XYgeom

  def solver_2d_to_3d(self):
    # --- switch to 3-D solver
    w3d.solvergeom = self.solvergeomcp
    top.fstype = self.fstypecp
    if self.l_rz:
      frz.basegrid=self.gridcp
      bg=frz.basegrid

  def switch2d3d(self,iz):
      # --- define shortcuts
      pg = top.pgroup
      bg = frz.basegrid
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
              brz.phi[1:-1,1:-1] += b3d.potentialarray[:,:,iz+1+enzl,0,0]
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
      # --- update guard cells of 2-D solver
      updateguardcells2d()
      if frz.l_get_fields_on_grid:getallfieldsfromphip()

  def push_electrons_2d(self):
       # --- define shortcuts
       pg = top.pgroup
       bg = frz.basegrid
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
            if self.l_verbose:print me,top.it,self.iz,'me = ',me,';fieldweight' # MR OK
            pg.ex[il:iu]=0.
            pg.ey[il:iu]=0.
            pg.ez[il:iu]=0.
            if self.l_rz:
              fieldweightr(pg.xp[il:iu],pg.yp[il:iu],pg.ex[il:iu],pg.ey[il:iu],np)
            else:
              fieldweightxz(pg.xp[il:iu],pg.yp[il:iu],pg.ex[il:iu],pg.ey[il:iu],np,0.,top.efetch[js])
#            fetche3d(pg,il+1,np,js+1)
            if self.l_verbose:print me,top.it,self.iz,'me = ',me,';epush'
            epush3d(np,pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],
                    pg.ex[il:iu],pg.ey[il:iu],pg.ez[il:iu],pg.sq[js],pg.sm[js],0.5*self.dt)
            gammaadv(np,pg.gaminv[il:iu],pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],
                     top.gamadv,top.lrelativ)
            if self.l_verbose:print me,top.it,self.iz,'me = ',me,';xpush'
            xpush3d(np,pg.xp[il:iu],pg.yp[il:iu],pg.zp[il:iu],
                    pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],pg.gaminv[il:iu],self.dt)
            pg.zp[il:iu]-=w3d.dz
#            pli(frz.basegrid.phi);refresh()
#            ppzx(js=1,color=red,msize=2);refresh()
#            window(3);ppzvx(js=1,msize=2);window(0)
            if self.l_verbose:print me,top.it,self.iz,'me = ',me,';stckxy3d'
            stckxy3d(np,pg.xp[il:iu],w3d.xmmax,w3d.xmmin,w3d.dx,
                     pg.yp[il:iu],w3d.ymmax,w3d.ymmin,w3d.dy,
                     pg.zp[il:iu],w3d.zmminlocal,w3d.dz,
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
            if self.l_verbose:print me,top.it,self.iz,'me = ',me,';fieldweight' # MR OK
            pg.ex[il:iu]=0.
            pg.ey[il:iu]=0.
            pg.ez[il:iu]=0.
            if self.l_rz:
              fieldweightr(pg.xp[il:iu],pg.yp[il:iu],pg.ex[il:iu],pg.ey[il:iu],np)
            else:
              fieldweightxz(pg.xp[il:iu],pg.yp[il:iu],pg.ex[il:iu],pg.ey[il:iu],np,top.zgrid,top.efetch[js])
#            fetche3d(pg,il+1,np,js+1)
            if self.l_verbose:print me,top.it,self.iz,'me = ',me,';epush'
            epush3d(np,pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],
                    pg.ex[il:iu],pg.ey[il:iu],pg.ez[il:iu],pg.sq[js],pg.sm[js],0.5*self.dt)
            gammaadv(np,pg.gaminv[il:iu],pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],
                     top.gamadv,top.lrelativ)

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
            enzu=b3d.extradimsupper[2]+w3d.nzlocal-w3d.nzp
          else:
            enzl=0
            enzu=w3d.nzlocal-w3d.nzp
#            enzu=1
          for iz in range(nz+enzl):
#            b3d.phi[:,:,iz] = b3d.phi[:,:,nz+enzl]
            b3d.potentialparray[:,:,iz,0,0] = b3d.potentialparray[:,:,nz+enzl,0,0]
          for iz in range(nz+enzu):
#            b3d.phi[:,:,-iz-1] = b3d.phi[:,:,-nz-1-enzu]
            b3d.potentialparray[:,:,-iz-1,0,0] = b3d.potentialparray[:,:,-nz-1-enzu,0,0]
      else:
        if self.l_rz:
          for iz in range(nz):
            bg.phi[:,iz] = bg.phi[:,nz]
            bg.phi[:,-iz-1] = bg.phi[:,-nz-1]
        else:
          for iz in range(nz):
            w3d.phi[:,:,iz] = w3d.phi[:,:,nz]
            w3d.phi[:,:,-iz-1] = w3d.phi[:,:,-nz-1]


