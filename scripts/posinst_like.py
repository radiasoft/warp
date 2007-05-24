from pos import *
from Secondaries import *
from appendablearray import *

class Posinst_Like:
  """
Class for generating photo-electrons
 - posinst_file: name of Posinst input file 
 - l_verbose   : sets verbosity (default=0). 
  """
  def __init__(self,posinst_file=None,l_verbose=0,l_3d=0,nx=64,ny=None,
               xmin=None,xmax=None,ymin=None,ymax=None,weight=None,
               conductors=None,l_secondaries=1,dtfact=1,l_switchyz=0, 
               electronnmax=None,l_posmgsolver=0,l_posscatter=0):
    top.lrelativ = true
    if posinst_file is not None:
      init_posinst_for_warp(posinst_file)
      self.nbuckets = pos.lastbkt+1   # nb of buckets
      self.nkicks   = pos.nkicks      # nb of kicks/bucket
      self.nsteps_gap=pos.nsteps
      self.beamnp   = pos.xnpnom      # nb of particles/bunch
      self.beamen   = pos.beamen      # nominal energy bunch
      self.sigz     = pos.sigz        # beam RMS length
      self.bucket_train = pos.rbktarr[:self.nbuckets]
      self.Lambda   = self.beamnp*echarge*pos.wght[:pos.nkicks]/(pos.beamvel*pos.dt[0])
#      self.Lambda/=10.
      self.sigx     = pos.sigx      
      self.sigy     = pos.sigy 
      top.prwall    = 2.*max(pos.ach,pos.bch,0.5)
      if pos.ispch==3:
        if nx is None:nx=pos.imax-1
        if ny is None:ny=pos.jmax-1
      if l_switchyz:
        top.bz0      = pos.bfield
      else:
        top.by0      = pos.bfield
      self.l_posmgsolver=l_posmgsolver
      self.l_posscatter =l_posscatter
    installothereuser(self.beam_kick)
#    installothereuser(self.external_field)
    if electronnmax is not None:
      top.wpid = nextpid()  # inserting particle weight as a variable in array pid
      self.ncull=0
      self.electronnmax=electronnmax
      self.wherecull=AppendableArray(typecode='d')
      self.qph=AppendableArray(typecode='d')
      self.qsec=AppendableArray(typecode='d')
      installbeforelr(self.randomcull)
    if weight is None:weight=pos.chmphelnom#/pos.slength#*self.sigz
    self.addkick=0
    self.phelectrons=Species(type=Electron,weight=weight)
    if l_secondaries:
      self.secelectrons=Species(type=Electron,weight=weight)
    self.PHEL=PhotoElectrons(l_verbose=0,l_switchyz=l_switchyz)
    self.PHEL.add(emitted_species=self.phelectrons)
#    top.lresetparticlee = false
#    top.lresetparticleb = false
    self.nparpgrp = 4096
    self.l_3d=l_3d
    self.dtfact=dtfact
    self.l_switchyz=l_switchyz
    if xmin is None:xmin=-pos.ach#*1.05
    if xmax is None:xmax= pos.ach#*1.05
    if ymin is None:ymin=-pos.bch#*1.05
    if ymax is None:ymax= pos.bch#*1.05
    w3d.nx=nx
    if ny is None:
      ny=nint(pos.bch*nx/pos.ach)
    w3d.xmmin=xmin
    w3d.xmmax=xmax
    if l_switchyz:
      w3d.zmmin=ymin
      w3d.zmmax=ymax
      w3d.nz=ny+0
    else:
      w3d.ymmin=ymin
      w3d.ymmax=ymax
      w3d.ny=ny+0
    if l_3d:
      pass
    else:
      if pos.ispch==0:
        frz.mgridrz_ncmax=0
        w3d.nx=1
        top.depos='none'
        if l_switchyz:
#          w3d.nz=1
          pass
        else:
          w3d.ny=1
      if l_switchyz:
        w3d.ymmin=-0.5#*pos.slength
        w3d.ymmax=-w3d.ymmin
        w3d.solvergeom=w3d.XZgeom
      else:
        w3d.zmmin=-0.5#*pos.slength
        w3d.zmmax=-w3d.zmmin
        w3d.solvergeom=w3d.XYgeom
      if self.l_posmgsolver:
        frz.mgridrz_ncmax=0
        installafterfs(self.pos_fieldsol)
      if self.l_posscatter:
        installafterfs(self.pos_scatter)
        installothereuser(self.pos_electronkick)
    package('w3d');generate()
    self.ibk=0
    if conductors is None:
      if l_switchyz:
        self.pipe = YCylinderEllipticOut(ellipticity = pos.bch/pos.ach,
                                                   radius      = pos.ach,
                                                   length      = (w3d.ymmax-w3d.ymmin)*10.,
                                                   ycent       = 0.5*(w3d.ymmin+w3d.ymmax),
                                                   condid      = 1)
        self.pipescraper = -Sphere(radius=pos.ach,condid=1)
        self.scrapegrid=Grid(nx=nx,ny=w3d.nz,nz=ny,nzfull=ny)
      else:
        self.pipe = ZCylinderEllipticOut(ellipticity = pos.bch/pos.ach,
                                                   radius      = pos.ach,
                                                   length      = w3d.zmmaxglobal-w3d.zmminglobal,
                                                   zcent       = 0.5*(w3d.zmmin+w3d.zmmax),
                                                   condid      = 1)
        self.scrapegrid=Grid(nx=nx,ny=ny,nz=w3d.nz,nzfull=w3d.nz)
        self.pipescraper = self.pipe
        import __main__
        __main__.sc=self.scrapegrid
      installconductors(self.pipe)
      self.scraper = ParticleScraper(self.pipescraper,
                                     lsaveintercept=1,
                                     lsavecondid=1,
                                     lcollectlpdata=1,
                                     grid=self.scrapegrid,
                                     lrefineallintercept=0)
#            self.scraper.l_print_timing=1
      self.conductors = [self.pipe]
    if l_secondaries:
      def set_params(maxsec,mat_num):
        pass
#        self.Sec.enpar=pos.enpar
#        self.Sec.pnpar=pos.pnpar
      self.Sec = Secondaries(min_age=None,set_params_user=set_params,vmode=1)
      self.set_params=set_params
      self.Sec.set_params_user=self.set_params
      self.Sec.enpar=pos.enpar
      self.Sec.pnpar=pos.pnpar
      self.Sec.add(incident_species = self.phelectrons,
                   emitted_species  = self.secelectrons,
                   conductor        = self.pipe,
                   interaction_type = 0)   
      self.Sec.add(incident_species = self.secelectrons,
                   emitted_species  = self.secelectrons,
                   conductor        = self.pipe,
                   interaction_type = 0)   

  def push_buckets(self):
    for ibk in range(self.nbuckets):
      self.push_bucket()

  def push_bucket(self):
    print '    *** bucket %g out of %g'%(self.ibk+1,self.nbuckets)
    for k in range(self.nkicks+self.nsteps_gap):
      self.ikick=k
      if k<self.nkicks:
        self.addkick=1
        self.PHEL.Lambda=self.Lambda[k]*self.bucket_train[self.ibk]
        top.dt=pos.dt[k]
      else:
        self.addkick=0
        self.PHEL.Lambda=0.
        top.dt=pos.deltat_g
      top.dt/=self.dtfact
      for i in range(self.dtfact):
        step()
        self.zero_z()
    self.ibk+=1
      
  def zero_z(self):
    pg = top.pgroup
    if sum(pg.nps)==0:return
    for js in range(pg.ns):
      if pg.nps[js]==0:continue
      ng = 1+pg.nps[js]/self.nparpgrp
      for ig in range(ng):
        il = pg.ins[js]-1+self.nparpgrp*ig
        iu = min(il+self.nparpgrp,pg.ins[js]-1+pg.nps[js])
        if self.l_switchyz:
          pg.yp[il:iu]=0.
        else:
          pg.zp[il:iu]=0.

  def beam_kick(self):
    if  self.ikick>=self.nkicks:return
#    print 'beam_kick',w3d.jmin,w3d.jmax
    fact=1.#/(1.*pi)
    # set linear charge density
    Lambda = self.Lambda[self.ikick]*self.bucket_train[self.ibk]
    pg = top.pgroup
    exim = zeros(w3d.jmax-w3d.jmin,'d')
    eyim = zeros(w3d.jmax-w3d.jmin,'d')
    x0=y0=0.
    il = w3d.jmin
    iu = w3d.jmax
    np = iu-il
    x = pg.xp[il:iu]
    ex = pg.ex[il:iu]
    if self.l_switchyz:
      y = pg.zp[il:iu]
      ey = pg.ez[il:iu]
    else:
      y = pg.yp[il:iu]
      ey = pg.ey[il:iu]
#    if(iden_xy==0) zdir=zef_direct_point(x,y,x0,y0)
    if(pos.iden_xy==1): 
      exbe,eybe = Bassetti_Erskine(fact*Lambda,x,y,x0,y0,self.sigx,self.sigy)
      ex[...]+=exbe
      ey[...]+=eybe
#    if(iden_xy==2) zdir=zef_direct_unifell(x,y,x0,y0,abeam,bbeam)
#    if(iden_xy==3) zdir=zef_direct_parabell(x,y,x0,y0,abeam,bbeam)
    if(pos.iim==1):
#    if(iim==1.and.ichsh==2) zim=zef_image_rect(x,y,x0,y0,ach,bch)
#    if(iim==1.and.ichsh==3) zim=zef_image_open_H(x,y,x0,y0,bch)
      if pos.ichsh==1: 
        ef_image_ell(np,exim,eyim,x,y,x0,y0,pos.ach,pos.bch)
      ex[...]+=fact*exim*Lambda/(4.*pi*eps0)
      ey[...]+=fact*eyim*Lambda/(4.*pi*eps0)

  def beam_kickold(self):
    fact=1.#/(1.*pi)
    if not self.addkick:
      top.lresetparticlee = true
      return
    else:
      top.lresetparticlee = false
    # set linear charge density
    Lambda = self.Lambda[self.ikick]*self.bucket_train[self.ibk]
    pg = top.pgroup
    if sum(pg.nps)==0:return
    exim0 = zeros(self.nparpgrp,'d')
    eyim0 = zeros(self.nparpgrp,'d')
    x0=y0=0.
    for js in range(pg.ns):
      if pg.nps[js]==0:continue
      ng = 1+pg.nps[js]/self.nparpgrp
      for ig in range(ng):
        il = pg.ins[js]-1+self.nparpgrp*ig
        iu = min(il+self.nparpgrp,pg.ins[js]-1+pg.nps[js])
        np = iu-il
        x = pg.xp[il:iu]
        ex = pg.ex[il:iu]
        if self.l_switchyz:
          y = pg.zp[il:iu]
          ey = pg.ez[il:iu]
        else:
          y = pg.yp[il:iu]
          ey = pg.ey[il:iu]
        exim = exim0[:np]
        eyim = eyim0[:np]
#        if(iden_xy==0) zdir=zef_direct_point(x,y,x0,y0)
        if(pos.iden_xy==1): 
          exbe,eybe = Bassetti_Erskine(fact*Lambda,x,y,x0,y0,self.sigx,self.sigy)
          ex[...]=exbe
          ey[...]=eybe
#        if(iden_xy==2) zdir=zef_direct_unifell(x,y,x0,y0,abeam,bbeam)
#        if(iden_xy==3) zdir=zef_direct_parabell(x,y,x0,y0,abeam,bbeam)
        if(pos.iim==1):
#        if(iim==1.and.ichsh==2) zim=zef_image_rect(x,y,x0,y0,ach,bch)
#        if(iim==1.and.ichsh==3) zim=zef_image_open_H(x,y,x0,y0,bch)
          if pos.ichsh==1: 
            ef_image_ell(np,exim,eyim,x,y,x0,y0,pos.ach,pos.bch)
          ex[...]+=fact*exim*Lambda/(4.*pi*eps0)
          ey[...]+=fact*eyim*Lambda/(4.*pi*eps0)

  def randomcull(self):
 #  Routine to cull every other particle (i.e., "random culling")
    pg = top.pgroup
    phe = self.phelectrons
    sece = self.secelectrons
#    print "Number of primaries=",phe.getn(),"Secondaries=",sece.getn()
    if (phe.getn()+sece.getn())>self.electronnmax:
     self.ncull = self.ncull+1 
     self.wherecull.append(top.it)
     for sp in [phe,sece]:
      for js in sp.jslist:
       if pg.nps[js]>1:
        il=pg.ins[js]-1
        iu=il+pg.nps[js]
        initial_charge = sum(pg.pid[il:iu,top.wpid-1])
        iunew=il+pg.nps[js]/2
        pg.gaminv[il:iu:2]=0
        processlostpart(top.pgroup,js+1,top.clearlostpart,top.time,top.zbeam) 
        il=pg.ins[js]-1
        iu=il+pg.nps[js]
        final_charge = sum(pg.pid[il:iu,top.wpid-1])
        newwtfac = initial_charge/final_charge
        pg.pid[il:iu,top.wpid-1]*=newwtfac
 #       print "for js=",js,"init chg=",initial_charge,"final chg=",final_charge
 #       print "After cull,","Number of primaries=",phe.getn(),"Secondaries=",sece.getn()
 #       print "Culled Electrons", "ncull=",self.ncull
    q=0.
    for js in phe.jslist:
     if phe.getn()>0:    
      il=pg.ins[js]-1
      iu=il+pg.nps[js]
      weights=phe.getpid(id=top.wpid-1)
      if me==0:q+=sum(weights)
    self.qph.append(q) 
    q=0.
    for js in sece.jslist:
     if sece.getn()>0:    
      il=pg.ins[js]-1
      iu=il+pg.nps[js]
      weights=sece.getpid(id=top.wpid-1)
      if me==0:q+=sum(weights)
    self.qsec.append(q) 

  def pos_fieldsol(self):
#     print 'pos_fieldsol'
    pel = self.phelectrons
    sel = self.secelectrons
    n=pel.getn()+sel.getn()
    if n==0:return
    wpel=top.pgroup.sw[pel.jslist[0]]
    wsel=top.pgroup.sw[sel.jslist[0]]
    self.rays=fzeros([3,n],'d')
    self.rays[0,:]=concatenate((pel.getx(),sel.getx(),))
    self.rays[1,:]=concatenate((pel.gety(),sel.gety(),))
    self.rays[2,:]=concatenate((wpel*pel.getpid(id=top.wpid-1),wsel*sel.getpid(id=top.wpid-1),))
    mgdeposit2dln(pos.xl,pos.yb,n,pos.imax,pos.jmax,pos.rhs0,self.rays)
#solve for the electric potential
#initm =2 starting from the coarsest grid iteration.
#initm =1 starting from the finest grid iteration.
#If we start from the finest grid iteration, the previous
#solution of phi can be used as a first guess.
    initm=1
    mgfieldsolve(initm,pos.imax,pos.jmax,pos.mngrids,pos.phi0,pos.rhs0)
#    print 'pos_fieldsol done'

  def pos_scatter(self):
    if self.l_posmgsolver:
      phi=pos.phi0
    else:
      phi=transpose(frz.basegrid.phi[1:-1,1:-1])
    mgexey(pos.imax,pos.jmax,phi,pos.ex,pos.ey)
    if not self.l_posmgsolver:
      frz.basegrid.phi[...]=0.

  def pos_electronkick(self):
#    print 'pos_electronkick'
    pg=top.pgroup
    il = w3d.jmin
    iu = w3d.jmax
    np = iu-il
    if np==0:return
    exeypt = fzeros([2,w3d.jmax-w3d.jmin],'d')
#	relec=ech**2/(fourpieps0*emass)	!class el. radius [m]
#	scale=4*pi/(cellszx*cellszy)
    x = pg.xp[il:iu]
    ex = pg.ex[il:iu]
    weights = pg.pid[il:iu,top.wpid-1]
    if self.l_switchyz:
      y = pg.zp[il:iu]
      ey = pg.ez[il:iu]
    else:
      y = pg.yp[il:iu]
      ey = pg.ey[il:iu]
    rays=fzeros([3,np],'d')
    rays[0,:]=x
    rays[1,:]=y
    rays[2,:]=weights
#    spcoeff=pos.relec*clight**2*emass/(echarge*pos.slength)	#[m**2/s]
#    coeffkick=-pos.scale*spcoeff
    if self.l_posmgsolver:
      coeffkick = echarge/(eps0*pos.cellszx*pos.cellszy*pos.slength)
    else:
      coeffkick = 1.
    mgscatter2dln(pos.xl,pos.yb,np,pos.imax,pos.jmax,pos.ex,pos.ey,rays,exeypt)
    ex[...]+=coeffkick*exeypt[0,:]
    ey[...]+=coeffkick*exeypt[1,:]
#    print 'pos_electronkick done'
