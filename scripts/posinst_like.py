from pos import *
from Secondaries import *
from appendablearray import *

class Posinst_Like:
  """
Class for generating photo-electrons
 - posinst_file: name of Posinst input file 
 - l_verbose   : sets verbosity (default=0). 
  """
  def __init__(self,posinst_file=None,l_verbose=0,nx=64,ny=None,
               xmin=None,xmax=None,ymin=None,ymax=None,weight=None,
               conductors=None,l_secondaries=1,dtfact=1,l_switchyz=0, 
               electronnmax=None,l_posmgsolver=0,l_posscatter=0,
               lrefineallintercept=0,aura=0.,
               l_3d=0,nz=4,nsz=1.):
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
      if l_3d:self.nz=nz
      if pos.ispch==3:
        if nx is None:nx=pos.imax-1+2
        if ny is None:ny=pos.jmax-1+2
#        if nx is None:nx=pos.imax-1
#        if ny is None:ny=pos.jmax-1
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
    self.addkick=0
    if pos.chmphelnom>0.:
      if weight is None:weight=pos.chmphelnom#/pos.slength#*self.sigz
      if l_3d:weight*=nsz*self.sigz#/self.nz
      self.phelectrons=Species(type=Electron,weight=weight)
      self.PHEL=PhotoElectrons(l_verbose=0,l_switchyz=l_switchyz)
      self.PHEL.add(emitted_species=self.phelectrons)
    if pos.chmionelnom>0.:
      if weight is None:weight=pos.chmionelnom#/pos.slength#*self.sigz
      if l_3d:weight*=nsz*self.sigz#/self.nz
      self.ioelectrons=Species(type=Electron,weight=weight)
      self.IOEL=IonizElectrons(l_verbose=0,l_switchyz=l_switchyz)
      self.IOEL.add(emitted_species=self.ioelectrons)
    if l_secondaries:
      self.secelectrons=Species(type=Electron,weight=weight)
#    top.lresetparticlee = false
#    top.lresetparticleb = false
    self.nparpgrp = 4096
    self.l_3d=l_3d
    self.dtfact=dtfact
    self.l_switchyz=l_switchyz
    dx = 2.*pos.ach/(pos.imax-1)
    dy = 2.*pos.bch/(pos.jmax-1)
    if xmin is None:xmin=-pos.ach-dx
    if xmax is None:xmax= pos.ach+dx
    if ymin is None:ymin=-pos.bch-dy
    if ymax is None:ymax= pos.bch+dy
#    if xmin is None:xmin=-pos.ach*1.05
#    if xmax is None:xmax= pos.ach*1.05
#    if ymin is None:ymin=-pos.bch*1.05
#    if ymax is None:ymax= pos.bch*1.05
    w3d.nx=nx
    if ny is None:
      ny=nint(pos.bch*nx/pos.ach)
    w3d.xmmin=xmin
    w3d.xmmax=xmax
    if l_3d:
      w3d.ymmin=ymin
      w3d.ymmax=ymax
      w3d.ny=ny
      w3d.nz = nz
      w3d.zmmin=0.#*pos.slength
      w3d.zmmax=nsz*self.sigz
      w3d.boundnz=w3d.bound0=periodic
      top.pboundnz=top.pbound0=periodic
      self.MRroot=MultiGrid()
      registersolver(self.MRroot)
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
        w3d.zmmin=ymin
        w3d.zmmax=ymax
        w3d.nz=ny
      else:
        w3d.ymmin=ymin
        w3d.ymmax=ymax
        w3d.ny=ny
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
        if pos.ichsh==1: # elliptical
          print w3d.nz,w3d.nzlocal
          self.pipe = YCylinderEllipticOut(ellipticity = pos.bch/pos.ach,
                                                   radius      = pos.ach,
                                                   length      = (w3d.ymmax-w3d.ymmin)*10.,
                                                   ycent       = 0.5*(w3d.ymmin+w3d.ymmax),
                                                   condid      = 1)
          self.pipescraper = self.pipe
#          self.pipescraper = -Sphere(radius=pos.ach,condid=1)
        if pos.ichsh==2: # rectangular
          self.pipe = Box(xsize=2.*pos.ach,
                          zsize=3.*pos.bch,
                          ysize=pos.slength,
                          xcent=2.*pos.ach) \
                    + Box(xsize=2.*pos.ach,
                          zsize=3.*pos.bch,
                          ysize=pos.slength,
                          xcent=-2.*pos.ach) \
                    + Box(xsize=2.*pos.ach,
                          zsize=2.*pos.bch,
                          ysize=pos.slength,
                          zcent=2.*pos.bch) \
                    + Box(xsize=2.*pos.ach,
                          zsize=2.*pos.bch,
                          ysize=pos.slength,
                          zcent=-2.*pos.bch)
          self.pipescraper = self.pipe
#        self.scrapegrid=Grid(nx=nx,ny=w3d.nz,nzlocal=ny,nz=ny)
#        self.scrapegrid=Grid(nx=nx,ny=ny,nzlocal=w3d.nzlocal,nz=w3d.nz)
      else:
        if pos.ichsh==1: # elliptical
          self.pipe = ZCylinderEllipticOut(ellipticity = pos.bch/pos.ach,
                                                   radius      = pos.ach,
                                                   length      = (w3d.zmmax-w3d.zmmin)*10.,
                                                   zcent       = 0.5*(w3d.zmmin+w3d.zmmax),
                                                   condid      = 1)
        if pos.ichsh==2: # rectangular
          self.pipe = Box(xsize=2.*pos.ach,
                          ysize=3.*pos.bch,
                          zsize=w3d.zmmaxglobal-w3d.zmminglobal,
                          xcent=2.*pos.ach) \
                    + Box(xsize=2.*pos.ach,
                          ysize=3.*pos.bch,
                          zsize=w3d.zmmaxglobal-w3d.zmminglobal,
                          xcent=-2.*pos.ach) \
                    + Box(xsize=2.*pos.ach,
                          ysize=2.*pos.bch,
                          zsize=w3d.zmmaxglobal-w3d.zmminglobal,
                          ycent=2.*pos.bch) \
                    + Box(xsize=2.*pos.ach,
                          ysize=2.*pos.bch,
                          zsize=w3d.zmmaxglobal-w3d.zmminglobal,
                          ycent=-2.*pos.bch)
        self.scrapegrid=Grid(nx=nx,ny=ny,nzlocal=w3d.nzlocal,nz=w3d.nz)
        self.pipescraper = self.pipe
        import __main__
        __main__.sc=self.scrapegrid
      installconductors(self.pipe)
      self.scraper = ParticleScraper(self.pipescraper,
                                     lsaveintercept=1,
                                     lsavecondid=1,
                                     lcollectlpdata=1,
 #                                    grid=self.scrapegrid,
                                     lrefineintercept=0,
                                     lrefineallintercept=lrefineallintercept,
                                     nstepsperorbit=8,
                                     aura=aura)
#            self.scraper.l_print_timing=1
      self.conductors = [self.pipe]
    if l_secondaries:
      def set_params(maxsec,mat_num):
        pass
#        self.Sec.enpar=pos.enpar
#        self.Sec.pnpar=pos.pnpar
      self.Sec = Secondaries(min_age=None,set_params_user=set_params,vmode=1,l_set_params_user_only=true)
      self.set_params=set_params
      self.Sec.set_params_user=self.set_params
      self.Sec.enpar=pos.enpar
      self.Sec.pnpar=pos.pnpar
      try:
        self.Sec.add(incident_species = self.phelectrons,
                     emitted_species  = self.secelectrons,
                     conductor        = self.pipe,
                     interaction_type = 0)   
      except:
        pass
      try:
        self.Sec.add(incident_species = self.ioelectrons,
                     emitted_species  = self.secelectrons,
                     conductor        = self.pipe,
                     interaction_type = 0)   
      except:
        pass
      self.Sec.add(incident_species = self.secelectrons,
                   emitted_species  = self.secelectrons,
                   conductor        = self.pipe,
                   interaction_type = 0)   

  def push_buckets(self):
    for ibk in range(self.nbuckets):
      self.push_bucket()

  def push_bucket(self):
    print '    *** bucket %g out of %g'%(self.ibk+1,self.nbuckets)
    if self.l_3d:
      self.ldz = pos.beamvel*pos.dt[0]
      self.lnz = int((w3d.zmmax-w3d.zmmin)/self.ldz)+1
      self.lzmin = w3d.zmmin
      self.lzmax = self.lzmin + self.lnz*self.ldz
      self.Lambdaz = zeros(self.lnz+1,'d')
      thislambda = zeros(self.nz+1,'d')
      z=w3d.zmmin+arange(w3d.nz+1)*w3d.dz
    for k in range(self.nkicks+self.nsteps_gap):
      self.ikick=k
      if self.l_3d:self.Lambdaz[1:] = self.Lambdaz[:-1]
      if k<self.nkicks:
        self.addkick=1
        top.dt=pos.dt[k]
        newlambda = self.Lambda[k]*self.bucket_train[self.ibk]
      else:
        self.addkick=0
        top.dt=pos.deltat_g
        newlambda = 0.
      if self.l_3d:
        self.Lambdaz[0] = newlambda
        getgrid1d(w3d.nz+1,z,thislambda,self.lnz,self.Lambdaz,self.lzmin,self.lzmax)
      else:
        thislambda = newlambda
      try:
        if self.l_3d:
          self.PHEL.nz   = w3d.nz
          self.PHEL.dz   = w3d.dz
          self.PHEL.zmin = w3d.zmmin            
        self.PHEL.Lambda=thislambda
      except:
        pass
      try:
        self.IOEL.Lambda=thislambda
      except:
        pass
      top.dt/=self.dtfact
#      window(3);fma();pla(thislambda);refresh();window(0)
#      window(4);fma();ppzy(color=blue,msize=2);ppzy(js=1,msize=2,color=red);limits(w3d.zmmin,w3d.zmmax,w3d.ymmin,w3d.ymmax);refresh();window(0)
      for i in range(self.dtfact):
        step()
        if not self.l_3d:
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
    if self.l_3d:
      self.beam_kick_3d()
    else:
      self.beam_kick_2d()

  def beam_kick_2d(self):
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

  def beam_kick_3d(self):
#    if  self.ikick>=self.nkicks:return
    fact=1.#/(1.*pi)
    pg = top.pgroup
    exim = zeros(w3d.jmax-w3d.jmin,'d')
    eyim = zeros(w3d.jmax-w3d.jmin,'d')
    x0=y0=0.
    il = w3d.jmin
    iu = w3d.jmax
    np = iu-il
    x = pg.xp[il:iu]
    ex = pg.ex[il:iu]
    y = pg.yp[il:iu]
    ey = pg.ey[il:iu]
    z = pg.zp[il:iu]
    Lambda = zeros(np,'d')
    getgrid1d(np,z,Lambda,self.lnz,self.Lambdaz,self.lzmin,self.lzmax)
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
    sps=[]
    ne=0
    try:
      sps.append(self.phelectrons)
      ne+=sps[-1].getn()
    except:
      pass
    try:
      sps.append(self.ioelectrons)
      ne+=sps[-1].getn()
    except:
      pass
    try:
      sps.append(self.secelectrons)
      ne+=sps[-1].getn()
    except:
      pass
#    print "Number of primaries=",phe.getn(),"Secondaries=",sece.getn()
    if ne>self.electronnmax:
     self.ncull = self.ncull+1 
     self.wherecull.append(top.it)
     for sp in sps:
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
    return
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
    if self.l_switchyz:
      self.rays[1,:]=concatenate((pel.getz(),sel.getz(),))
    else:
      self.rays[1,:]=concatenate((pel.gety(),sel.gety(),))
    if top.wpid>0:
      self.rays[2,:]=concatenate((wpel*pel.getpid(id=top.wpid-1),wsel*sel.getpid(id=top.wpid-1),))
    else:
      self.rays[2,:]=wpel
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
    if top.wpid>0:
      weights = pg.pid[il:iu,top.wpid-1]
    else:
      weights = ones(np,'d')
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

class IonizElectrons:
  """
Class for generating photo-electrons
 - posinst_file: name of Posinst input file 
 - xfloor      : photo-electrons generated by Posinst that have x<xfloor will be forced to x=xfloor
 - xceiling    : photo-electrons generated by Posinst that have x>xceiling will be forced to x=xceiling
 - yfloor      : photo-electrons generated by Posinst that have y<yfloor will be forced to y=yfloor
 - yceiling    : photo-electrons generated by Posinst that have y>xceiling will be forced to y=yceiling
 - nz          : number of longitudinal slices (default=100)
 - l_xmirror   : turns mirroring of emitted photo-electrons with regard to x-axis on/off
 - l_verbose   : sets verbosity (default=0). 
  """
  def __init__(self,posinst_file=None,xfloor=None,xceiling=None,yfloor=None,yceiling=None,
               nz=100,l_xmirror=0,l_switchyz=0,l_verbose=0):
     self.xfloor=xfloor
     self.xceiling=xceiling
     self.yfloor=yfloor
     self.yceiling=yceiling
     self.nz=nz
     self.l_xmirror=l_xmirror
     self.l_switchyz=l_switchyz
     self.l_verbose=l_verbose
     self.inter={}
     self.npmax=4096
     self.nps={}
     self.x={}
     self.y={}
     self.z={}
     self.vx={}
     self.vy={}
     self.vz={}
     self.pid={}
     self.Lambda=0.
     if posinst_file is not None:init_posinst_for_warp(posinst_file)
     self.install()
     
  def add(self,incident_species=None,emitted_species=None):
    isinc=incident_species
    issec=[]
    if not self.inter.has_key(isinc):
        self.inter[isinc]={}
        for key in ['incident_species','emitted_species']:
          self.inter[isinc][key]=[]
        self.inter[isinc]['incident_species']=incident_species
    self.inter[isinc]['emitted_species'] = emitted_species
    js=emitted_species.jslist[0]
    if not self.x.has_key(js):
      self.nps[js]=0
      self.x[js]=fzeros(self.npmax,'d')
      self.y[js]=fzeros(self.npmax,'d')
      self.z[js]=fzeros(self.npmax,'d')
      self.vx[js]=fzeros(self.npmax,'d')
      self.vy[js]=fzeros(self.npmax,'d')
      self.vz[js]=fzeros(self.npmax,'d')
      if top.wpid>0:
        self.pid[js]=fzeros([self.npmax,top.npid],'d')

  def install(self):
    if not isinstalledbeforeloadrho(self.generate):
      installbeforeloadrho(self.generate)

  def addpart(self,nn,x,y,z,vx,vy,vz,js,weight=None):
    if self.nps[js]+nn>self.npmax:self.flushpart(js)
    il=self.nps[js]
    iu=il+nn
    self.x[js][il:iu]=x
    self.y[js][il:iu]=y
    self.z[js][il:iu]=z
    self.vx[js][il:iu]=vx
    self.vy[js][il:iu]=vy
    self.vz[js][il:iu]=vz
    if weight is not None:self.pid[js][il:iu,top.wpid-1]=weight
    self.nps[js]+=nn
      
  def flushpart(self,js):
    if self.nps[js]>0:
       nn=self.nps[js]
       if top.wpid==0:
         addparticles(x=self.x[js][:nn],
                      y=self.y[js][:nn],
                      z=self.z[js][:nn],
                      vx=self.vx[js][:nn],
                      vy=self.vy[js][:nn],
                      vz=self.vz[js][:nn],
                      js=js,
                      lmomentum=true)
       else: 
         addparticles(x=self.x[js][:nn],
                      y=self.y[js][:nn],
                      z=self.z[js][:nn],
                      vx=self.vx[js][:nn],
                      vy=self.vy[js][:nn],
                      vz=self.vz[js][:nn],
                      pid=self.pid[js][:nn,:],
                      js=js,
                      lmomentum=true)
       self.nps[js]=0
         
  def generate(self):
    for ints in self.inter.keys():
     incident_species=self.inter[ints]['incident_species']
     emitted_species=self.inter[incident_species]['emitted_species']
     if incident_species is None:
       self.nz=1
       if self.l_switchyz:
         ymin=w3d.ymmin
         ymax=w3d.ymmax
         dy=(w3d.ymmax-w3d.ymmin)
       else:
         zmin=w3d.zmmin
         zmax=w3d.zmmax
         dz=(w3d.zmmax-w3d.zmmin)
     else:
       if self.l_switchyz:
         ymin=min(incident_species.gety())
         ymax=max(incident_species.gety())
         dy=(ymax-ymin)/self.nz
         self.Lambda = sum(sum(incident_species.get_density(nx=2, 
                                                            nz=2, 
                                                            ny=self.nz,
                                                            ymin=ymin,
                                                            ymax=ymax,
                                                            l_minmax_grid=false,
                                                            l_dividebyvolume=false,
                                                            charge=1),2),0)
       else:
         zmin=min(incident_species.getz())
         zmax=max(incident_species.getz())
         dz=(zmax-zmin)/self.nz
         self.Lambda = sum(sum(incident_species.get_density(nx=2, 
                                                            ny=2, 
                                                            nz=self.nz,
                                                            zmin=zmin,
                                                            zmax=zmax,
                                                            l_minmax_grid=false,
                                                            l_dividebyvolume=false,
                                                            charge=1),0),0)
     weightemit=top.pgroup.sw[emitted_species.jslist[0]]*abs(top.pgroup.sq[emitted_species.jslist[0]])
     torr_to_MKS=133.3224	#1 Torr=133.3224  N/m**2
     #ideal gas law: rho=p/(k*T)
     gas_prefactor=(torr_to_MKS/pos.boltzk)/294
     gasden=(gas_prefactor)*pos.pressure*(294/pos.temperature)
     rnipbppm=(1e-28)*pos.crossect*gasden
     for i in range(self.nz):
       if incident_species is None:
         rhel = self.Lambda*rnipbppm*clight*top.dt/weightemit
       else:
         rhel = self.Lambda[i]*rnipbppm*clight*top.dt/weightemit
       # rhel is the number of electrons created at each timestep
       n=int(rhel)
       if ranf()<rhel-n:n+=1  # randomly add one electrons based on rhel fractional part
       if self.l_verbose:print ' *** i,rhel,nemit= ',i,rhel,n
       if n==0:continue
       pos.nionel[0]=n   # tells Posinst to emit n photoelectrons
       gen_ionizelectrons(1) # number of beam slice in POSINST =1. Use only 1.
       if self.l_verbose:print 'nlast',pos.nlast,"nphel=",pos.nionel[0]

       if self.l_xmirror:
         # put photons on both sides of the vacuum chamber
         xran = ranf(pos.x[:pos.nlast])
         xran = where(xran>0.5,1.,-1.)
         pos.x[:pos.nlast] = pos.x[:pos.nlast]*xran
         pos.vgx[:pos.nlast] = pos.vgx[:pos.nlast]*xran

       if self.l_verbose:print "min and max of photoelectrons=",min((pos.z[:pos.nlast]/pos.slength-0.5)*dz+i*dz),\
                                                                max((pos.z[:pos.nlast]/pos.slength-0.5)*dz+i*dz)
       ns = pos.nlast
       js_new=emitted_species.jslist[0]
       usq = (pos.vgx[:pos.nlast]**2 + pos.vgy[:pos.nlast]**2 + pos.vgz[:pos.nlast]**2)/clight**2
       gaminv = 1./sqrt(1. + usq)
       dt=ranf(usq)*top.dt
       if self.xfloor is not None:
         pos.x[:pos.nlast]=where(pos.x[:pos.nlast]>self.xfloor,pos.x[:pos.nlast],self.xfloor)
       if self.xceiling is not None:
         pos.x[:pos.nlast]=where(pos.x[:pos.nlast]<self.xceiling,pos.x[:pos.nlast],self.xceiling)
       if self.yfloor is not None:
         pos.y[:pos.nlast]=where(pos.y[:pos.nlast]>self.yfloor,pos.y[:pos.nlast],self.yfloor)
       if self.yceiling is not None:
         pos.y[:pos.nlast]=where(pos.y[:pos.nlast]<self.yceiling,pos.y[:pos.nlast],self.yceiling)
       if top.wpid==0:
         weights = None
       else:
         weights = ones(pos.nlast,'d')
       if self.l_switchyz:
           self.addpart(ns,pos.x[:pos.nlast]+dt*pos.vgx[:pos.nlast]*gaminv,
                         pos.z[:pos.nlast]*0.,
                         pos.y[:pos.nlast]+dt*pos.vgy[:pos.nlast]*gaminv,
                         pos.vgx[:pos.nlast],
                         pos.vgz[:pos.nlast],
                         pos.vgy[:pos.nlast],
                         js_new,
                         weights)
       else:
           self.addpart(ns,pos.x[:pos.nlast]+dt*pos.vgx[:pos.nlast]*gaminv,
                         pos.y[:pos.nlast]+dt*pos.vgy[:pos.nlast]*gaminv,
                         (pos.z[:pos.nlast]/pos.slength)*dz+i*dz+zmin,
                         pos.vgx[:pos.nlast],
                         pos.vgy[:pos.nlast],
                         pos.vgz[:pos.nlast],
                         js_new,
                         weights)
       pos.nlast=0

    # --- make sure that all particles are added
    for js in self.x.keys():
      self.flushpart(js)
