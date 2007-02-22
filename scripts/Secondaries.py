"""
Secondaries: class for generating secondaries
"""
from warp import *
from pos import *
from species import *
try:
  import txphysics
  l_txphysics = 1
except:
  print 'WARNING: module txphysics is not accessible.'
  l_txphysics = 0
try:
  import desorb
  l_desorb = 1
except:
  print 'WARNING: module desorb is not accessible.'
  l_desorb = 0
import timing as t
import time

secondaries_version = "$Id: Secondaries.py,v 1.15 2007/02/22 23:55:17 dave Exp $"
def secondariesdoc():
  import Secondaries
  print Secondaries.__doc__

class Secondaries:
  """
Class for generating secondaries
 - isinc:      list of incident species
 - conductors: list of list of conductors 
 - issec:      list of list of secondary species
 - type:       list of list of types of interaction
               e- -> Cu              = 1 
               e- -> Stainless Steel = 2 
               H  -> Au              = 3 
               He -> Au              = 4 
               K  -> SS              = 5 
 - min_age: this sets a minimum age for a macroparticle in order to 
            qualify for emitting secondaries (must be old enough to procreate :-) ).
            This is in units of time steps. The default is None (no minimum age required).
 - set_params_user: function that sets SEY parameters provided by the user.
                    Default is None: default set_params routine is used.
 - vmode: 1 (default) -> uses instantaneous velocity to compute angle to normal of conductor
          2           -> uses difference between current and old positions to compute angle to normal of conductor
 - l_verbose: sets verbosity (default=0). 
  """
  def __init__(self,isinc=None,conductors=None,issec=None,set_params_user=None,material=None,
                    xoldpid=None,yoldpid=None,zoldpid=None,min_age=None,vmode=1,l_verbose=0):
    self.inter={}
    self.outparts=[]
#    self.isinc = isinc
#    self.conductors = conductors
#    self.issec = issec
#    self.type = type
    self.vmode=vmode
    self.l_verbose=l_verbose
    self.l_record_timing=0
#    self.condids={}
#    self.emitted={}
    self.set_params_user=set_params_user
    self.call_set_params_user(pos.maxsec,pos.mat_number)
    self.min_age=min_age
    if self.min_age is not None:
      w3d.l_inj_rec_inittime=true
      if top.tpid==0:
        top.tpid=nextpid()
        gchange('Particles')
    self.install()
    if xoldpid is None:
      self.xoldpid=top.npid-3
    else:
      self.xoldpid=xoldpid
    if yoldpid is None:
      self.yoldpid=top.npid-2
    else:
      self.yoldpid=yoldpid
    if zoldpid is None:
      self.zoldpid=top.npid-1
    else:
      self.zoldpid=zoldpid
    # set variables for secondary electrons routines
    self.secelec_ns = zeros(1)
    self.secelec_un = zeros(pos.maxsec,'d')
    self.secelec_ut = zeros(pos.maxsec,'d')
    self.secelec_uz = zeros(pos.maxsec,'d')
    self.secelec_ityps  = zeros(pos.maxsec)
    self.secelec_ekstot = zeros(pos.maxsec,'d')
    self.secelec_dele   = zeros(1,'d')
    self.secelec_delr   = zeros(1,'d')
    self.secelec_delts  = zeros(1,'d')
    # set work arrays for txphysics
    # --- This is done in a separate method since it will be called from
    # --- multiple places, here and in setstate.
    self.creattxphysicsarrays()
    # set arrays for emitted particles
    self.npmax=4096
    self.nps={}
    self.x={}
    self.y={}
    self.z={}
    self.vx={}
    self.vy={}
    self.vz={}
    self.pid={}
    if pos.nsteps==0:
      self.piditype=0
    else:
      self.piditype=nextpid()-1
    
    if isinc is None:return
    for iis,js in enumerate(isinc):
      for ics,cond in enumerate(conductors[iis]):
        if material is None: m = None
        else:                m = material[iis][ics]
        self.add(js,cond,issec[iis][ics],m)
    
  def creattxphysicsarrays(self):
    if l_txphysics:
      self.emitted_e=txphysics.doubleArray(1000)
      self.emitted_bn=txphysics.doubleArray(1000)
      self.emitted_bt=txphysics.doubleArray(1000)
      self.emitted_bz=txphysics.doubleArray(1000)

  def __getstate__(self):
    dict = self.__dict__.copy()

    # --- Functions cannot be pickled, so set_params_user has
    # --- to be specially handled. If it is defined, then only save
    # --- the name of the function so that it can be restored later.
    # --- This is dealt with in call_set_params_user so nothing needs
    # --- to happen in setstate. It is done there since this instance
    # --- may be restore before the function.
    if self.set_params_user is not None:
      try:
        dict['set_params_user'] = self.set_params_user.__name__
      except AttrbiuteError:
        print ("Warning: Secondaries set_params_user function '%s' is not "+
               "a proper function and will not be saved")%self.set_params_user
        dict['set_params_user'] = None

    # --- These arrays are swig arrays and cannot be pickled.
    # --- They are only temporary arrays and don't need to be
    # --- restore so they can just be deleted. They are recreated
    # --- in the setstate below.
    try:
      del dict['emitted_e']
      del dict['emitted_bn']
      del dict['emitted_bt']
      del dict['emitted_bz']
    except KeyError:
      pass

    return dict

  def __setstate__(self,dict):
    'Restore the dictionary and recreate the txphysics arrays'
    self.__dict__.update(dict)
    self.creattxphysicsarrays()

  def add(self,incident_species=None,conductor=None,emitted_species=None,material=None,interaction_type=None,scale_factor=None):
    if material is None: material = conductor.material
    if interaction_type is None:
      if material=='Cu' and incident_species.type is Electron: interaction_type=1
      if material=='SS' and incident_species.type is Electron: interaction_type=2
      if material=='Au' and incident_species.type is Hydrogen: interaction_type=3
      if material=='Au' and incident_species.type is Helium:   interaction_type=4
      if material=='SS' and incident_species.type is Potassium:interaction_type=5
    if interaction_type is None:raise('Error in Secondaries: invalid material or incident species')
    isinc=incident_species
    if type(emitted_species)<>type([]):emitted_species=[emitted_species]
    issec=[]
    for e in emitted_species:
      issec.append(e.jslist[0])
    if not self.inter.has_key(isinc):
        self.inter[isinc]={}
        for key in ['emitted','condids','issec','conductors','type','isneut','incident_species', \
                    'emitted_species','material','scale_factor']:
          self.inter[isinc][key]=[]
        self.inter[isinc]['incident_species']=incident_species
    self.inter[isinc]['condids']         += [conductor.condid]
    self.inter[isinc]['conductors']      += [conductor]
    self.inter[isinc]['emitted']         += [[]]
    for i in range(len(emitted_species)):
      self.inter[isinc]['emitted'][-1]   += [0.]
    self.inter[isinc]['issec']           += [issec]
    self.inter[isinc]['emitted_species'] += [emitted_species]
    self.inter[isinc]['material']        += [material]
    self.inter[isinc]['type']            += [interaction_type]
    if scale_factor is not None and scale_factor>1.:
      raise('Error in secondaries: scale_factor must be less than 1.')
    self.inter[isinc]['scale_factor']    += [scale_factor]
    for e in emitted_species:
      js=e.jslist[0]
      if not self.x.has_key(js):
        self.nps[js]=0
        self.x[js]=fzeros(self.npmax,'d')
        self.y[js]=fzeros(self.npmax,'d')
        self.z[js]=fzeros(self.npmax,'d')
        self.vx[js]=fzeros(self.npmax,'d')
        self.vy[js]=fzeros(self.npmax,'d')
        self.vz[js]=fzeros(self.npmax,'d')
        if top.wpid>0 or self.piditype>0:
          self.pid[js]=fzeros([self.npmax,top.npid],'d')

  def install(self):
    if not isinstalledbeforeloadrho(self.generate):
      installbeforeloadrho(self.generate)

  def addpart(self,nn,x,y,z,vx,vy,vz,js,weight=None,itype=None):
    if self.nps[js]+nn>self.npmax:self.flushpart(js)
    il=self.nps[js]
    iu=il+nn
    xx=fones(nn,'d')
    self.x[js][il:iu]=x*xx
    self.y[js][il:iu]=y*xx
    self.z[js][il:iu]=z*xx
    self.vx[js][il:iu]=vx
    self.vy[js][il:iu]=vy
    self.vz[js][il:iu]=vz
    if weight is not None:self.pid[js][il:iu,top.wpid-1]=weight
    if itype is not None:self.pid[js][il:iu,self.piditype]=itype.astype(Float)
    self.nps[js]+=nn
      
  def flushpart(self,js):
    if self.nps[js]>0:
       nn=self.nps[js]
       if top.wpid==0 and self.piditype==0:
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
         
  def printall(self,l_cgm=0):
    title='        *** Particle/wall interactions ***\n'
    textblock=''
    swidth=0
    cwidth=0
    ewidth={}
    for js in self.inter.keys():
      swidth=max(swidth,len(self.inter[js]['incident_species'].name))
      for ics,cond in enumerate(self.inter[js]['conductors']):
        cwidth=max(cwidth,len(cond.name))
        for ie,emitted_species in enumerate(self.inter[js]['emitted_species'][ics]):
          if not ewidth.has_key(ie):
            ewidth[ie]=len(emitted_species.name)
          else:
            ewidth[ie]=max(ewidth[ie],len(emitted_species.name))
    fs='%%-%gs'%swidth
    fc='%%-%gs'%cwidth
    fe={}
    for ie in ewidth.keys():
      fe[ie]='%%-%gs'%ewidth[ie]
    for js in self.inter.keys():
      sname=fs%self.inter[js]['incident_species'].name
      textblock+='\n'
      for ics,cond in enumerate(self.inter[js]['conductors']):
        cname=fc%cond.name
        for ie,emitted_species in enumerate(self.inter[js]['emitted_species'][ics]):
          if ie==0:
            ename=fe[ie]%emitted_species.name
          else:
            ename+=' + '+fe[ie]%emitted_species.name
        textblock += sname+' + '+cname+' => '+ename+'\n'
    if l_cgm:
      plt(title,0.22,0.905,justify="LT",height=14)
      plt(textblock,0.14,0.88,justify="LT",height=9,font='courier')
      fma()  
    else:
      print title
      print textblock

  def generate(self):
    # theta and phi are angles from normal to surface with regard to z and x axis respectively
    # psi is angle to rotate warp local frame to Posinst local frame around normal
    # eta is angle between incident velocity vector and normal to surface (called theta in Posinst)

    if self.l_record_timing:t1 = time.clock()

    # reset 'emitted' list to zero
    if self.l_verbose>1:print 'call secondaries'
    for js in self.inter.keys():
     for i in range(len(self.inter[js]['emitted'])):
      for j in range(len(self.inter[js]['emitted'][i])):  
       self.inter[js]['emitted'][i][j] = 0.

    # set computing box mins and maxs
    if w3d.l4symtry:
      xmin=-w3d.xmmax
    else:
      xmin=w3d.xmmin
    xmax=w3d.xmmax
    if w3d.l2symtry or w3d.l4symtry:
      ymin=-w3d.ymmax
    else:
      ymin=w3d.ymmin
    ymax=w3d.ymmax
    zmin=w3d.zmmin
    zmax=w3d.zmmax

    if self.l_record_timing:t2 = time.clock()
    tinit=tgen=tprepadd=tadd=0.
    # compute number of secondaries and create them
    for ints in self.inter.keys():
     incident_species=self.inter[ints]['incident_species']
     for js in incident_species.jslist:
      if self.l_verbose:print 'js',js
      if top.npslost[js]==0:continue
#      if top.npslost[js]==0 or top.it%top.pgroup.ndts[js]<>0:continue
      stride=top.pgroup.ndts[js]
      i1 = top.inslost[js] - 1 
      i2 = top.inslost[js] + top.npslost[js] - 1
      for ics,cond in enumerate(self.inter[incident_species]['conductors']):
        icond = cond.condid
        if self.l_verbose:print 'ics',ics
        iit = compress(top.pidlost[i1+top.it%stride:i2:stride,-1]==icond,arange(top.it%stride,top.npslost[js],stride))
        n = len(iit)
        if n==0:continue
        xplost = take(top.xplost[i1:i2],iit)
        yplost = take(top.yplost[i1:i2],iit)
        zplost = take(top.zplost[i1:i2],iit)
        # exclude particles out of computational box 
#        if w3d.solvergeom==w3d.RZgeom:
#          condition = (sqrt(xplost**2+yplost**2)<xmax) & \
#                      (zplost>zmin) & (zplost<zmax)
#        else:
#          condition = (xplost>xmin) & (xplost<xmax) & \
#                      (yplost>ymin) & (yplost<ymax) & \
#                      (zplost>zmin) & (zplost<zmax)
        # exclude particles recently created
        if self.min_age is not None:
          inittime = take(top.pidlost[i1:i2,top.tpid-1],iit)
          condition =  ((top.time-inittime)>self.min_age*top.dt)      
#          condition = condition & ((top.time-inittime)>self.min_age*top.dt)      
          iit2 = compress(condition,arange(n))
        else:
          iit2=arange(n)
        n = len(iit2)
        if self.l_verbose:print 'nlost=',n
        if n==0:continue
        xplost = take(xplost,iit2)
        yplost = take(yplost,iit2)
        zplost = take(zplost,iit2)
        iit    = take(iit,iit2)
        if self.l_record_timing:t.start()    
        uxplost = take(top.uxplost[i1:i2],iit)
        uyplost = take(top.uyplost[i1:i2],iit)
        uzplost = take(top.uzplost[i1:i2],iit)
        gaminvlost = take(top.gaminvlost[i1:i2],iit)
        if self.vmode==1:
          vxplost=uxplost*gaminvlost
          vyplost=uyplost*gaminvlost
          vzplost=uzplost*gaminvlost
        elif self.vmode==2:
          xplostold = take(top.pidlost[i1:i2,self.xoldpid],iit)
          yplostold = take(top.pidlost[i1:i2,self.yoldpid],iit)
          zplostold = take(top.pidlost[i1:i2,self.zoldpid],iit)
          vxplost = (xplost-xplostold)/top.dt
          vyplost = (yplost-yplostold)/top.dt
          vzplost = (zplost-zplostold)/top.dt
        else:
          raise('Error in Secondaries, one should have lmode=1 or 2, but have lmode=%g'%self.lmode)
        # set energy of incident particle in eV
        e0 = where(gaminvlost==1., \
                   0.5*top.pgroup.sm[js]*(uxplost**2+uyplost**2+uzplost**2)/top.echarge,
                   (1./gaminvlost-1.)*top.pgroup.sm[js]*clight**2/top.echarge)
        if self.l_verbose:
          print 'xplost',xplost
          print 'yplost',yplost
          print 'zplost',zplost
          print 'e0',e0,gaminvlost,uxplost,uyplost,uzplost
        v = array([vxplost,vyplost,vzplost])
#        u = array([uxplost,uyplost,uzplost])
        theta = take(top.pidlost[i1:i2,-3],iit)
        phi   = take(top.pidlost[i1:i2,-2],iit)
        if top.wpid>0: 
          weight = take(top.pidlost[i1:i2,top.wpid-1],iit)
        else:
          weight=1.
        costheta = cos(theta)
        sintheta = sin(theta)
        cosphi   = cos(phi)
        sinphi   = sin(phi)
        # theta is relative to the z axis, phi to the x axis in the x-y plane 
        n_unit0 = array([sintheta*cosphi,sintheta*sinphi,costheta])
        coseta = -sum(v*n_unit0)/sqrt(sum(v*v))
#        print 'coseta = ',coseta
#        return
#        coseta=0.5*ones(shape(e0)[0])
#        e0=30.*ones(shape(e0)[0])
        if self.piditype>0:
          pos.ek0av   += sum(weight*e0)	#cummulative collision kinetic energy [eV] this step
          pos.ek0max   = max(max(e0),pos.ek0max)	#maximum collision kinetic energy [eV]
          pos.costhav += sum(weight*abs(coseta))
          pos.rcoll   += sum(weight)
#        coseta = 0.*cosphi
        if 1:#cond.lcollectlpdata:
          if not cond.lostparticles_angles.has_key(js):
            cond.lostparticles_angles[js]=zeros(181,'d')
          if not cond.lostparticles_energies.has_key(js):
            cond.lostparticles_energies[js]=zeros(1001,'d')
          e0min = min(e0)
          e0max = max(e0)
          if not cond.lostparticles_minenergy.has_key(js):
            cond.lostparticles_minenergy[js]=e0min
          if not cond.lostparticles_maxenergy.has_key(js):
            cond.lostparticles_maxenergy[js]=e0max
          l_rescale_energy_array=0
          e0minnew=e0min
          e0maxnew=e0max
          if e0min<cond.lostparticles_minenergy[js]:
            e0minnew=0.9*e0min
            l_rescale_energy_array=1
          if e0max>cond.lostparticles_maxenergy[js]:
            e0maxnew=1.1*e0max
            l_rescale_energy_array=1
          if l_rescale_energy_array:
            newlostparticles_energies=zeros(1001,'d')
            tmpcount=zeros(1001,'d')
            e0minold = cond.lostparticles_minenergy[js]
            e0maxold = cond.lostparticles_minenergy[js]
            e0old = e0minold+arange(1001)*(e0maxold-e0minold)/1000
            deposgrid1d(1,1001,e0old,cond.lostparticles_energies[js],1000,newlostparticles_energies,tmpcount,e0minnew,e0maxnew)
            cond.lostparticles_minenergy[js]=e0minnew
            cond.lostparticles_maxenergy[js]=e0maxnew
            cond.lostparticles_energies[js]=newlostparticles_energies
          if top.wpid >0:
            pass
#            deposgrid1d(1,i2-i1,tl, ql,nt,qt,qtmp,tmin,tmax)
          else:
            setgrid1d(shape(coseta)[0],arccos(coseta),180,cond.lostparticles_angles[js],0.,pi)
            setgrid1d(shape(e0)[0],e0,1000,cond.lostparticles_energies[js],cond.lostparticles_minenergy[js],cond.lostparticles_maxenergy[js])
        for i in range(n):  
#          print 'v',v[0][i],v[1][i],v[2][i],i,iit[i],js
#          print 'x',[xplost[i],yplost[i],zplost[i]]
#          print 'xold',[xplostold[i],yplostold[i],zplostold[i]]
#          print 'u',[uxplost[i],uyplost[i],uzplost[i]]
#          print 'e0',e0[i]
#          coseta[i] = -sum(v[0][i]*n_unit0[0][i]+v[1][i]*n_unit0[1][i]+v[2][i]*n_unit0[2][i])/sqrt(sum(v[0][i]*v[0][i]+v[1][i]*v[1][i]+v[2][i]*v[2][i]))
#          print 'coseta',coseta[i]
          l_warning=0
          l_infinity=0
          if coseta[i]<0.:
            l_warning=1
            swarn = 'WARNING issued by Secondaries.generate: coseta<0.'
            coseta[i]=-coseta[i]
            n_unit0[0][i]=-n_unit0[0][i]
            n_unit0[1][i]=-n_unit0[1][i]
            n_unit0[2][i]=-n_unit0[2][i]
            costheta[i] = cos(pi+theta[i])
            sintheta[i] = sin(pi+theta[i])
#            print 'coseta 1,2 :',coseta[i],-sum(u[i]*n_unit0[i])/sqrt(sum(u[i]*u[i]))
#            print 'n 1, 2',n_unit0[0][i],n_unit0[1][i],n_unit0[2][i],sintheta[i]*cosphi[i],sintheta[i]*sinphi[i],costheta[i]
          if xplost[i]==largepos or yplost[i]==largepos or zplost[i]==largepos:
            l_warning=1
            l_infinity=1
            swarn = 'WARNING issued by Secondaries.generate: particle at infinity'
          if l_warning and self.l_verbose:
            print swarn
#          print 'phi, theta',phi[i],theta[i]
#          print 'n',n_unit0[0][i],n_unit0[1][i],n_unit0[2][i]
#          print 'u',u[0][i],u[1][i],u[2][i]
          if l_infinity:
            continue
          if self.l_record_timing:
            t.finish()
            tinit+=t.micro()
            t.start()
          if self.l_verbose:print 'e0, coseta',e0[i],coseta[i]
          for ie,emitted_species in enumerate(self.inter[incident_species]['emitted_species'][ics]):
           js_new=emitted_species.jslist[0]
           if emitted_species.type is Electron:
            if incident_species.type is Electron:
#             try: # need try since will generate error if no secondaries are created
               if top.wpid>0:
                 self.generate_secondaries(e0[i],coseta[i],weight[i],self.inter[incident_species]['type'][ics])
               else:
                 self.generate_secondaries(e0[i],coseta[i],weight,self.inter[incident_species]['type'][ics])
               ns=self.secelec_ns[0]
               ut=self.secelec_ut[:ns]
               un=self.secelec_un[:ns]
               uz=self.secelec_uz[:ns]
               if self.piditype==0:
                 itype=None
               else:
                 itype=self.secelec_ityps[:ns]
               if self.l_verbose:
                 print 'nb secondaries = ',self.secelec_ns,' from conductor ',icond, \
                 e0[i], coseta[i],i1,i2,iit[i],top.npslost,self.secelec_dele,self.secelec_delr,self.secelec_delts            
#             except:
#               ns=0
            else: # incidents are atoms or ions
             if incident_species.type.__class__ is Atom:
#              try:

               if self.inter[incident_species]['material'][ics]=='SS':target_num=10025
#              -- version 0.2.1
               ns=txphysics.ion_ind_elecs(e0[i]/(1.e6*incident_species.type.A),
                                          max(0.04,coseta[i]),
                                          float(incident_species.type.Z),
                                          incident_species.type.A,
                                          target_num,
                                          self.emitted_e,
                                          self.emitted_bn,
                                          self.emitted_bt,
                                          self.emitted_bz)
#              -- version 0.6.1
#               ion_ind_e0 = txphysics.doubleArray(1)
#               ion_ind_ct = txphysics.doubleArray(1)
#               ion_ind_e0[0] = e0[i]/(1.e6*incident_species.type.A)
#               ion_ind_ct[0] = max(0.04,coseta[i])
#               ns,self.emitted_e,self.emitted_bn,self.emitted_bt,self.emitted_bz = \
#               txphysics.ion_ind_elecs(1, # size of input array
#                                       ion_ind_e0,
#                                       ion_ind_ct,
#                                       float(incident_species.type.Z),
#                                       incident_species.type.A,
#                                       target_num,
#                                       0.)
               ns = min(ns,self.npmax)
               if self.inter[incident_species]['scale_factor'][ics] is not None:
                 # scale by user provided factor (default is not to scale)
                 fns=float(ns)*self.inter[incident_species]['scale_factor'][ics]
                 ns=int(fns)
                 if ranf()<(fns-ns):ns+=1
               self.coseta=coseta[i]
               itype=None
               nsemit=ns+0
               # for some reason, the energy is sometimes negative on OSX
               for iemit in range(ns):
                 if self.emitted_e[iemit]<0.:
                   print 'Warning: negative energy from txphysics.ion_ind_elecs:',self.emitted_e[iemit]
                   nsemit-=1
               un=zeros(nsemit,'d')
               ut=zeros(nsemit,'d')
               uz=zeros(nsemit,'d')
               iwemit=0
               for iemit in range(ns):
                 if self.emitted_e[iemit]>=0.:
                   gammac = clight/sqrt(1.-self.emitted_bn[iemit]**2 \
                                          +self.emitted_bt[iemit]**2 \
                                          +self.emitted_bz[iemit]**2)
                   un[iwemit]=gammac*self.emitted_bn[iemit] 
                   ut[iwemit]=gammac*self.emitted_bt[iemit] 
                   uz[iwemit]=gammac*self.emitted_bz[iemit] 
                   iwemit+=1
               ns=nsemit
               if self.l_verbose:
                 print 'nb secondaries = ',ns,' from conductor ',icond, e0[i], coseta[i],i1,i2,iit[i],top.npslost          
#              except:
#               ns=0
            if self.l_record_timing:
              t.finish()
              tgen+=t.micro()
              t.start()
            if ns>0:
             self.inter[incident_species]['emitted'][ics][ie] += ns*top.pgroup.sq[js_new]*top.pgroup.sw[js_new]
             if costheta[i]<1.-1.e-10:
              z_unit0 = array([-sinphi[i],cosphi[i],0.])
#              z       = -array([uzplost[i]*n_unit0[1][i]-uyplost[i]*n_unit0[2][i],
#                                uxplost[i]*n_unit0[2][i]-uzplost[i]*n_unit0[0][i],
#                                uyplost[i]*n_unit0[0][i]-uxplost[i]*n_unit0[1][i]])
              z       = -array([vzplost[i]*n_unit0[1][i]-vyplost[i]*n_unit0[2][i],
                                vxplost[i]*n_unit0[2][i]-vzplost[i]*n_unit0[0][i],
                                vyplost[i]*n_unit0[0][i]-vxplost[i]*n_unit0[1][i]])
              z_unit  = z/sqrt(sum(z*z))
              cospsi  = sum(z_unit*z_unit0)
              sinpsi  = sqrt(max(0.,1.-cospsi*cospsi))
#              bt0 = (cospsi*bt - sinpsi*bz)
#              bz0 = (sinpsi*bt + cospsi*bz)
              ut0 = -(cospsi*ut - sinpsi*uz)
              uz0 = -(sinpsi*ut + cospsi*uz)
              un0 = un
             else:
              ut0 = ut
              uz0 = uz
              un0 = un
             uxsec = cosphi[i]*costheta[i]*ut0 - sinphi[i]*uz0 + cosphi[i]*sintheta[i]*un0
             uysec = sinphi[i]*costheta[i]*ut0 + cosphi[i]*uz0 + sinphi[i]*sintheta[i]*un0
             uzsec =          -sintheta[i]*ut0                 +           costheta[i]*un0
             if self.l_record_timing:
               t.finish()
               tprepadd+=t.micro()
               t.start()
             xnew = xplost[i]+n_unit0[0][i]*1.e-10*w3d.dx
             ynew = yplost[i]+n_unit0[1][i]*1.e-10*w3d.dy 
             znew = zplost[i]+n_unit0[2][i]*1.e-10*w3d.dz
#            pid[:,self.xoldpid]=xnew-vx*top.dt
#            pid[:,self.yoldpid]=ynew-vy*top.dt
#            pid[:,self.zoldpid]=znew-vz*top.dt
             if w3d.solvergeom==w3d.RZgeom:
               condition = (sqrt(xnew**2+ynew**2)>xmax) or \
                           (znew<zmin) or (znew>zmax)
             else:
               condition = (xnew<xmin) or (xnew>xmax) or \
                           (ynew<ymin) or (ynew>ymax) or \
                           (znew<zmin) or (znew>zmax)
             if condition:
              print 'WARNING: new particle outside boundaries',xnew,ynew,znew
#              self.outparts+=[[xnew,ynew,znew,xplost[i],yplost[i],zplost[i], \
#              xplostold[i],yplostold[i],zplostold[i],n_unit0[0][i],n_unit0[1][i],n_unit0[2][i],icond]]
              self.outparts+=[[xnew,ynew,znew,xplost[i],yplost[i],zplost[i], \
              n_unit0[0][i],n_unit0[1][i],n_unit0[2][i],icond]]
             else:
               if top.wpid==0:
                self.addpart(ns,xnew,ynew,znew,uxsec,uysec,uzsec,js_new,itype=itype)
               else:
                self.addpart(ns,xnew,ynew,znew,uxsec,uysec,uzsec,js_new,ones(ns)*weight[i],itype)
              
           # emit neutrals
           if emitted_species.type.__class__ is not Particle and emitted_species.charge_state==0: 
            my_yield=1.+1.82e-4*exp(0.09*180./pi*arccos(coseta[i]))
            ns = int(my_yield)
            # --- The ns+1 is only a temporary fix to avoid an array out of
            # --- bounds errors. Once the desorb routine is fixed, this
            # --- code should be updated. Note that as the code is now,
            # --- in some cases, a desorbed particle will be thrown out.
            vx = desorb.floatArray(ns+1)
            vy = desorb.floatArray(ns+1)
            vz = desorb.floatArray(ns+1)
            vxnew = zeros(ns,'d')
            vynew = zeros(ns,'d')
            vznew = zeros(ns,'d')

            #compute the desorbed neutrals
            # note that the gamma0 (1.) and rel_weight (top.pgroup.sw[js_new]/top.pgroup.sw[js]) are actually not used
            desorb.desorb(my_yield,v[0][i],v[1][i],v[2][i],theta[i],phi[i],1.,0.4,top.pgroup.sm[js],top.pgroup.sw[js_new]/top.pgroup.sw[js],vx,vy,vz)
            for ivnew in range(ns):
              vxnew[ivnew]=vx[ivnew]
              vynew[ivnew]=vy[ivnew]
              vznew[ivnew]=vz[ivnew]
            xnew = xplost[i]+n_unit0[0][i]*1.e-10*w3d.dx
            ynew = yplost[i]+n_unit0[1][i]*1.e-10*w3d.dy
            znew = zplost[i]+n_unit0[2][i]*1.e-10*w3d.dz
#            pid[:,self.xoldpid]=xnew-vxnew*top.dt
#            pid[:,self.yoldpid]=ynew-vynew*top.dt
#            pid[:,self.zoldpid]=znew-vznew*top.dt
            if w3d.solvergeom==w3d.RZgeom:
               condition = (sqrt(xnew**2+ynew**2)>xmax) or \
                           (znew<zmin) or (znew>zmax)
            else:
               condition = (xnew<xmin) or (xnew>xmax) or \
                           (ynew<ymin) or (ynew>ymax) or \
                           (znew<zmin) or (znew>zmax)
            if condition:
              print 'WARNING: new particle outside boundaries',xnew,ynew,znew
              self.outparts+=[[xnew,ynew,znew,xplost[i],yplost[i],zplost[i], \
              xplostold[i],yplostold[i],zplostold[i],n_unit0[0][i],n_unit0[1][i],n_unit0[2][i],icond]]
            else:
              if top.wpid==0:
                self.addpart(ns,xnew,ynew,znew,vxnew,vynew,vznew,js_new,itype=None)
              else:
                self.addpart(ns,xnew,ynew,znew,vxnew,vynew,vznew,js_new,ones(ns)*weight[i],None)
            
          if self.l_record_timing:
            t.finish()
            tadd+=t.micro()

    # make sure that all particles are added
    for js in self.x.keys():
      self.flushpart(js)
      
    if self.l_record_timing:t3 = time.clock()
#    print "tinit,tgen,tadd:",tinit*1.e-6,tgen*1.e-6,tprepadd*1.e-6,tadd*1.e-6
    # --- append total emitted charge in conductors emitparticles_data arrays
    for js in self.inter.keys():
      for ics,c in enumerate(self.inter[js]['conductors']):
        for ie in range(len(self.inter[js]['emitted'][ics])):
          totemit = parallelsum(self.inter[js]['emitted'][ics][ie])
          if abs(totemit)>0.:
            c.emitparticles_data += [[top.time, 
                                      totemit,
                                      top.dt,
                                      self.inter[js]['emitted_species'][ics][ie].jslist[0]]]

#    w3d.lcallscraper=0
#    particleboundaries3d()
#    w3d.lcallscraper=1
##    top.npslost=0
    if self.l_record_timing:t4 = time.clock()
#    print 'time Secondaries = ',time.clock()-t1,'s',t2-t1,t3-t2,t4-t3
    
  def call_set_params_user(self,maxsec,mat_num=None):
    if self.set_params_user is None:
      self.set_params(maxsec,mat_num)
    else:
      if type(self.set_params_user) is StringType:
        # --- This is needed primarily since a user defined function cannot be
        # --- directly picklable and so only the name is saved.
        import __main__
        try:
          self.set_params_user = __main__.__dict__[self.set_params_user]
        except KeyError:
          # --- Maybe this should raise an error?
          print "Warning: Secondaries set_params_user function '%s' is not defined."%self.set_params_user
      if callable(self.set_params_user):
        self.set_params_user(maxsec,mat_num)
      else:
        # --- Maybe this should raise an error?
        print "Warning: Secondaries set_params_user function '%s' is not callable."%self.set_params_user

  def set_params(self,maxsec,mat_num=None):
# This routine sets the parameters for a given material.
# Written by Peter Stoltz.  All values taken from Furman's paper 
#
# Here is a list of ions & material numbers:
# e- -> Cu - mat_num=1 (default)
# e- -> Stainless Steel - mat_num=2 
# H  -> Au - mat_num=3 
# He -> Au - mat_num=4 
# K  -> SS - mat_num=5 
       if mat_num is None:
         raise("Error in set_params: mat_num has not been setup.")

# Set the values
# maxsec must be ten or greater
#       pos.th  = 0.
#       pos.pn = 0.
#       pos.pt = 0.
#       pos.pz = 0.
#       pos.ns = 0
       pos.matsurf =3
#       pos.ielswitch =1
#       pos.iredswitch =1
#       pos.itrueswitch =1
#       pos.irelk =0
#       pos.pmax =0.77939
#       pos.range = 0.
#       pos.freepath = 0.
       pos.iprob = 4
#       pos.ndelerm =0
#       pos.ndeltspm =0
#       pos.np0lt0 =0
#       pos.np1gt1 =0
#       pos.totsec1 =0.
#       pos.totsec2 =0.

# Here we set material-dependent parameters...Cu (mat_num=1) is default

       pos.enpar[0] = 1.5
       pos.enpar[1] = 1.75
       pos.enpar[2] = 1.
       pos.enpar[3] = 3.75
       pos.enpar[4] = 8.5
       pos.enpar[5] = 11.5
       pos.enpar[6] = 2.5
       pos.enpar[7] = 3.0
       pos.enpar[8] = 2.5
       pos.enpar[9] = 3.0

       pos.pnpar[0] = 2.5
       pos.pnpar[1] = 3.3
       pos.pnpar[2] = 2.5
       pos.pnpar[3] = 2.5
       pos.pnpar[4] = 2.8
       pos.pnpar[5] = 1.3
       pos.pnpar[6] = 1.5
       pos.pnpar[7] = 1.5
       pos.pnpar[8] = 1.5
       pos.pnpar[9] = 1.5

       pos.dtspk = 1.8848
       pos.dtotpk = 2.1
       pos.pangsec =1.
       pos.pr =0.5
       pos.sige =2.
       pos.Ecr =0.0409
       pos.E0tspk =276.812
       pos.E0epk = 0.
       pos.E0w =60.8614
       pos.rpar1 =0.26
       pos.rpar2 =2.
       pos.tpar1 =0.66
       pos.tpar2 =0.8
       pos.tpar3 =0.7
       pos.tpar4 =1.
       pos.tpar5 =0.
       pos.tpar6 =0.
       pos.epar1 =0.26
       pos.epar2 =2.
       pos.P1rinf =0.2
       pos.P1einf =0.02
       pos.P1epk =0.49623
       pos.powts =1.54033
       pos.powe =1.
       pos.qr=0.104045


# Now check for Stainless Steel
       if (mat_num == 2):

         pos.enpar[0] = 3.9
         pos.enpar[1] = 6.2
         pos.enpar[2] = 13.
         pos.enpar[3] = 8.8
         pos.enpar[4] = 6.25
         pos.enpar[5] = 2.25
         pos.enpar[6] = 9.2
         pos.enpar[7] = 5.3
         pos.enpar[8] = 17.8
         pos.enpar[9] = 10.
  
         pos.pnpar[0] = 1.6
         pos.pnpar[1] = 2.
         pos.pnpar[2] = 1.8
         pos.pnpar[3] = 4.7
         pos.pnpar[4] = 1.8
         pos.pnpar[5] = 2.4
         pos.pnpar[6] = 1.8
         pos.pnpar[7] = 1.8
         pos.pnpar[8] = 2.3
         pos.pnpar[9] = 1.8
  
         pos.pangsec =1.
         pos.pr =0.4
         pos.sige =1.9
         pos.dtspk = 1.22
         pos.dtotpk = 2.1
         pos.Ecr =40.0
         pos.E0tspk =310.
         pos.E0epk = 0.
         pos.E0w =100.
         pos.rpar1 =0.26
         pos.rpar2 =2.
         pos.tpar1 =0.66
         pos.tpar2 =0.8
         pos.tpar3 =0.7
         pos.tpar4 =1.
         pos.tpar5 =0.
         pos.tpar6 =0.
         pos.epar1 =0.26
         pos.epar2 =2.
         pos.P1rinf =0.74
         pos.P1einf =0.07
         pos.P1epk =0.5
         pos.powts =1.813
         pos.powe =0.9
         pos.qr=1.
   
# Now check for Ion->Au 
# Because it is Ion, we set elastic and rediffused to zero
# 
# First we take H ion (data from Eder, Rev Sci Inst 68 1(1997))
       if (mat_num == 3):
         pos.enpar[0] = 3.9
         pos.enpar[1] = 6.2
         pos.enpar[2] = 13.
         pos.enpar[3] = 8.8
         pos.enpar[4] = 6.25
         pos.enpar[5] = 2.25
         pos.enpar[6] = 9.2
         pos.enpar[7] = 5.3
         pos.enpar[8] = 17.8
         pos.enpar[9] = 10.
  
         pos.pnpar[0] = 1.6
         pos.pnpar[1] = 2.
         pos.pnpar[2] = 1.8
         pos.pnpar[3] = 4.7
         pos.pnpar[4] = 1.8
         pos.pnpar[5] = 2.4
         pos.pnpar[6] = 1.8
         pos.pnpar[7] = 1.8
         pos.pnpar[8] = 2.3
         pos.pnpar[9] = 1.8
  
         pos.pangsec =1.

         pos.dtspk = 2.382
         pos.dtotpk = 2.382
         pos.E0tspk = 104624.0
         pos.powts =1.44
         pos.tpar1 =0.66
         pos.tpar2 =0.8
         pos.tpar3 =0.7
         pos.tpar4 =1.
         pos.tpar5 =0.
         pos.tpar6 =0.

         pos.sige =0.0
         pos.pr =0.0
         pos.E0epk = 0.
         pos.E0w =0.
         pos.epar1 =0.0
         pos.epar2 =0.
         pos.P1einf =0.0
         pos.P1epk =0.0
         pos.powe =0.0

         pos.Ecr =0.0
         pos.rpar1 =0.0
         pos.rpar2 =0.
         pos.P1rinf =0.0
         pos.qr=0.

# Now we take He ion (data from Eder, Rev Sci Inst 68 1(1997))
       if (mat_num == 4):
         pos.enpar[0] = 3.9
         pos.enpar[1] = 6.2
         pos.enpar[2] = 13.
         pos.enpar[3] = 8.8
         pos.enpar[4] = 6.25
         pos.enpar[5] = 2.25
         pos.enpar[6] = 9.2
         pos.enpar[7] = 5.3
         pos.enpar[8] = 17.8
         pos.enpar[9] = 10.
  
         pos.pnpar[0] = 1.6
         pos.pnpar[1] = 2.
         pos.pnpar[2] = 1.8
         pos.pnpar[3] = 4.7
         pos.pnpar[4] = 1.8
         pos.pnpar[5] = 2.4
         pos.pnpar[6] = 1.8
         pos.pnpar[7] = 1.8
         pos.pnpar[8] = 2.3
         pos.pnpar[9] = 1.8
  
         pos.pangsec =1.

         pos.dtspk = 5.68
         pos.dtotpk = 5.68
         pos.E0tspk = 1410000.0
         pos.powts =1.044
         pos.tpar1 =0.66
         pos.tpar2 =0.8
         pos.tpar3 =0.7
         pos.tpar4 =1.
         pos.tpar5 =0.
         pos.tpar6 =0.

         pos.sige =0.0
         pos.pr =0.0
         pos.E0epk = 0.
         pos.E0w =0.
         pos.epar1 =0.0
         pos.epar2 =0.
         pos.P1einf =0.0
         pos.P1epk =0.0
         pos.powe =0.0

         pos.Ecr =0.0
         pos.rpar1 =0.0
         pos.rpar2 =0.
         pos.P1rinf =0.0
         pos.qr=0.

         # This is K -> SS (see Phys. Rev. ST Accel. Beams 6, 054701 [2003])
       if (mat_num == 5):
         pos.enpar[0] = 3.9
         pos.enpar[1] = 6.2
         pos.enpar[2] = 13.
         pos.enpar[3] = 8.8
         pos.enpar[4] = 6.25
         pos.enpar[5] = 2.25
         pos.enpar[6] = 9.2
         pos.enpar[7] = 5.3
         pos.enpar[8] = 17.8
         pos.enpar[9] = 10.
         if (maxsec > 10):
                pos.enpar[10:] = 5.
 
         pos.pnpar[0] = 1.6
         pos.pnpar[1] = 2.
         pos.pnpar[2] = 1.8
         pos.pnpar[3] = 4.7
         pos.pnpar[4] = 1.8
         pos.pnpar[5] = 2.4
         pos.pnpar[6] = 1.8
         pos.pnpar[7] = 1.8
         pos.pnpar[8] = 2.3
         pos.pnpar[9] = 1.8
         if (maxsec > 10):
                pos.pnpar[10:] = 2.
 
         pos.pangsec =1.

         pos.dtspk = 55.
         pos.dtotpk = 2.382
         pos.E0tspk = 5.9e7
         pos.powts =1.25
         pos.tpar1 =0.66
         pos.tpar2 =0.8
         pos.tpar3 =0.7
         pos.tpar4 =1.
         pos.tpar5 =0.
         pos.tpar6 =0.

         pos.sige =0.0
         pos.pr =0.0
         pos.E0epk = 0.
         pos.E0w =0.
         pos.epar1 =0.0
         pos.epar2 =0.
         pos.P1einf =0.0
         pos.P1epk =0.0
         pos.powe =0.0

         pos.Ecr =0.0
         pos.rpar1 =0.0
         pos.rpar2 =0.
         pos.P1rinf =0.0
         pos.qr=0.

         pos.tpar1 = -1.
         pos.tpar2 = -1.
         pos.tpar3 = 0.
         pos.tpar4 = 0.


  def generate_secondaries(self,Ek0,costheta,weight,itype,maxsec=10):
   """
Wrapper to secondary electrons routine secelec.
 - Ek0      # energy in eV
 - costheta # cosine of the angle; costheta = 1. is normal incidence
 - weight   # weight of macroparticle
 - itype    # type of interaction
            # = 1: e- -> Cu (default)
            # = 2: e- -> Stainless Steel 
            # = 3: H  -> Au  
            # = 4: He -> Au  
            # = 5: K  -> SS 
 - maxsec   # size of secondary arrays
 - set_params_user # alternate set_params routine provided by the user
 
Given the incident energy Ek0 and costheta, subroutine SECELEC first
decides how many secondary electrons (NS) are generated in the
collision (NS is constrained to lie in [0,MAXSEC]). It then generates the
kinetic energies (in [eV]) and the angles (in [rad]) of these NS
secondaries. From the energies and angles, it computes
the normalized velocity (bn,bt,bz)=(inward normal,CCW tangent,longitudinal)
components of the secondaries (dimensionless).
   """

   if(maxsec<>pos.maxsec):
    pos.maxsec = maxsec
    pos.gchange("bincoeff")
    init_pascal_triangle(pos.nbc,pos.maxsec)
 
   if  itype<>pos.mat_number:
    pos.mat_number=itype
    self.call_set_params_user(maxsec,pos.mat_number)

   ndelerm=0
   ndeltspm=0
   np0lt0=0
   np1gt1=0
   pos.secelec(Ek0,costheta,weight, #in
          self.secelec_ns,self.secelec_un,self.secelec_ut,self.secelec_uz,self.secelec_ityps,
          self.secelec_ekstot,self.secelec_dele,self.secelec_delr,self.secelec_delts,
          maxsec,pos.enpar,pos.pnpar,pos.matsurf,
          pos.pangsec,pos.pmax,pos.pr,pos.sige,
          pos.iprob,ndelerm,ndeltspm,np0lt0,np1gt1,
          pos.dtspk,pos.Ecr,pos.E0tspk,pos.E0epk,pos.E0w,
          pos.rpar1,pos.rpar2,pos.tpar1,pos.tpar2,pos.tpar3,
          pos.tpar4,pos.tpar5,pos.tpar6,pos.epar1,pos.epar2,
          pos.P1rinf,pos.P1einf,pos.P1epk,pos.powts,pos.powe,pos.qr,pos.nbc,
          pos.rp0lt0,pos.rp1gt1,pos.rdeltspm,pos.rdelerm)
#	call secelec(Ek0,costheta,chm(n),ns,vgns,vgts,vgzs,ityps,ens,
#     + dele,delr,delts,
#     + maxsec,enpar,pnpar,matsurf,
#     + pangsec,pmax,pr,sige,
#     + iprob,ndelerm,ndeltspm,np0lt0,np1gt1,
#     + dtspk,Ecr,E0tspk,E0epk,E0w,
#     + rpar1,rpar2,tpar1,tpar2,tpar3,tpar4,tpar5,tpar6,epar1,epar2,
#     + P1rinf,P1einf,P1epk,powts,powe,qr,nbc,
#     + rp0lt0,rp1gt1,rdeltspm,rdelerm)


  def getP1elast(self,E0,costheta,material,maxsec=10):
    self.set_params(maxsec,material)
    return P1elast(E0,costheta,pos.P1einf,pos.P1epk,
 	           pos.E0epk,pos.E0w,pos.powe,pos.epar1,pos.epar2)

  def getP1rediff(self,E0,costheta,material,maxsec=10):
    self.set_params(maxsec,material)
    return P1rediff(E0,costheta,pos.Ecr,pos.rpar1,pos.rpar2,pos.qr,pos.P1rinf)

  def getdeltats(self,E0,costheta,material,maxsec=10):
    self,set_params(maxsec,material)
    return deltats(E0,costheta,pos.dtspk,pos.E0tspk,
                   pos.powts,pos.tpar1,pos.tpar2,pos.tpar3,
                   pos.tpar4,pos.tpar5,pos.tpar6)

  def generate_probabilities(self,mye0,mycostheta,mymaterial,maxsec=10,iprob=4):

    if(maxsec<>pos.maxsec):
      pos.maxsec = maxsec
      pos.gchange("bincoeff")
      init_pascal_triangle(pos.nbc,pos.maxsec)
  
#    pos.enpar = zeros(maxsec,Float)
#    pos.pnpar = zeros(maxsec,Float)

 # Initialize all parameters  #
  # 1 = Cu (default)
  # 2 = Stainless Steel
  # 3 = H+ on Au
    mat_number = mymaterial
  
    self.set_params(maxsec,mat_number)
 
#  if pos.ielswitch or pos.iredswitch:   
#    pos.dtspk=pos.dtotpk-pos.P1einf-pos.P1rinf
  
  # Seed random number generator
#  semod.rnset(0)
  
# define inputs
    Ek0 = mye0 # energy in eV
    costheta=mycostheta # cosine of the angle; costheta = 1. is normal incidence
    pos.iprob = iprob
  
# define outputs
    dele   = zeros(1,'d')
    delr   = zeros(1,'d')
    delts  = zeros(1,'d')
    prob = zeros(maxsec+1,Float)
    probts = zeros(maxsec+1,Float)
  
    pos.gen_prob(Ek0,costheta,dele,delr,delts, #in
            maxsec,pos.iprob,prob,probts,
            pos.ndelerm,pos.ndeltspm,pos.pmax,pos.np0lt0,pos.np1gt1,
            pos.dtspk,pos.Ecr,pos.E0tspk,pos.E0epk,pos.E0w,
            pos.rpar1,pos.rpar2,pos.tpar1,pos.tpar2,pos.tpar3,
            pos.tpar4,pos.tpar5,pos.tpar6,pos.epar1,pos.epar2,
            pos.P1rinf,pos.P1einf,pos.P1epk,pos.powts,pos.powe,pos.qr,pos.nbc)
  # The result 'res' is a list of 4 things
  # Here I just assign them to more useful names
    return  prob,probts

  def sey2(self,energy):
    maxsec=10

    if(maxsec<>pos.maxsec):
      pos.maxsec = maxsec
      pos.gchange("bincoeff")
      init_pascal_triangle(pos.nbc,pos.maxsec)
  
#    pos.enpar = zeros(maxsec,Float)
#    pos.pnpar = zeros(maxsec,Float)
    n=shape(energy)[0]
    s1=zeros(n,Float)
    s2=zeros(n,Float)
    s3=zeros(n,Float)
    for i in range(n):
      s1[i]=getdeltats(E0=energy[i],costheta=1.,material=2,maxsec=10)
      s2[i]=getP1rediff(E0=energy[i],costheta=1.,material=2,maxsec=10)
      s3[i]=getP1elast(E0=energy[i],costheta=1.,material=2,maxsec=10)
    return s1+s2+s3
