"""
Secondaries: class for generating secondaries
"""
from warp import *
from secondaries import *
from species import *
import txphysics
try:
  import desorb
except:
  print 'WARNING: module desorb is not accessible.'
import timing as t
import time
import __main__

secondaries_version = "$Id: Secondaries.py,v 1.3 2005/10/31 23:14:08 jlvay Exp $"
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
 - set_params_user: function that sets CMEE SEY parameters provided by the user.
                    Default is None: CMEE set_params routine is used.
 - l_verbose: sets verbosity (default=0). 
  """
  def __init__(self,isinc=None,conductors=None,issec=None,set_params_user=None,material='SS',
                    xoldpid=None,yoldpid=None,zoldpid=None,min_age=None,l_verbose=0):
    if top.wpid>0:
      raise('Error in Secondaries: variable weights is not yet implemented.')
    self.inter={}
    __main__.outparts=[]
#    self.isinc = isinc
#    self.conductors = conductors
#    self.issec = issec
#    self.type = type
    self.l_verbose=l_verbose
#    self.condids={}
#    self.emitted={}
    self.set_params_user=set_params_user
    self.min_age=min_age
    if self.min_age is not None:
      w3d.l_inj_rec_inittime=true
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
    self.emitted_e=txphysics.doubleArray(1000)
    self.emitted_bn=txphysics.doubleArray(1000)
    self.emitted_bt=txphysics.doubleArray(1000)
    self.emitted_bz=txphysics.doubleArray(1000)
    self.npmax=4096
    self.nps={}
    self.x={}
    self.y={}
    self.z={}
    self.vx={}
    self.vy={}
    self.vz={}

    if isinc is None:return
    for iis,js in enumerate(isinc):
      for ics,cond in enumerate(conductors[iis]):
        self.add(js,cond,issec[iis][ics],material[iis][ics])
    
  def add(self,incident_species=None,conductor=None,emitted_species=None,material='SS'):
    interaction_type = 0
    if material=='Cu' and incident_species.type is Electron: interaction_type=1
    if material=='SS' and incident_species.type is Electron: interaction_type=2
    if material=='Au' and incident_species.type is Hydrogen: interaction_type=3
    if material=='Au' and incident_species.type is Helium:   interaction_type=4
    if material=='SS' and incident_species.type is Potassium:interaction_type=5
    if interaction_type==0:raise('Error in Secondaries: invalid material or incident species')
    isinc=incident_species.jslist[0]
    if type(emitted_species)<>type([]):emitted_species=[emitted_species]
    issec=[]
    for e in emitted_species:
      issec.append(e.jslist[0])
    if not self.inter.has_key(isinc):
        self.inter[isinc]={}
        for key in ['emitted','condids','issec','conductors','type','isneut','incident_species','emitted_species','material']:
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

  def install(self):
    if not isinstalledafterfs(self.generate):
      installafterfs(self.generate)

  def addpart(self,nn,x,y,z,vx,vy,vz,js):
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
    self.nps[js]+=nn
      
  def flushpart(self,js):
      if self.nps[js]>0:
         nn=self.nps[js]
         addparticles(x=self.x[js][:nn],
                      y=self.y[js][:nn],
                      z=self.z[js][:nn],
                      vx=self.vx[js][:nn],
                      vy=self.vy[js][:nn],
                      vz=self.vz[js][:nn],
                      js=js)
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
    # theta, phi are angles from normal to surface with regard to z and x axis rescpectively
    # psi is angle to rotate warp local frame to Posinst local frame around normal
    # eta is angle between incident velocity vector and normal to surface

    t1 = time.clock()

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

    t2 = time.clock()
    tinit=tgen=tprepadd=tadd=0.
    # compute number of secondaries and create them
    for ints in self.inter.keys():
     incident_species=self.inter[ints]['incident_species']
     for js in incident_species.jslist:
      if self.l_verbose:print 'js',js
      if top.npslost[js]==0:continue
#      if top.npslost[js]==0 or top.it%top.ndts[js]<>0:continue
      stride=top.ndts[js]
      i1 = top.inslost[js] - 1 
      i2 = top.inslost[js] + top.npslost[js] - 1
      for ics,cond in enumerate(self.inter[js]['conductors']):
        icond = cond.condid
        if self.l_verbose:print 'ics',ics
        iit = compress(top.pidlost[i1+top.it%stride:i2:stride,-1]==icond,arange(top.it%stride,top.npslost[js],stride))
        n = len(iit)
#        print '*1',n
        if n==0:continue
        xplost = take(top.xplost[i1:i2],iit)
        yplost = take(top.yplost[i1:i2],iit)
        zplost = take(top.zplost[i1:i2],iit)
        # exclude particles out of computational box 
        if w3d.solvergeom==w3d.RZgeom:
          condition = (sqrt(xplost**2+yplost**2)<xmax) & \
                      (zplost>zmin) & (zplost<zmax)
        else:
          condition = (xplost>xmin) & (xplost<xmax) & \
                      (yplost>ymin) & (yplost<ymax) & \
                      (zplost>zmin) & (zplost<zmax)
        # exclude particles recently created
        if self.min_age is not None:
          inittime = take(top.pidlost[i1:i2,top.tpid-1],iit)
          condition = condition & ((top.time-inittime)>self.min_age*top.dt)      
        iit2 = compress(condition,arange(n))
        n = len(iit2)
        if self.l_verbose:print 'nlost=',n
        if n==0:continue
        xplost = take(xplost,iit2)
        yplost = take(yplost,iit2)
        zplost = take(zplost,iit2)
        iit    = take(iit,iit2)
        t.start()    
        uxplost = take(top.uxplost[i1:i2],iit)
        uyplost = take(top.uyplost[i1:i2],iit)
        uzplost = take(top.uzplost[i1:i2],iit)
        xplostold = take(top.pidlost[i1:i2,self.xoldpid],iit)
        yplostold = take(top.pidlost[i1:i2,self.yoldpid],iit)
        zplostold = take(top.pidlost[i1:i2,self.zoldpid],iit)
        vxplost = (xplost-xplostold)/top.dt
        vyplost = (yplost-yplostold)/top.dt
        vzplost = (zplost-zplostold)/top.dt
        gaminvlost = take(top.gaminvlost[i1:i2],iit)
        e0 = where(gaminvlost==1., \
                   0.5*top.sm[js]*sqrt(uxplost**2+uyplost**2+uzplost**2)/top.echarge,
                   (1./gaminvlost-1.)*top.sm[js]*clight**2/top.echarge)
#        print 'e0',e0,gaminvlost,uxplost,uyplost,uzplost
        v = array([vxplost,vyplost,vzplost])
#        u = array([uxplost,uyplost,uzplost])
        theta = take(top.pidlost[i1:i2,-3],iit)
        phi   = take(top.pidlost[i1:i2,-2],iit)
        costheta = cos(theta)
        sintheta = sin(theta)
        cosphi   = cos(phi)
        sinphi   = sin(phi)
        # theta is relative to the z axis, phi to the x axis in the x-y plane 
        n_unit0 = array([sintheta*cosphi,sintheta*sinphi,costheta])
        coseta = -sum(v*n_unit0)/sqrt(sum(v*v))
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
          t.finish()
          tinit+=t.micro()
          t.start()
          if self.l_verbose:print 'e0, coseta',e0[i],coseta[i]
          for ie,emitted_species in enumerate(self.inter[js]['emitted_species'][ics]):
           js_new=emitted_species.jslist[0]
           if emitted_species.type is Electron:
            if incident_species.type is Electron:
             try: # need try since will generate error if no secondaries are created
               ns,bn,bt,bz,ekstot,dele,delr,delts = generate_secondaries(e0[i],
                                                                         coseta[i],
                                                                         self.inter[js]['type'][ics],
                                                                         set_params_user=self.set_params_user)
               if self.l_verbose:
                 print 'nb secondaries = ',ns,' from conductor ',icond, e0[i], coseta[i],i1,i2,iit[i],top.npslost,dele,delr,delts            
             except:
               ns=0
            else: # incidents are atoms or ions
             if incident_species.type.__class__ is Atom:
              try:
               if self.inter[js]['material'][ics]=='SS':target_num=10025
               ns=txphysics.ion_ind_elecs(e0[i]/(1.e6*incident_species.type.A),
                                          max(0.04,coseta[i]),
                                          incident_species.type.Z,
                                          incident_species.type.A,
                                          target_num,
                                          self.emitted_e,
                                          self.emitted_bn,
                                          self.emitted_bt,
                                          self.emitted_bz)
               ns = min(ns,self.npmax)
               ekstot=zeros(ns,'d')
               bn=zeros(ns,'d')
               bt=zeros(ns,'d')
               bz=zeros(ns,'d')
               for iemit in range(ns):
                 ekstot=self.emitted_e 
                 bn[iemit]=self.emitted_bn[iemit] 
                 bt[iemit]=self.emitted_bt[iemit] 
                 bz[iemit]=self.emitted_bz[iemit] 
               if self.l_verbose:
                 print 'nb secondaries = ',ns,' from conductor ',icond, e0[i], coseta[i],i1,i2,iit[i],top.npslost          
              except:
               ns=0
            t.finish()
            tgen+=t.micro()
            t.start()
            if ns>0:
             self.inter[js]['emitted'][ics][ie] += ns*top.sq[js_new]*top.sw[js_new]
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
              sinpsi  = sqrt(1.-cospsi*cospsi)
#              bt0 = (cospsi*bt - sinpsi*bz)
#              bz0 = (sinpsi*bt + cospsi*bz)
              bt0 = -(cospsi*bt - sinpsi*bz)
              bz0 = -(sinpsi*bt + cospsi*bz)
              bn0 = bn
             else:
              bt0 = bt
              bz0 = bz
              bn0 = bn
             bxsec = cosphi[i]*costheta[i]*bt0 - sinphi[i]*bz0 + cosphi[i]*sintheta[i]*bn0
             bysec = sinphi[i]*costheta[i]*bt0 + cosphi[i]*bz0 + sinphi[i]*sintheta[i]*bn0
             bzsec =          -sintheta[i]*bt0                 +           costheta[i]*bn0
             t.finish()
             tprepadd+=t.micro()
             t.start()
#            print 'b', bn,bt,bz
#            print 'b0',bt0,bz0,bn0
#            print 'bsec',bxsec,bysec,bzsec
#            print 'x,y,z',xplost[i],yplost[i],zplost[i],sqrt(xplost[i]**2+yplost[i]**2)
             vx=bxsec*clight
             vy=bysec*clight
             vz=bzsec*clight
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
               condition = (xnew<xmin) & (xnew>xmax) & \
                           (ynew<ymin) & (ynew>ymax) & \
                           (znew<zmin) & (znew>zmax)
             if condition:
              print 'WARNING: new particle outside boundaries',xnew,ynew,znew
              __main__.outparts+=[[xnew,ynew,znew,xplost[i],yplost[i],zplost[i], \
              xplostold[i],yplostold[i],zplostold[i],n_unit0[0][i],n_unit0[1][i],n_unit0[2][i],icond]]
             else:
              self.addpart(ns,xnew,ynew,znew,vx,vy,vz,js_new)
              
           # emit neutrals
           if emitted_species.type.__class__ is not Particle and emitted_species.charge_state==0: 
            my_yield=1.+1.82e-4*exp(0.09*180./pi*arccos(coseta[i]))
            ns = int(my_yield)
            vx = desorb.floatArray(ns)
            vy = desorb.floatArray(ns)
            vz = desorb.floatArray(ns)
            vxnew = zeros(ns,'d')
            vynew = zeros(ns,'d')
            vznew = zeros(ns,'d')

            #compute the desorbed neutrals
            desorb.desorb(my_yield,v[0][i],v[1][i],v[2][i],theta[i],phi[i],1./gaminvlost[i],0.4,top.sm[js],7000,vx,vy,vz)
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
            if xnew<=xmin or xnew>=xmax or ynew<=ymin or ynew>=ymax or znew<=zmin or znew>=zmax:
              print 'WARNING: new particle outside boundaries',xnew,ynew,znew
              __main__.outparts+=[[xnew,ynew,znew,xplost[i],yplost[i],zplost[i], \
              xplostold[i],yplostold[i],zplostold[i],n_unit0[0][i],n_unit0[1][i],n_unit0[2][i],icond]]
            else:
              addpart(ns,xnew,ynew,znew,vxnew,vynew,vznew,js_new)
            
          t.finish()
          tadd+=t.micro()

    # make sure that all particles are added
    for js in self.x.keys():
      self.flushpart(js)
      
    t3 = time.clock()
#    print "tinit,tgen,tadd:",tinit*1.e-6,tgen*1.e-6,tprepadd*1.e-6,tadd*1.e-6
    # --- append total emitted charge in conductors emitparticles_data arrays
    for js in self.inter.keys():
      for ics,c in enumerate(self.inter[js]['conductors']):
        for ie in range(len(self.inter[js]['emitted'][ics])):
          totemit = parallelsum(self.inter[js]['emitted'][ics][ie])
          if me==0 and abs(totemit)>0.:
            c.emitparticles_data += [[top.time, 
                                      totemit,
                                      top.dt,
                                      self.inter[js]['emitted_species'][ics][ie]]]

#    w3d.lcallscraper=0
#    particleboundaries3d()
#    w3d.lcallscraper=1
##    top.npslost=0
    t4 = time.clock()
    print 'time Secondaries = ',time.clock()-t1,'s',t2-t1,t3-t2,t4-t3