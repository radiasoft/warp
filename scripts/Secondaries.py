"""
Secondaries: class for generating secondaries
"""
from warp import *
from secondaries import *
try:
  import desorb
except:
  print 'WARNING: module desorb is not accessible.'
import timing as t
import time
import __main__

secondaries_version = "$Id: Secondaries.py,v 1.1 2005/06/21 22:00:30 jlvay Exp $"
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
  def __init__(self,isinc=None,conductors=None,issec=None,type=None,set_params_user=None,
                    xoldpid=None,yoldpid=None,zoldpid=None,min_age=None,l_verbose=0):
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
    if isinc is None:return
    for iis,js in enumerate(isinc):
#      self.condids[js]=[]
#      self.emitted[js]=[]
      for ics,cond in enumerate(conductors[iis]):
        self.add(js,cond,issec[iis][ics],type[iis][ics])
#        self.condids[js]+=[cond.condid]
#        self.emitted[js]+=[0.]
    
  def add(self,isinc,cond,issec,type,isneut=-1):
    if not self.inter.has_key(isinc):
        self.inter[isinc]={}
        for key in ['emitted','condids','issec','conductors','type','isneut']:
          self.inter[isinc][key]=[]
    self.inter[isinc]['condids']   +=[cond.condid]
    self.inter[isinc]['conductors']+=[cond]
    self.inter[isinc]['emitted']   +=[0.]
    self.inter[isinc]['issec']     +=[issec]
    self.inter[isinc]['isneut']    +=[isneut]
    self.inter[isinc]['type']      +=[type]

  def install(self):
    if not isinstalledafterfs(self.generate):
      installafterfs(self.generate)

  def generate(self):
    # theta, phi are angles from normal to surface with regard to z and x axis rescpectively
    # psi is angle to rotate warp local frame to Posinst local frame around normal
    # eta is angle between incident velocity vector and normal to surface

    t1 = time.clock()

    # reset 'emitted' list to zero
    if self.l_verbose>1:print 'call secondaries'
    for js in self.inter.keys():
     for i in range(len(self.inter[js]['emitted'])):
      self.inter[js]['emitted'][i] = 0.

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
    for js in self.inter.keys():
#      print 'js',js
      if top.npslost[js]==0 or top.it%top.ndts[js]<>0:continue
      i1 = top.inslost[js] - 1
      i2 = top.inslost[js] + top.npslost[js] - 1
      for ics,icond in enumerate(self.inter[js]['condids']):
#        print 'ics',ics
        iit = compress(top.pidlost[i1:i2,-1]==icond,arange(top.npslost[js]))
        n = len(iit)
#        print '*1',n
        if n==0:continue
        xplost = take(top.xplost[i1:i2],iit)
        yplost = take(top.yplost[i1:i2],iit)
        zplost = take(top.zplost[i1:i2],iit)
        # exclude particles out of computational box 
        condition = (xplost>xmin) & (xplost<xmax) & \
                    (yplost>ymin) & (yplost<ymax) & \
                    (zplost>zmin) & (zplost<zmax)
        # exclude particles recently created
        if self.min_age is not None:
          inittime = take(top.pidlost[i1:i2,top.tpid-1],iit)
          condition = condition & ((top.time-inittime)>2.*top.dt)      
        iit2 = compress(condition,arange(n))
        n = len(iit2)
#        print '*2',n
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
          js_new=self.inter[js]['issec'][ics]
          if js_new>-1:
           try:
            ns,bn,bt,bz,ekstot,dele,delr,delts = generate_secondaries(e0[i],
                                                                      coseta[i],
                                                                      self.inter[js]['type'][ics],
                                                                      set_params_user=self.set_params_user)
            if self.l_verbose:
              print 'nb secondaries = ',ns,' from conductor ',icond, e0[i], coseta[i],i1,i2,iit[i],top.npslost,dele,delr,delts            
           except:
             ns=0
           t.finish()
           tgen+=t.micro()
           t.start()
           if ns>0:
            self.inter[js]['emitted'][ics] += ns*top.sq[js_new]*top.sw[js_new]
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
            xnew = repeat(array([xplost[i]+n_unit0[0][i]*1.e-10*w3d.dx]),ns)
            ynew = repeat(array([yplost[i]+n_unit0[1][i]*1.e-10*w3d.dy]),ns)
            znew = repeat(array([zplost[i]+n_unit0[2][i]*1.e-10*w3d.dz]),ns)
            pid= zeros([ns,top.npidmax],Float)
            if w3d.l_inj_rec_inittime:
              pid[:,top.tpid-1]=top.time
#            pid[:,self.xoldpid]=xnew-vx*top.dt
#            pid[:,self.yoldpid]=ynew-vy*top.dt
#            pid[:,self.zoldpid]=znew-vz*top.dt
            if xnew[0]<=xmin or xnew[0]>=xmax or ynew[0]<=ymin or ynew[0]>=ymax or znew[0]<=zmin or znew[0]>=zmax:
              print 'WARNING: new particle outside boundaries',xnew[0],ynew[0],znew[0]
              __main__.outparts+=[[xnew[0],ynew[0],znew[0],xplost[i],yplost[i],zplost[i], \
              xplostold[i],yplostold[i],zplostold[i],n_unit0[0][i],n_unit0[1][i],n_unit0[2][i],icond]]
            else:
              addpart(ns,top.npid,
                    xnew,ynew,znew,
                    vx,vy,vz,
                    ones(ns,Float),
                    pid,
                    js_new+1,
                    false,
                    w3d.zmmin,
                    w3d.zmmax,
                    false)
          # emit neutrals
          if self.inter[js]['type'][ics]==5:
           js_new = self.inter[js]['isneut'][ics]
           if js_new>-1:
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
            xnew = repeat(array([xplost[i]+n_unit0[0][i]*1.e-10*w3d.dx]),ns)
            ynew = repeat(array([yplost[i]+n_unit0[1][i]*1.e-10*w3d.dy]),ns)
            znew = repeat(array([zplost[i]+n_unit0[2][i]*1.e-10*w3d.dz]),ns)
            pid= zeros([ns,top.npidmax],Float)
            if w3d.l_inj_rec_inittime:
              pid[:,top.tpid-1]=top.time
#            pid[:,self.xoldpid]=xnew-vxnew*top.dt
#            pid[:,self.yoldpid]=ynew-vynew*top.dt
#            pid[:,self.zoldpid]=znew-vznew*top.dt
            if xnew[0]<=xmin or xnew[0]>=xmax or ynew[0]<=ymin or ynew[0]>=ymax or znew[0]<=zmin or znew[0]>=zmax:
              print 'WARNING: new particle outside boundaries',xnew[0],ynew[0],znew[0]
              __main__.outparts+=[[xnew[0],ynew[0],znew[0],xplost[i],yplost[i],zplost[i], \
              xplostold[i],yplostold[i],zplostold[i],n_unit0[0][i],n_unit0[1][i],n_unit0[2][i],icond]]
            else:
              addpart(ns,top.npid,
                    xnew,ynew,znew,
                    vxnew,vynew,vznew,
                    ones(ns,Float),
                    pid,
                    js_new+1,
                    false,
                    w3d.zmmin,
                    w3d.zmmax,
                    false)
            
          t.finish()
          tadd+=t.micro()

    t3 = time.clock()
#    print "tinit,tgen,tadd:",tinit*1.e-6,tgen*1.e-6,tprepadd*1.e-6,tadd*1.e-6
    # --- append total emitted charge in conductors emitparticles_data arrays
    for js in self.inter.keys():
      for ics,c in enumerate(self.inter[js]['conductors']):
        totemit = parallelsum(self.inter[js]['emitted'][ics])
        if me==0 and abs(totemit)>0.:
          c.emitparticles_data += [[top.time, 
                                    totemit,
                                    top.dt,
                                    self.inter[js]['issec'][ics]]]

#    w3d.lcallscraper=0
#    particleboundaries3d()
#    w3d.lcallscraper=1
##    top.npslost=0
    t4 = time.clock()
    print 'time Secondaries = ',time.clock()-t1,'s',t2-t1,t3-t2,t4-t3