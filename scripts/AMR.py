from warp import *
from copy import *
from MeshRefinement import *
from pyOpenDX import Visualizable,DXCollection,viewboundingbox,DXImage
import time

class AMRtree(Visualizable):
    """
  Adaptive Mesh Refinement class.
    """
    def __init__(self):
      self.baseblock = None
      if w3d.solvergeom==w3d.XYZgeomMR:
        self.blocks = MRBlock()
      else:
        self.blocks = frz.basegrid
      self.colors=['yellow','red','blue','green','cyan','magenta','white']
      self.ntransit=2
        
    def getabsgrad(self,f,dx,dy,dz):
      """
    get average of absolute value of grad(f).
      """
      if rank(f)==2:
        s=shape(f)
        nx = s[0]
        ny = s[1]
        g = zeros([nx+2,ny+2],Float)
        # fill interior
        g[1:-1,1:-1] = f
        # set boundaries
        g[0,   1:-1] = 2.*f[0, :]-f[1, :]
        g[-1,  1:-1] = 2.*f[-1,:]-f[-2,:]
        g[1:-1,0   ] = 2.*f[:, 0]-f[:, 1]
        g[1:-1,-1  ] = 2.*f[:,-1]-f[:,-2]
        # computes average of gradients
        gr = 0.5*abs((g[1:-1,1:-1]-g[ :-2,1:-1])/dx) \
           + 0.5*abs((g[2:,  1:-1]-g[1:-1,1:-1])/dx) \
           + 0.5*abs((g[1:-1,1:-1]-g[1:-1, :-2])/dy) \
           + 0.5*abs((g[1:-1,2:  ]-g[1:-1,1:-1])/dy)
      else:
        s=shape(f)
        nx = s[0]
        ny = s[1]
        nz = s[2]
        g = zeros([nx+2,ny+2,nz+2],Float)
        # fill interior
        g[1:-1,1:-1,1:-1] = f
        # set boundaries
        g[0   ,1:-1,1:-1]=2.*f[0 ,: ,: ]-f[1 ,: ,: ]
        g[-1  ,1:-1,1:-1]=2.*f[-1,: ,: ]-f[-2,: ,: ]
        g[1:-1,0   ,1:-1]=2.*f[: ,0 ,: ]-f[: ,1 ,: ]
        g[1:-1,-1  ,1:-1]=2.*f[: ,-1,: ]-f[: ,-2,: ]
        g[1:-1,1:-1,0   ]=2.*f[: ,: ,0 ]-f[: ,: ,1 ]
        g[1:-1,1:-1,-1  ]=2.*f[: ,: ,-1]-f[: ,: ,-2]
        # computes average of gradients
        gr = 0.5*abs((g[1:-1,1:-1,1:-1]-g[ :-2,1:-1,1:-1])/dx) \
           + 0.5*abs((g[2:,  1:-1,1:-1]-g[1:-1,1:-1,1:-1])/dx) \
           + 0.5*abs((g[1:-1,1:-1,1:-1]-g[1:-1, :-2,1:-1])/dy) \
           + 0.5*abs((g[1:-1,2:  ,1:-1]-g[1:-1,1:-1,1:-1])/dy) \
           + 0.5*abs((g[1:-1,1:-1,1:-1]-g[1:-1,1:-1, :-2])/dz) \
           + 0.5*abs((g[1:-1,1:-1,2:  ]-g[1:-1,1:-1,1:-1])/dz)
      return gr

    def getedges_byslice(self,f,dx,dy,dz,m):
      """
    Returns array with non-zero values only at edges.
    For each line (horizontals and verticals), the values which
    are above m*max(values of f in line) are selected as edges.
      """
      # get average of absolute value of gradient of f
      fg = self.getabsgrad(f,dx,dy,dz)

      dim = rank(f)
      # get edges using vertical lines
      g1 = zeros(shape(fg),Float)
      if dim==2:
        for j in range(shape(fg)[1]):
          g1[:,j] = where(fg[:,j]>m*max(fg[:,j]),fg[:,j],0.)
      else: #dim=3
        for k in range(shape(fg)[2]):
         for j in range(shape(fg)[1]):
          g1[:,j,k] = where(fg[:,j,k]>m*max(fg[:,j,k]),fg[:,j,k],0.)
        
      # get edges using horizontal lines
      g2 = zeros(shape(fg),Float)
      if dim==2:
        for i in range(shape(fg)[0]):
          g2[i,:] = where(fg[i,:]>m*max(fg[i,:]),fg[i,:],0.)
      else: #dim=3
        for k in range(shape(fg)[2]):
         for j in range(shape(fg)[0]):
          g2[j,:,k] = where(fg[j,:,k]>m*max(fg[j,:,k]),fg[j,:,k],0.)

      # take max of g1 and g2
      g = where(g1>g2,g1,g2)

      if dim==3:
        g3 = zeros(shape(fg),Float)
        for k in range(shape(fg)[1]):
         for j in range(shape(fg)[0]):
          g3[j,k,:] = where(fg[j,k,:]>m*max(fg[j,k,:]),fg[j,k,:],0.)
        # returns max of g and g3
        g = where(g>g3,g,g3)
      return g

    def getnbcell_edges(self,f,dx,dy,dz,m,RMR):
      """
    Returns array with non-zero value RMR at edges, zero elsewhere.
    For each line (horizontals and verticals), the values which
    are above m*max(values of f in line) are selected as edges.
      """
      fg = self.getedges_byslice(f,dx,dy,dz,m)
      return where(fg>0.,RMR,0.)

    def getnbcell_rho(self,f,nmax,b=2):
      """
    returns nb cells proportional to density f, with mesh refinement factor of b
    (default=2) and maximum number of refined cells per coarse cell nmax (along
    one dimension).
      """
      # get dimension (2-D or 3-D)
      dim = rank(f)
      # get number of refinement levels
      n = nint(log(nmax)/log(b)) 
      # get nb cells proportional to f
      fg=b**(dim*(n+1))*f/maxnd(f)
      fg=where(fg>1.,fg,1.)
      return b**int(log(fg)/log(dim**b))
#      return b**int(log(fg)/log(b**dim))

    def getnbcells(self,f,dx,dy,dz,m,nmax,rmr,b=2,l_removesinglecells=1,lmax=4):
      self.MRfact=MRfact
      fg1 = self.getnbcell_edges(f,dx,dy,dz,m,rmr)
      fg2 = self.getnbcell_rho(f,nmax,b)
      f = int(where(fg1>fg2,fg1,fg2))
      # next loop removes isolated blocks of lmax cells or less
      # this needs improvements and is only 2-D for now
      if(l_removesinglecells and rank(f)==2):
        nr=shape(f)[0]
        nz=shape(f)[1]
        t = zeros([nr,nz])
        r = maxnd(f)
        while r>=1:
          sum_neighbors(where(f==r,0,1),t,nr-1,nz-1)
          fl = where(f==1 and t<=lmax,1,0)
          f = where(fl==1,r,f)
          r=r/b
      return f

    def setlist(self,f,rl,b,progressive=true,nooverlap=true):
      p = progressive
      nlevels = nint(log(maxnd(f))/log(b))+1
      dim = rank(f)
      nx = shape(f)[0]
      ny = shape(f)[1]
      if dim==3:
        nz = shape(f)[2]
      else:
        nz=0; iz=0; izm=0; l=0
      self.listblocks = [0]

      # loop all refinement levels
      for i in range(nlevels-1,0,-1):
        r = rl[i]
        ib = b**i
        if(progressive and i<nlevels-1):
          f0 = where(self.sumpatch(listpatches,nx,ny,nz,dim)>0,ib,f)
        else:
          f0 = f.copy()
        f1 = ravel(f0)
        if(progressive):
          listnodes = nonzero(f1>=ib)
        else:
          listnodes = nonzero(f1==ib)
        listpatches = []
        if(nooverlap):
          if dim==2:
            fno = zeros([nx,ny])
          else:
            fno = zeros([nx,ny,nz])
        # loop all nodes where refinement is needed
        for n in listnodes:
          ix  = 1
          iy  = 1
          ixm = 0
          iym = 0
          if dim==2:
            j = n/ny
            k = n%ny
            f0t = f0[j,k]
          else: #dim=3
            j = n/(ny*nz)
            k = n%(ny*nz)/nz
            l = n%nz
            iz  = 1
            izm = 0
            f0t = f0[j,k,l]
          # check if a patch is present
          if(progressive):
            cond = f0t>=ib
          else:
            cond = f0t==ib
          if(cond):
            tryit = true
            # try to expand the patch
            while tryit:
              ix0  = ix+0
              iy0  = iy+0
              ixm0 = ixm+0
              iym0 = iym+0
              if dim==3:
                iz0  = iz+0
                izm0 = izm+0
              # x up
              if(j+ix<nx):
                if dim==2:
                  f0t = f0[j+ix,k-iym:k+iy]
                else:
                  f0t = f0[j+ix,k-iym:k+iy,l-izm:l+iz]
                if(progressive):
                  cond = sum(f0t>=ib)>0
                else:
                  cond = sum(f0t==ib)>0
                if(nooverlap):
                  if dim==2:
                    cond2 = sum(fno[j+ix,k-iym:k+iy])==0
                  else:
                    cond2 = sum(fno[j+ix,k-iym:k+iy,l-izm:l+iz])==0
                else:
                  cond2 = true
                if(cond and cond2):
                  if(self.get_area_fraction(f0,j,k,l,ix+1,iy,iz,ixm,iym,izm,ib,p)>=r):ix+=1
              # y up
              if(k+iy<ny):
                if dim==2:
                  f0t = f0[j-ixm:j+ix,k+iy]
                else:
                  f0t = f0[j-ixm:j+ix,k+iy,l-izm:l+iz]
                if(progressive):
                  cond = sum(f0t>=ib)>0
                else:
                  cond = sum(f0t==ib)>0
                if(nooverlap):
                  if dim==2:
                    cond2 = sum(fno[j-ixm:j+ix,k+iy])==0
                  else:
                    cond2 = sum(fno[j-ixm:j+ix,k+iy,l-izm:l+iz])==0
                else:
                  cond2 = true
                if(cond and cond2):
                  if(self.get_area_fraction(f0,j,k,l,ix,iy+1,iz,ixm,iym,izm,ib,p)>=r):iy+=1
              # z up
              if dim==3:
               if(l+iz<nz):
                f0t = f0[j-ixm:j+ix,k-iym:k+iy,l+iz]
                if(progressive):
                  cond = sum(f0t>=ib)>0
                else:
                  cond = sum(f0t==ib)>0
                if(nooverlap):
                  cond2 = sum(fno[j-ixm:j+ix,k-iym:k+iy,l+iz])==0
                else:
                  cond2 = true
                if(cond and cond2):
                  if(self.get_area_fraction(f0,j,k,l,ix,iy,iz+1,ixm,iym,izm,ib,p)>=r):iz+=1
              # x down
              if(j-ixm-1>-1):
                if dim==2:
                  f0t = f0[j-ixm-1,k-iym:k+iy]
                else:
                  f0t = f0[j-ixm-1,k-iym:k+iy,l-izm:l+iz]
                if(progressive):
                  cond = sum(f0t>=ib)>0
                else:
                  cond = sum(f0t==ib)>0
                if(nooverlap):
                  if dim==2:
                    cond2 = sum(fno[j-ixm-1,k-iym:k+iy])==0
                  else:
                    cond2 = sum(fno[j-ixm-1,k-iym:k+iy,l-izm:l+iz])==0
                else:
                  cond2 = true
                if(cond and cond2):
                  if(self.get_area_fraction(f0,j,k,l,ix,iy,iz,ixm+1,iym,izm,ib,p)>=r):ixm+=1
              # y down
              if(k-iym-1>-1):
                if dim==2:
                  f0t = f0[j-ixm:j+ix,k-iym-1]
                else:
                  f0t = f0[j-ixm:j+ix,k-iym-1,l-izm:l+iz]
                if(progressive):
                  cond = sum(f0t>=ib)>0
                else:
                  cond = sum(f0t==ib)>0
                if(nooverlap):
                  if dim==2:
                    cond2 = sum(fno[j-ixm:j+ix,k-iym-1])==0
                  else:
                    cond2 = sum(fno[j-ixm:j+ix,k-iym-1,l-izm:l+iz])==0
                else:
                  cond2 = true
                if(cond and cond2):
                  if(self.get_area_fraction(f0,j,k,l,ix,iy,iz,ixm,iym+1,izm,ib,p)>=r):iym+=1
              # z down
              if(l-izm-1>-1):
                f0t = f0[j-ixm:j+ix,k-iym:k+iy,l-izm-1]
                if(progressive):
                  cond = sum(f0t>=ib)>0
                else:
                  cond = sum(f0t==ib)>0
                if(nooverlap):
                  cond2 = sum(fno[j-ixm:j+ix,k-iym:k+iy,l-izm-1])==0
                else:
                  cond2 = true
                if(cond and cond2):
                  if(self.get_area_fraction(f0,j,k,l,ix,iy,iz,ixm,iym,izm+1,ib,p)>=r):izm+=1
              if dim==2:
                if(ix==ix0 and iy==iy0 and ixm==ixm0 and iym==iym0):tryit=false
              else:
                if(ix==ix0 and iy==iy0 and iz==iz0 and 
                   ixm==ixm0 and iym==iym0 and izm==izm0):tryit=false
            if dim==2:
              f0[j-ixm:j+ix,k-iym:k+iy] = 0
              if(nooverlap):fno[j-ixm:j+ix,k-iym:k+iy] = 1
              listpatches.append([j-ixm,k-iym,ixm+ix,iym+iy])
            else:
              f0[j-ixm:j+ix,k-iym:k+iy,l-izm:l+iz] = 0
              if(nooverlap):fno[j-ixm:j+ix,k-iym:k+iy,l-izm:l+iz] = 1
              listpatches.append([j-ixm,k-iym,l-izm,ixm+ix,iym+iy,izm+iz])
        self.listblocks.insert(1,listpatches)

    def get_area_fraction(self,f,j,k,l,ix,iy,iz,ixm,iym,izm,ib,progressive=true):
      if rank(f)==2:
#        print ix,ixm,iy,iym,((ix+ixm)*(iy+iym))
        if(progressive):
          return float(sum(sum(f[j-ixm:j+ix,k-iym:k+iy]>=ib)))/((ix+ixm)*(iy+iym))
        else:
          return float(sum(sum(f[j-ixm:j+ix,k-iym:k+iy]==ib)))/((ix+ixm)*(iy+iym))
      else:
        if(progressive):
          return float(sum(sum(sum(f[j-ixm:j+ix,k-iym:k+iy,l-izm:l+iz]>=ib))))/((ix+ixm)*(iy+iym)*(iz+izm))
        else:
          return float(sum(sum(sum(f[j-ixm:j+ix,k-iym:k+iy,l-izm:l+iz]==ib))))/((ix+ixm)*(iy+iym)*(iz+izm))

    def sumpatch(self,listpatches,nx,ny,nz,dim):
      if dim==2:
        f = zeros([nx,ny])
        for patch in listpatches:
          j=patch[0]
          k=patch[1]
          ix = patch[2]
          iy = patch[3]
          f[j:j+ix,k:k+iy] += 1
      else:
        f = zeros([nx,ny,nz])
        for patch in listpatches:
          j=patch[0]
          k=patch[1]
          l=patch[2]
          ix = patch[3]
          iy = patch[4]
          iz = patch[5]
          f[j:j+ix,k:k+iy,l:l+iz] += 1
      return f


    def setblocks2d(self,rmin0, rmax0, zmin0, zmax0, dr, dz, MRfact=2, transit_rmin=2, transit_rmax=2, transit_zmin=2, transit_zmax=2):
      mothergrid = frz.basegrid
      for ii,blocks in enumerate(self.listblocks[1:]):
        i = ii+1
        if i>1:mothergrid = mothergrid.down  
        r=MRfact**i
        for patch in blocks:
            nr = nint(patch[2]*r)
            nz = nint(patch[3]*r)
            drnew = dr/r
            dznew = dz/r
            drmother = drnew*MRfact
            dzmother = dznew*MRfact
            rmin = rmin0+patch[0]*dr
            rmax = rmin+nr*drnew
            zmin = zmin0+patch[1]*dz
            zmax = zmin+nz*dznew
            
            rmin_try = rmin-transit_rmin*drmother
            t_rmin = min(transit_rmin, max(0,transit_rmin-int((rmin0-rmin_try)/drmother))) 
            rmin = rmin-t_rmin*drmother
            nr = nr+MRfact*t_rmin
            
            rmax_try = rmax+transit_rmin*drmother
            t_rmax = min(transit_rmax, max(0,transit_rmax-int((rmax_try-rmax0)/drmother))) 
            rmax = rmax+t_rmax*drmother
            nr = nr+MRfact*t_rmax
            
            zmin_try = zmin-transit_zmin*dzmother
            t_zmin = min(transit_zmin, max(0,transit_zmin-int((zmin0-zmin_try)/dzmother))) 
            zmin = zmin-t_zmin*dzmother
            nz = nz+MRfact*t_zmin
            
            zmax_try = zmax+transit_zmin*dzmother
            t_zmax = min(transit_zmax, max(0,transit_zmax-nint((zmax_try-zmax0)/dzmother))) 
            zmax = zmax+t_zmax*dzmother
            nz = nz+MRfact*t_zmax
            add_subgrid(mothergrid.gid[0],nr,nz,drnew,dznew,rmin,zmin,t_rmin*MRfact,t_rmax*MRfact,t_zmin*MRfact,t_zmax*MRfact,false)
      for i in range(frz.ngrids-1):
        self.blocks += [frz.basegrid]
      g = frz.basegrid
      self.blocks[0]=g
      for i in range(1,frz.ngrids):
        try:
          g = g.next
        except:
          g=g.down
        self.blocks[g.gid[0]-1] = g

    def setblocks(self):
      mothergrid = self.blocks
      xmin0 = w3d.xmmin; xmax0 = w3d.xmmax; dx = w3d.dx
      if w3d.solvergeom == w3d.XYZgeomMR:
        ymin0 = w3d.ymmin; ymax0 = w3d.ymmax; dy = w3d.dy
        zmin0 = w3d.zmmin; zmax0 = w3d.zmmax; dz = w3d.dz
      if w3d.solvergeom == w3d.RZgeom or w3d.solvergeom == w3d.XZgeom:
        ymin0 = w3d.zmmin; ymax0 = w3d.zmmax; dy = w3d.dz
      if w3d.solvergeom == w3d.XYgeom:
        ymin0 = w3d.ymmin; ymax0 = w3d.ymmax; dy = w3d.dy
      for ii,blocks in enumerate(self.listblocks[1:]):
        i = ii+1
        if i>1:
          if w3d.solvergeom == w3d.XYZgeomMR:
            mothergrid = mothergrid.children[0]
          else:
            mothergrid = mothergrid.down
        r=self.MRfact**i
        for patch in blocks:
          if w3d.solvergeom == w3d.XYZgeomMR:
            nx = nint(patch[3]*r)
            ny = nint(patch[4]*r)
            nz = nint(patch[5]*r)
            dxnew = dx/r
            dynew = dy/r
            dznew = dz/r
            dxmother = dxnew*self.MRfact
            dymother = dynew*self.MRfact
            dzmother = dznew*self.MRfact
            xmin = xmin0 + patch[0]*dx
            xmax = xmin  + nx*dxnew
            ymin = ymin0 + patch[1]*dy
            ymax = ymin  + ny*dynew
            zmin = zmin0 + patch[2]*dz
            zmax = zmin  + nz*dznew
            nx, xmin, xmax = self.add_transit(nx, xmin, xmax, dxmother, xmin0, xmax0)
            ny, ymin, ymax = self.add_transit(ny, ymin, ymax, dymother, ymin0, ymax0)
            nz, zmin, zmax = self.add_transit(nz, zmin, zmax, dzmother, zmin0, zmax0)
            mothergrid.addchild(None,None,[xmin,ymin,zmin,],[xmax,ymax,zmax,])
          else:
            nx = nint(patch[2]*r)
            ny = nint(patch[3]*r)
            dxnew = dx/r
            dynew = dy/r
            dxmother = dxnew*MRfact
            dymother = dynew*MRfact
            xmin = xmin0 + patch[0]*dx
            xmax = xmin  + nx*dxnew
            ymin = ymin0 + patch[1]*dy
            ymax = ymin  + ny*dynew
            nx, xmin, xmax = self.add_transit(nx, xmin, xmax, dxmother, xmin0, xmax0)
            ny, ymin, ymax = self.add_transit(ny, ymin, ymax, dymother, ymin0, ymax0)
            add_subgrid(mothergrid.gid[0],nx,ny,dxnew,dynew,xmin,ymin,
                        t_xmin*self.MRfact,t_xmax*self.MRfact,
                        t_zmin*self.MRfact,t_zmax*self.MRfact,false)

      if w3d.solvergeom == w3d.XYZgeomMR:
        self.blocks.finalize()
        registersolver(self.blocks)
      else:
        for i in range(frz.ngrids-1):
          self.blocks += [frz.basegrid]
        g = frz.basegrid
        self.blocks[0]=g
        for i in range(1,frz.ngrids):
          try:
            g = g.next
          except:
            g=g.down
          self.blocks[g.gid[0]-1] = g

    def add_transit(self, nx, xmin, xmax, dxmother, xmin0, xmax0):
      nt = self.ntransit
      n = self.MRfact
            
      xmin_try = xmin-nt*dxmother
      xmax_try = xmax+nt*dxmother
      t_xmin   = min(nt, max(0,nt-int((xmin0-xmin_try)/dxmother))) 
      t_xmax   = min(nt, max(0,nt-int((xmax_try-xmax0)/dxmother))) 
      xmin    -= t_xmin*dxmother
      xmax    += t_xmax*dxmother
      nx      += n*(t_xmin+t_xmax)
      
      return nx, xmin, xmax

    def del_blocks2d(self,g=None):
      if g==None: g=frz.basegrid
      try:
        del_blocks(g.next)
      except:
        try:
          del_blocks(g.down)
        except:
          pass    
      if g is not frz.basegrid:
        del_subgrid(g.gid[0])
      else:
        frz.ngrids=1
        g=[frz.basegrid]
        g[0].loc_part=g[0].gid[0]
        g[0].loc_part_fd=g[0].gid[0]

    def generate(self,MRfact,ref_rho,ref_edge=None,r1=0.8,r2=0.8,ntransit=2,redge=0.5):
      if w3d.solvergeom==w3d.XYZgeomMR:
        pass
      else:
        pass
      self.base = MRfact
      self.ref_rho=ref_rho
      if ref_edge is None:
        self.ref_edge = ref_rho
      else:
        self.ref_edge = ref_edge
      if(maxnd(abs(frz.basegrid.rho))==0.):return
      self.nbcells=self.getnbcells(frz.basegrid.rho,w3d.dz,w3d.dx,0.,redge,ref_rho,ref_edge,b=base)
      self.setlist(self.nbcells[:-1,:-1],[0.,r2,r2,r2,r1],2,true)
      self.setblocks2d(rmin0=w3d.xmmin, rmax0=w3d.xmmax, 
                       zmin0=w3d.zmmin, zmax0=w3d.zmmax, dr=w3d.dx, dz=w3d.dz, 
                       MRfact=MRfact, transit_rmin=ntransit, transit_rmax=ntransit, 
                       transit_zmin=ntransit, transit_zmax=ntransit)
      g = frz.basegrid
      adjust_lpfd(self.nbcells,g.nr,g.nz,g.rmin,g.rmax,g.zmin,g.zmax)
      loadrho()

    def draw_blocks2d(self,level=None,color='black',width=1.,allmesh=0,f=1):
      for i,blocks in enumerate(self.listblocks[1:]):
        if level is None or level==i:
          for patch in blocks:
            j=patch[0]
            k=patch[1]
            l = patch[2]
            h = patch[3]
            r=(float(self.MRfact)**i)/f
            nx = nint(l*r)
            ny = nint(h*r)
            if w3d.solvergeom<>w3d.RZgeom:
              xmin=w3d.xmmin
              ymin=w3d.ymmin
              dx=w3d.dx
              dy=w3d.dy
              if(allmesh):
                self.draw_mesh(nx,ny,xmin+j*dx,ymin+k*dy,dx/r,dy/r,color=self.colors[i],width=width)
              else:
                self.draw_box(ymin+k*dy, ymin+k*dy+h*dy, xmin+j*dx, xmin+j*dx+l*dx, color=self.colors[i],width=width)
            else:
              xmin=w3d.xmmin
              ymin=w3d.zmmin
              dx=w3d.dx
              dy=w3d.dz
              if(allmesh):
                self.draw_mesh(ny,nx,ymin+k*dy,xmin+j*dx,dy/r,dx/r,color=self.colors[i],width=width)
              else:
                self.draw_box(xmin+j*dx, xmin+j*dx+l*dx, ymin+k*dy, ymin+k*dy+h*dy, color=self.colors[i],width=width)
               
    def draw_mesh(self,nx,ny,xmin,ymin,dx,dy,color='black',width=1):
      x = xmin+arange(nx+1)*dx
      y = ymin+arange(ny+1)*dy
      xxmin = xmin*ones(ny+1)
      yymin = ymin*ones(nx+1)
      pldj(x,yymin,x,yymin+ny*dy,color=color,width=width)
      pldj(xxmin,y,xxmin+nx*dx,y,color=color,width=width)


    def draw_box(self,rmin, rmax, zmin, zmax, color='blue',width=1):
      pldj([zmin,zmin,zmin,zmax],
           [rmin,rmax,rmin,rmin],
           [zmax,zmax,zmin,zmax],
           [rmin,rmax,rmax,rmax],color=color,width=width)

    def createdxobject(self,kwdict={},**kw):
      """
    Create DX object drawing the object.
    - withguards=1: when true, the guard cells are included in the bounding box
    
      """
      kw.update(kwdict)
      withguards = kw.get('withguards',1)
      level = kw.get('level',None)
      dxlist = []
      for i,blocks in enumerate(self.listblocks[1:]):
       if level is None or level==i:
        for patch in blocks:
            xmin = w3d.xmmin + patch[0]*w3d.dx
            xmax =     xmin  + patch[3]*w3d.dx
            ymin = w3d.ymmin + patch[1]*w3d.dy
            ymax =     ymin  + patch[4]*w3d.dy
            zmin = w3d.zmmin + patch[2]*w3d.dz
            zmax =     zmin  + patch[5]*w3d.dz
            dxlist.append(viewboundingbox(xmin,xmax,ymin,ymax,zmin,zmax,self.colors[i]))
      self.dxobject = DXCollection(*dxlist)
    def draw(self,level=None):
        if w3d.solvergeom==w3d.XYZgeomMR:
          self.createdxobject(level=level)
          DXImage(self)
        else:
          self.draw_blocks2d(level=level)
            
            
def draw_mesh(nx,ny,xmin,ymin,dx,dy,color='black',width=1):
      x = xmin+arange(nx+1)*dx
      y = ymin+arange(ny+1)*dy
      xxmin = xmin*ones(ny+1)
      yymin = ymin*ones(nx+1)
      pldj(x,yymin,x,yymin+ny*dy,color=color,width=width)
      pldj(xxmin,y,xxmin+nx*dx,y,color=color,width=width)
   
def plphirz(grid=None,which='phi',cmin=None,cmax=None,
                 border=1,bordercolor='yellow',borderwidth=1,
                 mesh=0,meshcolor='white',meshwidth=1,meshr=1,
                 siblings=1,children=1,firstcall=1,level=1,maxlevel=0,delay=0):
    if grid is None:
        g = frz.basegrid
    else:
        g = grid

    zmin = g.zmin-0.5*g.dz
    rmin = g.rmin-0.5*g.dr
    zmax = zmin+(g.nz+1)*g.dz
    rmax = rmin+(g.nr+1)*g.dr
    if(which=='phi'):
      f = g.phi[g.nguardx:-g.nguardx,g.nguardz:-g.nguardz]
#      zmin = zmin-g.nguardz*g.dz
#      zmax = zmax+g.nguardz*g.dz
#      rmin = rmin-g.nguardx*g.dr
#      rmax = rmax+g.nguardx*g.dr
    if(which=='rho'):
      f = g.rho
    if(which=='lp'):
      f = g.loc_part
    if(which=='lpfd'):
      f = g.loc_part_fd
    if(firstcall):
      if(which=='phi'):
        if cmin is None:cmin = minnd(frz.basegrid.phi)
        if cmax is None:cmax = maxnd(frz.basegrid.phi)
      elif(which=='rho'):
        if cmin is None:cmin = minnd(frz.basegrid.rho)
        if cmax is None:cmax = maxnd(frz.basegrid.rho)
      else:
        cmin=0
        cmax=frz.ngrids
#    pli(f[g.nguardx:-g.nguardx,g.nguardz:-g.nguardz],
    pli(f,zmin,rmin,zmax,rmax,cmin=cmin,cmax=cmax)
    if(mesh):
        nr = nint(float(g.nr)/meshr)
        nz = nint(float(g.nz)/meshr)
        dr = g.dr*meshr
        dz = g.dz*meshr
        draw_mesh(nz,nr,g.zmin,g.rmin,dz,dr,color=meshcolor,width=meshwidth)
    if(border):
        draw_box(rmin, rmax, zmin, zmax, color=bordercolor,width=borderwidth)
    time.sleep(delay)
    pyg_pending()
    pyg_idler()
    if(siblings):
      try:
         plphirz(g.next,which,cmin,cmax,border,bordercolor,borderwidth,mesh,meshcolor,meshwidth,meshr,
                    siblings,children=0,firstcall=0,level=level,maxlevel=maxlevel,delay=delay)
      except:
         pass
    if(children):
      if maxlevel==0 or level<maxlevel:
        try:
          plphirz(g.down,which,cmin,cmax,border,bordercolor,borderwidth,mesh,meshcolor,meshwidth,meshr,
                     siblings,children,firstcall=0,level=level+1,maxlevel=maxlevel,delay=delay)
        except:
          pass
    if(firstcall):
      colorbar(cmin,cmax,view=plsys())
   
def plrhorz(**args):
    plphirz(which='rho',**args)

def plcondrz(grid=None,border=1,bordercolor='yellow',mesh=0,meshcolor='white',meshr=1,
                 siblings=1,children=1,firstcall=1,level=1,maxlevel=0,delay=0):
    if grid is None:
        g = frz.basegrid
    else:
        g = grid

    zmin = g.zmin-g.dz
    rmin = g.rmin-g.dr
    zmax = g.zmax+g.dz
    rmax = g.rmax+g.dr
    b = g.bndfirst
    for i in range(b.nb_conductors):
        if(i==0):
            c = b.cndfirst
        else:
            c = c.next
        color=red
        for ic in range(c.nbbnd):
            if ic>=c.nbbndred:color=green
            z=zmin+c.kk[ic]*g.dz
            x=rmin+c.jj[ic]*g.dr
            if(c.dxm[ic]<g.dr):pldj([z],[x],[z],[x-c.dxm[ic]],color=color)
            if(c.dxp[ic]<g.dr):pldj([z],[x],[z],[x+c.dxp[ic]],color=color)
            if(c.dzm[ic]<g.dz):pldj([z],[x],[z-c.dzm[ic]],[x],color=color)
            if(c.dzp[ic]<g.dz):pldj([z],[x],[z+c.dzp[ic]],[x],color=color)

    if(mesh):
        nr = nint(float(g.nr)/meshr)
        nz = nint(float(g.nz)/meshr)
        dr = g.dr*meshr
        dz = g.dz*meshr
        draw_mesh(nz,nr,g.zmin,g.rmin,dz,dr,color=meshcolor)
    if(border):
        draw_box(rmin, rmax, zmin, zmax, color=bordercolor)
    time.sleep(delay)
    pyg_pending()
    pyg_idler()
    if(siblings):
      try:
         plcondrz(g.next,border,bordercolor,mesh,meshcolor,meshr,
                    siblings,children=0,firstcall=0,level=level,maxlevel=maxlevel,delay=delay)
      except:
         pass
    if(children):
      if maxlevel==0 or level<maxlevel:
        try:
          plcondrz(g.down,border,bordercolor,mesh,meshcolor,meshr,
                     siblings,children,firstcall=0,level=level+1,maxlevel=maxlevel,delay=delay)
        except:
          pass
   
def draw_box(rmin, rmax, zmin, zmax, color='blue',width=1):
       pldj([zmin,zmin,zmin,zmax],
             [rmin,rmax,rmin,rmin],
             [zmax,zmax,zmin,zmax],
             [rmin,rmax,rmax,rmax],color=color,width=width)
 
