from warp import *

#=================================================================
#setbgrd:
#   generates array of grad B^2 data
# 
#  Arguments:
#    nx,ny,nz: number of grid cells in x,y,z. Input ignored if there is 
#     gridded B data; will calculate instead from grid data.
#    xmin,xmax,ymin,ymax,zmin,zmax: min, max values of coordinates
#     over which grad B is to be calculated.
#  Options:
#    griddedBonly: if true, just passes gridded B data to bx,by,bz.
#     if false, finds B from call to geteb.
#    symmetry: assume quadrupole symmetry of B field data if == 2.
#    zonly: if true, set top.bgrd to have only dB^2/dz (for use
#     when the other components of grad B^2 are to be calculated from
#     multipole approximations)
#getbarray:
#    fills b arrays consisting of sum of all B field data
#fillbgrd:
#  Calculates grad B^2 from Bfield data on grid and loads into top.brgd
#=================================================================
true = 1; false = 0
def setbgrd(nx=0,ny=0,nz=0,xmin=0,xmax=0,ymin=0,ymax=0,zmin=0,
            zmax=0,griddedBOnly=false,symmetry=2,zonly=false):
    print "Setting grad b**2 array"
    # check to see if there is gridded B data.  Use top.bgrdbx as test
    #  If so, use its array to define dx, dy, dz.
    global npuse,x,y,z,uzd,gaminv,bx,by,bz,bendres,bendradi,gaminv,dtl,dtr, \
           ex,ey,ez,dt
    try:
        dx=top.bgrddx;dy=top.bgrddy;dz=top.bgrddz
        nx=top.bgrdnx;ny=top.bgrdny;nz=top.bgrdnz
        # tentative: assume only one kind of element
        zmin=top.bgrdzs[1];xmin=top.bgrdxs[1];ymin=top.bgrdys[1]
    except:
       print "no allocated gridded B data; defining a grid now"
       # Note in this case, griddedBOnly shouldn't have been set true
       if griddedBOnly == true:
           print "resetting griddedBOnly to false"
           griddedBOnly == false
       # calculate dx,dy,dz
       if nx == 0 or ny == 0 or nz == 0:
           raise "Error setting grad Bsq: nx, ny, or nz = 0"
       if xmax-xmin == 0 or ymax-ymin == 0 or zmax-zmin == 0:
           raise "Error setting grad Bsq: xmin=xmax or ymin=ymax"
       dx=(xmax-xmin)/nx
       dy=(ymax-ymin)/ny
       dz=(zmax-zmin)/nz
       top.bgrdnx=nx;top.bgrdny=ny;top.bgrdnz=nz
       top.bgrdns=1
       top.bgrdnc=3
       top.nbgrd=1
       gchange("Lattice")
       gchange("BGRDdata")
       top.bgrds=1
       top.bgrddx=dx
       top.bgrddy=dy
       top.bgrddz=dz
       top.bgrdxs[1]=xmin
       top.bgrdys[1]=ymin
       top.bgrdzs[1]=zmin
       top.bgrdze[1]=zmax
       top.bgrdid[1]=1
       top.bgrdsy[1]=2
    if griddedBOnly == true:
        # set bx,by,bz to top.bgrdbx, etc.
        bx=top.bgrdbx[:,:,:,0]
        by=top.bgrdby[:,:,:,0]
        bz=top.bgrdbz[:,:,:,0]
    else:
        x,y,z=getmesh3d(xmin,dx,nx,ymin,dy,ny,zmin,dz,nz)
        print "finished getmesh3d"
        dt=top.dt;dtr=.5*dt;dtl=-dtr
        nx1=nx+1;ny1=ny+1;nz1=nz+1
        npuse=nx1*ny1*nz1
        bx=zeros((nx1,ny1,nz1),"d")
        by=zeros((nx1,ny1,nz1),"d")
        bz=zeros((nx1,ny1,nz1),"d")
        ex=zeros((nx1,ny1,nz1),"d")
        ey=zeros((nx1,ny1,nz1),"d")
        ez=zeros((nx1,ny1,nz1),"d")
        gaminv=ones((nx1,ny1,nz1),"d")
        uzd=ones((nx1,ny1,nz1),"d")
        bendres=zeros((nx1,ny1,nz1),"d")
        bendradi=zeros((nx1,ny1,nz1),"d")
        # now fetch the B array 
        print "about to call exteb3d"
        print "shapes", shape(x),shape(uzd),shape(gaminv)
        print shape(bx),shape(ex)
        w3d.exteb3d(npuse,x,y,z,uzd,gaminv,dtl,dtr,
               bx,by,bz,ex,ey,ez,top.sm[0],
                    top.sq[0],bendres,bendradi,gaminv,dt)
        print "called exteb3d"
    fillbgrd(bx,by,bz,dx,dy,dz,symmetry,zonly)
    print "Done setting grad B^2 array"
def fillbgrd(bx,by,bz,dx,dy,dz,symmetry=0,zonly=false):
    global bsq,twodxi,dbsqdx,nx,dbsqdz,nz1,dzi,nz
    #number of grid points
    nx1=shape(bx)[0]
    ny1=shape(bx)[1]
    nz1=shape(bx)[2]
    #number of cells:
    nx=nx1-1;ny=ny1-1;nz=nz1-1
    # check to see that top.bgrd has appropriate dimensions; if
    # not, print warnings
    needsgchange=0
    if zonly == false:
      if top.bgrdnc < 3:
        print "top.bgrdnc too small; fixing"
        needsgchange=1
        top.bgrdnc=3
    else:
      if top.bgrdnc != 1:
        print "top.brdnc is not 1; fixing"
        top.bgrdnc=1
    if top.bgrdnx != nx:
        print "top.bgrdnx < nx; fixing"
        needsgchange=1
        top.bgrdnx=nx
    if top.bgrdny != ny:
        print "top.bgrdny < ny; fixing"
        needsgchange=1
        top.bgrdny=ny
    if top.bgrdnz != nz:
        print "top.bgrdnz < nz; fixing"
        needsgchange=1
        top.bgrdnz=nz
    if needsgchange == 1:
        top.bgrdns=1
        gchange("BGRDdata")
    dxi=1./dx
    dyi=1./dy
    dzi=1./dz
    twodxi=.5*dxi
    twodyi=.5*dyi
    twodzi=.5*dzi
    bsq=bx*bx+by*by+bz*bz
    # equivalence dbsqdx to the first entry in top.bgrd, and similarly
    # for dbsqdy, dbsqdz
    dbsqdz=top.bgrd[0,:,:,:,0]
    if zonly == false:
      dbsqdx=top.bgrd[1,:,:,:,0]
      dbsqdy=top.bgrd[2,:,:,:,0]
      dbsqdx[1:nx,:,:]=(bsq[2:nx1,:,:]-bsq[0:nx-1,:,:])*twodxi
      dbsqdy[:,1:ny,:]=(bsq[:,2:ny1,:]-bsq[:,0:ny-1,:])*twodyi
    dbsqdz[:,:,1:nz]=(bsq[:,:,2:nz1]-bsq[:,:,0:nz-1])*twodzi
    # This takes care of interior points.  Now take care of boundaries
    # First upper boundaries, one-sided differences
    dbsqdz[:,:,nz]=(bsq[:,:,nz]-bsq[:,:,nz-1])*dzi
    if zonly == false:
      dbsqdx[nx,:,:]=(bsq[nx,:,:]-bsq[nx-1,:,:])*dxi
      dbsqdy[:,ny,:]=(bsq[:,ny,:]-bsq[:,ny-1,:])*dyi
    # now lower boundaries.  Treat differently depending on if symmetry.
    # NOTE this is still within the "if zonly == false" and so indented
      if symmetry == 0:
          # no symmetry
          dbsqdx[0,:,:]=(bsq[1,:,:]-bsq[0,:,:])*dxi
          dbsqdy[:,0,:]=(bsq[:,1,:]-bsq[:,0,:])*dyi
      else:
        if symmetry == 2:
          # quadrupole symmetry
          dbsqdx[0,:,:]=0.
          dbsqdy[:,0,:]=0.
        else:
          raise("Unimplemented data symmetry parameter")
    dbsqdz[:,:,0]=(bsq[:,:,1]-bsq[:,:,0])*dzi
