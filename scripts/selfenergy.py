from warp import *
selfenergy_version = "$Id: selfenergy.py,v 1.1.1.1 2000/10/16 18:34:19 dave Exp $"
# These routines calculate the self-electric field energy. There is also an
# expression for the estimated emittance growth. See the paper Wrangler, IEEE
# Trans. in Nuc. Sci., 32, 2196, (1985) for more information.


def selfenergyr():
  er = zeros(w3d.nx+1,'d')
  er[1:-1] = (w3d.phi[:-2,w3d.iy_axis,0]-w3d.phi[2:,w3d.iy_axis,0])/(2.*w3d.dx)
  return (sum(w3d.xmesh*er**2)*w3d.dx/
         ((echarge*top.npmax*top.sw[0])**2/16.))*(pi*eps0)**2
def selfenergy():
  ex=zeros((w3d.nx+1,w3d.ny+1),'d')
  ey=zeros((w3d.nx+1,w3d.ny+1),'d')
  ex[1:-1,:] = (w3d.phi[:-2,:,0]-w3d.phi[2:,:,0])/(2.*w3d.dx)
  ey[:,1:-1] = (w3d.phi[:,:-2,0]-w3d.phi[:,2:,0])/(2.*w3d.dy)
  exsq = 4.*sum(sum(ex[1:-1,1:-1]**2)) + 2.*(sum(ex[0,:]**2)+sum(ex[:,0]**2))+\
         ex[0,0]**2
  eysq = 4.*sum(sum(ey[1:-1,1:-1]**2)) + 2.*(sum(ey[0,:]**2)+sum(ey[:,0]**2))+\
         ey[0,0]**2
  return 0.5*(((exsq+eysq)*w3d.dx*w3d.dy)/
          ((echarge*top.npmax*top.sw[0])**2/16.))*pi*eps0**2

hselfenergy = zeros(1201,'d')
def gethdiags():
  global hselfenergy
  if top.it > len(hselfenergy)-1:
    h = zeros(len(hselfenergy)+1200,'d')
    h[:-1200] = hselfenergy
    hselfenergy = h
  hselfenergy[top.it] = selfenergy()
installafterstep(gethdiags)

def getdeleps():
  return sqrt(1.-(hselfenergy-hselfenergy[0])/2.*
              ((top.sigma0x/top.sigmax)**2-1.))

