"""Smoothing and filtering utilities
class Smoother - applies [(1-a)/2, a, (1-a)/2] filter
"""
__all__ = ['Smoother']

from warp import *


class Smoother(object):
    """
  Creates a smoother object which can be used to smooth a an array in x, y and/or z.

  A number of smoothing sequences are applied, each being a succesion of n passes
  of a three points stencil with coefficients [(1-a)/2, a, (1-a)/2]. The default, a
  coefficient of 0.5 results in the stencil [0.25,0.5,0.25] which totally suppresses
  signals at the wavelength of twice the cell length in the dimension in which it is
  applied.

  For suppressing signals at longer wavelength, a stride is being used.
  For example, for suppressing exactly signals at four, six, or N times the cell
  size along a given dimension, respectively use strides of 2, 3 or N/2.

  Wide band filtering is obtained using a succession of sequences with different
  strides.

  Optionally, a compensation step may be added to a sequence, reducing smoothing
  at long wavelengths. The compensation coefficient is calculated automatically.

  Use the 'add' function to add a smoothing sequence, and the 'apply' function for
  application to an array.

  Available methods:
  add: Adds smoothing sequence
  apply: applies smoother to array f

  pltx: plots stencil of smoothing function in x
  plty: plots stencil of smoothing function in y
  pltz: plots stencil of smoothing function in z
  pltxfft: plots spectrum of smoothing function in x
  pltyfft: plots spectrum of smoothing function in y
  pltzfft: plots spectrum of smoothing function in z
    """

    def __init__(self):
        self.npass  = []
        self.alpha  = []
        self.stride  = []
        self.n      = 0

    def add(self,npass=[1,1,1],alpha=[0.5,0.5,0.5],stride=[1,1,1],l_compensate=False):
        """
    Adds smoothing sequence:
        - npass  = [ 1 , 1 , 1 ]: number of passes in x, y and z,
        - alpha  = [0.5,0.5,0.5]: smoothing coefficients in x ,y and z,
        - stride = [ 1 , 1 , 1 ]: stride in x, y and z,
        - l_compensate = False  : if True, a compensation step will be added to the sequence.
        """
        self.npass.append(array(npass))
        self.alpha.append(array(alpha))
        self.stride.append(array(stride))
        self.n+=1

        if l_compensate:
            alphac=[0,0,0]
            npassc=[0,0,0]
            for i in range(3):
                if npass[i]>0:
                    alphac[i] = npass[i]*(1.-alpha[i])+1.
#          alphac[i] = self.get_binomial_filter_factors(npass[i],3)
#          print 'test : ',alphac[i],npass[i]*(1.-alpha[i])+1.
                    npassc[i] = 1
            self.npass.append(array(npassc))
            self.alpha.append(array(alphac))
            self.stride.append(array(stride))
            self.n+=1

    def apply(self,f):
        """
        applies smoother to array f.
        """
        # --- determine rank of array
        r = rank(f)
        # 1D
        if r == 1:
            nx = shape(f)[0]
            ny=nz=1
            f = reshape(f,(nx,ny,nz))
        # 2D
        if r == 2:
            nx,ny = shape(f)
            nz=1
            f = reshape(f,(nx,ny,nz))
        # 3D
        if r == 3:
            nx,ny,nz = shape(f)

        # --- apply smoothing
        for i in range(self.n):
            if all(self.stride[i]==1):
                npass = array(self.npass)
                alpha = array(self.alpha)
                smooth3d_121(f,nx-1,ny-1,nz-1,
                             npass[i,:],
                             alpha[i,:])
            else:
                npass = array(self.npass)
                alpha = array(self.alpha)
                stride = array(self.stride)
                smooth3d_121_stride(f,nx-1,ny-1,nz-1,
                             npass[i,:],
                             alpha[i,:],
                             stride[i,:])
#                            self.npass[i,:],
#                            self.alpha[i,:],
#                            self.stride[i,:])

        # --- restore array shape for 1D and 2D arrays
        if r == 1:
            f = reshape(f,(nx))
        if r == 2:
            f = reshape(f,(nx,ny))

    def get_binomial_filter_factors(self,n,o):
        """
        set costheta as order n expansion of cos(theta)
        """
        costheta = zeros(o,'d')
        costheta[0] = 1.
        s = -1.
        f = 2.
        for i in range(2,o,2):
            costheta[i] = s / f
            f = f * (f+1.) * (f+2.)
            s*=-1.

        # get ((1 + cos(theta))/2)^n
        cf = self.binomial_expansion(1.,costheta,n)/(2.**n)

        # get coeffs for compensation
        a = cf[2]
        wcomp1 = -a/(2*a-1.)
        coefsm1 = 1./(1.+2.*wcomp1)

#    return cf,wcomp1,coefsm1
        return coefsm1

    def binomial_expansion(self,x,y,n):
        """
        performs binomial expansion
        """
        # --- we assume that x is real and y is an array of coefficients
        # initializes f
        f = 0.*y
        # initializes Pascal triangle
        pt = zeros(n)
        pt[0] = 1.
        for i in range(n):
            f += pt[i]*x**(n-i)*self.getpower(y,i)
            if i<n-1:pt[1:] = pt[1:]+pt[:-1]
        return f

    def getpower(self,y,n):
        """
        computes coefficients of y**n where y is a list of coefficients
        """
        o = shape(y)[0]
        f = y.copy()
        for k in range(n):
            x = f.copy()
            for i in range(o):
                for j in range(o):
                    if (i+j)<o:
                        f[i+j] += y[i]*x[j]
        return f

    def getstencilx(self,nx=None):
        """
        returns stencil of smoothing function in x
        """
        ny = nz = 1
        if nx is None:
            npass = array(self.npass)
            stride = array(self.stride)
            nx = 2*sum(npass[:,0]*stride[:,0]+1)+1
        a = zeros([nx,ny,nz],'d')
        a[nx/2,:,:]=1.
        self.apply(a)
        return a[:,0,0]

    def getstencily(self,ny=100):
        """
        returns stencil of smoothing function in y
        """
        nx = nz = 1
        if ny is None:
            npass = array(self.npass)
            stride = array(self.stride)
            ny = 2*sum((npass[:,1]*stride[:,1]+1))+1
        a = zeros([nx,ny,nz],'d')
        a[:,ny/2,:]=1.
        self.apply(a)
        return a[0,:,0]

    def getstencilz(self,nz=100):
        """
        returns stencil of smoothing function in z
        """
        nx = ny = 1
        if nz is None:
            npass = array(self.npass)
            stride = array(self.stride)
            nz = 2*sum((npass[:,2]*stride[:,2]+1))+1
        a = zeros([nx,ny,nz],'d')
        a[:,:,nz/2]=1.
        self.apply(a)
        return a[0,0,:]

    def pltx(self,nx=None,color='black',width=1.):
        """
        plots stencil of smoothing function in x
        """
        s = self.getstencilx(nx)
        n = shape(s)[0]
        c = arange(n)-n/2
        plg(s,c,color=color,width=width)
        ptitles('stencil','cell position','s')

    def plty(self,ny=None,color='black',width=1.):
        """
        plots stencil of smoothing function in y
        """
        s = self.getstencily(ny)
        n = shape(s)[0]
        c = arange(n)-n/2
        plg(s,c,color=color,width=width)
        ptitles('stencil','cell position','s')

    def pltz(self,nz=None,color='black',width=1.):
        """
        plots stencil of smoothing function in z
        """
        s = self.getstencilz(nz)
        n = shape(s)[0]
        c = arange(n)-n/2
        plg(s,c,color=color,width=width)
        ptitles('stencil','cell position','s')

    def pltxfft(self,nx=200,color='black',width=1.):
        """
        plots spectrum of smoothing function in x
        """
        f=abs(fft.fft(self.getstencilx(nx)))
        theta=(arange(shape(f)[0]))*2.*pi/shape(f)[0]
        plg(f,2.*pi/theta,color=color,width=width)
        limits(2.,shape(f)[0],1.e-2,1.)
        logxy(1,1)
        ptitles('spectrum','lambda/dx','s')

    def pltyfft(self,ny=200,color='black',width=1.):
        """
        plots spectrum of smoothing function in y
        """
        f=abs(fft.fft(self.getstencily(ny)))
        theta=(arange(shape(f)[0]))*2.*pi/shape(f)[0]
        plg(f,2.*pi/theta,color=color,width=width)
        limits(2.,shape(f)[0],1.e-2,1.)
        logxy(1,1)
        ptitles('spectrum','lambda/dy','s')

    def pltzfft(self,nz=200,color='black',width=1.):
        """
        plots spectrum of smoothing function in z
        """
        f=abs(fft.fft(self.getstencilz(nz)))
        theta=(arange(shape(f)[0]))*2.*pi/shape(f)[0]
        plg(f,2.*pi/theta,color=color,width=width)
        limits(2.,shape(f)[0],1.e-2,1.)
        logxy(1,1)
        ptitles('spectrum','lambda/dz','s')


if __name__ == "__main__":
    l_compensate = 1

    sm = Smoother()
    sm.add(l_compensate=l_compensate)
    sm.add(stride=[2,2,2],l_compensate=l_compensate)

    n = 100

    env = exp(-((arange(n)-n/2)*5./n)**2)
    base = env * cos(arange(n)*2*pi*5./n)

    a = base.copy()
    a += 0.2*cos(arange(n)*2*pi*0.5*n/n) * env
    a += 0.2*cos(arange(n)*2*pi*0.25*n/n) * env

    f=ones([n,n,n],'d')

    for i in range(n):
        f[0,:,i]*=a

    for i in range(n):
        f[0,i,:]*=a

    for i in range(1,n):
        f[i,:,:]=f[0,:,:]*a[i]

    f[0,:,:]*=a[i]

    winon()
    plsys(3)
    pla(base,color=red,width=2)
    pla(a)
    ptitles('before filtering','# cells','f')
    plsys(4)
    theta=(1.+arange(n))*2.*pi/n
    pla(abs(fft.fft(a)),theta,color=blue)
    limits(0.,pi)
    ptitles('before filtering','2*pi*dx/lambda','FFT(f)')

    window(1);
    ppgeneric(f[:,:,n/2])
    ptitles('f before filtering','X','Y')

    sm.apply(f)

    window(2)
    ppgeneric(f[:,:,n/2])
    ptitles('f after filtering','X','Y')

    window(0);
    plsys(5);
    pla(base,color=red,width=4)
    pla(f[:,n/2,n/2])
    ptitles('after filtering','# cells','f')
    plsys(6);
    theta=(1.+arange(n))*2.*pi/n
    pla(abs(fft.fft(f[:,n/2,n/2])),theta,color=blue)
    limits(0.,pi)
    ptitles('after filtering','2*pi*dx/lambda','FFT(f)')

    window(3)
    plsys(9)
    sm.pltx()
    plsys(10)
    sm.pltxfft(n)
