from warp import *
import RandomArray
spsa_version = "$Id: spsa.py,v 1.1.1.1 2000/10/16 18:34:19 dave Exp $"

class Spsa:
  """
Implements the Simultaneaous Perturbation Stochastic Approximation
minimization algorithm. To use, create an instance of the Spsa class
and then call the iter method.
  """
  def __init__(self,nparams,params,func,lossfunc,c1,a1,a2=100.,
               paramsmin=None,paramsmax=None):
    """
Creates an instance of the Spsa class.
  - nparams Number of params to vary
  - params Initial and current values of the parameters
  - func Method which does any calculations given the varying parameters
  - lossfunc Method which calculates the loss
  - c1 Amount by which params are varied (initially) to calculate gradient
  - a1 Amount params are changed by (initially), scaled by the gradient (which
       is approximately the change in loss when params are changed by c1
       divided by c1)
  - a2=100. Scale (in iteration numbers) over which c and a are decreased
  - paramsmin=-1.e+36 Min value of the parameters
  - paramsmax=+1.e+36 Max value of the parameters
    """
    self.nparams = nparams
    self.params = params
    self.func = func
    self.loss = lossfunc
    self.k = 1
    self.a1 = a1
    self.a2 = a2
    self.c1 = c1
    self.hloss = []
    if not paramsmin:
      self.paramsmin = -ones(nparams)*1.e+36
    else:
      self.paramsmin = paramsmin
    if not paramsmax:
      self.paramsmax = +ones(nparams)*1.e+36
    else:
      self.paramsmax = paramsmax
  def ak(self):
    return self.a1/(self.k+self.a2)**0.602
  def ck(self):
    return self.c1/self.k**.101
  def gradloss(self):
  # --- Calculated the approximated gradient of the loss function.
    deltak = (2*RandomArray.random(self.nparams)).astype(Int) - .5
    self.func(self.params + self.ck()*deltak)
    lplus = self.loss()
    self.func(self.params - self.ck()*deltak)
    lminus = self.loss()
    return (lplus - lminus)/(2.*self.ck()*deltak)

  def iter(self,err=1.e-9,imax=10000,kprint=10):
    """
Function to do iterations.
  - err=err=1.e-9 Convergence criteria
  - imax=10000 Maximum number of iterations
  - kprint=10 exponential scale for printing out current status
    """
    i = 0
    while (self.loss()>err and i < imax):
      i = i + 1
      # --- calculate new parameters
      #self.params = self.params - self.ak()*self.gradloss()
      dp = self.gradloss()
      #print "ak = %f" % self.ak()
      #print "gradloss = " + repr(dp)
      #print "params = " + repr(self.params)
      self.params = self.params - self.ak()*dp
      print "new params = " + repr(self.params)
      # --- Makes sure all params are within bounds
      self.params = minimum(maximum(self.params,self.paramsmin),self.paramsmax)
      #print "new params = " + repr(self.params)
      # --- Calculate function with new params
      self.func(self.params)
      # --- Save the latest value of the loss function.
      self.hloss.append(self.loss())
      # --- Increment the counter
      self.k = self.k + 1
      # --- Print out loss function
      if (self.k <= kprint):
        print "iterations = %d  error = %e %e" %(self.k,self.loss(),self.loss()/err)
      elif ((self.k>kprint) and (self.k<=kprint**2) and ((self.k%kprint)==0)):
        print "iterations = %d  error = %e %e" %(self.k,self.loss(),self.loss()/err)
      elif ( (self.k%(kprint**2)) == 0):
        print "iterations = %d  error = %e %e" %(self.k,self.loss(),self.loss()/err)
    # --- Print out the resulting params
    for i in xrange(self.nparams):
      print '%15.12e'%self.params[i]
  
