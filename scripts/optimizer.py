from warp import *
import RandomArray
optimizer_version = "$Id: optimizer.py,v 1.4 2002/02/28 21:17:37 dave Exp $"
"""
This file contains several optimizers, including:
  Simultaneaous Perturbation Stochastic Approximation
  Evolution
  simple search
"""

class Spsa:
  """
Implements the Simultaneaous Perturbation Stochastic Approximation
minimization algorithm. To use, create an instance of the Spsa class
and then call the iter method.
  """
  def __init__(self,nparams,params,func,lossfunc,c1,a1,a2=100.,
               paramsmin=None,paramsmax=None,verbose=0):
    """
Creates an instance of the Spsa class.
  - nparams: Number of params to vary
  - params: Initial and current values of the parameters
  - func: Method which does any calculations given the varying parameters
  - lossfunc: Method which calculates the loss
  - c1: Amount by which params are varied (initially) to calculate gradient
  - a1: Amount params are changed by (initially), scaled by the gradient (which
       is approximately the change in loss when params are changed by c1
       divided by c1)
  - a2=100.: Scale (in iteration numbers) over which c and a are decreased
  - paramsmin=-1.e+36: Min value of the parameters. This can be a function
                       which takes self as its single argument.
  - paramsmax=+1.e+36: Max value of the parameters. This can be a function
                       which takes self as its single argument.
  - verbose=0: when true, print diagnostics
    """
    self.nparams = nparams
    self.params = params
    self.func = func
    self.loss = lossfunc
    self.k = 1
    self.a1 = a1
    self.a2 = a2
    self.c1 = c1
    self.verbose = verbose
    self.hloss = []
    if paramsmin is None:
      self.paramsmin = -ones(nparams)*1.e+36
    else:
      self.paramsmin = paramsmin
    if paramsmax is None:
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
    self.func(self.constrainparams(self.params + self.ck()*deltak))
    lplus = self.loss()
    self.func(self.constrainparams(self.params - self.ck()*deltak))
    lminus = self.loss()
    return (lplus - lminus)/(2.*self.ck()*deltak)
  def printerror(self,err):
    print "iterations = %d  error = %e %e" %(self.k,self.loss(),self.loss()/err)
  def printparams(self):
    for i in xrange(self.nparams): print '%15.12e'%self.params[i]
  def getparamsmin(self):
    """Returns the min limit of parameters."""
    if type(self.paramsmin) is FunctionType: return self.paramsmin(self)
    return self.paramsmin
  def getparamsmax(self):
    """Returns the max limit of parameters."""
    if type(self.paramsmax) is FunctionType: return self.paramsmax(self)
    return self.paramsmax
  def constrainparams(self,params):
    "Makes sure all params are within bounds"
    params = maximum(params,self.getparamsmin())
    params = minimum(params,self.getparamsmax())
    return params

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
      if self.verbose:
        print "ak = %f" % self.ak()
        print "gradloss = " + repr(dp)
        print "params = " + repr(self.params)
      self.params = self.constrainparams(self.params - self.ak()*dp)
      if self.verbose: print "new params = " + repr(self.params)
      # --- Calculate function with new params
      self.func(self.params)
      # --- Save the latest value of the loss function.
      self.hloss.append(self.loss())
      # --- Increment the counter
      self.k = self.k + 1
      # --- Print out loss function
      if (self.k <= kprint):
        self.printerror(err)
      elif ((self.k>kprint) and (self.k<=kprint**2) and ((self.k%kprint)==0)):
        self.printerror(err)
      elif ( (self.k%(kprint**2)) == 0):
        self.printerror(err)
        self.printparams()
    # --- Print out the resulting params
    self.printerror(err)
    self.printparams()



###########################################################################
###########################################################################
###########################################################################
# Performs global optimization using genetic evolution algorithms.
# Algorithm taken from Dr Dobb's Journal, April 1997, K Price, R. Storn
###########################################################################

class Evolution:
  """
Differential Evolution
  Performs global optimization using genetic evolution algorithms.
  Algorithm taken from Dr Dobb's Journal, April 1997, K Price, R. Storn
  """
  def __init__(self,npop,nparams,params,evaluate,crossover=.5,f=.7):
    """
Differential Evolution
Input:
  npop = size of population (must be greater than 3)
  nparams = number of parameters
  crossover = fraction of crossovers, in range [0,1)
  f = differential factor, in range (0,1.2]
  evaluate(params) is function, given a set of parameters, returns a score
Output:
  best_params = array hold parameters which give the lowest score
The function evolve_init can be called before evolve to initialize
the population or it can be done by hand.
    """
    self.npop = npop
    self.nparams = params
    self.crossover = crossover
    self.evaluate = evaluate
    if (crossover < 0. or crossover > 1.):
      print "Warning: crossover outside of the range [0,1)"
    if (f < 0. or f > 1.2):
      print "Warning: differential factor f outside of the range (0,1.2]"
    if (npop < 4):
      raise "Error: number of populations, npop, must be greater than 3"
    self.trial = zeros(nparams,'d')
    self.x1 = zeros((npop,nparams),'d')
    self.x2 = zeros((npop,nparams),'d')
    self.cost = zeros(npop,'d')
  def best_params(self):
    "Function to return best set of parameters so far"
    imin = 0
    costmin = cost[imin]
    for i in xrange(1,self.npop):
      if cost[i] <  costmin:
        imin = i
        costmin = cost[i]
    return self.x1[imin,:]

  def evolve_init(sample,deltas=None):
    """
Function to initialize the population.
Picks parameters randomly distributed by deltas about a base sample set of
parameters.
  - sample is the initial set of parameters
  - delta=0.01 is the fractional variation about the sample
    It can either be a scalar or an array the same size as sample.
    """
    if not deltas:
      deltas = ones(shape(sample),'d')/100.
    elif type(delta) == type(1.):
      deltas = ones(shape(sample),'d')*deltas
    self.x1[1,:] = sample
    self.cost[1] = self.evaluate(sample)
    for i in xrange(1,npop):
      self.x1[i,:] = sample*(1.+2.*(RandomArray.random(self.nparams)-.5)*deltas)
      self.cost[i] = self.evaluate(self.x1[i,:])

  def evolve_reset(self):
    "Reset cost function.  Used for example if cost function is changed."
    for i in xrange(npop):
      cost[i] = evaluate(x1[i,:])

  def evolve(self,gen_max=1):
    """
    Do the optimization
      - gen_max=1 maximum number of generations to run through
    """
    self.score = self.cost[0]

    # --- Loop over the generations
    for count in xrange(gen_max):

      # --- loop through population
      for i in xrange(self.npop):

        # Mutate/Recombine

        # --- Randomly pick three vectors different from each other and 'i'.
        a=ranf()*self.npop
        b=ranf()*self.npop
        c=ranf()*self.npop
        while (a == i): a=ranf()*self.npop
        while (b == i or b == a): b=ranf()*self.npop
        while (c == i or c == a or c == b): c=ranf()*self.npop

        # --- Randomly pick the first parameter
        j = ranf()*self.nparams

        # --- Load parameters into trial, performing binomial trials
        for k in xrange(self.nparams):
          if (ranf() < self.crossover or k == self.nparams):
            # --- Source for trial is a random vector plus weighted differential
            # --- The last parameter always comes from noisy vector
            self.trial[j] = self.x1[c,j] + self.f*(self.x1[a,j] - self.x1[b,j])
          else:
            # --- Trial parameter come from x1 itself
            self.trial[j] = self.x1[i,j]
          # --- get next parameter
          j = (j+1)%self.nparams

        # Evaluate/Select

        # --- Evaluate trial function
        self.score = self.evaluate(self.trial)

        if (self.score <= self.cost[i]):
          # --- If trial improves on x1, move trial to secondary array
          # --- and save the new score
          self.x2[i,:] = self.trial
          self.cost[i] = self.score
        else:
          # --- otherwise move x1 to secondary array
          self.x2[i,:] = self.x1[i,:]
  
      # --- End of population loop, so copy new parameters into x1
      self.x1[...] = self.x2[...]


##############################################################################
# Simple optimization over each parameter. Simple means simplistic algorithm
# not ease of use. This is probably not very robust.
class Simpleoptimizer:
  def pxpone(self,id):
    return self.params[id] + self.x[id,+1]*abs(self.params[id])
  def pxmone(self,id):
    return self.params[id] + self.x[id,-1]*abs(self.params[id])
  def paramspone(self,id):
    result = self.params + 0.
    result[id] = self.pxpone(id)
    return result
  def paramsmone(self,id):
    result = self.params + 0.
    result[id] = self.pxmone(id)
    return result
  def checklimits(self,id):
    if self.pxmone(id) < self.paramsmin[id]:
      self.x[id,-1] = self.paramsmin[id]/self.params[id] - 1. + 1.e-14
    if self.pxpone(id) > self.paramsmax[id]:
      self.x[id,+1] = self.paramsmax[id]/self.params[id] - 1. - 1.e-14
    self.x[id,+1] = min(self.x[id,+1],+10.*self.vary)
    self.x[id,-1] = max(self.x[id,-1],-10.*self.vary)
    self.x[id,+1] = max(self.x[id,+1],+1.e-14)
    self.x[id,-1] = min(self.x[id,-1],-1.e-14)
    #if self.params[id] < self.paramsmin[id] or self.params[id] > self.paramsmax[id]:
      #print id
      #raise "Params out of bounds"
  def __init__(self,params,func,loss,vary=0.01,paramsmin=None,paramsmax=None,
               maxxdisparity=1.e5 ):
    self.nparams = len(params)
    self.params = params
    self.func = func
    self.loss = loss
    self.vary = vary
    self.maxxdisparity = maxxdisparity
    if not paramsmin:
      self.paramsmin = -ones(self.nparams)*1.e+36
    else:
      if type(paramsmin)==type(array([])) or type(paramsmin)==type(list([])):
        self.paramsmin = array(paramsmin)
      else:
        self.paramsmin = ones(self.nparams)*paramsmin
    if not paramsmax:
      self.paramsmax = +ones(self.nparams)*1.e+36
    else:
      if type(paramsmax)==type(array([])) or type(paramsmax)==type(list([])):
        self.paramsmax = array(paramsmax)
      else:
        self.paramsmax = ones(self.nparams)*paramsmax
    for i in range(self.nparams):
      if not (self.paramsmin[i] < self.params[i] < self.paramsmax[i]):
        raise "ERROR: Starting value is outside the parameter limits"
        return
    self.x = zeros((self.nparams,3),'d')
    self.f = zeros((self.nparams,3),'d')
    self.func(self.params)
    self.f[:,0] = self.loss()
    for i in range(self.nparams):
      self.x[i,-1] =  - vary
      self.checklimits(i)
      self.func(self.paramsmone(i))
      self.f[i,-1] = self.loss()
      self.x[i,+1] = + vary
      self.checklimits(i)
      self.func(self.paramspone(i))
      self.f[i,+1] = self.loss()
  def minimize1d(self,id,niters=1):
    for ii in range(niters):
      self.func(self.paramspone(id))
      self.f[id,+1] = self.loss()
      self.func(self.params)
      self.f[id,0] = self.loss()
      self.func(self.paramsmone(id))
      self.f[id,-1] = self.loss()
      #print self.params
      print "loss = %e"%(self.f[id,0])
      if self.f[id,-1] < self.f[id,0] < self.f[id,1]:
        self.x[id,+1] = - self.x[id,-1] + self.x[id,+1]
        self.params[id] = self.pxmone(id)
        self.x[id,-1] = 2.*self.x[id,-1]
        self.checklimits(id)
        self.f[id,0] = self.f[id,-1]
        self.func(self.paramsmone(id))
        self.f[id,-1] = self.loss()
      elif self.f[id,-1] > self.f[id,0] > self.f[id,1]:
        self.x[id,-1] = - self.x[id,+1] + self.x[id,-1]
        self.params[id] = self.params[id]*(1.+self.x[id,+1])
        self.x[id,+1] = 2.*self.x[id,+1]
        self.checklimits(id)
        self.f[id,0] = self.f[id,+1]
        self.func(self.paramspone(id))
        self.f[id,+1] = self.loss()
      else:
        slopem1 = (self.f[id,-1] - self.f[id,0])/abs(self.x[id,-1])
        slopep1 = (self.f[id,+1] - self.f[id,0])/abs(self.x[id,+1])
        if slopem1 < slopep1:
          self.x[id,-1] = self.x[id,-1]/2.
          self.params[id] = self.pxmone(id)
          self.func(self.params)
          self.f[id,0] = self.loss()
          self.func(self.paramspone(id))
          self.f[id,+1] = self.loss()
        else:
          self.x[id,+1] = self.x[id,+1]/2.
          self.params[id] = self.pxpone(id)
          self.func(self.params)
          self.f[id,0] = self.loss()
          self.func(self.paramsmone(id))
          self.f[id,-1] = self.loss()
  def rebalancex(self):
    # --- Make sure that the difference in the sizes of x doesn't become
    # --- too large.
    xmin = min(min(abs(self.x[:,+1])),min(abs(self.x[:,-1])))
    xmax = max(max(abs(self.x[:,+1])),max(abs(self.x[:,-1])))
    if xmax/xmin > self.maxxdisparity:
      newx = 0.5*(ave(abs(self.x[:,-1]))+ave(abs(self.x[:,+1])))
      self.x[:,-1] = - newx
      self.x[:,+1] = + newx
      for i in range(self.nparams):
        self.checklimits(i)
        self.func(self.paramsmone(i))
        self.f[i,-1] = self.loss()
        self.func(self.paramspone(i))
        self.f[i,+1] = self.loss()
  def minimize(self,niters=1,nsubiters=1,tol=1.e-10):
    self.ii = 0
    while self.ii < niters and max(self.f[:,0]) > tol:
      self.rebalancex()
      self.ii = self.ii + 1
      for i in range(self.nparams):
        self.minimize1d(i,nsubiters)
      print "Iteration number %d" %(self.ii)
      print "Current values = "
      print self.params
      print "loss = %e %e"%(self.f[-1,0],self.f[-1,0]/tol)


