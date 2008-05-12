from warp import *
optimizer_version = "$Id: optimizer.py,v 1.10 2008/05/12 16:34:17 dave Exp $"
"""
This file contains several optimizers, including:
  Spsa: Simultaneaous Perturbation Stochastic Approximation
  Evolution: Genetic based algorithm
  Simpleoptimizer: simple search
"""

class Spsa:
  """
Implements the Simultaneaous Perturbation Stochastic Approximation
minimization algorithm. To use, create an instance of the Spsa class
and then call the iter method.

opt = Spsa(...)
opt.iter(...)

Note that the parameters are all varied by approximately the same amount and so
should be normalized to be of the same order of magnitude.

Constructor arguments:
  - nparams: Number of params to vary
  - params: Initial values of the parameters
  - func: Method which does any calculations given the varying parameters.
          It takes a single argument, the list of parameters.
  - lossfunc: Method which calculates the loss and returns it.
  - c1: Amount by which params are initially varied to calculate gradient
  - a1: Amount by which the params are changed, scaled by the gradient (the gradient
        is the change in loss when params are changed by c1, divided by c1).
        To check the scaling, call the method gradloss, which returns the gradient.
  - a2=100.: Scale (in iteration numbers) over which c and a are decreased
  - paramsmin=-1.e+36: Min value of the parameters. This can be a function
                       which takes the params as its single argument and returns
                       the parameter mins.
  - paramsmax=+1.e+36: Max value of the parameters. This can be a function
                       which takes the params as its single argument and returns
                       the parameter maxes.
  - paramsave=None: average value of each parameter, used to scale params
                    internally to improve performance
  - paramsrms=None: RMS values of each parameter, used to scale params
                    internally to improve performance
  - verbose=0: when true, print diagnostics
  - errmax=+1.e36: Maximum acceptable value of the error. If the error
                   is greater, the iteration is skipped.
  - saveparamhist=false: When true, saves the history of the parameters in the
                         attribute hparam. Note that the history of the loss is
                         always saved in hloss.
  """

  def __init__(self,nparams,params,func,lossfunc,c1,a1,a2=100.,
               paramsmin=None,paramsmax=None,paramsave=None,paramsrms=None,
               verbose=0,errmax=top.largepos,
               saveparamhist=false):
    """
Creates an instance of the Spsa class.
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
    self.errmax = errmax
    self.saveparamhist = saveparamhist
    if paramsmin is None: self.paramsmin = -ones(nparams)*1.e+36
    else:                 self.paramsmin = paramsmin
    if paramsmax is None: self.paramsmax = +ones(nparams)*1.e+36
    else:                 self.paramsmax = paramsmax
    self.paramsave = paramsave
    self.paramsrms = paramsrms
    self.params = self.scaledparams(self.params)
    if self.saveparamhist: self.hparam = []
  def ak(self):
    return self.a1/(self.k+self.a2)**0.602
  def ck(self):
    return self.c1/self.k**.101
  def scaledparams(self,params=None):
    if params is None: params = self.params
    if self.paramsave is not None: params = params - self.paramsave
    if self.paramsrms is not None: params = params/self.paramsrms
    return params
  def unscaledparams(self,params=None):
    if params is None: params = self.params
    if self.paramsrms is not None: params = params*self.paramsrms
    if self.paramsave is not None: params = params + self.paramsave
    return params
  def gradloss(self):
  # --- Calculated the approximated gradient of the loss function.
    deltak = (2*random.random(self.nparams)).astype(Int) - .5
    nextparams = self.constrainparams(self.params + self.ck()*deltak)
    nextparams = self.unscaledparams(nextparams)
    self.func(nextparams)
    lplus = self.loss()
    nextparams = self.constrainparams(self.params - self.ck()*deltak)
    nextparams = self.unscaledparams(nextparams)
    self.func(nextparams)
    lminus = self.loss()
    return (lplus - lminus)/(2.*self.ck()*deltak)
  def printerror(self,err):
    print "iterations = %d  error = %e %e" %(self.k,self.loss(),self.loss()/err)
  def printparams(self):
    pp = self.unscaledparams(self.params)
    for i in xrange(self.nparams): print '%15.12e'%pp[i]
  def getparamsmin(self,params):
    """Returns the min limit of parameters."""
    if type(self.paramsmin) is FunctionType:
      return self.paramsmin(self.unscaledparams(params))
    return self.paramsmin
  def getparamsmax(self,params):
    """Returns the max limit of parameters."""
    if type(self.paramsmax) is FunctionType:
      return self.paramsmax(self.unscaledparams(params))
    return self.paramsmax
  def constrainparams(self,params):
    "Makes sure all params are within bounds"
    params = self.unscaledparams(params)
    params = maximum(params,self.getparamsmin(params))
    params = minimum(params,self.getparamsmax(params))
    params = self.scaledparams(params)
    return params

  def iter(self,err=1.e-9,imax=100,verbose=None,kprint=10,kprintlogmax=3):
    """
Function to do iterations.
  - err=err=1.e-9: Convergence criteria
  - imax=10000: Maximum number of iterations
  - verbose=None: can override main value of verbose
  - kprint=10: exponential scale for printing out current status
  - kprintlogmax=3: max value of exponential scale
    """
    if verbose is None: verbose = self.verbose
    i = 0
    while (self.loss()>err and i < imax):
      i = i + 1
      # --- calculate new parameters
      #self.params = self.params - self.ak()*self.gradloss()
      dp = self.gradloss()
      if verbose:
        print "ak = %f" % self.ak()
        print "gradloss = " + repr(dp)
        print "params = " + repr(self.unscaledparams(self.params))
      oldparams = self.params + 0.
      self.params = self.constrainparams(self.params - self.ak()*dp)
      if verbose: print "new params = " + repr(self.unscaledparams(self.params))
      # --- Calculate function with new params
      if self.saveparamhist:
        self.hparam.append(self.unscaledparams(self.params))
      self.func(self.unscaledparams(self.params))
      # --- Check if loss it too great.
      if self.loss() > self.errmax:
        self.params = oldparams
        if verbose: print "Skipping"
        continue
      # --- Save the latest value of the loss function.
      self.hloss.append(self.loss())
      # --- Increment the counter
      self.k = self.k + 1
      # --- Print out loss value
      klog = int(log(self.k)/log(kprint))
      klog = min(klog,kprintlogmax)
      if ((self.k%(kprint**klog)) == 0):
        self.printerror(err)
        if klog == kprintlogmax: self.printparams()
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
  def __init__(self,npop,nparams,evaluate,params,deltas=None,
               crossover=.5,f=.7,paramsmin=None,paramsmax=None):
    """
Differential Evolution
Input:
  - npop: size of population (must be greater than 3)
  - nparams: number of parameters
  - evaluate(params): is function, given a set of parameters, returns a score
  - params: intial set of parameters
  - deltas=0.01: fractional variation of parameters to fill initial population
  - crossover=0.5: fraction of crossovers, in range [0,1)
  - f=0.7: differential factor, in range (0,1.2]
  - paramsmin=-1.e+36: Min value of the parameters. This can be a function
                       which takes the params as its single argument.
  - paramsmax=+1.e+36: Max value of the parameters. This can be a function
                       which takes the params as its single argument.
Output:
  best_params: returns parameters which give the lowest score
    """
    self.npop = npop
    self.nparams = nparams
    self.crossover = crossover
    self.f = f
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
    self.count = 0
    if paramsmin is None:
      self.paramsmin = -ones(nparams)*1.e+36
    else:
      self.paramsmin = paramsmin
    if paramsmax is None:
      self.paramsmax = +ones(nparams)*1.e+36
    else:
      self.paramsmax = paramsmax
    self.evolve_init(params,deltas)
  def best_params(self):
    "Function to return best set of parameters so far"
    imin = 0
    costmin = self.cost[imin]
    for i in xrange(1,self.npop):
      if self.cost[i] <  costmin:
        imin = i
        costmin = self.cost[i]
    return self.x1[imin,:]
  def printbestcost(self):
    print "Generation %d, best cost %f worst cost %f"% \
          (self.count,min(self.cost),max(self.cost))
  def getparamsmin(self,params):
    """Returns the min limit of parameters."""
    if type(self.paramsmin) is FunctionType: return self.paramsmin(params)
    return self.paramsmin
  def getparamsmax(self,params):
    """Returns the max limit of parameters."""
    if type(self.paramsmax) is FunctionType: return self.paramsmax(params)
    return self.paramsmax
  def constrainparams(self,params):
    "Makes sure all params are within bounds"
    params = maximum(params,self.getparamsmin(params))
    params = minimum(params,self.getparamsmax(params))
    return params

  def evolve_init(self,sample,deltas=None):
    """
Function to initialize the population.
Picks parameters randomly distributed by deltas about a base sample set of
parameters.
  - sample: is the initial set of parameters
  - delta=0.01: is the fractional variation about the sample
                It can either be a scalar or an array the same size as sample.
    """
    if deltas is None:             deltas = ones(shape(sample),'d')*0.01
    elif type(deltas) == type(1.): deltas = ones(shape(sample),'d')*deltas
    elif type(deltas) == ListType: deltas = array(deltas)
    self.x1[0,:] = sample
    self.cost[0] = self.evaluate(sample)
    for i in xrange(1,self.npop):
      trial = sample*(1.+2.*(ranf(self.x1[i,:])-.5)*deltas)
      self.x1[i,:] = self.constrainparams(trial)
      self.cost[i] = self.evaluate(self.x1[i,:])

  def evolve_reset(self):
    "Reset cost function.  Used for example if cost function is changed."
    for i in xrange(npop):
      cost[i] = evaluate(x1[i,:])

  def evolve(self,gen_max=1,nprint=100):
    """
    Do the optimization
      - gen_max=1: number of generations to run through
      - nprint=100: base frequency to print cost
    """
    self.score = self.cost[0]

    # --- Loop over the generations
    for count in xrange(gen_max):
      self.count = self.count + 1

      # --- loop through population
      for i in xrange(self.npop):

        # Mutate/Recombine

        # --- Randomly pick three vectors different from each other and 'i'.
        a = i
        b = i
        c = i
        while (a == i):                     a = int(ranf()*self.npop)
        while (b == i or b == a):           b = int(ranf()*self.npop)
        while (c == i or c == a or c == b): c = int(ranf()*self.npop)

        # --- Randomly pick the first parameter
        j = int(ranf()*self.nparams)

        # --- Load parameters into trial, performing binomial trials
        for k in xrange(self.nparams):
          if (ranf() < self.crossover or k == self.nparams-1):
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
        self.trial = self.constrainparams(self.trial)
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

      # --- Print out loss function
      if (self.count <= nprint):
        self.printbestcost()
      elif ((self.count>nprint) and (self.count<=nprint**2) and
            ((self.count%nprint)==0)):
        self.printbestcost()
      elif ( (self.count%(nprint**2)) == 0):
        self.printbestcost()
        print self.best_params()

    self.printbestcost()
    print self.best_params()


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


