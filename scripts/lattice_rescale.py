"""Function to scale lattice focusing strength to yield a target value of x-plane or y-plane undepressed phase advance.

See Sven Chilton's Master's Report for full documentation: "Implementation of an iterative matching scheme for the Kapchinskij-Vladimirskij equations in WARP," UC Berkeley Nuclear Engineering Department, March 2008; LBL Report Number ?? 

The functions listed below are included in this module.  To access more detailed documentation on a given function, enter either help(function_name) or doc(function_name) into the Python command line.

  - rescalefunc(): Overarching function.  Determines the scaling factor by which to multiply the lattice focusing strength to yield the target x- or y-plane undepressed phase advance, then resets lattice quantities appropriately.  Valid for target phase advances between 0 and 180 degrees.
  - rescale(): Multiplies all allocated top variables having to do with lattice focusing strength by scale factor a.
  - phasediff(): Temporarily resets lattice focusing strength, determines the corresponding value of x-plane undepressed phase advance, and returns the difference between this and the target x-plane phase advance.
"""

########################################################################
# Import necessary packages from warp and scipy                        #
########################################################################

from warp import *
from scipy import optimize

########################################################################
# Load in lattice-related functions from matching program              #
########################################################################

import envmatch_KVinvariant as matching

########################################################################
# Add script name version and script documentation function            #
########################################################################

lattice_rescale_version = "$Id: lattice_rescale.py,v 1.2 2008/02/26 17:33:23 sven Exp $"
def lattice_rescaledoc():
  import lattice_rescale
  print  lattice_rescale.__doc__

########################################################################
# Overarching lattice rescaling function                               #
########################################################################

def rescalefunc(sigma0target,plane,da=0.1,rftol=1.e-16,
                steps=1000,error_stop=True):
  """
Determines by bisection root-finding the scale factor by which to multiply the lattice focusing strength so as to achieve a target value (sigma0target: in degrees per period) of x- or y-plane phase advance (sigma0), then resets appropriate lattice parameters.  The x- or y-plane choice (plane) is user-set.  If the x-plane is chosen, the y-plane phase advance will be consistently modified and vice-versa.  Valid for target sigma0 values between 0 and 180 degrees.

Arguments:
  - sigma0target: Desired value of undepressed phase advance in degrees
  - plane = 'x': Specifies sigma0target as an x-plane quantity
          = 'y': Specifies sigma0target as a  y-plane quantity
  - da = 0.1: Size by which to increase or decrease lattice strength scale to establish a bracket for root-finding
  - rftol = 1.e-16: Root find tolerance for lattice rescale 
  - steps = 1000: Number of evenly spaced intervals per lattice period
  - error_stop = True: Stops the program if an error trap is tripped

Notes:

See Sven Chilton's Master's Report for full documentation: "Implementation of an iterative matching scheme for the Kapchinskij-Vladimirskij equations in WARP," UC Berkeley Nuclear Engineering Department, March 2008; LBL Report Number ?? 

This package works by integrating principal orbit functions and root finding on a simple phase advance constraint equation.  The principal orbits are integrated using functions in the envelope matching package envmatch_KVinvariant.py.   

If the bisection root-finding scheme ventures into an area of parameter space corresponding to a complex undepressed phase advance, the program will print warning messages.  As long as the target undepressed phase advance is between 0 and 180 degrees, simple tests indicate that the program will still converge to the physical answer in spite of these warning messages.  However, the user is advised to check answers carefully in such exceptional situations.  
  """
  # Error Traps
  # --- Ensure plane is set correctly
  assert plane == 'x' or plane == 'y', 'Reset plane variable to x or y, in quotes' 
  # --- Check that period is physical 
  if error_stop:
    assert top.zlatperi > 0., 'Lattice period length must be strictly positive.'

  #  --- Ensure that the parameter steps is a postive integer
  steps = int(abs(steps))
  
  # Numerics
  
  # Parameters for numerical integration of the principal orbit functions.
  #   Uses the envmatch_KVinvariant package 
  matching.lperiod = top.zlatperi
  
  matching.steps = steps
  
  matching.si   = top.zlatstrt
  
  matching.cxi  = 1.;  matching.cxpi = 0.;
  matching.sxi  = 0.;  matching.sxpi = 1.;
  matching.cyi  = 1.;  matching.cypi = 0.; 
  matching.syi  = 0.;  matching.sypi = 1.;

  # Establish a bracket of lattice strength scale for the root find.
  #   This works by first calculating the phase advances for the present
  #   lattice and then taking steps up or down in strength to establish
  #   a bracket.
  
  # --- Calculate phase difference of orbit relative to target
  sigma0diff = phasediff(1.,plane,steps,error_stop,sigma0target)
  
  # --- Phase advance of lattice is within tolerance rftol of target 
  #     phase advance
  if abs(sigma0diff)/sigma0target <= rftol:
    # --- x-plane case
    if plane == 'x':
      print 'Lattice x-plane phase advance is within fractional tolerance'
      print rftol,' of target phase advance', sigma0target
      print 'Lattice strength variables (Voltages, fields, etc.)'
      print 'have not been modified.'
      return
    # --- y-plane case
    if plane == 'y':
      print 'Lattice y-plane phase advance is within fractional tolerance'
      print rftol,' of target phase advance', sigma0target
      print 'Lattice strength variables (Voltages, fields, etc.)'
      print 'have not been modified.'
      return
  # --- Phase advance of lattice is not within tolerance rftol of target 
  #     phase advance
  else:
    # --- Phase advance of lattice is higher than target: scale down  
    if sigma0diff >= 0.:
      aupper = 1.
      alower = 1.-da
      dupper = sigma0diff
      dlower = phasediff(1.-da,plane,steps,error_stop,sigma0target)
      while sign(dlower) == sign(dupper) and alower >= 0.:
        aupper = alower
        alower = alower-da
        dupper = dlower
        dlower = phasediff(alower,plane,steps,error_stop,sigma0target)
    # --- Phase advance of lattice is lower than target: scale up 
    else:
      alower = 1.
      aupper = 1.+da
      dlower = sigma0diff
      dupper = phasediff(1.+da,plane,steps,error_stop,sigma0target)
      while sign(dupper) == sign(dlower):
        alower = aupper
        aupper = aupper+da
        dlower = dupper
        dupper = phasediff(aupper,plane,steps,error_stop,sigma0target)
    
    # Use bracketed root-finding to compute the scale factor ac such that 
    # ac*kappax yields the target value of sigma0x (kappax is the x-plane 
    # lattice focusing function).    
    
    ac = optimize.bisect(phasediff,alower,aupper,
                         args=(plane,steps,error_stop,sigma0target),
                         rtol=rftol)
    
    # Reset WARP lattice with the scale found and output message to user 
    rescale(ac)
    matching.latfunc(steps,error_stop)
    
    print 'Lattice focusing strength has been scaled by a factor ',ac
    print 'to achieve a target ',plane,'-plane phase advance sigma0 ='
    print sigma0target,' deg/period.'
    print 'Lattice strength variables (Voltages, fields, etc.) have all been'
    print 'rescaled consistently.'
    
    return


########################################################################
# Auxiliary functions                                                  #
########################################################################

def rescale(a):
  """
Scale lattice variables associated with the applied focusing forces by a factor a.

Comments:
  Voltages associated with electric quadrupole elements are also rescaled.  These voltages are NOT used by the lattice force accumulations in the envelope package that are called when calculating the principal orbits.   However, they are consistently rescaled in case multipole moment strengths etc were extracted from more detailed element descriptions (e.g., plate and rod formulations of a quadrupole) and rescaled.
  """
  # Continuous focusing forces 
  top.dedr  = a*top.dedr
  top.dexdx = a*top.dexdx
  top.deydy = a*top.deydy
  top.dbdr  = a*top.dbdr
  # Hard edge quadrupoles
  if top.quads:
    top.quaddb = a*top.quaddb
    top.quadde = a*top.quadde
    top.quadvx = a*top.quadvx
    top.quadvy = a*top.quadvy
  # Hard edge multipole elements
  if top.heles:
    top.heleae = a*top.heleae
    top.heleam = a*top.heleam
  # Electric multipole moments 
  if top.emlts:
    top.emltsc = a*top.emltsc
  # Magnetic multipole moments 
  if top.mmlts:
    top.mmltsc = a*top.mmltsc
  # Gridded electric field 
  if top.egrds:
    top.egrdsc = a*top.egrdsc
  # Gridded magnetic field 
  if top.bgrds:
    top.bgrdsc = a*top.bgrdsc
  # Gridded electrostatic potential 
  if top.pgrds:
    top.pgrdsc = a*top.pgrdsc
  return

########################################################################

def phasediff(a,plane,steps,error_stop,sigma0target):
  """
Temporarily rescales the lattice focusing strength by a factor a and computes the difference between the corresponding x- or y-plane undepressed phase advance after the rescale and the target undepressed phase advance

Input parameters:

  - a: Scale factor by which the lattice focusing strength is multiplied
  - plane = 'x': Selects the x-plane
          = 'y': Selects the y-plane
  - steps: Number of evenly spaced intervals per lattice period
  - error_stop = True: Stops the program if an error trap is tripped
  - sigma0target: Target value of the undepressed phase advance sigma0

Output: Difference between the undepressed phase advance of the scaled lattice and the target undepressed phase advance.
  """
  # Rescale the lattice by a 
  rescale(a)

  # Calculate phase advance of undepressed particle orbits in
  # the x- and y-planes with the envmatch_KVinvariant package 
  matching.latfunc(steps,False)
  sigma0xlocal = top.sigma0x
  sigma0ylocal = top.sigma0y

  # Remove the lattice rescale and recalculate phase advances consistently
  rescale(1./a)
  matching.latfunc(steps,error_stop)

  # Return phase advance difference 
  if plane == 'x':
    return sigma0xlocal-sigma0target
  if plane == 'y':
    return sigma0ylocal-sigma0target

########################################################################


  
