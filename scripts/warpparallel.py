"""Functions used in parallel version.
Most important ones are the paralleldump and parallelrestore functions.
"""
from warp import *
import mpi
import __main__
import copy
warpparallel_version = "$Id: warpparallel.py,v 1.78 2008/11/19 18:30:00 dave Exp $"

def warpparalleldoc():
  import warpparallel
  print warpparallel.__doc__

top.my_index = me
top.nprocs = npes
top.nslaves = top.nprocs

# ---------------------------------------------------------------------------
def gatherallzarray(a,zaxis=0):
  """Gathers and broadcasts the data in a z-array which is decomposed
the same way as the particle domains. Each processor contributes the
data from within the particle decomposition region it owns. This works
with any array from the groups Z_arrays and Z_Moments.
 - first argument is the z-array
 - zaxis: axis which is decomposed in z
  """
  if not lparallel: return a
  # --- Get start and end of particle decomposition region
  iz1 = 0
  if me < npes-1: iz2 = top.izpslave[me+1] - 1 - top.izpslave[me]
  else:           iz2 = w3d.nz - top.izpslave[me]
  # --- Rearrange array to put the decomposed axis first
  if zaxis != 0: a = swapaxes(a,0,zaxis)
  # --- Gather and broadcast it
  result = gatherarray(a[iz1:iz2+1,...],bcast=1)
  # --- Rearrange array to put the decomposed axis back where it started
  if zaxis != 0: result = swapaxes(result,0,zaxis)
  return result
 
# ---------------------------------------------------------------------------
def scatterallzarray(a,zaxis=0):
  """Scatters the data in a z-array which is decomposed the same way as
the particle domains. Each processor contributes the data from within
the particle decomposition region it owns. This works with any array
from the groups Z_arrays and Z_Moments.
 - first argument is the z-array
 - zaxis: axis which is decomposed in z
  """
  if not lparallel: return a
  # --- Rearrange array to put the decomposed axis first
  if zaxis != 0: a = swapaxes(a,0,zaxis)
  # --- Get the appropriate subsection
  result = a[top.izpslave[me]:top.izpslave[me]+top.nzpslave[me] + 1,...]
  # --- Rearrange array to put the decomposed axis back where it started
  if zaxis != 0: result = swapaxes(result,0,zaxis)
  return result

# ---------------------------------------------------------------------------
def gatherallzfsarray(a,zaxis=0):
  """Gathers and broadcasts the data in a z-array decomposed in the same
way as the field grid. Each processor contributes the data from within
the field-solve decomposition region it owns.
 - first argument is the z-array
 - zaxis: axis which is decomposed in z
  """
  if not lparallel: return a
  # --- Get start and end of field-solve decomposition region
  iz1 = 0
  if me < npes-1: iz2 = top.izfsslave[me+1] - 1 - top.izfsslave[me]
  else:           iz2 = w3d.nz - top.izfsslave[me]
  # --- Rearrange array to put the decomposed axis first
  if zaxis != 0: a = swapaxes(a,0,zaxis)
  # --- Gather and broadcast it
  result = gatherarray(a[iz1:iz2+1,...],bcast=1)
  # --- Rearrange array to put the decomposed axis back where it started
  if zaxis != 0: result = swapaxes(result,0,zaxis)
  return result
 
# ---------------------------------------------------------------------------
def scatterallzfsarray(a,zaxis=0):
  """Scatters the data in a z-array decomposed in the same way as the
field grid. Each processor contributes the data from within the
field-solve decomposition region it owns.
 - first argument is the z-array
 - zaxis: axis which is decomposed in z
  """
  if not lparallel: return a
  # --- Rearrange array to put the decomposed axis first
  if zaxis != 0: a = swapaxes(a,0,zaxis)
  # --- Get the appropriate subsection
  result = a[top.izfsslave[me]:top.izfsslave[me]+top.nzfsslave[me] + 1,...]
  # --- Rearrange array to put the decomposed axis back where it started
  if zaxis != 0: result = swapaxes(result,0,zaxis)
  return result

#-------------------------------------------------------------------------
def convertiztope(iz):
  """Given an iz value, returns the processor number whose particle region
contains that value."""
  if 0 <= iz <= w3d.nz:
    # --- This finds all of the processors for which have iz within their
    # --- domain. The last one is selected since in the regions which
    # --- overlap, the standard is that the processor to the right has
    # --- priority for that region, i.e. the processor which has the
    # --- overlapping on it left hand edge.
    pe = compress(logical_and(less_equal(top.izpslave,iz),
                    less_equal(iz,top.izpslave+top.nzpslave)),arange(npes))[-1]
  else:
    pe = None
  return pe
convertizptope = convertiztope

def convertizfstope(iz):
  """Given an iz value, returns the processor number whose field solve region
contains that value."""
  if 0 <= iz <= w3d.nz:
    # --- This finds all of the processors for which have iz within their
    # --- domain. The last one is selected since in the regions which
    # --- overlap, the standard is that the processor to the right has
    # --- priority for that region, i.e. the processor which has the
    # --- overlapping on it left hand edge.
    pe = compress(logical_and(less_equal(top.izfsslave,iz),
                  less_equal(iz,top.izfsslave+top.nzfsslave)),arange(npes))[-1]
  else:
    pe = None
  return pe

###########################################################################
###########################################################################
###########################################################################
###########################################################################
###########################################################################
# Make a restart dump for a parallel simulation. There are three basic steps.
#  1 PE0 writes out all non-parallel data (data which is the same on all
#    processors, such as the lattice)
#  2 For parallel data, PE0 creates entries in the dump file big enough to
#    hold the data from all of the processors
#  3 All of the processors, in parallel, write out the parallel data
# An attempt was made to make the parallel restart file the same format as
# a serial restart dump. This allows a serial job to use a parallel restart
# dump without modification. It also would make is easier for a restart
# dump to be use by jobs with differing numbers of processors (not implemented
# yet). To do that, global values of scalars are written out and arrays are
# written out in the same format as serial dump.
def paralleldump(fname,attr='dump',vars=[],serial=0,histz=2,varsuffix=None,
                 verbose=false,hdf=0):

  # --- Convert attr into a list if needed
  if not (type(attr) == type([])): attr = [attr]

  # --- Gather nps from all processors
  nps_p = gatherarray(top.pgroup.nps,bcast=1)
  nps_p.shape = (top.nzprocs,top.ns)
  top.npmax = top.pgroup.npmax

  # --- Same for npslost
  npslost_p = gatherarray(top.npslost,bcast=1)
  npslost_p.shape = (top.nzprocs,top.ns)
  npslost_p0 = zeros((top.nzprocs+1,top.ns+1),'l')
  npslost_p0[1:,1:] = npslost_p

  # --- Need boundnz from the right most processor.
  boundnz_p = gatherarray(w3d.boundnz)

  # --- PE0 first writes out all of the non-parallel stuff and sets up
  # --- space for parallel stuff.
  if me==0:
    # --- PE0 creates the file
    ff = PW.PW(fname)

    # --- Write out all of the non-parallel data using the standard dump
    # --- routine. Note that the serial flag is switched on so that
    # --- pydump skips parallel variables. This also writes out any
    # --- python variables.
    fobjlist = pydump(fname=None,attr=attr,vars=vars,serial=1,ff=ff,
                      varsuffix=varsuffix,verbose=verbose,hdf=hdf,
                      returnfobjlist=1)

    # --- Might as well write out *_p data right now while I'm thinking of it.
    # --- But only if they have the attribute attr.
    for a in attr:
      if re.search(a,top.getvarattr("npslost")):
        ff.write('npslost_p@parallel',npslost_p)

    # --- Loop through all variables, getting the ones with attribute attr
    packagelist = copy.copy(package())
    packagelist.reverse()
    for p in packagelist:
      pkg = __main__.__dict__[p]
      # --- Get list of variables in the package p with attribute attr
      vlist = []
      for a in attr: vlist = vlist + pkg.varlist(a)
      for vname in vlist:
        # --- Get attributes of the variable
        a = pkg.getvarattr(vname)
        # --- Check if variable has the attribute 'parallel'
        parallelvar = re.search('parallel',a)
        if not parallelvar: continue

        # --- Check if can get a python object for the variable.
        # --- This check is for dynamic arrays - if array is unallocated,
        # --- the getpyobject routine returns None.
        v = pkg.getpyobject(vname)

        if v is None:
          # --- If v is None, then replace it with a dummy array. In this
          # --- part of the routine, the data is not actually written out,
          # --- but only the type of the variable is needed. The variable
          # --- can't just be skipped since it may be allocated on some
          # --- other processors, so space must still be made in the dump file.
          vartype = pkg.getvartype(vname)
          if   vartype == 'double':  v = zeros(1,'d')
          elif vartype == 'integer': v = zeros(1,'l')

        # --- Skip complex since they can't currently be written out to
        # --- the file.
        if type(v) == ArrayType and gettypecode(v) == Complex: continue

        # --- Set name of variable in pdb file
        if varsuffix is None:
          pdbname = vname+'@'+p
        else:
          pdbname = vname+varsuffix

        if verbose:
          print "writing "+p+"."+vname+" as "+pdbname+" -first pass"

        # --- First, deal with scalars. The values written for the scalars
        # --- are the global values. For many of these, the value on PE0
        # --- is the correct value, for example, zmmin. Exceptions are
        # --- those that differ, for example zmmaxlocal and nz. Also, space
        # --- is created in the pdb file so that all of the processors
        # --- can write there own values of each of the parallel scalars.
        # --- That is done with the 'defent' call below.
        # --- Note that parallel scalars are written out during a serial
        # --- dump.
        #if type(v) != type(array([])):
        if len(shape(v)) == 0:
          # --- First, deal with exceptions
          if p == 'w3d' and vname in ['zmmaxlocal','zmmaxp']:
            ff.write(pdbname,w3d.zmmax)
          elif p == 'w3d' and vname in ['zmminp']:
            ff.write(pdbname,w3d.zmmin)
          elif p == 'w3d' and vname in ['nzlocal','izfsmax','nz_selfe','nzp']:
            ff.write(pdbname,w3d.nz)
          elif ((p=='top' and vname in ['np','nplive','npmax']) or
                (p=='wxy' and vname in ['npmaxxy'])):
            ff.write(pdbname,sum(nps_p)[0])
          elif p=='top' and vname in ['npmaxlost']:
            ff.write(pdbname,sum(npslost_p)[0])
          else:
            # --- Otherwise, write out value as is on PE0
            ff.write(pdbname,v)
          # --- For all parallel scalars, create entry to write data too
          if not serial:
            ff.defent(vname+'@'+p+'@parallel',array([v]),(top.nzprocs,))

        elif not serial:
          # --- Now arrays...
          # --- No parallel arrays are written out during a serial dump.
          # --- The majority of the arrays are 1-D arrays which are
          # --- domain decomposed in z. For these, an entry is made in
          # --- the file big enough to fit the data from all of the
          # --- processors. Other arrays have special requirements.
          # --- First deal with the exceptions
          if vname == 'inslost' and p == 'top':
            # --- This is set to be correct globally
            iii = array([1])
            if top.ns > 1:
              iii = array([1]+list(cumsum(sum(npslost_p[:,:-1]))+array([1])))
            ff.write(pdbname,iii)
            ff.defent(vname+'@'+p+'@parallel',v,(top.nzprocs,top.ns))
          elif vname == 'npslost' and p == 'top':
            # --- This is set to be correct globally
            ff.write(pdbname,sum(npslost_p))
            ff.defent(vname+'@'+p+'@parallel',v,(top.nzprocs,top.ns))
          elif p == 'top' and vname in ['xplost','yplost','zplost',
                                        'uxplost','uyplost','uzplost',
                                        'gaminvlost','tplost']:
            # --- For the particle data, a space big enough to hold
            # --- all of the data is created.
            if sum(sum(npslost_p)) > 0:
              ff.defent(pdbname,v,(sum(sum(npslost_p)),))
          elif p == 'top' and vname == 'pidlost':
            # --- For the particle data, a space big enough to hold
            # --- all of the data is created.
            if sum(sum(npslost_p)) > 0:
              ff.defent(pdbname,v,(sum(sum(npslost_p)),top.npidlost))
          elif p == 'f3d' and vname == 'conductors':
            # --- This is written out below
            pass
          elif p == 'top' and vname == 'pgroup':
            # --- This is written out below
            pass
          elif p == 'w3d' and vname in ['rho']:
            # --- Be prepared to dump out rho in case it is needed.
            # --- For example the egun script wants rho saved.
            ff.defent(pdbname,v,(w3d.nx+1,w3d.ny+1,w3d.nz+1))
          elif p == 'w3d' and vname in ['phi']:
            # --- Be prepared to dump out phi in case it is needed.
            ff.defent(pdbname,v,
                             (w3d.nx+1,w3d.ny+1,w3d.nz+2+w3d.izextra))
          else:
            # --- The rest are domain decomposed Z arrays
            ff.defent(pdbname,v,(w3d.nz+1,))

    # --- PE0 closes the file at this point
    ff.close()

  # --- None of the processors can procede past this point until PE0 has
  # --- completed the above.
  mpi.barrier()

  # --- If only writing non-parallel data, then return here
  if serial: return

  # --- Now that that is all done, the parallel data can actually be written
  # --- out now.

  # --- This is ugly...
  # --- Each processor makes space for its conductors and pgroup data.
  # --- First check if f3d.conductors or top.pgroup is being written out.
  cattr = f3d.getvarattr('conductors')
  pattr = top.getvarattr('pgroup')
  if (max(map(lambda a:re.search(a,cattr),attr)) or
      max(map(lambda a:re.search(a,pattr),attr))):
    if me > 0: mpirecv(me-1)
    ff = PW.PW(fname,'a')
    if max(map(lambda a:re.search(a,cattr),attr)):
      if verbose: print "Writing out the conductors"
      pydumpforthonobject(ff,[''],'conductors',f3d.conductors,
                          '@conductors%d@parallel'%me,[],[],0,
                          verbose=verbose,lonlymakespace=1)
    if max(map(lambda a:re.search(a,pattr),attr)):
      if verbose: print "Writing out the pgroups"
      pydumpforthonobject(ff,[''],'pgroup',top.pgroup,
                          '@pgroup%d@parallel'%me,[],[],0,
                          verbose=verbose,lonlymakespace=1)
    ff.close()
    if me < npes-1: mpi.send(1,me+1)
    mpi.barrier()

  # --- All of the processors open the file for appending
  ff = PW.PW(fname,'a')

  # --- Now we gotta go through all of this mess again!
  # --- Loop through all variables, getting the ones with attribute attr
  # --- See comments above for more details.
  for p in package():
    pkg = packageobject(p)
    vlist = []
    for a in attr: vlist = vlist + pkg.varlist(a)
    for vname in vlist:
      v = pkg.getpyobject(vname)
      if v is None: continue
      a = pkg.getvarattr(vname)
      parallelvar = re.search('parallel',a)
      if not parallelvar: continue
      if varsuffix is None:
        pdbname = vname+'@'+p
      else:
        pdbname = vname+varsuffix
      if verbose:
        print "writing "+p+"."+vname+" as "+pdbname+" -second pass"

      # --- First, deal with scalars. Each processor writes out it's
      # --- own value of the scalars into a special array.
      if type(v) != type(array([])):

        if p == 'f3d' and vname == 'conductors':
          # --- Each process writes out its conductor object with the
          # --- process number and '@parallel' appended.
          vs = '@conductors%d@parallel'%me
          pydumpforthonobject(ff,[''],'conductors',v,vs,[],[],0,verbose,0)

        elif p == 'top' and vname == 'pgroup':
          # --- Each process writes out its pgroup object with the
          # --- process number and '@parallel' appended.
          vs = '@pgroup%d@parallel'%me
          pydumpforthonobject(ff,[''],'pgroup',v,vs,[],[],0,verbose,0)

        else:
          ff.write(vname+'@'+p+'@parallel',array([v]),indx=(me,))
          # --- That was easy.

      else:
        # --- Now arrays...
        # --- The data is now written into the space which was set aside
        # --- above by PE0.
        # --- First the exceptions
        if p == 'top' and vname in ['ins','nps','inslost','npslost']:
          # --- Write out to parallel space
          ff.write(vname+'@'+p+'@parallel',array([v]),indx=(me,0))
        elif p == 'top' and vname in ['xplost','yplost','zplost',
                                      'uxplost','uyplost','uzplost',
                                      'gaminvlost','tplost']:
          # --- Write out each species seperately.
          for js in range(top.ns):
            if top.npslost[js] > 0:
              ipmin = sum(sum(npslost_p0[:,0:js+1])) + sum(npslost_p0[:me+1,js+1])
              ff.write(pdbname,v[top.inslost[js]-1:top.inslost[js]+top.npslost[js]-1],
                       indx=(ipmin,))
        elif p == 'top' and vname == 'pidlost':
          # --- Write out each species seperately.
          for js in range(top.ns):
            if top.npslost[js] > 0:
              ipmin = sum(sum(npslost_p0[:,0:js+1])) + sum(npslost_p0[:me+1,js+1])
              ff.write(pdbname,v[top.inslost[js]-1:top.inslost[js]+top.npslost[js]-1,:],
                       indx=(ipmin,0))
        elif p == 'w3d' and vname in ['rho']:
          iz1 = 0
          if me < npes-1: iz2 = top.izfsslave[me+1] - top.izfsslave[me]
          else:           iz2 = iz1 + top.nzfsslave[me] + 1
          ppp = w3d.rho[:,:,iz1:iz2]
          ff.write(pdbname,ppp,indx=(0,0,top.izfsslave[me]))
        elif p == 'w3d' and vname in ['phi']:
          iz1 = 0
          if me == 0: iz1 = iz1 - 1
          if me < npes-1: iz2 = top.izfsslave[me+1] - top.izfsslave[me]
          else:           iz2 = iz1 + top.nzfsslave[me] + 1
          ppp = w3d.phi[:,:,iz1+1:iz2+1]
          if me == 0: izmin = 0
          else:       izmin = top.izfsslave[me]+1
          ff.write(pdbname,ppp,indx=(0,0,izmin))
        elif p == 'w3d' and vname in ['zmeshlocal']:
          # --- This has the same decomposition as the fields
          iz1 = 0
          iz2 = top.nzfsslave[me] + 1
          ff.write(pdbname,v[iz1:iz2],indx=(top.izfsslave[me],))
        else:
          # --- The rest are domain decomposed Z arrays
          # --- Assuming they have the same decomposition as the particles.
          iz1 = 0
          iz2 = top.nzpslave[me] + 1
          ff.write(pdbname,v[iz1:iz2],indx=(top.izpslave[me],))

  # --- Everybody closes the file at this point
  ff.close()

  # --- The parallel dump is done.

############################################################################
# This assumes that the dump was made using the paralleldump routine
# with the default attribute, 'dump'.
# There are two main steps
#  1 Read in some initial data which is needed to setup array sizes or
#    needed in the second step
#  2 Read the rest of the data in
#
def parallelrestore(fname,verbose=false,skip=[],varsuffix=None,ls=0,lreturnff=0):
  # --- All PE's open the file for reading.
  ff = PR.PR(fname)

  # --- Long list of parallel variables that are to be skipped in the serial
  # --- restore.
  skipparallel = ['pgroup',
    'xplost','yplost','zplost','uxplost','uyplost','uzplost','gaminvlost',
    'tplost','pidlost','conductors','rho','phi','zmeshlocal']

  # --- Read in all of the serial data. This is only really needed
  # --- to deal with fortran derived types properly.
  fobjdict = pyrestore(verbose=verbose,ff=ff,varsuffix=varsuffix,ls=ls,
                       skip=skip+skipparallel,lreturnfobjdict=1)

  # --- Get list of all of the variables in the restart file
  vlist = ff.inquire_names()

  # --- Get a list of all of the variables with the parallel suffix.
  groups = sortrestorevarsbysuffix(vlist,skip)
  vlistparallel = groups['parallel']

  # --- Remove skipped variables from vlistparallel
  vlistparallelcopy = 1*vlistparallel
  for v in vlistparallelcopy:
    if v in skip:
      vlistparallel.remove(v)
      continue
    if len(v) > 4 and v[-4] == '@':
      if v[:-4] in skip or v[-3:]+'.'+v[:-4] in skip: vlistparallel.remove(v)
  del vlistparallelcopy

  # --- First, setup some arrays that need special handling.
  # --- The particles arrays are returned to there original size as in
  # --- the simulation which made the restart dump. The following
  # --- variables are needed to get that setup correctly.
  # --- In all cases, check if the variable was written out first.
  if 'npmax@top' in vlistparallel:
    itriple = array([me,me,1])
    top.npmax = ff.read_part('npmax@top@parallel',itriple)[0]

  if 'npslost_p' in vlistparallel:
    npslost_p = ff.read('npslost_p@parallel')
    npslost_p0 = zeros((top.nzprocs+1,top.ns+1),'l')
    npslost_p0[1:,1:] = npslost_p
  itriple = array([me,me,1,0,top.ns-1,1])
  if 'inslost@top' in vlistparallel:
    top.inslost[:] = ff.read_part('inslost@top@parallel',itriple)[0,...]
  if 'npslost@top' in vlistparallel:
    top.npslost[:] = ff.read_part('npslost@top@parallel',itriple)[0,...]

  # --- Loop over the list of all of the variables in the restart file.
  # --- Read in all of the scalars first - this ensures that all of the
  # --- integers which describe the size of dynamics arrays are read in
  # --- before the arrays, a requirement of the f90 version.
  for v in vlistparallel:
    if len(v) > 4 and v[-4]=='@':
      # --- Variable is a fortran variable
      vname = v[:-4]
      p = v[-3:]
      pkg = packageobject(p)
      pname = p+'.'+vname
      # --- The shape is used determine whether the variable is an array
      # --- or not. When the length of the shape is zero, then the
      # --- variable is a scalar.
      if len(ff.inquire_shape(v)) != 0: continue
      # --- Make sure that the variable is still valid. If not
      # --- (e.g. it has been deleted) then don't try to restore it.
      try:
        a = pkg.getvarattr(vname)
      except:
        print "Warning: There was a problem %s - it can't be found."%(pname)
        continue
#     parallelvar = re.search('parallel',a)
#     if not parallelvar: continue
      if verbose: print "reading "+p+"."+vname
      # --- Get data saved with parallel suffix.
      try:
        data = ff.read_part(vname+"@"+p+"@parallel",array([me,me,1]))[0]
        setattr(pkg,vname,data)
      except:
        print "Warning: There was a problem restoring %s"%(pname)

  # --- Get a list of all of the conductor variables
  groups = sortrestorevarsbysuffix(vlistparallel,[])
  if 'conductors%d'%me in groups.keys():
    vlistconductors = groups['conductors%d'%me]
    pyrestoreforthonobject(ff,'f3d.conductors',vlistconductors,fobjdict,
                           varsuffix,verbose,doarrays=0,
                           gpdbname='conductors%d@parallel'%me)
    pyrestoreforthonobject(ff,'f3d.conductors',vlistconductors,fobjdict,
                           varsuffix,verbose,doarrays=1,
                           gpdbname='conductors%d@parallel'%me)

  if 'pgroup%d'%me in groups.keys():
    vlistpgroup = groups['pgroup%d'%me]
    pyrestoreforthonobject(ff,'top.pgroup',vlistpgroup,fobjdict,
                           varsuffix,verbose,doarrays=0,
                           gpdbname='pgroup%d@parallel'%me)
    pyrestoreforthonobject(ff,'top.pgroup',vlistpgroup,fobjdict,
                           varsuffix,verbose,doarrays=1,
                           gpdbname='pgroup%d@parallel'%me)

  # --- Allocate any groups with parallel arrays
  gchange("LostParticles")

  # --- Now read in the arrays.
  for v in vlist:
    if len(v) > 4 and v[-4]=='@':
      # --- Variable is a fortran variable
      vname = v[:-4]
      p = v[-3:]
      pkg = packageobject(p)
      pname = p+'.'+vname
      # --- The shape is used determine whether the variable is an array
      # --- or not. When the length of the shape is zero, then the
      # --- variable is a scalar.
      if len(ff.inquire_shape(v)) == 0: continue
      if verbose: "reading "+p+"."+vname
      "reading "+p+"."+vname
      # --- Make sure that the variable is still valid. If not
      # --- (e.g. it has been deleted) then don't try to restore it.
      try:
        a = pkg.getvarattr(vname)
      except:
        print "Warning: %s no longer a variable."%(pname)
        continue
      parallelvar = re.search('parallel',a)
      if not parallelvar: continue
      # --- Many arrays need special handling. These are dealt with first.
      if p == 'top' and vname in ['ins','nps']:
        # --- These have already been restored above since they are
        # --- needed to read in the particles.
        continue
      elif p == 'top' and vname in ['inslost','npslost']:
        # --- These have already been restored above since they are
        # --- needed to read in the lost particles.
        continue
      elif p == 'top' and vname in ['xplost','yplost','zplost',
                                    'uxplost','uyplost','uzplost',
                                    'gaminvlost','tplost']:
        # --- Read in each species seperately.
        for js in range(top.ns):
          if top.npslost[js] > 0:
            ipmin = sum(sum(npslost_p0[:,0:js+1]))+sum(npslost_p0[:me+1,js+1])
            itriple = array([ipmin,ipmin+top.npslost[js]-1,1])
            lhs = getattr(pkg,vname)
            rhs = ff.read_part(v,itriple)
            lhs[top.inslost[js]-1:top.inslost[js]+top.npslost[js]-1] = rhs
      elif p == 'top' and vname == 'pidlost':
        # --- Read in each species seperately.
        for js in range(top.ns):
          if top.npslost[js] > 0:
            ipmin = sum(sum(npslost_p0[:,0:js+1]))+sum(npslost_p0[:me+1,js+1])
            itriple = array([ipmin,ipmin+top.npslost[js]-1,1,
                             0,top.npidlost-1,1])
            lhs = getattr(pkg,vname)
            rhs = ff.read_part(v,itriple)
            lhs[top.inslost[js]-1:top.inslost[js]+top.npslost[js]-1,:] = rhs
      elif p == 'w3d' and vname in ['rho']:
        itriple = array([0,w3d.nx,1,0,w3d.ny,1,
                    top.izpslave[me],top.izpslave[me]+top.nzpslave[me],1])
        setattr(pkg,vname,ff.read_part(v,itriple))
      elif p == 'w3d' and vname in ['phi']:
        itriple = array([0,w3d.nx,1,0,w3d.ny,1,
            top.izfsslave[me]-1+1,top.izfsslave[me]+top.nzfsslave[me]+2,1])
        setattr(pkg,vname,ff.read_part(v,itriple))
      elif p == 'w3d' and vname in ['zmeshlocal']:
        # --- The rest are domain decomposed Z arrays
        itriple = array([top.izfsslave[me],
                         top.izfsslave[me]+top.nzfsslave[me],1])
        setattr(pkg,vname,ff.read_part(v,itriple))
      else:
        # --- The rest are domain decomposed Z arrays
        itriple = array([top.izpslave[me],
                         top.izpslave[me]+top.nzpslave[me],1])
        setattr(pkg,vname,ff.read_part(v,itriple))

  if not lreturnff:
    ff.close()
  else:
    return ff

##############################################################################
