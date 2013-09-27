"""
Utility and convenience routines used in Warp
"""
from __future__ import generators # needed for yield statement for P2.2
from warp import *
import struct # needed for makefortranordered
import appendablearray
import re

warputils_version = "$Id: warputils.py,v 1.42 2011/11/04 21:39:11 grote Exp $"

def warputilsdoc():
    import warputils
    print warputils.__doc__

if sys.hexversion >= 0x03000000:
    # --- Hacky replacement of execfile which is not included in python3.
    # --- I'll probably regret having this.
    def execfile(filename, globals=None, locals=None):
        if globals is None:
            globals = sys._getframe(1).f_globals
        if locals is None:
            locals = sys._getframe(1).f_locals
        with open(filename, "r") as fh:
            exec(fh.read()+"\n", globals, locals)

# --- Make a wrapper around the Forthon doc function, which cleans up
# --- reStructuredText markup.
def doc(name,printit=1):
    import Forthon
    result = Forthon.doc(name,printit=0)
    # --- This changes function references into doc(funcname).
    result = re.sub(r':py:func:`(([~\w]*\.)*)([\w]+)(( \<[\w]+\>)*)`',r'\3',result)
    #result = re.sub(r'\A\|',r'',result)
    if printit: print result
    else: return result

# --- Convenience function modeled after the iota of basis
def iota(low,high=None,step=1):
    """
  Returns a array of sequential integers, with the low end defaulting to 1.
    """
    if high is None:
        if step > 0:
            return arange(1,low+1,step)
        else:
            return arange(low,0,step)
    else:
        if step > 0:
            return arange(low,high+1,step)
        else:
            return arange(low,high-1,step)

# --- Convenience function to do printing
def remark(s):
    """
  Same as print
    """
    print s

numpysign = sign
def sign(x,y=None):
    """
  Replicate the sign function with two arguments. If only one is
  given, return the value from the numpy sign function.
  This should really be removed.
    """
    if y is None: return numpysign(x)
    if isinstance(x,ndarray):
        result = where(greater(y,0.),abs(x),-abs(x))
        result = where(equal(y,0.),0.,result)
        return result
    else:
        if y > 0:
            return abs(x)
        elif y < 0:
            return -abs(x)
        else:
            return 0

# --- Convenience function which returns meshes filled with the coordinates.
def getmeshcoordinates(mins,dds,nns):
    """
  getmeshcoordinates(mins,dds,nns)
  Returns arrays holding the coordinates of the mesh points.
  Length of list of inputs determines number of dimensions.
    """
    nns = tuple(array(nns) + 1)
    cc = indices(nns,'d')
    for i in range(len(mins)): cc[i] = mins[i] + cc[i]*dds[i]
    clist = []
    for i in range(len(mins)): clist.append(cc[i])
    return tuple(clist)

def getmesh2d(xmin,dx,nx,ymin,dy,ny):
    """
  Returns a tuple of 2 2-d arrays holding the coordinates of the mesh points
  in two dimensions.
   - xmin,dx,nx,ymin,dy,ny: mesh description
    """
    return getmeshcoordinates([xmin,ymin],[dx,dy],[nx,ny])

def getmesh3d(xmin,dx,nx,ymin,dy,ny,zmin,dz,nz):
    """
  getmesh3d(xmin,dx,nx,ymin,dy,ny,zmin,dx,nz)
  Returns a tuple of 3 3-d arrays holding the coordinates of the mesh points
  in three dimensions.
   - xmin,dx,nx,ymin,dy,ny,zmin,dz,nz: mesh description
    """
    return getmeshcoordinates([xmin,ymin,zmin],[dx,dy,dz],[nx,ny,nz])

def arrayappend(x,a):
    """
  Appends a new element to the end of an array.
  It is not very efficient since it creates a whole new array each time.
  It is better to use a standard python list, and then convert to an array
  afterwards.
    """
    xshape = list(shape(x))
    if isinstance(a,ndarray):
        pass
    elif isinstance(a,list):
        a = array(a)
    else:
        a = array([a])
    ashape = list(shape(a))
    if len(xshape)==1 and len(ashape)==1:
        xshape[0] = xshape[0] + ashape[0]
        y = fzeros(xshape,gettypecode(x))
        y[0:xshape[0]-ashape[0]] = x
        y[xshape[0]-ashape[0]:] = a
    elif len(xshape) == len(ashape)+1 and xshape[:-1] == ashape:
        xshape[-1] = xshape[-1] + 1
        y = fzeros(xshape,gettypecode(x))
        y[...,0:-1] = x
        y[...,-1] = a
    return y

# Convenience function which returns true if variable exists
def exists(x):
    """
  Checks whether or not the variable whose name is specified exists in the
  main dictionary.
   - x: Name of variable - must be a string.
    """
    import __main__
    if x in __main__.__dict__: return true
    if x in locals(): return true
    if x in globals(): return true
    return false

def ave(x,index=0):
    """
  Returns the average of the input array
    """
    if len(shape(x)) == 0: return x
    if shape(x)[index] > 0:
        return sum(x,index)/shape(x)[index]
    else:
        return 0.

def averagezdata(qty,navg=0,nlines=100,n1=None,n2=None,istep=None,
                 includezeros=false):
    """
  Averages data over local region. It also can down select data in the other
  dimension.
    - qty: Data to be smoothed. Can be either a 1-D or 2-D array.
    - navg=0: number of data points to average over
    - nlines=100: number of lines from second dimension to choose.
    - n1=shape(qty)[0]-1:
    - n2=shape(qty)[1]-1:
    - istep=max(1,n2/nlines):
    - includezeros=false: by default, only non-zero data is averaged over.
    """
    if shape(qty)[0] < navg+navg+1: return qty
    if navg == 0 or nlines == 0: return qty
    if len(shape(qty)) == 1:
        fixqty = 1
        qty.shape = (len(qty),1)
    else:
        fixqty = 0
    if not n1: n1 = shape(qty)[0] - 1
    if not n2: n2 = shape(qty)[1] - 1
    if istep is None: istep = max(1,n2/nlines)
    hl = qty[:,::istep] + 0.
    hl[navg,:] = sum(qty[navg-navg:navg+navg+1,::istep])
    nn = 2*navg+1 + zeros(shape(hl),'l')
    if not includezeros:
        nn[navg,:] = sum(where(qty[navg-navg:navg+navg+1,::istep]==0.,0,1),0)
    for j in range(navg+1,n1-navg-1):
        hl[j,:] = hl[j-1,:] + (qty[j+navg,::istep] - qty[j-navg-1,::istep])
        nn[j,:] = nn[j-1,:] + (+ where(qty[j+navg,::istep]==0,0,1)
                               - where(qty[j-navg-1,::istep]==0,0,1))
    nn = where(nn==0,1,nn)
    hl = where(qty[:,::istep]==0.,0.,hl)
    hl[navg:n1-navg-1,:] = hl[navg:n1-navg-1,:]/nn[navg:n1-navg-1,:]
    if fixqty: qty.shape = (shape(qty)[0],)
    if n2 > 1: return hl
    else: return hl[:,0]

# --- Returns the max of the multiarray
def maxnd(x,defaultval=-1.e36):
    """
  Return the max element of an array of any dimension
    - x: a scalar or any sequence that can be converted to an array
    - defaultval=-1.e36: value returned if the input has zero length
    """
    if isinstance(x,ndarray) and x.size > 0: return x.max()
    xtemp = reshape(x,tuple([product(array(shape(x)))]))
    if len(xtemp) == 0: return defaultval
    return max(xtemp)
# --- Returns the min of the multiarray
def minnd(x,defaultval=+1.e36):
    """
  Return the min element of an array of any dimension
    - x: a scalar or any sequence that can be converted to an array
    - defaultval=+1.e36: value returned if the input has zero length
    """
    if isinstance(x,ndarray) and x.size > 0: return x.min()
    xtemp = reshape(x,tuple([product(array(shape(x)))]))
    if len(xtemp) == 0: return defaultval
    return min(xtemp)
# --- Returns the sum of the multiarray
def sumnd(x,defaultval=0.):
    """
  Return the total sum of an array of any dimension
    - x: a scalar or any sequence that can be converted to an array
    - defaultval=0.: value returned if the input has zero length
    """
    if isinstance(x,ndarray) and x.size > 0: return x.sum()
    xtemp = reshape(x,tuple([product(array(shape(x)))]))
    if len(xtemp) == 0: return defaultval
    return sum(xtemp)
# --- Returns the sum of the multiarray
def avend(x,defaultval=0.):
    """
  Return the average of an array of any dimension
    - x: a scalar or any sequence that can be converted to an array
    - defaultval=0.: value returned if the input has zero length
    """
    if isinstance(x,ndarray) and x.size > 0: return x.sum()/x.size
    xtemp = reshape(x,tuple([product(array(shape(x)))]))
    if len(xtemp) == 0: return defaultval
    return sum(xtemp)/len(xtemp)

def span(lo, hi, num):
    """
  Returns an array of num equally spaced numbers starting with lo and
  ending with hi.
    """
    if num == 1: return array([lo])
    return lo + (hi - lo)*arange(num)/(num-1.)

def makefortranordered(x):
    """
  Given an array, returns the same data but with fortran ordering.
  If the array already has the correct ordering, the array is just
  returned as is. Otherwise, a new array is created and the data copied.
    """
    assert isinstance(x,ndarray),"Input value must be an array."
    # --- Pick any package, since all have the getstrides method
    pkg = packageobject(getcurrpkg())
    # --- An array is fortran ordered if the strides are increasing,
    # --- starting at the elemental size. This means that the first
    # --- index varies the fastest in memory.
    strides = pkg.getstrides(x)
    fordered = 1
    ss = struct.calcsize(gettypecode(x))
    for id in range(len(x.shape)):
        if strides[id] != ss: fordered = 0
        ss = ss*x.shape[id]

    if fordered:
        # --- No change needed.
        return x
    else:
        xf = fzeros(x.shape,gettypecode(x))
        xf[...] = x
        return xf

# Gets next available filename with the format 'root.nnn.suffix'.
def getnextfilename(root,suffix):
    """
  Finds next available file name in a numeric sequence. The filename is assumed
  to have the format root.XXX.suffix, where XXX is three digits.
    """
    dir = ' '.join(os.listdir('.'))
    i = 0
    name = root+('.%03d.'%i)+suffix
    # --- escape adds a backslash to all non-alphanumeric characters. This is
    # --- needed in case special re search characters are used as part of the
    # --- file name, for example '+' or '.'.
    while re.search(re.escape(name),dir):
        i = i + 1
        name = root+('.%03d.'%i)+suffix
    return name

# --- Another Example profiler - this prints out the function name when the
# --- function is called and when it returns. Indentation is used to signify
# --- the call level.
# --- To use, execute the command sys.setprofile(warpprofile).
def warpprofilesample(frame,event,arg):
    try:    warpprofile.level
    except: warpprofile.level = 0
    if event == 'return': warpprofile.level = warpprofile.level - 1
    print "%s %s %s"%(warpprofile.level*'  ',event,frame.f_code.co_name)
    if event == 'call': warpprofile.level = warpprofile.level + 1

# --- Create enumerate function, which is defined in python2.3 but not 2.2
try:
    enumerate
except:
    def enumerate(ll):
        tt = []
        for i in range(len(ll)):
            tt.append((i,ll[i]))
        return tt

# --- Convenience function to read in data from a text file
def getdatafromtextfile(filename,nskip=0,dims=None,dtype='d',fortranordering=True,
                        converter=float,mode='r',get_header=False,l_checkdims=False,
                        separator=None,readuntilerror=False):
    """
  Reads data in from a text file. The data is assumed to be laid out on a
  logically Cartesian mesh.
   - filename: must be supplied

   - converter=float: Function which converts the strings into numbers. This should
                      be the type of result desired. It should only be float or int,
                      unless you know what you are doing.
   - dims=[]: The size of each dimension.
              The last dimension can be None - it will be determined
              automatically based on the amount of data in the file.
              If not supplied, a 2-D array will be returned, with the first dimension
              determined by the number of values in the first line of data.
   - fortranordering=True: when true, the data will be in fortran ordering, where
                           the index that varies that fastest in the file will be
                           the first index. Otherwise C ordering will be used.
   - get_header=False: When true, the header is returned as a list of strings along with
                       the data, in the tuple (data,header).
   - if l_checkdims=False: when true, check that the dimensions of the array in the ascii
                           file match the ones provided in dims.
   - mode='r': mode used to open the file.
   - nskip=0: number of lines at the beginning of the file to skip
              e.g. lines with file info and comments.
              If get_header is true, the skipped lines will be returned as the header.
   - separator=' ': separating character between data values. If not specified,
                    either a space or a comma will be used depending on what is in
                    the file.
   - readuntilerror=False: If True, when reading in open ended data, an error in the
                           number of data points in a line or in the data conversion is
                           taken to be a flag signifying the end of the data.

  Here's an example data file called 'testdata'::

    this line is skipped
    1 0.0379
    2 0.0583
    3 0.0768
    4 0.1201

  This can be read in with::

    >>> dd = getdatafromtextfile('testdata',nskip=1,dims=[2,4])

  to produce an array of shape (2,4) that looks like::

    >>> print dd
    [[ 1.    , 2.    , 3.    , 4.    ,]
     [ 0.0379, 0.0583, 0.0768, 0.1201,]]

  The second dimension of dims can be None, meaning that all of the data
  will be read in and the size of that dimension will be determined
  automatically.::

    >>> dd = getdatafromtextfile('testdata',nskip=1,dims=[2,None])
    >>> print dd
    [[ 1.    , 2.    , 3.    , 4.    ,]
     [ 0.0379, 0.0583, 0.0768, 0.1201,]]

    """
    ff = open(filename,mode)

    # --- Skip the number of lines at the top of the file as specified
    header = []
    # --- This attempts to scan through the top lines of a file and skip the
    # --- lines not containing numbers.
    # --- This is probably too complicated to do, since the number of possible
    # --- formats for numbers of very large.
    #if nskip is None:
    #  notnumbers = re.compile('[^0-9edED\.\-+ ,\n%s]'%seperator)
    #  while 1:
    #    lastpos = ff.tell()
    #    line = ff.readline()
    #    if notnumbers.search(line):
    #      header.append(line)
    #    else:
    #      ff.seek(lastpos)
    #      break
    #else:
    for i in range(nskip):
        header.append(ff.readline())

    # --- Read the first line to interrogate it if needed.
    if dims is None or separator is None:
        lastpos = ff.tell()
        firstline = ff.readline()
        ff.seek(lastpos)

    # --- Check which separator to use. Look for commas, and if none found,
    # --- use a space.
    if separator is None:
        if re.search(',',firstline):
            separator = ','

    # --- If dims is not supplied, set it so that a 2-D array is returned.
    # --- The first dimension is set to the number of values on the first line
    # --- of the file.
    if dims is None:
        nvalues = len(firstline.split(separator))
        dims = [nvalues,None]

    # --- If last dimension is given, calculate the total amount of data to
    # --- be read in to give an exit condition for the loop below.
    # --- This allows less than the whole file to be read in.
    if dims[-1] is not None:
        ntot = product(dims)
    else:
        ntot = None

    # --- Loop over the file, reading in one line at a time.
    # --- Each whole line is put into data at once.
    # --- For the conversion of the strings into numbers, it is faster to use
    # --- float or int rather than eval.
    data = appendablearray.AppendableArray(typecode=dtype)
    lines=ff.readlines()
    if l_checkdims:assert len(lines)==dims[1],"ERROR reading data from "+filename+": dimensions 2nd axis are incompatible. %g passed as argument, %g found in file."%(dims[1],len(lines))
    for line in lines:
        words = line.split(separator)
        if l_checkdims:assert len(words)==dims[0],"ERROR reading data from "+filename+": dimensions 1st axis are incompatible. %g passed as argument, %g found in file."%(dims[0],len(words))
        if readuntilerror and len(words) != dims[0]:
            break
        try:
            dataline = map(converter,words)
        except ValueError:
            if readuntilerror:
                break
            else:
                raise sys.exc_info()
        else:
            data.append(dataline)
        if ntot is not None and len(data) >= ntot: break

    # --- Get the data out of the appendable array.
    data = data[...]

    # --- If last dimension is None, calculate it from the amount of data.
    if dims[-1] is None:
        ndata = len(data)
        ndims = product(dims[:-1])
        nlast = int(ndata/ndims)
        assert nlast*ndims == ndata,"Amount of data does not conform to dims"
        dims[-1] = nlast

    if fortranordering:
        # --- reverse the order of the dims
        dims.reverse()

    # --- Set array to have proper shape.
    data.shape = tuple(dims)
    if fortranordering:
        data = transpose(data)
    else:
        data = transpose(data,[len(dims)-1]+range(len(dims)-1))

    if get_header:
        return data,header
    else:
        return data

class RandomStream(object):
    """Used to create independent streams of random numbers. ::

    >>> r1 = RandomStream()
    >>> r2 = RandomStream()
    >>> r1.seed([1,2,3])
    >>> r2.seed([1,2,3])
    >>> print r1.random()
    0.609861272287
    >>> print r2.random()
    0.609861272287

    """
    def __init__(self):
        # --- Save the initial random number state
        self._state = random.get_state()
    def __getattr__(self,name):
        # --- For all attribute access, save the attribute name and return the
        # --- _wrapper method.
        self._savedname = name
        return self._wrapper
    def _wrapper(self,*args,**kw):
        # --- This method calls the function from the random module with the saved
        # --- name, though first it switches to the stream's state, and then
        # --- switches back afterward.
        savedstate = random.get_state()
        random.set_state(self._state)
        result = getattr(random,self._savedname)(*args,**kw)
        self._state = random.get_state()
        random.set_state(savedstate)
        return result
    def __getnewargs__(self):
        # --- This method needs to be given explicitly, otherwise it will go
        # --- through getattr and _wrapper, causing problems.
        return tuple()
    def __getstate__(self):
        # --- create the thing to be pickled
        return {'_state':self._state}
    def __setstate__(self,dict):
        # --- set the state when being unpickled
        self._state = dict['_state']

def asciipart(s,name='particles.dat',format='%10.3e',mode='w'):
    """
  Writes the particles data from a species to a text file.
   - s: species instance or species number
   - name='particles.dat': file name
   - format='%10.3e': format used to write the data
   - mode='w': mode used when opening the file
    """
    f = open(name,mode)
    if isinstance(s,Species):
        n=s.getn()
        x=s.getx()
        y=s.gety()
        z=s.getz()
        ux=s.getux()
        uy=s.getuy()
        uz=s.getuz()
    else:
        n=getn(js=s)
        x=getx(js=s)
        y=gety(js=s)
        z=getz(js=s)
        ux=getux(js=s)
        uy=getuy(js=s)
        uz=getuz(js=s)
    npad = len(format%1.)
    if mode=='w':
        f.write('x'+' '*npad)
        f.write('y'+' '*npad)
        f.write('z'+' '*npad)
        f.write('ux'+' '*(npad-1))
        f.write('uy'+' '*(npad-1))
        f.write('uz'+' '*(npad-1))
        f.write('\n')
    for i in range(n):
        f.write(format%x[i]+' ')
        f.write(format%y[i]+' ')
        f.write(format%z[i]+' ')
        f.write(format%ux[i]+' ')
        f.write(format%uy[i]+' ')
        f.write(format%uz[i]+' ')
        f.write('\n')
    f.close()

def bisection(f,a,b,f0=0.,xtol=1.e-10,maxiter=100):
    """
  Given a function that is monotonic on an interval [a,b], finds x such that
  f(x)=f0 using the bisection method.
    """
    if f(a)-f0 <= 0.:
        lo = a
        hi = b
    else:
        lo = b
        hi = a

    mid = lo + (hi-lo)/2
    iter = 0
    while (abs(hi-lo)>xtol) or ((mid != lo) and (mid != hi)) and iter < maxiter:
        iter += 1
        if f(mid)-f0 <= 0.:
            lo = mid
        else:
            hi = mid
        mid = lo + (hi-lo)/2

    return mid

def getruntimememory(diag='vsz'):
    """Return the memory being used during the run time. It returns the
virtual size in bytes."""
    try:
        memstring = os.popen('ps -p %d -o %s=""'%(os.getpid(),diag)).read().strip()
        memlist = memstring.split('\n')
        if len(memlist) > 1:
                # --- Might by busybox style ouput. The results needs to be
                # --- searched through to find the data for the process.
            for pp in memlist:
                dd = pp.strip().split()
                if len(dd) > 0 and dd[0] == repr(os.getpid()):
                    memstring = dd[2]
                    break
        if memstring[-1] in ['M','m']:
            m = float(memstring[:-1])*1.e6
        elif memstring[-1] in ['K','k']:
            m = float(memstring[:-1])*1.e3
        else:
            m = float(memstring)*1.e3
    except:
        m = 0
    return m

def getallobjectsizes(minsize=1.e20,withpackages=1,withwarpglobals=0):
    """Makes a printout of all of the user created python objects and their
  sizes, in bytes.
   - minsize=1: only prints objects with a size greater than minsize
   - withpackages=1: when true, also includes the sizes of the fortran packages
    """
    import __main__
    import warp
    import types
    dd = {}
    for k,v in __main__.__dict__.iteritems():
        if not withwarpglobals and k in warp.initial_global_dict_keys: continue
        if isinstance(v,types.ModuleType): continue
        if k == 'controllerfunctioncontainer': continue
        if k == 'registeredsolverscontainer': continue
        dd[k] = v
    totals = 0
    for k,v in [[kk,__main__.__dict__[kk]] for kk in package()]:
        if k in dd: del dd[k]
        s = getobjectsize(v)*8
        totals += s
        if s > minsize:
            print k,s,'bytes'
    for k,v in dd.iteritems():
        s = getobjectsize(v)*8
        totals += s
        if s > minsize:
            print k,s,'bytes'
    return totals
