import string as st
#from numpy  import *
#from gist import *
from warp import *
try:
  import Image as Im
except:
  pass
  
def getpalhrgb(pname):
    f = file(pname+'.gp','r')
    lines = f.readlines()
    f.close()
    i = 0
    line = lines[i]
    words = st.split(line)
    while len(words)==0:
      i += 1
      line = lines[i]
      words = st.split(line)
    while st.find(line,'=')>=0 or st.find(line,'#')>=0:
        i += 1
        line = lines[i]
        words = st.split(line)
        while len(words)==0:
          i += 1
          line = lines[i]
          words = st.split(line)
    h = lines[:i]
    r = []
    g = []
    b = []
    for line in lines[i:]:
        words = st.split(line)
        if len(words)==3:
          r += [int(words[0])]
          g += [int(words[1])]
          b += [int(words[2])]
    return h,array(r),array(g),array(b)

def pltpalrgb(pname):
    h,r,g,b = getpalhrgb(pname)
    plg(r,color=red)
    plg(g,color=green)
    plg(b,color=blue)
    plg(r+g+b)
    ptitles('palette '+pname+'.gp','n','value')
    
    
def gfpmaker(pname,br=None):
    h,r,g,b = getpalhrgb(pname)
    f = file(pname+'_gfp.gp','w')
    f.writelines(h)
    n = shape(r)[0]
    s = r+g+b
    r = where(s==0,1,r)
    g = where(s==0,1,g)
    b = where(s==0,1,b)
    s = where(s==0,3,s)
    if br is None:
      t = 3*arange(n,typecode=Float)
    else:
      br[1][0]  = 0
      br[1][-1] = n-1
      t = array([])
      for i in range(shape(br[1])[0]-1):
        dy = float(br[0][i+1]-br[0][i])
        dx = float(br[1][i+1]-br[1][i])
        if dy==0.:
          t=concatenate([t,br[0][i]+zeros(br[1][i+1]-br[1][i],Float)])
        else:
          t=concatenate([t,arange(float(br[0][i]),float(br[0][i+1])-0.5*dy/dx,dy/dx,Float)])
      t = concatenate([t,array([float(br[0][-1])])])
    x = t/s
    rn = nint(r*x)
    gn = nint(g*x)
    bn = t-rn-gn
    bn = where(bn>=0,bn,0)
    fn = 255./float(max(max(rn),max(gn),max(bn)))
    rn = nint(rn*fn)
    gn = nint(gn*fn)
    bn = aint(t*fn)-rn-gn
    bn = where(bn>=0,bn,0)
    for i in range(n):
      f.write('%8d%8d%8d\n'%(rn[i],gn[i],bn[i]))
    f.close()

def image2rgb(filename):
  im = Im.open(filename)
  nx,ny = im.size
  rgbh = array(list(im.getdata()))
  rgbh.resize((ny,nx,4,))
  return ny,nx,rgbh[:,:,0],rgbh[:,:,1],rgbh[:,:,2]
  
def rgb2palette(r,g,b,filename='newpalette',header=None,l_inverse=0,l_plot=0):
    f = file(filename+'.gp','w')
    f.write('# Palette produced by image2rgb\n')
    if header is not None:f.write(header+'\n')
    ncolors = shape(r)[0]
    f.write('ncolors= %g\n'%ncolors)
    f.write('#  r   g   b\n')
    if not l_inverse:
     for i in range(ncolors):
      f.write('%8d%8d%8d\n'%(r[i],g[i],b[i]))
     if l_plot:
       plg(r,color='red')
       plg(g,color='green')
       plg(b,color='blue')
    else:
     for i in range(ncolors-1,-1,-1):
      f.write('%8d%8d%8d\n'%(r[i],g[i],b[i]))
     if l_plot:
       plg(r[::-1],color='red')
       plg(g[::-1],color='green')
       plg(b[::-1],color='blue')
    f.close()
      
