"""
# Module mphoto.py
#
# by: Agust Valfells
# created: Sep. 22, 2000
#
#	Last Modified: 4/20/2005
#
# Set of functions to generate tif photos and other goodies from
# 'rho' arrays in WARP, processed to look like experiment pictures.
# This module can be imported anytime, and the function take_photo
# will take a snapshot of the density at the particular time step it is
# invoked.
#
# ==================
# take_photo  ... generates tif file from rho using NumPy
# * write_info_file   ... Saves resolution and other info into text file
#
#   The following funcs are automatically invoked from within take_photo,
#   but are available for use in partial processing:
# median_filter ... Provides median filtering capability of the photo
# unfold    ... unfolds the symmetry of the 'rho' information
# ==================
"""
#===========
# GLOBALS + INITIALIZATION

yes = 1; no = 0
from warp import *
import Numeric, os, string
from tifrw import *

mphoto_version = "$Id: mphoto.py,v 1.3 2005/05/13 06:02:44 ramiak Exp $"
def mphotodoc():
  import mphoto
  print mphoto.__doc__

################ Median Filtering #######################################

def V_median(vector):
    """V_median(vector)
    """
    V = Numeric.sort(vector)
    if len(V) == 1:
        return V[len(V)/2]
    else:
        if len(V) % 2 == 0:
            return V[len(V)/2 - 1]
        else:
            return V[len(V)/2]

def median_filter(matrix,D_x = 1, D_y = 1):
    """median_filter(matrix,D_x = 1, D_y = 1)
    Applies a median filter to the 2-D array "matrix".
The pixel in the center of a moving window of size (2*D_x+1, 2*D_y+1)
is replaced with the median value of all the pixels in the window.
    """
    M, N = Numeric.shape(matrix)
    outmatrix = Numeric.zeros([M,N])*1.0
    D_x = int(Numeric.absolute(D_x))
    D_y = int(Numeric.absolute(D_y))

    if D_y > M/2:
        D_y = M/2
        print("D_y to large; decreased to D_y =", D_y)
    if D_x > N/2:
        D_x = N/2
        print("D_x to large; decreased to D_x =", D_x)

    for i in range(0,M):

        min_h = -D_y
        max_h = D_y + 1

        if i < D_y:             min_h = -i
        if i > (M - 1 - D_y):   max_h = M - i

        Ly = max_h - min_h

        for j in range(0,N):

            min_k = -D_x
            max_k = D_x + 1

            if j < D_x:         min_k = -j
            if j > (N-1-D_x):   max_k = N - j

            Lx = max_k - min_k
            Axy = Lx * Ly
            V = Numeric.zeros([Axy,])*1.0	
				
            for h in range(min_h, max_h):
                for k in range (min_k, max_k):
                    V[(h - min_h)*Lx  + k - min_k] = matrix[i+h, j+k]
            outmatrix[i,j] = V_median(V)
    return outmatrix

########################  nfold Image for Given Symmetry ################

def unfold(matrix, symmetry=0):
	"""unfold(matrix, symmetry=0)
    Takes a 2_D array "matrix" and unfolds it depending on the value of
'symmetry' (0 for none -> no unfolding; 2 for 2-fold -> unfold over y;
and 4 for 4-fold -> unfold over x and y)
	"""
	if symmetry == 4:
		M, N = Numeric.shape(matrix)
		F_1 = Numeric.zeros([M,N-1])*1.0
		F_2 = Numeric.zeros([M-1,2*N-1])*1.0
		n = range(N-1)
		for i in n:
			F_1[:,i] = matrix[:,(N-1) - i]
		F_3 = Numeric.concatenate((F_1,matrix),1)
		del(F_1)
		del(n)

		F_2 = Numeric.zeros([M-1,2*N-1])*1.0
		m = range(M-1)
		for i in m:
			F_2[i,:] = F_3[(M-1) - i,:]
		photo_array = Numeric.concatenate((F_2,F_3),0)

		del(F_2)
		del(F_3)
		del(m)

	elif symmetry == 2:
		M, N = Numeric.shape(matrix)
		F_2 = Numeric.zeros([M, 2*N-1])*1.0	                            	
       	
		for i in range(2*N-1):
			if i < N:
				F_2[:,i] = matrix[:, N-i-1]
			else:
				F_2[:,i] = matrix[:, i-N+1]


		photo_array = F_2
                      	
		del(F_2)
		
	else:   photo_array = matrix
	
	return photo_array


# ################## Save Array to Tif ##########################################
# 
# def save_tif(matrix, filename = None):
#     """ save_tif(matrix, filename = "temp.tif")
#     Saves a 2-D array "matrix" into a tif picture file.
#     """
#     S = Numeric.ravel(matrix)
#     M,N = Numeric.shape(matrix)
# 
#     min_val = float(Numeric.minimum.reduce(S))
#     max_val = float(Numeric.maximum.reduce(S))
#     if max_val <> min_val:
#         matrix = (matrix - min_val) / (max_val - min_val) *255		#Preprocessor
#     matrix = matrix.astype('b')					#Convert to binary
#     matrix = Numeric.transpose(matrix)			#Preprocess for tif-ization
# 
#     if filename is None:    filename = "temp.tif"
# 
#     tif = "P5\n#TIF version of array\n%d %d\n255\n%s" % (M, N,
#                                 Numeric.ravel(matrix).tostring())
#     f_tif = open(filename,'wb')
#     f_tif.write(tif)
#     f_tif.close()
# 
# ################## Read Array from Tif ##########################################
# 
# def read_tif(phpath):
#     """ read_tif(phpath): read tif photo speicified by phpath and
#         return as a 2-D array
#     """
#     f_tif = open(phpath,'rb')
#     tif = f_tif.read()
#     f_tif.close()
#     #
#     header, matrix = tif.split('\n255\n')
#     dims = tuple([int(s) for s in header.split('\n')[-1].split()])
#     #
#     matrix = Numeric.array(list(matrix))
#     matrix = Numeric.reshape(matrix, dims)
#     dummy  = matrix.astype('l')
#     matrix = Numeric.where(dummy<0, dummy+255, dummy)
#     return Numeric.transpose(matrix), dims


########################### Info File ################################################

def write_info_file(empty = 1):
    """ write_info_file(empty = 1)
    <Currently unimplemented ?!>
    """
    tempstring = arraytostr(top.runid) + ".info"	

    if ~empty:
        f = open(tempstring,"w")
        f.write("Information file for density photos.  Runi ID: %s.\n"+
                "Format is: normal; position; peak_coord_1; peak_coord_2;"+
                " meshsize_1; meshsize_2 \n" % arraytostr(top.runid))
        f.close()	
    else:
        f = open(tempstring,"a")
        f.write("%s; %d; %d; %d; %E; %E \n" % (normal, position, peak1, peak2,
                                                mesh1, mesh2))


################################################################

def take_photo(numphoto='loc', peakx=None, peaky=None, peakz=None,
				iz=0, iy=None, ix=None,
				Filt1=0, Filt2=None, lprofiles=no, rho=None):
	""" take_photo(numphoto='loc', peakx=None, peaky=None, peakz=None,
                    iz=0, iy=None, ix=None,
                    Filt1=0, Filt2=None, lprofiles=no, rho=None)

    Extracts beam density information (using getrho()) and saves into a
tif file, including options for unfolding and filtering.

        'numphoto' is the photo number; sequentially if numphoto = 'seq';
          and zbeam in filename (in mm) if numphoto = 'loc'.
        'peakx', 'peaky', 'peakz' determin limits at which picture needs
          be cropped
        'iz', 'iy', 'ix' (only one of which should be defined)
          determines location at which beam should be sliced.
        'Filt1' and 'Filt2' determine the size of the filter to be used in
          each of the 2 directions (default is no filtering).
        'lprofiles' determines if profiles are also to be plotted.
	"""
	runid = arraytostr(top.runid)

# --- Calculate photo name and number
	num = '00'; lfirst = yes
	if numphoto == 'seq':
		for file in os.listdir(os.curdir):
			fname = string.split(file, '.')
			try:
				if (fname[0]==runid) and (fname[2]=='tif') and (fname[1][1:] > num):
					num = fname[1][1:]
					lfirst = no
			except IndexError:  continue
		num = `(int(num)+1)`
	elif numphoto == 'loc':
		num = "%06d" % (1000.0 * top.zbeam)
	else:
		num = `numphoto`
	if len(num) < 2:    num = '0'+num

    ### Add error correction ###
	if ix <> None:      label = "x"; iz = None; iy = None
	elif iy <> None:    label = "y"; iz = None; ix = None
	else:               label = "z"; iy = None; ix = None

	filename = runid+'.'+label+num+".tif"
	#write_info_file(lfirst, runid+'.info')

# --- Determine borders
	if peakx is None: peakx=w3d.xmmax
	if peaky is None: peaky=w3d.ymmax
	if peakz is None: peakz=w3d.zmmax

	
	nxh=nint(w3d.nx*(peakx/w3d.xmmax))+1; xb=nint((w3d.nx+1-nxh)/2); xt=nint((w3d.nx+1+nxh)/2)
	nyh=nint(w3d.ny*(peaky/w3d.ymmax))+1; yb=nint((w3d.ny+1-nyh)/2); yt=nint((w3d.ny+1+nyh)/2)
	nzh=nint(w3d.nz*(peakz/w3d.zmmax))+1; zb=nint((w3d.nz+1-nzh)/2); zt=nint((w3d.nz+1+nzh)/2)

	if w3d.l2symtry:
		if iy <> None:
			iy = iy
		elif ix <> None:
			ix = ix + nint(w3d.nx/2)
	elif w3d.l4symtry:
		if iy <> None:
			iy = iy
		elif ix <> None:
			ix = ix
	else:
		if iy <> None:
			iy = iy + nint(w3d.ny/2)
		elif ix <> None:
			ix = ix + nint(w3d.nx/2)

# --- Extract slice
	if rho is None: slice = getrho(ix, iy, iz)
	else: slice = rho

	if w3d.l2symtry:        # --- Take care of symmetry
		if   iz<>None: F=unfold(slice[xb:xt, 0:nyh], 2)
		elif iy<>None: F=unfold(Numeric.transpose(slice[xb:xt, zb:zt]), 0)
		elif ix<>None: F=unfold(Numeric.transpose(slice[0:nyh, zb:zt]), 2)
	elif w3d.l4symtry:
		if   iz<>None: F=unfold(slice[0:nxh, 0:nyh], 4)
		elif iy<>None: F=unfold(Numeric.transpose(slice[0:nxh, zb:zt]), 2)
		elif ix<>None: F=unfold(Numeric.transpose(slice[0:nyh, zb:zt]), 2)
	else:
		if   iz<>None: F=unfold(slice[xb:xt, yb:yt], 0)
		elif iy<>None: F=unfold(Numeric.transpose(slice[xb:xt, zb:zt]), 0)
		elif ix<>None: F=unfold(Numeric.transpose(slice[yb:yt, zb:zt]), 0)

	if Filt1:               # --- Filtering
		if not Filt2: Filt2 = Filt1
		F = median_filter(F, Filt1, Filt2)

	save_tif(abs(F), filename)   # --- Write to disk

	if lprofiles:           # --- Profiles
		pass

