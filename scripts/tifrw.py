"""
# Module tifrw.py
#
# by: Rami Kishek and A. Valfells
# created: April 20, 2005
#
#	Last Modified: 4/20/2005
#
# Set of functions to read and write tif photos from Numpy arrays in Python
# save_tif ...  Saves photo array into 'tif' file
# read_tif ...  Reads photo array from 'tif' file
# ==================
"""
# GLOBALS + INITIALIZATION

yes = 1; no = 0
import Numeric

tifrw_version = "$Id: tifrw.py,v 1.2 2008/01/04 00:10:47 dave Exp $"
def tifrwdoc():
  import tifrw
  print tifrw.__doc__

################## Save Array to Tif ##########################################

def save_tif(matrix, filename = None):
    """ save_tif(matrix, filename = "temp.tif")
    Saves a 2-D array "matrix" into a tif picture file.
    """
    S = Numeric.ravel(matrix)
    M,N = Numeric.shape(matrix)

    min_val = float(Numeric.minimum.reduce(S))
    max_val = float(Numeric.maximum.reduce(S))
    if max_val <> min_val:
        matrix = (matrix - min_val) / (max_val - min_val) *255		#Preprocessor
    matrix = matrix.astype(ubyte)					#Convert to binary
    matrix = Numeric.transpose(matrix)			#Preprocess for tif-ization

    if filename is None:    filename = "temp.tif"

    tif = "P5\n#TIF version of array\n%d %d\n255\n%s" % (M, N,
                                Numeric.ravel(matrix).tostring())
    f_tif = open(filename,'wb')
    f_tif.write(tif)
    f_tif.close()

################## Read Array from Tif ##########################################

def read_tif(phpath):
    """ read_tif(phpath): read tif photo speicified by phpath and
        return as a 2-D array
    """
    f_tif = open(phpath,'rb')
    tif = f_tif.read()
    f_tif.close()
    #
    header, matrix = tif.split('\n255\n')
    dims = tuple([int(s) for s in header.split('\n')[-1].split()])
    #
    matrix = Numeric.array(list(matrix))
    matrix = Numeric.reshape(matrix, dims)
    dummy  = matrix.astype('l')
    matrix = Numeric.where(dummy<0, dummy+255, dummy)
    return Numeric.transpose(matrix), dims


